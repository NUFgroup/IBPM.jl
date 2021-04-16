"""
Linear operations associated with the governing flow equations

Note: all matrices in this file are "static" and do not vary in time. They may
therefore be pre-computed as a pre-processing step.
"""

"""
Transpose of discrete curl (R matrix)

Γ = rot( q )
"""
function rot!( Γ, q, grid )
    for j=1:grid.ny-1
        for i=1:grid.nx-1
            Γ[grid.ω_idx(i, j)] = q[grid.v_idx(i+1, j+1)] - q[grid.v_idx(i, j+1)] - q[grid.u_idx(i+1, j+1)] + q[grid.u_idx(i+1, j)]
        end
    end
end


"""
!***************************************************************!
!*   returns curl(x) given x and the values of s'fun on bdy    *!
!***************************************************************!
"""
function curl!( q, ψ, grid )
   nx, ny = grid.nx, grid.ny
   for j=2:ny-1
      for i=2:nx
         q[grid.u_idx(i, j)] = ψ[grid.ω_idx(i-1, j)] - ψ[grid.ω_idx(i-1, j-1)]
      end
   end

   for i=2:nx
      q[grid.u_idx(i, 1)] = ψ[grid.ω_idx(i-1, 1)]
      q[grid.u_idx(i, ny)] = ψ[grid.ω_idx(i-1, ny-1)]
   end

   for j=2:ny
       q[grid.v_idx(1,j)] = -ψ[grid.ω_idx(1, j-1)]
       for i=2:nx-1
           q[grid.v_idx(i,j)] = ψ[grid.ω_idx(i-1,j-1)] - ψ[grid.ω_idx(i,j-1)]
       end
       q[grid.v_idx(nx,j)] = ψ[grid.ω_idx(nx-1,j-1)]
    end
    return q
end

function vort2flux!( ψ, q, Γ, model::IBModel{UniformGrid, <:Body} )
    mul!(ψ, model.mats.Δinv, Γ)     # Solve Poisson problem for streamfunction
    mul!(q, model.mats.C, ψ)          # q = ∇ x ψ
end

"""
Solving linear systems involving C^TC comes up in multiple places:
    -backing out vel flux from circulation, done multiple times in a time step
    -Solving a modified Poisson system (I + dt/2 * Beta * RC) that arises from
        the implicit treatment of the Laplacian

Both of these systems can be solved efficiently (on uniform grids or nested
    uniform grids) via FFTs. The three functions below enable the use of FFTs
    to solve these linear systems (Lap_eigs gives the eigenvalues of the matrix
    and Λinv_fn! provides a diagonalization via sine transforms),
    and to save the action of the inverse operator efficiently as a LinearMap
    (get_RCinv)
"""
function lap_eigs( grid::T ) where T <: Grid
    # eigenvalues of RC (negative of the evals of the 5point stencil Lap)
    nx = grid.nx; ny = grid.ny;
    Λ = -2*( cos.( π*(1:(nx-1))/nx ) .+ cos.( π*(1:(ny-1))/ny )' .- 2);
    return Λ
end

function get_Ainv(grid, dst_plan, Re, dt::Float64; lev::Int=1)
    """
    Construct LinearMap to solve
        (I + dt/2 * Beta * RC) * x = b for x
    where Beta = 1/(Re * h^2)

    Compare to get_RCinv in "ib_mats.jl"
    """

    hc = grid.h * 2^( lev - 1);  # Grid size at this level
    # Solve by transforming to and from Fourier space and scaling by evals
    Λ = lap_eigs( grid )
    Λ̃ = 1 .+ Λ * dt/( 2 * Re * hc^2 );

    # give output in same size as input b (before being reshaped)
    return get_lap_inv(grid, Λ̃, dst_plan)
end



function get_lap_inv( grid::T,
                      Λ::AbstractArray,
                      dst_plan::Tuple{Any, Array{Float64, 2}}) where T <: Grid
    """
    Construct LinearMap to solve inverse problem with Laplacian
    """
    nx, ny = grid.nx, grid.ny
    # give output in same size as input b (before being reshaped)
    return LinearMap(grid.nΓ; issymmetric=true, ismutating=true) do x, b
        # reshape for inversion in fourier space
        b = reshape( b, nx-1, ny-1)
        x = reshape( x, nx-1, ny-1)
        dst_inv!(x, b, Λ, dst_plan; scale=1/(4*nx*ny)); # Include scale to make fwd/inv transforms equal
        x = reshape( x, (nx-1)*(ny-1), 1 )
    end
end


function get_AB(model::IBModel, dt::Float64)
    A, Ainv = get_A(model, dt)
    Binv = get_B(model, Ainv)
    return A, Ainv, Binv
end

"""
Compute the matrix (as a LinearMap) that represents the modified Poisson
operator (I + dt/2 * Beta * RC) arising from the implicit treatment of the
Laplacian. A system involving this matrix is solved to compute a trial
circulation that doesn't satisfy the BCs, and then again to use the surface
stresses to update the trial circulation so that it satisfies the BCs
"""
function get_A(model::IBModel{UniformGrid, <:Body}, dt::Float64)
    Δ = (model.mats.C'*model.mats.C)/model.Re
    A = I - (dt/2 / (model.grid.h^2))*Δ
    Ainv = get_Ainv(model, dt)
    return A, Ainv
end

function get_A(model::IBModel{MultiGrid, <:Body}, dt::Float64)
    hc = [model.grid.h * 2^(lev-1) for lev=1:model.grid.mg]
    Δ = (model.mats.C'*model.mats.C)/model.Re
    A = [I - (dt/2 / (hc[lev]^2))*Δ for lev=1:model.grid.mg]
    Ainv = [get_Ainv(model, dt; lev=lev) for lev=1:model.grid.mg]
    return A, Ainv
end

function get_B(model::IBModel{UniformGrid, RigidBody{T}} where T <: Motion, Ainv::LinearMap)
    """
    Precompute 'B' matrix by evaluating mat-vec products for unit vectors

    This is a big speedup when the interpolation operator E isn't going to
    change (no FSI, for instance)
    """
    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    # need to build and store surface stress matrix and its inverse if at first time step
    B = zeros( nftot, nftot );
    # Pre-allocate arrays
    e = zeros( nftot );         # Unit vector

    # TODO: Alternative... could create a dummy state to operate on here
    b = zeros( nftot );         # Working array
    Γ = zeros(model.grid.nΓ, 1)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, 1)    # Working array for streamfunction
    q = zeros(model.grid.nq, 1)    # Working array for velocity flux

    for j = 1 : nftot
        e .*= 0.0
        e[j] = 1.0;

        B_times!( b, e, Ainv, model, Γ, ψ, q );
        B[:, j] = b
    end
    #Binv = inv(B)
    #println(B[1, 50:52]*model.grid.h)

    return B
end

function B_times!(x::AbstractArray,
                  z::AbstractArray,
                  Ainv::LinearMap,
                  model::IBModel{UniformGrid, RigidBody{T}} where T<:Motion,
                  Γ::Array{Float64, 2},
                  ψ::Array{Float64, 2},
                  q::Array{Float64, 2})
    """
    %Performs one matrix multiply of B*z, where B is the matrix used to solve
    %for the surface stresses that enforce the no-slip boundary condition.
    % (B arises from an LU factorization of the full system)
    Note ψ is just a dummy work array for circ2_st_vflx
    """
    E, C = model.mats.E, model.mats.C

    # Get circulation from surface stress Γ = Ainv * (EC)' * z
    #  We don't include BCs for Ainv because ET*z is compact
    mul!(Γ, Ainv, (E*C)'*z)  # Γ = ∇ x (E'*fb)

    #-- get vel flux q from circulation
    vort2flux!( ψ, q, Γ, model );  # THIS SHOULD BE THE RIGHT ONE

    #--Interpolate onto the body
    @views mul!(x, E, q[:, 1])
end
