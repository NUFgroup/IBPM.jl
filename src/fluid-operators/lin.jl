"""
Linear operations associated with the governing flow equations

Note: all matrices in this file are "static" and do not vary in time. They may
therefore be pre-computed as a pre-processing step.
"""

"""
    rot!( Γ, q, grid )

Transpose of discrete curl (R matrix)

Γ = rot( q )
"""
function rot!( Γ, q, grid )
    u, v, ω = grid.u_idx, grid.v_idx, grid.ω_idx
    i = 1:grid.nx-1; j = (1:grid.ny-1)'
    @. Γ[ω(i, j)] = q[v(i+1, j+1)] - q[v(i, j+1)] - q[u(i+1, j+1)] + q[u(i+1, j)]
end

"""
    curl!( q, ψ, grid )

Discrete curl operator

q = curl( ψ )
"""
function curl!( q, ψ, grid::UniformGrid )
   nx, ny = grid.nx, grid.ny
   u, v, ω = grid.u_idx, grid.v_idx, grid.ω_idx  # Indices to x-flux, y-flux, vorticity/streamfunction

   # x-fluxes
   i=2:nx; j=(2:ny-1)'
   @. q[u(i, j)] = ψ[ω(i-1,j)] - ψ[ω(i-1, j-1)]
   j=1;  @. q[u(i, j)] =  ψ[ω(i-1, j)]          # Top boundary
   j=ny; @. q[u(i, j)] = -ψ[ω(i-1, j-1)]        # Bottom boundary

   # y-fluxes
   i=2:nx-1;  j=(2:ny)'
   @. q[v(i,j)] = ψ[ω(i-1,j-1)] - ψ[ω(i,j-1)]
   i=1;  @. q[v(i,j)] = -ψ[ω(i, j-1)]           # Left boundary
   i=nx; @. q[v(i,j)] =  ψ[ω(i-1, j-1)]         # Right boundary
end


"""
    vort2flux!( ψ, q, Γ, model::IBModel{UniformGrid, <:Body} )
"""
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

"""
Construct LinearMap to solve inverse problem with Laplacian
"""
function get_lap_inv( grid::T,
                      Λ::AbstractArray,
                      dst_plan::Tuple{Any, Array{Float64, 2}}) where T <: Grid
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

"""
Compute the matrix (as a LinearMap) that represents the modified Poisson
operator (I + dt/2 * Beta * RC) arising from the implicit treatment of the
Laplacian. A system involving this matrix is solved to compute a trial
circulation that doesn't satisfy the BCs, and then again to use the surface
stresses to update the trial circulation so that it satisfies the BCs
"""
function get_Ainv(model::IBModel, dt::Float64, h::Float64)
    """
    Construct LinearMap to solve
        (I + dt/2 * Beta * RC) * x = b for x
    where Beta = 1/(Re * h^2)

    Solve by transforming to and from Fourier space and scaling by evals
    """
    Λ̃ = 1 .+ model.mats.Λ * dt/( 2*model.Re*h^2 );  # Implicit eigenvalues
    return get_lap_inv(model.grid, Λ̃, model.mats.dst_plan)
end

function get_A(model::IBModel, dt::Float64, h::Float64)
    Λ̃ = 1 .- model.mats.Λ * dt/( 2*model.Re*h^2 );  # Explicit eigenvalues
    return get_lap_inv(model.grid, 1 ./ Λ̃, model.mats.dst_plan)
end

function get_AB(model::IBModel{UniformGrid, <:Body}, dt::Float64)
    A = get_A(model, dt, model.grid.h)
    Ainv = get_Ainv(model, dt, model.grid.h)
    Binv = get_B(model, Ainv)
    return A, Ainv, Binv
end

function get_B(model::IBModel{<:Grid, RigidBody{T}} where T <: Motion, Ainv::LinearMap)
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
    Γ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for streamfunction
    q = zeros(model.grid.nq, model.grid.mg)    # Working array for velocity flux

    for j = 1 : nftot
        e .*= 0.0
        e[j] = 1.0;

        B_times!( b, e, Ainv, model, Γ, ψ, q );
        B[:, j] = b
    end

    return inv(B)
    #return cholesky( 0.5*(B + B'))
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
    vort2flux!( ψ, q, Γ, model );

    #--Interpolate onto the body
    @views mul!(x, E, q[:, 1])
end

function get_AB(model::IBModel{MultiGrid, <:Body}, dt::Float64)
    hc = [model.grid.h * 2^(lev-1) for lev=1:model.grid.mg]
    A = [get_A(model, dt, h) for h in hc]
    Ainv = [get_Ainv(model, dt, h) for h in hc]
    Binv = get_B(model, Ainv[1])
    return A, Ainv, Binv
end
