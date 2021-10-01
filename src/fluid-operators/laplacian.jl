"""
Solving linear systems involving the Laplacian C^TC comes up in multiple places:
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
                      Λ::Array{Float64, 2},
                      dst_plan::Tuple{Any, Array{Float64, 2}}) where T <: Grid
    nx, ny = grid.nx, grid.ny

    # To avoid alignment issues with FFT plan, copy data to/from a working array
    #notaligned(x) = mod(UInt(pointer(x)), 16) == 0  # Test alignment
    b_temp = zeros(Float64, grid.nx-1, grid.ny-1)  # Input
    x_temp = zeros(Float64, grid.nx-1, grid.ny-1)  # Output

    # Solve Ax = b for x
    return LinearMap(grid.nΓ; issymmetric=true, ismutating=true) do x, b
        # reshape for inversion in fourier space
        b_temp .= reshape( b, grid.nx-1, grid.ny-1)
        dst_inv!(x_temp, b_temp, Λ, dst_plan; scale=1/(4*nx*ny)); # Include scale to make fwd/inv transforms equal
        x .= reshape( x_temp, size(x) )
        return nothing
    end

    """# WITH FFTW.UNALIGNED
    return LinearMap(grid.nΓ; issymmetric=true, ismutating=true) do x, b
        # reshape for inversion in fourier space
        b = reshape( b, grid.nx-1, grid.ny-1)
        x = reshape( x, grid.nx-1, grid.ny-1)
        dst_inv!(x, b_temp, Λ, dst_plan; scale=1/(4*nx*ny)); # Include scale to make fwd/inv transforms equal
        x = reshape( x, grid.nΓ, 1 )
        return nothing
    end"""
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

"""function get_A(model::IBModel, dt::Float64, h::Float64)
    C = model.mats.C
    return I - (dt/(2*model.Re*h^2))*(C'C)
end"""

function get_AB(model::IBModel{MultiGrid, <:Body}, dt::Float64)
    hc = [model.grid.h * 2^(lev-1) for lev=1:model.grid.mg]
    A = [get_A(model, dt, h) for h in hc]
    Ainv = [get_Ainv(model, dt, h) for h in hc]
    Binv = get_Binv(model, Ainv[1])  # Only need this on the finest grid
    return A, Ainv, Binv
end


"""
Precompute 'Binv' matrix by evaluating mat-vec products for unit vectors

This is a big speedup when the interpolation operator E isn't going to
change (no FSI, for instance)
"""
function get_Binv(model::IBModel{<:Grid, RigidBody{T}} where T <: Motion, Ainv::LinearMap)
    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    B = zeros( nftot, nftot );
    # Pre-allocate arrays
    e = zeros( nftot );         # Unit vector

    b = zeros( nftot );         # Working array

    Γ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for streamfunction
    q = zeros(model.grid.nq, model.grid.mg)    # Working array for velocity flux

    for j = 1 : nftot
        e .*= 0.0; e[j] = 1.0;  # Construct unit vector
        B_times!( b, e, Ainv, model, Γ, ψ, q );
        B[:, j] = b
    end
    #sleep(100)
    return inv(B)
end



function get_Binv(model::IBModel{MultiGrid, RigidBody{MotionFunction}}, Ainv::LinearMap)
    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    # TODO: Alternative... could create a dummy state to operate on here
    Γ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for streamfunction
    q = zeros(model.grid.nq, model.grid.mg)    # Working array for velocity flux

    # f = B*g
    B = LinearMap((f, g) -> B_times!(f, g, Ainv, model, Γ, ψ, q),
                  nftot; issymmetric=true, ismutating=true)

    # solves B*f = g for f... so f = Binv * g
    Binv = LinearMap((f, g) -> cg!(f, B, g, maxiter=5000, reltol=1e-12),
                     nftot; issymmetric=true, ismutating=true)

    return Binv
end


function get_Binv(model::IBModel{MultiGrid, DeformingBody{<:Motion}},
    Ainv::LinearMap, Khatmat::Array{Float64,2}, QWmat::Array{Float64,2})

    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    # TODO: Alternative... could create a dummy state to operate on here
    Γ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for streamfunction
    q = zeros(model.grid.nq, model.grid.mg)    # Working array for velocity flux

    # f = B*g
    B = LinearMap((f, g) -> B_times!(f, g, Ainv, model, Γ, ψ, q),
                  nftot; issymmetric=true, ismutating=true)

    # solves B*f = g for f... so f = Binv * g
    Binv = LinearMap((f, g) -> bicgstabl!(f, B+Khatmat*QWmat,
                    g, maxiter=5000, reltol=1e-12),
                    nftot; issymmetric=false, ismutating=true)

    return Binv
end

function B_times!(x::AbstractVector,
                  z::AbstractVector,
                  Ainv::LinearMap,
                  model::IBModel{MultiGrid, RigidBody{T}} where T<:Motion,
                  Γ::Array{Float64, 2},
                  ψ::Array{Float64, 2},
                  q::Array{Float64, 2})
    """
    %Performs one matrix multiply of B*z, where B is the matrix used to solve
    %for the surface stresses that enforce the no-slip boundary condition.
    % (B arises from an LU factorization of the full system)
    MultiGrid version:
    Note ψ is just a dummy variable for computing velocity flux
        Also this only uses Ainv on the first level
    """
    E, C = model.mats.E, model.mats.C
    Γ .*= 0.0

    # Get circulation from surface stress
    Γ[:, 1] = Ainv * ( (E*C)'*z )    # Γ = ∇ x (E'*fb)

    #-- get vel flux from circulation
    vort2flux!( ψ, q, Γ, model, model.grid.mg );

    #--Interpolate onto the body
    @views mul!(x, E, q[:, 1])
end
