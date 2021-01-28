mutable struct WorkingMemory
    q1::AbstractArray
    q2::AbstractArray
    q3::AbstractArray
    q4::AbstractArray
    Γ1::AbstractArray
    Γ2::AbstractArray
    Γ3::AbstractArray
    bc::AbstractArray
end


function init_memory(grid::UniformGrid)
    return WorkingMemory(
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nq, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ, 1),
        zeros(grid.nΓ)
    )
end

function init_memory(grid::MultiGrid)
    return WorkingMemory(
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nq, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ, grid.mg),
        zeros(grid.nΓ)
    )
end



function get_Ainv(model::IBModel{<:Grid, <:Body}, dt::Float64;
                  lev::Int=1)
    """
    Construct LinearMap to solve
        (I + dt/2 * Beta * RC) * x = b for x
    where Beta = 1/(Re * h^2)

    Compare to get_RCinv in "ib_mats.jl"

    Benchmarks with nΓ = 39601
    Original function (non-mutating DST):
        1.611 ms (22 allocations: 2.12 MiB)
    Original function (mutating DST):
        1.526 ms (14 allocations: 1.51 MiB)
    LinearMap (non-mutating)
        1.413 ms (9 allocations: 619.25 KiB)
    LinearMap (mutating)
        1.281 ms (5 allocations: 192 bytes)
    """

    hc = model.grid.h * 2^( lev - 1);  # Grid size at this level
    # Solve by transforming to and from Fourier space and scaling by evals
    Λ̃ = 1 .+ model.mats.Λ * dt/( 2 * model.Re * hc^2 );

    # give output in same size as input b (before being reshaped)
    return LinearMap((x, b) -> Λinv_fn!(x, b, model.grid.nx, model.grid.ny, Λ̃, model.mats.dst_plan),
                     model.grid.nΓ; issymmetric=true, ismutating=true)
end


function B_times!(x::Array{Float64, 2},
                  z::Array{Float64, 2},
                  Ainv::LinearMap,
                  model::IBModel{UniformGrid, RigidBody{Static}},
                  Γ::Array{Float64, 2},
                  ψ::Array{Float64, 2},
                  q::Array{Float64, 2})
    """
    %Performs one matrix multiply of B*z, where B is the matrix used to solve
    %for the surface stresses that enforce the no-slip boundary condition.

    % (B arises from an LU factorization of the full system)

    Note ψ is just a dummy work array for circ2_st_vflx
    """
    # -- get circulation from surface stress  circ = Ainv * R * E' * z
    #     We don't include BCs for Ainv because ET*z is compact
    #circ = Ainv( mats.R*(mats.ET*z), dt, model );
    mul!(Γ, Ainv, model.mats.RET*z)

    #-- get vel flux from circulation
    #vflx, _ = circ2_st_vflx( circ, model );
    circ2_st_vflx!( ψ, q, Γ, model );

    #--Interpolate onto the body and scale by h
    #x = (model.mats.E*vflx) / model.grid.h;
    mul!(x, model.mats.E, q)
    rmul!(x, 1/model.grid.h)
end


function get_A(model::IBModel{UniformGrid, <:Body}, dt::Float64)
    A = I - (dt/2 / (model.grid.h^2))*model.mats.Lap
    Ainv = get_Ainv(model, dt)
    return A, Ainv
end

function get_A(model::IBModel{MultiGrid, <:Body}, dt::Float64)
    hc = [model.grid.h * 2^(lev-1) for lev=1:model.grid.mg]
    A = [I - (dt/2 / (hc[lev]^2)) *model.mats.Lap for lev=1:model.grid.mg]
    Ainv = [get_Ainv(model, dt; lev=lev) for lev=1:model.grid.mg]
    return A, Ainv
end

function get_B(model::IBModel{UniformGrid, RigidBody{Static}}, Ainv::LinearMap)
    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    # need to build and store surface stress matrix and its inverse if at first time step
    B = zeros( nftot, nftot );
    # Pre-allocate arrays
    e = zeros( nftot, 1 );         # Unit vector

    # TODO: Alternative... could create a dummy state to operate on here
    b = zeros( nftot, 1 );         # Working array
    Γ = zeros(model.grid.nΓ, 1)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, 1)    # Working array for streamfunction
    q = zeros(model.grid.nq, 1)    # Working array for velocity flux

    for j = 1 : nftot
        if j>1
            e[j-1] = 0.0
        end
        e[j] = 1.0;

        B_times!( b, e, Ainv, model, Γ, ψ, q );
        B[:, j] = b
    end
    Binv = inv(B)
    return Binv
end



function get_B(model::IBModel{MultiGrid, RigidBody{Static}}, Ainv)
    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    # need to build and store surface stress matrix and its inverse if at first time step
    B = zeros( nftot, nftot );
    # Pre-allocate arrays
    e = zeros( nftot, 1 );         # Unit vector

    # TODO: Alternative... could create a dummy state to operate on here
    b = zeros( nftot, 1 );         # Working array
    Γ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for streamfunction
    q = zeros(model.grid.nq, model.grid.mg)    # Working array for velocity flux

    for j = 1 : nftot
        if j>1
            e[j-1] = 0.0
        end
        e[j] = 1.0;

        # Only fine-grid Ainv is used here
        B_times!( b, e, Ainv[1], model, Γ, ψ, q );
        B[:, j] = b
    end
    Binv = inv(B)
    return Binv
end




function get_AB(model::IBModel, dt::Float64)
    A, Ainv = get_A(model, dt)
    Binv = get_B(model, Ainv)
    return A, Ainv, Binv
end


"""
RotatingCyl functions
"""


function B_times!(x::AbstractArray,
                  z::AbstractArray,
                  Ainv::LinearMap,
                  model::IBModel{MultiGrid, RigidBody{T}} where T<:Motion,
                  Γ::Array{Float64, 2},
                  ψ::Array{Float64, 2},
                  q::Array{Float64, 2})
    """
    %Performs one matrix multiply of B*z, where B is the matrix used to solve
    %for the surface stresses that enforce the no-slip boundary condition.

    % (B arises from an LU factorization of the full system)

    Note ψ is just a dummy variable for computing velocity flux
        Also this only uses Ainv on the first level

    7.080 ms (50 allocations: 3.34 MiB
    3.971 ms (30 allocations: 627.05 KiB)
    """

    # --Initialize
    grid = model.grid
    mats = model.mats

    # -- get circulation from surface stress
    Γ[:, 1] = Array( Ainv * (mats.RET*z) );
    #@views mul!( Γ[:, 2], mats.RET, z )  # Using Γ[:, 2] as dummy array for multiplication
    #@views mul!( Γ[:, 1], Ainv, Γ[:, 2] )

    # Coarsify circulation to second grid level to get BCs for stfn
    @views coarsify!( Γ[:,1], Γ[:,2], grid );

    #-- get vel flux from circulation
    circ2_st_vflx!( ψ, q, Γ, model, 2 );  # THIS IS THE MOST EXPENSIVE THING

    #--Interpolate onto the body and scale by h
    @views mul!(x, mats.E, q[:, 1])
    rmul!(x, 1/grid.h)
end


function get_B(model::IBModel{MultiGrid, RigidBody{T}} where T <: Motion, Ainv)
    nb, nf = get_body_info(model.bodies)
    nftot = sum(nf)

    # TODO: Alternative... could create a dummy state to operate on here
    Γ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for circulation
    ψ = zeros(model.grid.nΓ, model.grid.mg)    # Working array for streamfunction
    q = zeros(model.grid.nq, model.grid.mg)    # Working array for velocity flux

    # f = B*g
    B = LinearMap((f, g) -> B_times!(f, g, Ainv[1], model, Γ, ψ, q),
                  nftot; issymmetric=true, ismutating=true)

    # solves f = B*g for g... so g = Binv * f
    Binv = LinearMap((f, g) -> cg!(f, B, g, maxiter=5000, reltol=1e-12),
                     nftot; issymmetric=true, ismutating=true)

    return Binv
end
