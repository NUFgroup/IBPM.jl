"""
Linear operations associated with the governing flow equations

Note: all matrices in this file are "static" and do not vary in time. They may
therefore be pre-computed as a pre-processing step.
"""

include("utilities/dst-inversion.jl")
include("utilities/indexing-utils.jl")
include("utilities/multigrid-utils.jl")

"""
Discrete curl operator -- relates streamfunction to velocity flux, q = Cs

Note that R := C^T is also a discrete curl operator that relates vel flux to
    circulation, gamma = R q
As such, the vel flux can be backed out from the circulation via
    q = C(C'C)^-1 gamma
This relationship is implemented in circ2_st_vflux
"""
function get_C( grid::T ) where T <: Grid
    m = grid.nx;
    n = grid.ny;

    nrows = get_velx_ind( m-1, n, grid ) +
            get_vely_ind( m, n-1, grid ) ;
    ncols = get_vort_ind( m-1, n-1, grid ) ;
    C = spzeros(nrows, ncols)

    # --First build block corresponding to x-velocities

    # vorticity points above x-velocity point
    velx_ind = get_velx_ind(x1_ind(m, n), y1_ind(m, n), grid);

    vort_ind = get_vort_ind(x1_ind(m, n), y1_ind(m, n), grid);

    C .+= sparse( velx_ind, vort_ind,
              ones(size(vort_ind)), nrows, ncols);

    #vorticity points below x-velocity point (only the velx_ind chages)
    velx_ind = get_velx_ind(x1_ind(m, n), y2_ind(m, n), grid);

    C .-= sparse( velx_ind, vort_ind,
              ones(size(vort_ind)), nrows, ncols);

    #--Now build y-velocity block

    # rows start at end of x-velocity block:
    n_add = get_velx_ind( m-1, n, grid );

    # vorticity points to the right of y-velocity point
    vely_ind = get_vely_ind(x1_ind(m, n), y1_ind(m, n), grid);

    C .-= sparse(n_add .+ vely_ind, vort_ind,
                ones(size(vort_ind)), nrows, ncols);

    # vorticity points to the left of y-velocity point
    vely_ind = get_vely_ind(x2_ind(m, n), y1_ind(m, n), grid);

    C .+= sparse(n_add .+ vely_ind, vort_ind,
                ones(size(vort_ind)), nrows, ncols);

    return C
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

function Lap_eigs( grid::T ) where T <: Grid
    # eigenvalues of RC (negative of the evals of the 5point stencil Lap)
    m = grid.nx;
    n = grid.ny;

    ii, jj = meshgrid( 1:(m-1), 1:(n-1) );  # Should replace with broadcasting
    Λ = -2*( cos.( π*ii/m ) .+ cos.( π*jj/n ) .- 2);
    return Λ
end


function Λinv_fn!(x::AbstractArray,
                  b::AbstractArray,
                  m::Int,
                  n::Int,
                  Λ::AbstractArray,
                  dst_plan::Tuple{Any, Array{Float64, 2}})
    """
    Solve inverse Laplacian RC (or similar for A matrix)
    Used to construct both Ainv and RCinv operators
    """
    # reshape for inversion in fourier space
    b = reshape( b, m-1, n-1)
    x = reshape( x, m-1, n-1)
    dst_inv!(x, b, Λ, dst_plan);
    # Include scale to make fwd/inv transforms equal
    rmul!(x, 4.0/( m*n ))
    x = reshape( x, (m-1)*(n-1), 1 )
end


function get_RCinv( grid::T,
                    Λ::AbstractArray,
                    dst_plan::Tuple{Any, Array{Float64, 2}}) where T <: Grid
    """
    Construct LinearMap to solve
        (I + dt/2 * Beta * RC) * x = b for x
    where Beta = 1/(Re * h^2)

    Benchmarks with nΓ = 39601
    Original function (non-mutating DST):
        1.506 ms (10 allocations: 928.73 KiB)
    LinearMap (mutating)
        1.251 ms (4 allocations: 160 bytes)
    """

    # give output in same size as input b (before being reshaped)
    return LinearMap((x, b) -> Λinv_fn!(x, b, grid.nx, grid.ny, Λ, dst_plan),
                     grid.nΓ; issymmetric=true, ismutating=true)
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


"""
Solve Poisson problem to go from circulation to vorticity flux.
Used to:
    - back out a trial velocity flux after the trial circulation,
        to get the surface stress
    - update the velocity flux at the end of a time step, once the circulation
        that satisfies the no-slip BCs has been computed
"""
function circ2_st_vflx!( ψ::AbstractArray,
                         q::AbstractArray,
                         Γ::AbstractArray,
                         model::IBModel{UniformGrid, <:Body} )
    """
    Take vorticity and return velocity flux and streamfunction
    """
    # Solve for streamfcn  ψ = RCinv * Γ
    mul!(ψ, model.mats.RCinv, Γ)

    #--Get velocity flux from curl of stream function
    mul!(q, model.mats.C, ψ)
end


function circ2_st_vflx!( ψ::AbstractArray,
                         q::AbstractArray,
                         Γ::AbstractArray,
                         model::IBModel{MultiGrid, <:Body},
                         ngrids::Int )
         grid = model.grid

         for lev = ngrids:-1:1

             grid.stbc .*= 0.0   # Reset boundary conditions in pre-allocated memory
             #--Solve Poisson problem for ψ

             #BCs for Poisson problem (will be overwritten if mg > 1)
             if ( lev < ngrids )
                get_stfn_BCs!( grid.stbc, @view(ψ[:, lev+1]), model );

                # Solve for streamfcn
                # TODO: in place addition with @view macro
                ψ[:, lev] = Array( model.mats.RCinv * (Γ[:, lev] .+ sum(grid.stbc, dims=2)) )
             else  # don't need bcs for largest grid
             # TODO: Figure out how to use view for in-place
                ψ[:, lev] = model.mats.RCinv * Γ[:, lev]
                #@views mul!(ψ[:, lev], model.mats.RCinv, Γ[:, lev])
             end

             #--Get velocity on first grid from stream function
             @views curl!( q[:,lev], ψ[:,lev], grid.stbc, model );
         end
end
