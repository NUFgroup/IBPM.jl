"""
Linear operations associated with the governing flow equations

Note: all matrices in this file are "static" and do not vary in time. They may
therefore be pre-computed as a pre-processing step.
"""

"""
Transpose of discrete curl (R matrix)

Γ = rot( q )
"""
function rot( q, grid )
    Γ = zeros(grid.nΓ)
    for j=1:grid.ny-1
        for i=1:grid.nx-1
            Γ[grid.ω(i, j)] = q[grid.v(i+1, j+1)] - q[grid.v(i, j+1)] - q[grid.u(i+1, j+1)] + q[grid.u(i+1, j)]
        end
    end

    return Γ
end


"""
!***************************************************************!
!*   returns curl(x) given x and the values of s'fun on bdy    *!
!***************************************************************!
"""
function curl!( q, ψ, grid )
   nx, ny = grid.nx, grid.ny
   for j=2:ny-1
      for i=2:nx-1
         q[grid.u(i, j)] = ψ[grid.ω(i, j+1)] - ψ[grid.ω(i, j)]
      end
   end

   for i=1:nx
      j=1
      q[grid.u(i, j)] = ψ[grid.ω(i, j+1)]
      j=ny
      q[grid.u(i, j)] = ψ[grid.ω(i, j)]
   end

   for j=2:ny
       i=1
       q[grid.v(i,j)] = - ψ[grid.ω(i+1, j)]
       for i=2:nx-1
          q[grid.v(i,j)] = ψ[grid.ω(i,j)] - ψ[grid.ω(i+1,j)]
      end
       i=nx
       q[grid.v(i,j)] = ψ[grid.ω(i,j)]
   end
   return q
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
    nx = grid.nx; ny = grid.ny;
    Λ = -2*( cos.( π*(1:(ny-1))/ny ) .+ cos.( π*(1:(nx-1))/nx )' .- 2);
    return Λ
end


function Λinv_fn!(x::AbstractArray,
                  b::AbstractArray,
                  nx::Int,
                  ny::Int,
                  Λ::AbstractArray,
                  dst_plan::Tuple{Any, Array{Float64, 2}})
    """
    Solve inverse Laplacian RC (or similar for A matrix)
    Used to construct both Ainv and RCinv operators
    """
    # reshape for inversion in fourier space
    b = reshape( b, nx-1, ny-1)
    x = reshape( x, nx-1, ny-1)
    dst_inv!(x, b, Λ, dst_plan);
    # Include scale to make fwd/inv transforms equal
    rmul!(x, 4.0/( nx*ny ))
    x = reshape( x, (nx-1)*(ny-1), 1 )
end


function get_lap_inv( grid::T,
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



function get_Ainv(grid, dst_plan, Re, dt::Float64; lev::Int=1)
    """
    Construct LinearMap to solve
        (I + dt/2 * Beta * RC) * x = b for x
    where Beta = 1/(Re * h^2)

    Compare to get_RCinv in "ib_mats.jl"
    """

    hc = grid.h * 2^( lev - 1);  # Grid size at this level
    # Solve by transforming to and from Fourier space and scaling by evals
    Λ = Lap_eigs( grid )
    Λ̃ = 1 .+ Λ * dt/( 2 * Re * hc^2 );

    # give output in same size as input b (before being reshaped)
    return LinearMap((x, b) -> Λinv_fn!(x, b, grid.nx, grid.ny, Λ̃, dst_plan),
                     grid.nΓ; issymmetric=true, ismutating=true)
end
