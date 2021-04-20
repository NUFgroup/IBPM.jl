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
    #i = 1:grid.nx-1; j = (1:grid.ny-1)'
    #@. Γ[ω(i, j)] = q[v(i+1, j+1)] - q[v(i, j+1)] - q[u(i+1, j+1)] + q[u(i+1, j)]
    i = 2:grid.nx; j = (2:grid.ny)'
    @. Γ[ω(i-1, j-1)] = q[v(i, j)] - q[v(i-1, j)] - q[u(i, j)] + q[u(i, j-1)]
end

"""
    curl!( q, ψ, grid )

Discrete curl operator

q = curl( ψ )
"""
function curl!( q, ψ, grid )
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


function curl!( q, ψ, ψbc, grid::MultiGrid )
    """
    Compute velocity flux from streamfunction.
       Note: requires streamfunction from coarser grid on edge
             of current domain (stored in ψbc)

    ! for multigrid boundary conditions
    DIMENSION(2*(m+1)+2*(n+1)) :: ψbc
    """
   nx, ny = grid.nx, grid.ny
   u, v, ω = grid.u_idx, grid.v_idx, grid.ω_idx  # Indices to x-flux, y-flux, vorticity/streamfunction
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs

   ### x-fluxes
   i=2:nx; j=(2:ny-1)'
   @. q[u(i, j)] = ψ[ω(i-1,j)] - ψ[ω(i-1, j-1)]

   # Top/bottom boundaries
   j=1;  @. q[u(i, j)] = ψ[ω(i-1, j)] - ψbc[B+i-1 ]
   j=ny; @. q[u(i, j)] = ψbc[T+i] - ψ[ω(i-1, j-1)]

   # Left/right boundaries
   j=(1:ny)';
   i=1;     @. q[u(i, j)] = ψbc[L+j+1] - ψbc[L+j]
   i=nx+1;  @. q[u(i, j)] = ψbc[R+j+1] - ψbc[R+j]

   #### y-fluxes
   i=2:nx-1;  j=(2:ny)'
   @. q[v(i,j)] = ψ[ω(i-1,j-1)] - ψ[ω(i,j-1)]

   # Left/right boundaries
   i=1;  @. q[v(i,j)] = -ψ[ω(i, j-1)] + ψbc[L+j]
   i=nx; @. q[v(i,j)] = -ψbc[R+j] + ψ[ω(i-1, j-1)]

   # Top/bottom boundaries
   i=1:nx
   j=1;    @. q[v(i,j)] = ψbc[B+i] - ψbc[B+i+1]
   j=ny+1; @. q[v(i,j)] = ψbc[T+i] - ψbc[T+i+1]
end



function curl!( q, ψ, ψc, grid::MultiGrid )
    """
    Compute velocity flux from streamfunction.
       Note: requires streamfunction from coarser grid on edge
             of current domain (stored in ψc)

    ! for multigrid boundary conditions
    DIMENSION(2*(m+1)+2*(n+1)) :: ψbc
    """
   nx, ny = grid.nx, grid.ny
   u, v, ω = grid.u_idx, grid.v_idx, grid.ω_idx  # Indices to x-flux, y-flux, vorticity/streamfunction
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs

   ### x-fluxes
   i=2:nx; j=(2:ny-1)'
   @. q[u(i, j)] = ψ[ω(i-1,j)] - ψ[ω(i-1, j-1)]

   # Top/bottom boundaries
   j=1;  @. q[u(i, j)] = ψ[ω(i-1, j)] - ψbc[B+i-1 ]
   j=ny; @. q[u(i, j)] = ψbc[T+i] - ψ[ω(i-1, j-1)]

   # Left/right boundaries
   j=(1:ny)';
   i=1;     @. q[u(i, j)] = ψbc[L+j+1] - ψbc[L+j]
   i=nx+1;  @. q[u(i, j)] = ψbc[R+j+1] - ψbc[R+j]

   #### y-fluxes
   i=2:nx-1;  j=(2:ny)'
   @. q[v(i,j)] = ψ[ω(i-1,j-1)] - ψ[ω(i,j-1)]

   # Left/right boundaries
   i=1;  @. q[v(i,j)] = -ψ[ω(i, j-1)] + ψbc[L+j]
   i=nx; @. q[v(i,j)] = -ψbc[R+j] + ψ[ω(i-1, j-1)]

   # Top/bottom boundaries
   i=1:nx
   j=1;    @. q[v(i,j)] = ψbc[B+i] - ψbc[B+i+1]
   j=ny+1; @. q[v(i,j)] = ψbc[T+i] - ψbc[T+i+1]
end

"""
    vort2flux!( ψ, q, Γ, model::IBModel{UniformGrid, <:Body} )
"""
function vort2flux!( ψ::AbstractArray,
                     q::AbstractArray,
                     Γ::AbstractArray,
                     model::IBModel{UniformGrid, <:Body} )
    mul!(ψ, model.mats.Δinv, Γ)     # Solve Poisson problem for streamfunction
    mul!(q, model.mats.C, ψ)          # q = ∇ x ψ
end


"""
!***************************************************************!
!*    Multiscale method to solve C^T C s = omega               *!
!*    and return the velocity, C s.                            *!
!*    Results are returned in vel on each of the first nlev    *!
!*    grids.                                                   *!
!*    Warning: the vorticity field on all but the finest mesh  *!
!*     is modified by the routine in the following way:        *!
!*     the value in the center of the domain is interpolated   *!
!*     from the next finer mesh (the value near the edge is    *!
!*     not changed.                                            *!
!***************************************************************!
"""
function vort2flux!( ψ::AbstractArray,
                     q::AbstractArray,
                     Γ::AbstractArray,
                     model::IBModel{MultiGrid, <:Body},
                     ngrids::Int )
   grid = model.grid
   nx, ny = grid.nx, grid.ny

   #println("=== VORT2FLUX ===")
   # Interpolate values from finer grid to center region of coarse grid
   for lev=2:ngrids
      @views coarsify!( Γ[:, lev], Γ[:, lev-1], grid);
   end

   # Invert Laplacian on largest grid
   ψ .*= 0.0

   # TODO: In place??
   ψ[:, ngrids] = model.mats.Δinv * Γ[:, ngrids]    # Δψ = Γ
   @views curl!(q[:, ngrids], ψ[:, ngrids], grid )  # q = ∇×ψ

   # Telescope in to finer grids
   for lev=(ngrids-1):-1:1
      # TODO REDO WITHOUT ALLOCATION
      Γ̃ = copy(Γ[:, lev])

      @views get_bc!(ψbc, ψ[:, lev+1], 1.0, grid)
      @views apply_bc!(Γ̃, ψbc, 1.0, grid)

      # TODO: In place??
      ψ[:, lev] = model.mats.Δinv * Γ̃             # Δψ = Γ

      if ngrids >= lev
         @views curl!(q[:, lev], ψ[:, lev], ψbc, grid )   # q = ∇×ψ
      end
   end

   #println(sum(q.^2))
   #sleep(100)
end




function vort2flux_OLD!( ψ::AbstractArray,
                     q::AbstractArray,
                     Γ::AbstractArray,
                     model::IBModel{MultiGrid, <:Body},
                     ngrids::Int )
   grid = model.grid
   nx, ny = grid.nx, grid.ny

   # TODO: PREALLOCATE
   ψbc = zeros(2*(nx+1) + 2*(ny+1))   # Streamfunction boundary conditions

   #println("=== VORT2FLUX ===")
   # Interpolate values from finer grid to center region of coarse grid
   for lev=2:ngrids
      @views coarsify!( Γ[:, lev], Γ[:, lev-1], grid);
   end

   # Invert Laplacian on largest grid
   ψ .*= 0.0

   # TODO: In place??
   ψ[:, ngrids] = model.mats.Δinv * Γ[:, ngrids]             # Δψ = Γ
   @views curl!(q[:, ngrids], ψ[:, ngrids], ψbc, grid )      # q = ∇×ψ

   # Telescope in to finer grids
   for lev=(ngrids-1):-1:1
      # TODO REDO WITHOUT ALLOCATION
      Γ̃ = copy(Γ[:, lev])

      @views get_bc!(ψbc, ψ[:, lev+1], 1.0, grid)
      @views apply_bc!(Γ̃, ψbc, 1.0, grid)

      # TODO: In place??
      ψ[:, lev] = model.mats.Δinv * Γ̃             # Δψ = Γ

      if ngrids >= lev
         @views curl!(q[:, lev], ψ[:, lev], ψbc, grid )   # q = ∇×ψ
      end
   end

   #println(sum(q.^2))
   #sleep(100)
end
