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
    u, v = grid.u_idx, grid.v_idx
    i = 2:grid.nx; j = 2:grid.ny
    Γ = reshape(Γ, grid.nx-1, grid.ny-1)
    qx, qy = grid.split_flux(q)

    @views Γ[i.-1, j.-1] = qy[i, j] - qy[i.-1, j] - qx[i, j] + qx[i, j.-1]
    Γ = reshape(Γ, grid.nΓ, 1)
end

"""
    curl!( q, ψ, grid )

Discrete curl operator

q = curl( ψ )
"""
function curl!( q, ψ, grid )
   nx, ny = grid.nx, grid.ny
   u, v = grid.u_idx, grid.v_idx  # Indices to x-flux, y-flux, vorticity/streamfunction
   ψ = reshape(ψ, grid.nx-1, grid.ny-1)
   qx, qy = grid.split_flux(q)

   # x-fluxes
   i=2:nx; j=2:ny-1
   @views qx[i, j] = ψ[i.-1,j] - ψ[i.-1, j.-1]
   j=1;  @views qx[i, j] =  ψ[i.-1, j]      # Top boundary
   j=ny; @views qx[i, j] = -ψ[i.-1, j.-1]   # Bottom boundary

   # y-fluxes
   i=2:nx-1;  j=2:ny
   @views qy[i,j] = ψ[i.-1,j.-1] - ψ[i,j.-1]
   i=1;  @views qy[i,j] = -ψ[i, j.-1]           # Left boundary
   i=nx; @views qy[i,j] =  ψ[i.-1, j.-1]        # Right boundary

   ψ = reshape(ψ, grid.nΓ, 1)
end

function curl!( q, ψ, ψbc, grid::MultiGrid )
    """
    Compute velocity flux from streamfunction.
       Note: requires streamfunction from coarser grid on edge
             of current domain (stored in ψbc)
    """
   nx, ny = grid.nx, grid.ny
   u, v = grid.u_idx, grid.v_idx  # Indices to x-flux, y-flux, vorticity/streamfunction
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
   ψ = reshape(ψ, grid.nx-1, grid.ny-1)
   qx, qy = grid.split_flux(q)

   ### x-fluxes
   i=2:nx; j=2:ny-1
   @views qx[i, j] = ψ[i.-1,j] - ψ[i.-1, j.-1]

   # Top/bottom boundaries
   j=1;  @views qx[i, j] = ψ[i.-1, j] - ψbc[(B-1).+i]
   j=ny; @views qx[i, j] = ψbc[i.+T] - ψ[i.-1, j.-1]

   # Left/right boundaries
   j=1:ny;
   i=1;     @views qx[i, j] = ψbc[(L+1).+j] - ψbc[L.+j]
   i=nx+1;  @views qx[i, j] = ψbc[(R+1).+j] - ψbc[R.+j]

   #### y-fluxes
   i=2:nx-1;  j=2:ny
   @views qy[i,j] = ψ[i.-1,j.-1] - ψ[i,j.-1]

   # Left/right boundaries
   @views i=1;  qy[i,j] = -ψ[i, j.-1][:] + ψbc[L.+j]
   @views i=nx; qy[i,j] = -ψbc[R.+j] + ψ[i.-1, j.-1][:]

   # Top/bottom boundaries
   i=1:nx
   j=1;    @views qy[i,j] = ψbc[B.+i] - ψbc[(B+1).+i]
   j=ny+1; @views qy[i,j] = ψbc[T.+i] - ψbc[(T+1).+i]

   ψ = reshape(ψ, grid.nΓ, 1)
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
   ψbc = model.work.Γbc  # Same shape for both boundary conditions
   Γwork = @view(model.work.Γ3[:, 1])  # Working memory

   #println("=== VORT2FLUX ===")
   # Interpolate values from finer grid to center region of coarse grid
   for lev=2:ngrids
      @views coarsify!( Γ[:, lev], Γ[:, lev-1], grid);
   end

   # Invert Laplacian on largest grid
   ψ .*= 0.0

   @views mul!(ψ[:, ngrids], model.mats.Δinv, Γ[:, ngrids])  # Δψ = Γ
   @views curl!(q[:, ngrids], ψ[:, ngrids], ψbc, grid )      # q = ∇×ψ

   # Telescope in to finer grids
   for lev=(ngrids-1):-1:1
      Γwork .= Γ[:, lev]
      @views get_bc!(ψbc, ψ[:, lev+1], 1.0, grid)
      @views apply_bc!(Γwork, ψbc, 1.0, grid)

      @views mul!(ψ[:, lev], model.mats.Δinv, Γwork)   # Δψ = Γ
      if lev < ngrids
         @views curl!(q[:, lev], ψ[:, lev], ψbc, grid )   # q = ∇×ψ
      end
   end

   #println(sum(q.^2))
   #sleep(100)
end
