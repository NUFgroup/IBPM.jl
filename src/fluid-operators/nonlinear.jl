
"""
Compute average fluxes across cells

Return views into working memory prob.model.work.q2
"""
function avg_flux(state, prob)
   grid = prob.model.grid
   nx, ny = grid.nx, grid.ny

   # Get reshaped 2D views into flux
   qx, qy = grid.split_flux(state.q);
   q0x, q0y = grid.split_flux(state.q0);  # Base flux

   # Same for working memory
   Q = prob.model.work.q2;  Q.*=0.0;
   Qx, Qy = grid.split_flux(Q);

   # Index into Qx from (1:nx+1)×(2:ny)
   i=2:nx+1; j=2:ny
   broadcast!(+, @view(Qx[i,j]),
              @view(qx[i,j]), @view(q0x[i,j]),
              @view(qx[i,j.-1]), @view(q0x[i,j.-1]) );
   rmul!(Qx, 0.5);

   # Index into Qy from (2:nx)×(1:ny+1)
   i=2:nx; j=1:ny+1
   broadcast!(+, @view(Qy[i,j]),
              @view(qy[i,j]), @view(q0y[i,j]),
              @view(qy[i.-1,j]), @view(q0y[i.-1,j]) );
   rmul!(Qy, 0.5);

   return Qx, Qy
end

"""
!***************************************************************!
!*   nonlinear terms in rotational form                        *!
!***************************************************************!
"""
function nonlinear!( nonlin::AbstractArray,
                     state::IBState{UniformGrid},
                     prob::IBProblem )
   """
   Note factor of 1/2 in boundary conditions absorbed to eliminate
   "lastbc" from Fortran code

   Original:
      9.909 ms (236 allocations: 5.47 MiB)
   Views and precompute avg flux:
      7.538 ms (268 allocations: 4.27 MiB)


   """
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   C = prob.model.mats.C

   # Alias working memory
   fq = prob.model.work.q3; fq.*=0.0;
   fqx, fqy = grid.split_flux(fq);
   Γ = reshape(state.Γ, nx-1, ny-1)

   Qx, Qy = avg_flux(state, prob)

   i=2:nx; j=2:ny-1
   @views fqx[i,j] = 0.5*( Qy[i,j.+1].*Γ[i.-1,j] .+ Qy[i,j].*Γ[i.-1,j.-1] )
   j=1;  @views fqx[i,j] = 0.5*( Qy[i,j.+1].*Γ[i.-1,j] )
   j=ny; @views fqx[i,j] = 0.5*( Qy[i,j].*Γ[i.-1,j.-1] )

   i=2:nx-1; j=2:ny
   @views fqy[i,j] = -0.5*( Qx[i.+1,j].*Γ[i,j.-1] .+ Qx[i,j].*Γ[i.-1,j.-1] )
   i=1;  @views fqy[i,j]  = -0.5*( Qx[i.+1,j].*Γ[i,j.-1] )
   i=nx; @views fqy[i,j]  = -0.5*( Qx[i,j].*Γ[i.-1,j.-1] )

   mul!( nonlin, C', fq )    # OK TO HERE
   nonlin .*= 1/grid.h^2  # Scaling factor for differencing

   return nothing
end


"""
Compute average fluxes across cells

Return views into working memory prob.model.work.q2
"""
function avg_flux(state, prob, lev)
   grid = prob.model.grid
   nx, ny = grid.nx, grid.ny

   # Get reshaped 2D views into flux
   qx, qy = grid.split_flux(state.q, lev);
   q0x, q0y = grid.split_flux(state.q0, lev);  # Base flux

   # Same for working memory
   Q = prob.model.work.q2;  Q.*=0.0;
   Qx, Qy = grid.split_flux(Q, lev);

   # Index into Qx from (1:nx+1)×(2:ny)
   i=1:nx+1; j=2:ny
   @views broadcast!(+, Qx[i,j], qx[i,j], q0x[i,j], qx[i,j.-1], q0x[i,j.-1] )
   Qx .*= 0.5;

   # Index into Qy from (2:nx)×(1:ny+1)
   i=2:nx; j=1:ny+1
   @views broadcast!(+, Qy[i,j], qy[i,j], q0y[i,j], qy[i.-1,j], q0y[i.-1,j] )
   Qy .*= 0.5;

   return Qx, Qy
end


function nonlinear!( nonlin::AbstractArray,
                     state::IBState{MultiGrid},
                     Γbc::AbstractArray,
                     lev::Int,
                     prob::IBProblem )
   """
   Note factor of 1/2 in boundary conditions absorbed to eliminate
   "lastbc" from Fortran code

   11.066 ms (323 allocations: 4.28 MiB)

   9.984 ms (400 allocations: 1.84 MiB)
   """
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
   C = prob.model.mats.C

   hc = grid.h * 2^( lev - 1);  # Coarse grid spacing

   # Alias working memory
   fq = @view(prob.model.work.q3[:, lev]); fq.*=0.0;
   fqx, fqy = grid.split_flux(fq);

   Γ = reshape(@view(state.Γ[:, lev]), nx-1, ny-1)
   Qx, Qy = avg_flux(state, prob, lev)

   # Don't need bc's for anything after this, so we can rescale in place
   Γbc .*= 0.5

   i=2:nx; j=2:ny-1
   @views @. fqx[i,j] = 0.5*( Qy[i,j+1]*Γ[i-1,j] + Qy[i,j]*Γ[i-1,j-1] )
   j=1;  @views @. fqx[i,j] = 0.5*( Qy[i,j+1]*Γ[i-1,j] + Qy[i,j]*Γbc[B+i] )
   j=ny; @views @. fqx[i,j] = 0.5*( Qy[i,j]*Γ[i-1,j-1] + Qy[i,j+1]*Γbc[T+i] )

   i=2:nx-1; j=2:ny
   @views @. fqy[i,j] = -0.5*( Qx[i+1,j]*Γ[i,j-1] + Qx[i,j]*Γ[i-1,j-1] )
   i=1;  @views @. fqy[i,j] = -0.5*( Qx[i+1,j]*Γ[i,j-1] + Qx[i,j]*Γbc[L+j] )
   i=nx; @views @. fqy[i,j] = -0.5*( Qx[i,j]*Γ[i-1,j-1] + Qx[i+1,j]*Γbc[R+j] )

   mul!( nonlin, C', fq )
   nonlin .*= 1/hc^2   # Scaling factor for differencing

   return nothing
end
