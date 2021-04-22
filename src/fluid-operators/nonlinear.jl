" fq = scale*(Q1*Γ1 + Q2*Γ2) "
function nl_avg!(fq, Q1, Γ1, Q2, Γ2, fq1, fq2; scale=nothing)
   broadcast!(*, fq1, Q1, Γ1)
   broadcast!(*, fq2, Q2, Γ2)
   broadcast!(+, fq, fq1, fq2)
   if scale != nothing
      fq .*= scale
   end
   return nothing
end

"""
Compute average fluxes across cells

Return views into working memory prob.model.work.q2
"""
function avg_flux(state, prob; lev=1)
   grid = prob.model.grid
   nx, ny = grid.nx, grid.ny

   # Get reshaped 2D views into flux
   qx, qy = grid.split_flux(state.q; lev=lev);
   q0x, q0y = grid.split_flux(state.q0; lev=lev);  # Base flux

   # Same for working memory
   Q = prob.model.work.q2;  Q.*=0.0;
   Qx, Qy = grid.split_flux(Q);

   # Index into Qx from (1:nx+1)×(2:ny)
   i=1:nx+1; j=2:ny
   @views broadcast!(+, Qx[i,j], qx[i,j], q0x[i,j], qx[i,j.-1], q0x[i,j.-1] )

   # Index into Qy from (2:nx)×(1:ny+1)
   i=2:nx; j=1:ny+1
   @views broadcast!(+, Qy[i,j], qy[i,j], q0y[i,j], qy[i.-1,j], q0y[i.-1,j] )

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
   Broadcasting:
      6.902 ms (389 allocations: 1.23 MiB)
   Updated rot:
      5.033 ms (387 allocations: 2.44 MiB)
   Eliminated multiplication
      1.790 ms (383 allocations: 2.44 MiB)
   """
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   C = prob.model.mats.C

   # Alias working memory  (q1 is qs in advance, q2 becomes average Q)
   fq = prob.model.work.q3; fq.*=0.0;
   fqx, fqy = grid.split_flux(fq);

   fq1 = prob.model.work.q4; fq1.*=0.0;
   fq1x, fq1y = grid.split_flux(fq1);

   fq2 = prob.model.work.q5; fq2.*=0.0;
   fq2x, fq2y = grid.split_flux(fq2);

   Γ = reshape(state.Γ, nx-1, ny-1)

   Qx, Qy = avg_flux(state, prob)

   i=2:nx; j=2:ny-1
   @views nl_avg!(fqx[i,j], Qy[i,j.+1], Γ[i.-1,j],
         Qy[i,j], Γ[i.-1,j.-1], fq1x[i,j], fq2x[i,j])
   """
   @views broadcast!(*, fq1x[i,j], Qy[i,j.+1], Γ[i.-1,j])
   @views broadcast!(*, fq2x[i,j], Qy[i,j], Γ[i.-1,j.-1])
   @views broadcast!(+, fqx[i,j], fq1x[i,j], fq2x[i,j])
   """
   j=1;  @views broadcast!(*, fqx[i,j], Qy[i,j.+1], Γ[i.-1,j]);
   j=ny; @views broadcast!(*, fqx[i,j], Qy[i,j], Γ[i.-1,j.-1]);

   i=2:nx-1; j=2:ny
   @views nl_avg!(fqy[i,j], Qx[i.+1,j], Γ[i,j.-1],
         Qx[i,j], Γ[i.-1,j.-1], fq1y[i,j], fq2y[i,j]; scale=-1.0)
   """
   @views broadcast!(*, fq1y[i,j], -Qx[i.+1,j], Γ[i,j.-1])
   @views broadcast!(*, fq2y[i,j], -Qx[i,j], Γ[i.-1,j.-1])
   @views broadcast!(+, fqy[i,j], fq1y[i,j], fq2y[i,j])
   """
   i=1;  @views broadcast!(*, fqy[i,j], -Qx[i.+1,j], Γ[i,j.-1])
   i=nx; @views broadcast!(*, fqy[i,j], -Qx[i,j], Γ[i.-1,j.-1])

   mul!( nonlin, C', fq )
   nonlin .*= 1/(4*grid.h^2)  # Scaling factor for Γ->ω and flux averaging

   return nothing
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
   4.049 ms (505 allocations: 2.46 MiB)
   1.945 ms (427 allocations: 1.24 MiB)
   """
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
   C = prob.model.mats.C

   hc = grid.h * 2^( lev - 1);  # Coarse grid spacing

   # Alias working memory
   #fq = @view(prob.model.work.q3[:, lev]); fq.*=0.0;
   fq = prob.model.work.q3; fq.*=0.0;
   fqx, fqy = grid.split_flux(fq);

   #fq1 = @view(prob.model.work.q4[:, lev]); fq1.*=0.0;
   fq1 = prob.model.work.q4; fq1.*=0.0;
   fq1x, fq1y = grid.split_flux(fq1);

   #fq2 = @view(prob.model.work.q5[:, lev]); fq2.*=0.0;
   fq2 = prob.model.work.q5; fq2.*=0.0;
   fq2x, fq2y = grid.split_flux(fq2);

   Γ = reshape(@view(state.Γ[:, lev]), nx-1, ny-1)
   Qx, Qy = avg_flux(state, prob; lev=lev)

   # Don't need bc's for anything after this, so we can rescale in place
   Γbc .*= 0.25  # Account for scaling between grids

   ### x-fluxes
   i=2:nx; j=2:ny-1
   @views nl_avg!(fqx[i,j], Qy[i,j.+1], Γ[i.-1,j],
         Qy[i,j], Γ[i.-1,j.-1], fq1x[i,j], fq2x[i,j])
   j=1;  # Bottom boundary
   @views nl_avg!(fqx[i,j], Qy[i,j.+1], Γ[i.-1,j],
         Qy[i,j], Γbc[B.+i], fq1x[i,j], fq2x[i,j])
   j=ny;  # Top boundary
   @views nl_avg!(fqx[i,j], Qy[i,j], Γ[i.-1,j.-1],
         Qy[i,j.+1], Γbc[T.+i], fq1x[i,j], fq2x[i,j])

   ### y-fluxes
   i=2:nx-1; j=2:ny
   @views nl_avg!(fqy[i,j], Qx[i.+1,j], Γ[i,j.-1],
         Qx[i,j], Γ[i.-1,j.-1], fq1y[i,j], fq2y[i,j]; scale=-1.0)
   i=1;  # Left boundary
   @views nl_avg!(fqy[i,j], Qx[i.+1,j], Γ[i,j.-1],
         Qx[i,j],  Γbc[L.+j], fq1y[i,j], fq2y[i,j]; scale=-1.0)
   i=nx;
   @views nl_avg!(fqy[i,j], Qx[i,j], Γ[i.-1,j.-1],
         Qx[i.+1,j],  Γbc[R.+j], fq1y[i,j], fq2y[i,j]; scale=-1.0)

   mul!( nonlin, C', fq )

   # Scaling: 1/hc^2 to convert circulation to vorticity
   #   plus factor of 4 from averaging across cells above
   #  (applies to flux field, but this avoids a second operation)
   nonlin .*= 1/(4*hc^2)

   return nothing
end
