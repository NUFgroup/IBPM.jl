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
function avg_flux(state::IBState, prob::IBProblem; lev=1)
   grid = prob.model.grid
   nx, ny = grid.nx, grid.ny

   # Get reshaped 2D views into flux
   qx, qy = grid.split_flux(state.q; lev=lev);
   q0x, q0y = grid.split_flux(state.q0; lev=lev);  # Base flux from background flow

   # Same for working memory
   Q = prob.model.work.q2;  Q.*=0.0;
   Qx, Qy = grid.split_flux(Q);

   # Index into Qx from (1:nx+1)×(2:ny)
   i=1:nx+1; j=2:ny
   @views broadcast!(+, Qx[i,j], qx[i,j], q0x[i,j], qx[i,j.-1], q0x[i,j.-1] )

   # Index into Qy from (2:nx)×(1:ny+1)
   i=2:nx; j=1:ny+1
   @views broadcast!(+, Qy[i,j], qy[i,j], q0y[i,j], qy[i.-1,j], q0y[i.-1,j] )

   return Q
end

"""
Compute average fluxes across cells in linearized problem (no freestream)

Return views into working memory prob.model.work.q2
"""
function avg_flux(state::IBState, prob::LinearizedIBProblem; lev=1)
   grid = prob.model.grid
   nx, ny = grid.nx, grid.ny

   # Get reshaped 2D views into flux (note no background flux here)
   qx, qy = grid.split_flux(state.q; lev=lev);

   # Same for working memory
   Q = prob.model.work.q2;  Q.*=0.0;
   Qx, Qy = grid.split_flux(Q);

   # Index into Qx from (1:nx+1)×(2:ny)
   i=1:nx+1; j=2:ny
   @views broadcast!(+, Qx[i,j], qx[i,j], qx[i,j.-1] )

   # Index into Qy from (2:nx)×(1:ny+1)
   i=2:nx; j=1:ny+1
   @views broadcast!(+, Qy[i,j], qy[i,j], qy[i.-1,j] )

   return Q
end

"""
   direct_product_loops!(fq, )
"""
function direct_product_loops!(fq, Q, Γ, Γbc, grid, fq_tmp1, fq_tmp2)
   nx, ny = grid.nx, grid.ny
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs

   Qx, Qy = grid.split_flux(Q)   # These are views into Q (though won't be changed here)
   fqx, fqy = grid.split_flux(fq);  # Views into fq - changing them will change fq

   #-- Alias working memory
      fq_tmp1.*=0.0;
      fq1x, fq1y = grid.split_flux(fq_tmp1);  # Views of fq1 (which is not itself needed anywhere)

      fq_tmp2.*=0.0;
      fq2x, fq2y = grid.split_flux(fq_tmp2);  # Views of fq2
   #--

   #-- Product of flux and circulation
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
   #--

end

"""
   direct_product!(fq, Q, Γ, prob::IBProblem)

Compute the direct product of Q and Γ for nonlinear term.

Needs three flux-shaped matrices from working memory
"""
function direct_product!(fq, Q, Γ, Γbc, prob::IBProblem)
   grid = prob.model.grid; work = prob.model.work

   # fq is the output array: the product of flux and circulation such that the nonlinear
   #   term is C'*fq (or ∇⋅fq)
   fq.*=0.0;  # Zero out in case some locations aren't indexed

   # Call helper function to loop over the arrays and store product in fq
   direct_product_loops!(fq, Q, Γ, Γbc, grid, work.q4, work.q5)

   return nothing
end

"""
   direct_product!(fq, Q, Γ, prob::LinearizedIBProblem)

Compute the base flow advection terms for the linear problem.

Needs FOUR flux-shaped arrays from working memory (one extra compared to nonlinear solver)
"""
function direct_product!(fq, Q, Γ, Γbc, prob::LinearizedIBProblem)
   grid = prob.model.grid; work = prob.model.work

   # fq is the output array: the product of flux and circulation such that the nonlinear
   #   term is C'*fq (or ∇⋅fq)
   fq.*=0.0;  # Zero out in case some locations aren't indexed
   fq_tmp = work.q6; fq_tmp.*=0.0  # Extra array for second product (saves allocation)

   # Call helper function to loop over the arrays and store product in fq
   direct_product_loops!(fq, Q, ΓB, Γbc, grid, work.q4, work.q5)
   direct_product_loops!(fq_tmp, QB, Γ, Γbc, grid, work.q4, work.q5)

   fq .+= fq_tmp  # Add results
   return nothing
end

"""
    nonlinear!(nonlin, state, Γbc, lev, prob)
"""
function nonlinear!( nonlin::AbstractArray,
                     state::IBState{MultiGrid},
                     Γbc::AbstractArray,
                     lev::Int,
                     prob::AbstractIBProblem )
   """
   Note factor of 1/2 in boundary conditions absorbed to eliminate
   "lastbc" from Fortran code
   """
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   C = prob.model.mats.C

   hc = grid.h * 2^( lev - 1);  # Coarse grid spacing

   # Don't need bc's for anything after this, so we can rescale in place
   Γbc .*= 0.25  # Account for scaling between grids

   Γ = reshape(@view(state.Γ[:, lev]), nx-1, ny-1)  # Circulation at this grid level
   Q = avg_flux(state, prob; lev=lev)  # Compute average fluxes across cells

   fq = prob.model.work.q3; # Alias working memory
   direct_product!(fq, Q, Γ, Γbc, prob)  # Product of flux and circulation

   # Divergence of flux-circulation product
   mul!( nonlin, C', fq )

   # Scaling: 1/hc^2 to convert circulation to vorticity
   #   plus factor of 4 from averaging across cells above
   #  (applies to flux field, but this avoids a second operation)
   nonlin .*= 1/(4*hc^2)

   return nothing
end

"""
    nonlinear!(nonlin, state, Γbc, lev, prob)
"""
function nonlinear_OLD!( nonlin::AbstractArray,
                     state::IBState{MultiGrid},
                     Γbc::AbstractArray,
                     lev::Int,
                     prob::IBProblem )
   """
   Note factor of 1/2 in boundary conditions absorbed to eliminate
   "lastbc" from Fortran code
   """
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
   C = prob.model.mats.C

   hc = grid.h * 2^( lev - 1);  # Coarse grid spacing

   # Don't need bc's for anything after this, so we can rescale in place
   Γbc .*= 0.25  # Account for scaling between grids

   #-- Alias working memory
      fq = prob.model.work.q3; fq.*=0.0;  # Zero out in case some locations aren't indexed
      fqx, fqy = grid.split_flux(fq);  # These are views into fq - changing them will change fq

      fq1 = prob.model.work.q4; fq1.*=0.0;
      fq1x, fq1y = grid.split_flux(fq1);  # Views of fq1 (which is not itself needed anywhere)

      fq2 = prob.model.work.q5; fq2.*=0.0;
      fq2x, fq2y = grid.split_flux(fq2);  # Views of fq2
   #--

   #-- Compute average fluxes across cells
      Γ = reshape(@view(state.Γ[:, lev]), nx-1, ny-1)
      Qx, Qy = avg_flux(state, prob; lev=lev)
   #--

   #-- Product of flux and circulation
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
   #--

   # Divergence of flux-circulation product
   mul!( nonlin, C', fq )

   # Scaling: 1/hc^2 to convert circulation to vorticity
   #   plus factor of 4 from averaging across cells above
   #  (applies to flux field, but this avoids a second operation)
   nonlin .*= 1/(4*hc^2)

   return nothing
end
