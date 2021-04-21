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
   """
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   u, v = grid.u_idx, grid.v_idx
   C = prob.model.mats.C

   # Alias working memory
   fq = prob.model.work.q3; fq.*=0.0;
   fqx, fqy = grid.split_flux(fq);
   qx, qy = grid.split_flux(state.q);
   q0x, q0y = grid.split_flux(state.q0);

   # Average net fluxes across cell
   U(i, j) = 0.5*( qx[i,j] .+ q0x[i,j] .+ qx[i,j.-1] .+ q0x[i,j.-1]  )
   V(i, j) = 0.5*( qy[i,j] .+ q0y[i,j] .+ qy[i.-1,j] .+ q0y[i.-1,j]  )
   Γ = reshape(state.Γ, nx-1, ny-1)

   i=2:nx; j=2:ny-1
   @views fqx[i,j] = 0.5*( V(i,j.+1).*Γ[i.-1,j] .+ V(i,j).*Γ[i.-1,j.-1] )
   j=1;  @views fqx[i,j] = 0.5*( V(i,j.+1).*Γ[i.-1,j] )
   j=ny; @views fqx[i,j] = 0.5*( V(i,j).*Γ[i.-1,j.-1] )

   i=2:nx-1; j=2:ny
   @views fqy[i,j] = -0.5*( U(i.+1,j).*Γ[i,j.-1] .+ U(i,j).*Γ[i.-1,j.-1] )
   i=1;  @views fqy[i,j]  = -0.5*( U(i.+1,j).*Γ[i,j.-1] )
   i=nx; @views fqy[i,j]  = -0.5*( U(i.+1,j).*Γ[i.-1,j.-1] )

   mul!( nonlin, C', fq )    # OK TO HERE
   rmul!(nonlin, 1/grid.h^2)  # Scaling factor for differencing
end

function nonlinear!( nonlin::AbstractArray,
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
   u, v = grid.u_idx, grid.v_idx  # Shift by -1 if using ω_idx
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
   C = prob.model.mats.C

    # Coarse grid spacing
    hc = grid.h * 2^( lev - 1);

   # Alias working memory
   fq = @view(prob.model.work.q3[:, lev]); fq.*=0.0;
   Q(idx) = state.q[idx, lev] + state.q0[idx, lev]  # Net flux

   # Average fluxes across cell
   U(i, j) = 0.5*( Q(u(i, j)) + Q(u(i,j-1)) )  # j=2:ny, i=1:nx+1
   V(i, j) = 0.5*( Q(v(i, j)) + Q(v(i-1,j)) )  # j=1:ny+1, i=2:nx
   #Γ(i, j) = state.Γ[grid.ω_idx(i, j), lev]   # Indexing function into vorticity
   Γ = reshape(@view(state.Γ[:, lev]), nx-1, ny-1)

   i=2:nx; j=2:ny-1
   @. fq[u(i,j')] = 0.5*( V(i,j'+1)*Γ[i-1,j] + V(i,j')*Γ[i-1,j-1] )
   j=1;  @. fq[u(i,j)] = 0.5*( V(i,j+1)*Γ[i-1,j] + 0.5*V(i,j)*Γbc[B+i] )
   j=ny; @. fq[u(i,j)] = 0.5*( V(i,j)*Γ[i-1,j-1] + 0.5*V(i,j+1)*Γbc[T+i] )

   i=2:nx-1; j=2:ny
   @. fq[v(i,j')] = -0.5*( U(i+1,j')*Γ[i,j-1] + U(i,j')*Γ[i-1,j-1] )
   i=1;  @. fq[v(i,j')]  = -0.5*( U(i+1,j')*Γ[i,j'-1] + 0.5*U(i,j')*Γbc[L+j'] )
   i=nx; @. fq[v(i,j')]  = -0.5*( U(i+1,j')*Γ[i-1,j'-1] + 0.5*U(i+1,j')*Γbc[R+j'] )

   mul!( nonlin, C', fq )    # OK TO HERE
   rmul!(nonlin, 1/hc^2)  # Scaling factor for differencing
end



function nonlinear_NEW!( nonlin::AbstractArray,
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
   u, v = grid.u_idx, grid.v_idx  # Shift by -1 if using ω_idx
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
   C = prob.model.mats.C

   # Coarse grid spacing
   hc = grid.h * 2^( lev - 1);

   # Alias working memory
   fq = @view(prob.model.work.q3[:, lev]); fq.*=0.0;
   fqx, fqy = grid.split_flux(fq);
   qx, qy = grid.split_flux(state.q, lev);
   q0x, q0y = grid.split_flux(state.q0, lev);

   # Average net fluxes across cell
   U(i, j) = 0.5*( qx[i,j] .+ q0x[i,j] .+ qx[i,j.-1] .+ q0x[i,j.-1]  )
   V(i, j) = 0.5*( qy[i,j] .+ q0y[i,j] .+ qy[i.-1,j] .+ q0y[i.-1,j]  )
   Γ = reshape(@view(state.Γ[:, lev]), nx-1, ny-1)

   i=2:nx; j=2:ny-1
   @views fqx[i,j] = 0.5*( V(i,j.+1).*Γ[i.-1,j] .+ V(i,j).*Γ[i.-1,j.-1] )
   j=1;  @views fqx[i,j] = 0.5*( V(i,j.+1).*Γ[i.-1,j] .+ 0.5*V(i,j).*Γbc[B.+i] )
   j=ny; @views fqx[i,j] = 0.5*( V(i,j).*Γ[i.-1,j.-1] .+ 0.5*V(i,j.+1).*Γbc[T.+i] )

   i=2:nx-1; j=2:ny
   @views fqy[i,j] = -0.5*( U(i.+1,j).*Γ[i,j.-1] .+ U(i,j).*Γ[i.-1,j.-1] )
   i=1;  @views fqy[i,j]  = -0.5*( U(i.+1,j).*Γ[i,j.-1] + 0.5.*U(i,j).*Γbc[L.+j] )
   i=nx; @views fqy[i,j]  = -0.5*( U(i.+1,j).*Γ[i.-1,j.-1] + 0.5*U(i.+1,j).*Γbc[R.+j] )

   mul!( nonlin, C', fq )    # OK TO HERE
   rmul!(nonlin, 1/hc^2)  # Scaling factor for differencing
end
