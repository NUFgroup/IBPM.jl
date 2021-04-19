"""
!***************************************************************!
!*   nonlinear terms in rotational form                        *!
!***************************************************************!

#function nonlinear!(fq, Γ, q, q0, grid) # INCLUDE BC HERE
"""
function nonlinear!( nonlin::AbstractArray,
                      state::IBState{UniformGrid},
                      prob::IBProblem )
   grid = prob.model.grid
   nx, ny = grid.nx, grid.ny
   u, v, ω = grid.u_idx, grid.v_idx, grid.ω_idx  # Shift by -1 if using ω_idx
   C = prob.model.mats.C

   # Alias working memory
   fq = prob.work.q3
   Γ = state.Γ

   Q(idx) = state.q[idx] + state.q0[idx]  # Net flux
   # Average fluxes across cell
   uavg(i, j) = 0.5*( Q(u(i, j)) + Q(u(i,j-1)) )  # j=2:ny, i=1:nx+1
   vavg(i, j) = 0.5*( Q(v(i, j)) + Q(v(i-1,j)) )  # j=1:ny+1, i=2:nx

   for j=2:ny-1
      for i=2:nx
         fq[u(i, j)] = 0.5*(vavg(i,j+1)*Γ[ω(i-1, j)] +
                            vavg(i,j)*Γ[ω(i-1,j-1)])
      end
   end

   for i=2:nx
      fq[u(i,1)]  = 0.5*( vavg(i,2)*Γ[ω(i-1,1)] )     # + vavg(i,j)*bc(bottom+i) )
      fq[u(i,ny)] = 0.5*( vavg(i,ny)*Γ[ω(i-1,ny-1)])  # + vavg(i,j+1)*bc(top+i)
   end
   # note...we don't need result for i=1 or 1=m+1 since it isn't needed by rot

   for j=2:ny
      fq[v(1, j)]  = -0.5*( uavg(2,j)*Γ[ω(1,j-1)] )  # + uavg(i,j)*bc(left+j) )
      for i=2:nx-1
         fq[v(i, j)] = -0.5*( uavg(i+1,j)*Γ[ω(i,j-1)] +
                              uavg(i,j)*Γ[ω(i-1,j-1)] )
      end
      fq[v(nx,j)]  = -0.5*( uavg(nx+1,j)*Γ[ω(nx-1,j-1)] )  # + uavg(i+1,j)*bc(right+j)
   end
   # note...we don't need result for j=1 or j=n+1 since it isn't needed by rot

   mul!( nonlin, C', fq )
end
