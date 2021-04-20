"""
!***************************************************************!
!*   nonlinear terms in rotational form                        *!
!***************************************************************!
"""
function nonlinear!( nonlin::AbstractArray,
                     state::IBState{UniformGrid},
                     prob::IBProblem )
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   u, v, ω = grid.u_idx, grid.v_idx, grid.ω_idx  # Shift by -1 if using ω_idx
   C = prob.model.mats.C

   # Alias working memory
   fq = prob.work.q3; fq.*=0.0;
   Q(idx) = state.q[idx] + state.q0[idx]  # Net flux

   # Average fluxes across cell
   U(i, j) = 0.5*( Q(u(i, j)) + Q(u(i,j-1)) )  # j=2:ny, i=1:nx+1
   V(i, j) = 0.5*( Q(v(i, j)) + Q(v(i-1,j)) )  # j=1:ny+1, i=2:nx
   Γ(i, j) = state.Γ[grid.ω_idx(i, j)]   # Indexing function into vorticity

   i=2:nx; j=(2:ny-1)'
   @. fq[u(i,j)] = 0.5*( V(i,j+1)*Γ(i-1,j) + V(i,j)*Γ(i-1,j-1) )
   j=1;  @. fq[u(i,j)] = 0.5*V(i,j+1)*Γ(i-1,j)  # + vavg(i,j)*bc(bottom+i) )
   j=ny; @. fq[v(i,j)] = 0.5*V(i,j)*Γ(i-1,j-1)  # + vavg(i,j+1)*bc(top+i)

   i=2:nx-1; j=(2:ny)'
   @. fq[v(i,j)] = -0.5*( U(i+1,j)*Γ(i,j-1) + U(i,j)*Γ(i-1,j-1) )
   i=1;  @. fq[v(i,j)]  = -0.5*( U(i+1,j)*Γ(i,j-1) )    # + uavg(i,j)*bc(left+j) )
   i=nx; @. fq[v(i,j)]  = -0.5*( U(i+1,j)*Γ(i-1,j-1) )  # + uavg(i+1,j)*bc(right+j)

   mul!( nonlin, C', fq )
   rmul!(nonlin, 1/grid.h^2)  # Scaling factor for differencing
end


function nonlinear!( nonlin::AbstractArray,
                     state::IBState{MultiGrid},
                     Γbc::AbstractArray,
                     prob::IBProblem,
                     lev::Int )
   grid = prob.model.grid;
   nx, ny = grid.nx, grid.ny
   u, v, ω = grid.u_idx, grid.v_idx, grid.ω_idx  # Shift by -1 if using ω_idx
   T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
   C = prob.model.mats.C

    # Coarse grid spacing
    hc = grid.h * 2^( lev - 1);

   # Alias working memory
   fq = @view(prob.work.q3[:, lev]); fq.=0.0;
   Q(idx) = state.q[idx, lev] + state.q0[idx, lev]  # Net flux

   # Average fluxes across cell
   U(i, j) = 0.5*( Q(u(i, j)) + Q(u(i,j-1)) )  # j=2:ny, i=1:nx+1
   V(i, j) = 0.5*( Q(v(i, j)) + Q(v(i-1,j)) )  # j=1:ny+1, i=2:nx
   Γ(i, j) = state.Γ[grid.ω_idx(i, j), lev]   # Indexing function into vorticity

   i=2:nx; j=(2:ny-1)'
   @. fq[u(i,j)] = 0.5*( V(i,j+1)*Γ(i-1,j) + V(i,j)*Γ(i-1,j-1) )
   j=1;  @. fq[u(i,j)] = 0.5*( V(i,j+1)*Γ(i-1,j) + V(i,j)*Γbc[B+i] )
   j=ny; @. fq[v(i,j)] = 0.5*( V(i,j)*Γ(i-1,j-1) + V(i,j+1)*Γbc[T+i] )

   i=2:nx-1; j=(2:ny)'
   @. fq[v(i,j)] = -0.5*( U(i+1,j)*Γ(i,j-1) + U(i,j)*Γ(i-1,j-1) )
   i=1;  @. fq[v(i,j)]  = -0.5*( U(i+1,j)*Γ(i,j-1) + U(i,j)*Γbc[L+j] )
   i=nx; @. fq[v(i,j)]  = -0.5*( U(i+1,j)*Γ(i-1,j-1) + U(i+1,j)*Γbc[R+j] )

   mul!( nonlin, C', fq )    # OK TO HERE

   #println(sum(Γbc.^2))
   #println(sum(state.q[:, lev].^2))
   #println(sum(fq.^2))
   #println(sum(nonlin.^2))

   #rmul!(nonlin, 1/hc^2)  # Scaling factor for differencing
end
