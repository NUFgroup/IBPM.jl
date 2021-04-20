"""
!***************************************************************!
!*   given vorticity on a larger, coarser mesh, interpolate    *!
!*   it's values to the edge of a smaller, finer mesh          *!
!***************************************************************!
"""
function get_bc!(rbc, r, fac, grid)
    T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
    nx, ny = grid.nx, grid.ny
    r = reshape(r, nx-1, ny-1)

    # Get interpolated boundary conditions on finer grid
    i=0:2:nx
    @. rbc[B+i+1] = r[nx÷4+i÷2,   ny÷4]
    @. rbc[T+i+1] = r[nx÷4+i÷2, 3*ny÷4]

    i=1:2:nx-1
    @. rbc[B+i+1] = 0.5*( r[nx÷4+(i+1)÷2,   ny÷4] + r[nx÷4-1+(i+1)÷2,   ny÷4] )
    @. rbc[T+i+1] = 0.5*( r[nx÷4+(i+1)÷2, 3*ny÷4] + r[nx÷4-1+(i+1)÷2, 3*ny÷4] )

    j=0:2:ny
    @. rbc[L+j+1] = r[  nx÷4, ny÷4+j÷2]
    @. rbc[R+j+1] = r[3*nx÷4, ny÷4+j÷2]

    j=1:2:ny-1
    @. rbc[L+j+1] = 0.5*( r[  nx÷4, ny÷4+(j+1)÷2] + r[  nx÷4, ny÷4-1+(j+1)÷2] )
    @. rbc[R+j+1] = 0.5*( r[3*nx÷4, ny÷4+(j+1)÷2] + r[3*nx÷4, ny÷4-1+(j+1)÷2] )

    rbc .*= fac

    r = reshape(r, grid.nΓ, 1)
end

"""
!***************************************************************!
!*   given vorticity at edges of domain, rbc, (from larger,    *!
!*   coarser mesh), add values to correct laplacian of         *!
!*   vorticity  on the (smaller, finer) domain, r.             *!
!***************************************************************!

r is a vorticity-like array of size (nx-1)×(ny-1)
"""
function apply_bc!(r, rbc, fac, grid)
    T, B, L, R = grid.TOP, grid.BOT, grid.LEFT, grid.RIGHT  # Constant offsets for indexing BCs
    nx, ny = grid.nx, grid.ny
    r = reshape(r, nx-1, ny-1)

    # add bc's from coarser grid
    i=1:nx-1
    j=1;    @. r[i,j] += fac*rbc[B+i+1]
    j=ny-1; @. r[i,j] += fac*rbc[T+i+1]

    j=1:ny-1
    i=1;    @. r[i,j] += fac*rbc[L+j+1]
    i=nx-1; @. r[i,j] += fac*rbc[R+j+1]

    r = reshape(r, grid.nΓ, 1)
end

"""
!***************************************************************!
!*   given vorticity on a smaller, fine mesh, (rhs) interp.    *!
!*   values to the center region of a larger, coarser mesh     *!
!*   (crhs).  The values outside the center region are         *!
!*   not unmodified. Result is placed in arhs                  *!
!***************************************************************!
"""
function coarsify!(Γc::AbstractVector,
                   Γ::AbstractVector,
                   grid::MultiGrid)
    nx, ny = grid.nx, grid.ny
    Γc = reshape(Γc, nx-1, ny-1)
    Γ  = reshape(Γ,  nx-1, ny-1)

    ic = nx÷2 .+ ( (-nx÷4+1):(nx÷4-1) )  # Coarse cell centers
    i = 2 .* ic .- nx÷2
    jc = ny÷2 .+ ( (-ny÷4+1):(ny÷4-1) )'
    j = 2 .* jc .- ny÷2

    @. Γc[ic, jc] = Γ[i, j] +
                 0.5*( Γ[i+1,j] + Γ[i,j+1] + Γ[i-1,j] + Γ[i,j-1] ) +
                0.25*( Γ[i+1,j+1] + Γ[i+1,j-1] + Γ[i-1,j-1] + Γ[i-1,j+1] )

    Γc = reshape(Γc, grid.nΓ, 1)
    Γ  = reshape(Γ,  grid.nΓ, 1)
end


"""
FUNCTION coarsify( crhs, rhs ) RESULT( arhs )

  !***************************************************************!
  !*   given vorticity on a smaller, fine mesh, (rhs) interp.    *!
  !*   values to the center region of a larger, coarser mesh     *!
  !*   (crhs).  The values outside the center region are         *!
  !*   not(? JC) unmodified. Result is placed in arhs                  *!
  !***************************************************************!
  REAL(KIND(0.D0)), DIMENSION(:,:)                     :: crhs, rhs
  REAL(KIND(0.D0)), DIMENSION(SIZE(rhs,1),SIZE(rhs,2)) :: arhs
  INTEGER                                              :: i,j,indi,indj

  arhs = crhs
  DO j=-n/4+1,n/4-1
     indj = n/2+2*j
     DO i=-m/4+1,m/4-1
        indi = m/2+2*i
        arhs(m/2+i,n/2+j) = rhs(indi  ,indj)   + &
                    0.5d0*( rhs(indi+1,indj)   + rhs(indi  ,indj+1)   + &
                            rhs(indi-1,indj)   + rhs(indi  ,indj-1) ) + &
                   0.25d0*( rhs(indi+1,indj+1) + rhs(indi+1,indj-1)   + &
                            rhs(indi-1,indj-1) + rhs(indi-1,indj+1) )
      ENDDO
   ENDDO

END FUNCTION coarsify
"""
