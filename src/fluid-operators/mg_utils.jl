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
    i=(nx÷4).+(0:nx÷2); ibc=1:2:nx+1;
    @views rbc[B.+ibc] = r[i,   ny÷4]
    @views rbc[T.+ibc] = r[i, 3*ny÷4]

    i=(nx÷4).+(1:nx÷2); ibc=2:2:nx;
    @views rbc[B.+ibc] = 0.5*( r[i,   ny÷4] + r[i.-1,   ny÷4] )
    @views rbc[T.+ibc] = 0.5*( r[i, 3*ny÷4] + r[i.-1, 3*ny÷4] )

    j=(ny÷4).+(0:ny÷2); jbc=1:2:ny+1
    @views rbc[L.+jbc] = r[  nx÷4, j]
    @views rbc[R.+jbc] = r[3*nx÷4, j]

    j=(ny÷4).+(1:ny÷2); jbc=2:2:ny
    @views rbc[L.+jbc] = 0.5*( r[  nx÷4, j] + r[  nx÷4, j.-1] )
    @views rbc[R.+jbc] = 0.5*( r[3*nx÷4, j] + r[3*nx÷4, j.-1] )

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
    j=1;    @views r[i,j] += fac*rbc[(B+1).+i]
    j=ny-1; @views r[i,j] += fac*rbc[(T+1).+i]

    j=1:ny-1
    i=1;    @views r[i,j] += fac*rbc[(L+1).+j]
    i=nx-1; @views r[i,j] += fac*rbc[(R+1).+j]

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

    @views Γc[ic, jc] = Γ[i, j] +
                 0.5*( Γ[i.+1,j] + Γ[i,j.+1] + Γ[i.-1,j] + Γ[i,j.-1] ) +
                0.25*( Γ[i.+1,j.+1] + Γ[i.+1,j.-1] + Γ[i.-1,j.-1] + Γ[i.-1,j.+1] )

    Γc = reshape(Γc, grid.nΓ, 1)
    Γ  = reshape(Γ,  grid.nΓ, 1)
end
