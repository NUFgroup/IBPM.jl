using Plots


function plot_state(state::IBState, grid::T) where T <:Grid
    ω = reshape( state.Γ[:, 1], grid.nx-1, grid.ny-1 ) / grid.h^2;
    heatmap(ω', c=:RdBu, size=(4*grid.nx, 4*grid.ny), colorbar=:false)

    #ψ = reshape( state.ψ[:, 1], grid.nx-1, grid.ny-1 ) / grid.h^2;
    #contour!(ψ', color=:black)
end
