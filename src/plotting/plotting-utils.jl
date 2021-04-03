using Plots
# using Colors


function plot_state(state::IBState, grid::T) where T <:Grid
    h = grid.h
    xlen = grid.len
    ylen = grid.len*(grid.ny/grid.nx)
    x = range(-grid.offx+h, xlen-grid.offx-h, length=grid.nx-1)
    y = range(-grid.offy+h, ylen-grid.offy-h, length=grid.ny-1)

    ω = reshape( state.Γ[:, 1], grid.nx-1, grid.ny-1 ) / grid.h^2

    display(contourf(x,y, ω', c=cgrad(:seaborn_icefire_gradient),
        size=(4*grid.nx, 4*grid.ny),
        colorbar=:false, lw=0, levels=30, aspect_ratio=:equal,
        framestyle=:box))

end

function plot_body(body)
    #anoying hack: find interior point of body for fill
    interior = sum(body.xb[:,2])/length(body.xb[:,2])
    
    display(plot!(body.xb[:, 1], body.xb[:, 2],
        linecolor=:white, lw=2,
        fillrange=0.2, fillcolor=:white, legend=false))
end
