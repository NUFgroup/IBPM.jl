using Plots
# using Colors


function plot_state(state::IBState, grid::T) where T <:Grid
    h = grid.h
    xlen = grid.len
    ylen = grid.len*(grid.ny/grid.nx)
    x = -grid.offx+h:h:xlen-grid.offx-h
    y = -grid.offy+h:h:ylen-grid.offy-h

    ω = reshape( state.Γ[:, 1], grid.nx-1, grid.ny-1 ) / grid.h^2;
    display(contourf(x,y, ω', c=cgrad(:seaborn_icefire_gradient),
        size=(4*grid.nx, 4*grid.ny),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
        colorbar=:false, lw=0, levels=30, aspect_ratio=:equal,
        framestyle=:box,
        clim=(-5, 5)))

end

function plot_body(body)
    display(plot!(body.xb[:, 1], body.xb[:, 2], lw=0, fillrange=0.0,
        fillcolor=:gray, legend=false))
end
