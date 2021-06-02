using Plots
# using Colors


function plot_state(state::IBState, grid::T; clims=nothing, lev=1) where T <:Grid
    h = grid.h
    xlen = grid.len
    ylen = grid.len*(grid.ny/grid.nx)
    x = -grid.offx+h:h:xlen-grid.offx-h
    y = -grid.offy+h:h:ylen-grid.offy-h

    ω = reshape( state.Γ[:, lev], grid.nx-1, grid.ny-1 ) / grid.h^2;
    display(contourf(x,y, ω', c=cgrad(:seaborn_icefire_gradient),
        size=(4*grid.nx, 4*grid.ny),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
        colorbar=:false, lw=0, levels=30, aspect_ratio=:equal,
        framestyle=:box, clims=clims))
end

function plot_cyl(body)
    display(plot!(body.xb[:, 1], body.xb[:, 2], lw=0, fillrange=0.0,
        fillcolor=:gray, legend=false))
end


function plot_u(state::IBState, grid::T; lev=1) where T <:Grid
    nx, ny = grid.nx, grid.ny
    nu = ny * (nx-1); nv = nx * (ny-1);
    hc = grid.h*2^(lev-1)
    x = hc*(1:nx-1) .- grid.offx
    y = hc*(1:ny) .- grid.offy

    qx = state.q[1:nu, lev] .+ state.q0[1:nu, lev]
    u = reshape( qx/hc, length(x), length(y) )
    display(contourf(x, y, u', c=cgrad(:seaborn_icefire_gradient),
        size=(4*grid.nx, 4*grid.ny),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
        colorbar=:true, lw=0, levels=30, aspect_ratio=:equal,
        framestyle=:box))
end


function plot_ψ(state::IBState, grid::T) where T <:Grid
    nx, ny = grid.nx, grid.ny
    x = grid.h*(1:nx-1) .- grid.offx
    y = grid.h*(1:ny-1) .- grid.offy

    ψ = reshape( state.ψ, length(x), length(y) )
    display(contourf(x, y, ψ', c=cgrad(:seaborn_icefire_gradient),
        size=(4*grid.nx, 4*grid.ny),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
        colorbar=:true, lw=1, levels=30, aspect_ratio=:equal,
        framestyle=:box))
end
