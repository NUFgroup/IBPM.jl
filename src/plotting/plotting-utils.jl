function plot_state( prob, state, t;
    var=:omega, xlims=:auto, ylims=:auto, clims=:auto, clevs=30)

    if (xlims==:auto) & (ylims==:auto)
        #Would be nice to add vel mag
        if var == :omega
            plot_ω(state, prob.model.grid, lev=1, clims=clims, xlims=xlims,
                ylims=ylims, clevs=clevs )
            plot_bodies(state.xb)
        elseif var == :vel
            plt = plot(layout=(2,1))
            plot_u(state, prob.model.grid, lev=1, clims=clims, xlims=xlims,
                ylims=ylims, plt=plt, clevs=clevs )
            plot_bodies(prob.model.bodies, plt=plt )
        elseif var == :psi
            plot_ψ(state, prob.model.grid, lev=glev, clims=clims, xlims=xlims,
                ylims=ylims, clevs=clevs )
            plot_bodies(prob.model.bodies)
        end
    else
        if var == :vel
            plt = plot(layout=(2,1))
        end

        for glev = prob.model.grid.mg : -1 : 1

            #Would be nice to add vel mag
            if var == :omega
                plot_ω(state, prob.model.grid, lev=glev, clims=clims, xlims=xlims,
                    ylims=ylims, clevs=clevs )
                if glev==1
                    plot_bodies(state.xb)
                end
            elseif var == :vel
                plot_u(state, prob.model.grid, lev=glev, clims=clims, xlims=xlims,
                    ylims=ylims, plt=plt, clevs=clevs )
                if glev==1
                    plot_bodies(prob.model.bodies, plt=plt )
                end
            elseif var == :psi
                plot_ψ(state, prob.model.grid, lev=glev, clims=clims, xlims=xlims,
                    ylims=ylims, clevs=clevs )
                if glev==1
                    plot_bodies(prob.model.bodies)
                end
            end

        end
    end


end


function plot_ω(state::IBState, grid::T; lev=1,
    xlims=:auto, ylims=:auto, clims=:auto,
    colorbar=:false, framestyle=:box, clevs=30) where T <:Grid

    h = grid.h
    len = grid.len

    fac = 2.0^(Float64(lev-1))
    δ = h * fac
    xlen = len*fac

    ylen = xlen*(grid.ny/grid.nx)
    offx = fac * len/2.0 - len/2.0 + grid.offx
    offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

    x = range(-offx+δ, xlen-offx-δ, length=grid.nx-1)
    y = range(-offy+δ, ylen-offy-δ, length=grid.ny-1)

    ω = reshape( state.Γ[:, lev], grid.nx-1, grid.ny-1 )' / δ^2.0

    if xlims==:auto
        xlims=(x[1], x[end])
    end

    if ylims==:auto
        ylims=(y[1], y[end])
    end

    if clims ≠ :auto
        ω[ ω .>= clims[2] ] .= clims[2]
        ω[ ω .<= clims[1] ] .= clims[1]
    end

    display(contourf!(x,y, ω, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:RdBu_11),#c=cgrad(:diverging_gkr_60_10_c40_n256), #
        colorbar=colorbar, lw=0, levels=clevs, aspect_ratio=:equal,
        clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle))
end


function plot_bodies(xbv; plt=missing)

    if ismissing(plt)
        for j = 1:length(xbv)
            xb = xbv[j]

            #anoying hack: find interior point of body for fill
            interior = sum(xb[:,2])/length(xb[:,2])

            display(plot!(xb[:, 1], xb[:, 2], lw=3, fillrange=interior,
            fillcolor=:gray, legend=false, color=:gray))
        end
    else
        for j = 1:length(bodies)
            xb = xbv[j]
            for jj = 1 : length(plt)
                #anoying hack: find interior point of body for fill
                interior = sum(xb[:,2])/length(xb[:,2])

                display(plot!(plt[jj],xb[:, 1], xb[:, 2], lw=3,
                fillrange=interior, fillcolor=:gray, legend=false, color=:gray))
            end
        end
    end
end


function plot_u(state::IBState, grid::T; lev=1,
    xlims=:auto, ylims=:auto, clims=:auto,
    colorbar=:false, framestyle=:box, plt=missing, clevs=30) where T <:Grid
    nx, ny = grid.nx, grid.ny
    nu = ny * (nx+1); nv = nx * (ny+1);
    # hc = grid.h*2^(lev-1)
    # x = hc*(1:nx-1) .- grid.offx
    # y = hc*(1:ny) .- grid.offy

    h = grid.h
    len = grid.len

    fac = 2.0^(Float64(lev-1))
    δ = h * fac
    xlen = len*fac

    ylen = xlen*(grid.ny/grid.nx)
    offx = fac * len/2.0 - len/2.0 + grid.offx
    offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

    #-- x vels
    x = range(-offx, xlen-offx, length=grid.nx+1)
    y = range(-offy+δ/2.0, ylen-offy-δ/2.0, length=grid.ny)
    qx = state.q[1:nu, lev] .+ state.q0[1:nu, lev]

    #for some reason the top left velocity on grid 1 is erroneous.
    #extract interior velocities....
    #TODO: look into this...
    u = (reshape( qx, length(x), length(y) ) / δ)'
    u = u[2:end-1, 2:end-1]
    x = x[2:end-1]
    y = y[2:end-1]

    if xlims==:auto
        xlims=(x[1], x[end])
    end

    if ylims==:auto
        ylims=(y[1], y[end])
    end

    if ismissing(plt)
        display(
            contourf!(x, y, u, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
            colorbar=:true, lw=0, levels=clevs, aspect_ratio=:equal,
            clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
            legend=:false, widen=false)
            )
    else

        display(
            contourf!(plt[1], x, y, u, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
            colorbar=:true, lw=0, levels=clevs, aspect_ratio=:equal,
            clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
            legend=:false, widen=false)
            )
    end

    #-- y vels
    x = range(-offx+δ/2.0, xlen-offx-δ/2.0, length=grid.nx)
    y = range(-offy, ylen-offy, length=grid.ny+1)
    qy = state.q[nu.+(1:nv), lev] .+ state.q0[nu.+(1:nv), lev]
    v = (reshape( qy, length(x), length(y) ) / δ)'

    # vplot = contourf!(x, y, v, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
    #     colorbar=:true, lw=0, levels=30, aspect_ratio=:equal,
    #     clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
    #     legend=:false)

    if xlims==:auto
        xlims=(x[1], x[end])
    end

    if ylims==:auto
        ylims=(y[1], y[end])
    end

    if ismissing(plt)
        display(
            contourf!(x, y, v, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
            colorbar=:true, lw=0, levels=clevs, aspect_ratio=:equal,
            clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
            legend=:false, widen=false)
            )
    else
        display(
            contourf!(plt[2], x, y, v, c=cgrad(:seaborn_icefire_gradient),#c=cgrad(:temperaturemap), #c=cgrad(:diverging_gkr_60_10_c40_n256),#c=cgrad(:blackbody),#c=cgrad(:temperaturemap),#,#
            colorbar=:true, lw=0, levels=clevs, aspect_ratio=:equal,
            clims=clims, xlims=xlims, ylims=ylims, framestyle=framestyle,
            legend=:false, widen=false)
            )
    end

end


function plot_ψ(state::IBState, grid::T; lev=1,
    xlims=:auto, ylims=:auto, clims=:auto,
    colorbar=:false, framestyle=:box, clevs=30) where T <:Grid

    h = grid.h
    len = grid.len

    fac = 2.0^(Float64(lev-1))
    δ = h * fac
    xlen = len*fac

    ylen = xlen*(grid.ny/grid.nx)
    offx = fac * len/2.0 - len/2.0 + grid.offx
    offy = fac * (grid.ny*h)/2.0 - (grid.ny*h)/2.0 + grid.offy

    x = range(-offx+δ, xlen-offx-δ, length=grid.nx-1)
    y = range(-offy+δ, ylen-offy-δ, length=grid.ny-1)

    ψ = reshape( state.ψ[:, lev], length(x), length(y) )'
    display(contourf!(x, y, ψ, c=cgrad(:seaborn_icefire_gradient),
        colorbar=:true, lw=1, levels=clevs, aspect_ratio=:equal,
        framestyle=framestyle, clims=clims, xlims=xlims, ylims=ylims,
        background_color=:transparent, foreground_color=:transparent))
end
