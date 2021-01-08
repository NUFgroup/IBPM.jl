include("src/ibpm.jl")


# Define grid
nx = 200  # num of x points on finest domain
ny = 200  # num of y points on finest domain (length in y-dirn is len/m * n )
mg = 5   # num domains
offx = 1.0; # offset in x dirn (on fine domain, x-grid runs from -offx to len-offx.
offy = 2.0; # offset in y dirn (same as offx but in y-dirn)

len = 4.0  # length of domain in x-direction

# MultiGrid
grid = ibpm.make_grid(nx, ny, offx, offy, len; mg=mg)

# Create a NACA airfoil
bodies = [ibpm.make_naca(0.25, 50, "2412", ibpm.Static())]

# Create full IBPM problem
Re = 40.0
dt = 1e-2

prob = ibpm.init_prob(grid, bodies, Re, dt);
state = ibpm.init_state(prob.model);

Uinf = 1.0;   # Free-stream flow
α = 10.0 * π/180.0;      # Angle of attack
ibpm.base_flux!(state, grid, Uinf, α)  # Initialize irrotational base flux

function run_sim(it_stop)
        for it=1:it_stop
                ibpm.advance!(it, state, prob)
                println([it, state.CD, state.CL, state.cfl])
        end
end


run_sim(1)  # First step to compile
println(@elapsed run_sim(350))
#ibpm.plot_state(state, prob.model.grid)


function plot_state(state::ibpm.IBPMState, grid::T) where T <:ibpm.Grid
    h = grid.h
    xlen = grid.len
    ylen = grid.len*(grid.ny/grid.nx)
    x = -grid.offx+h:h:xlen-grid.offx-h
    y = -grid.offy+h:h:ylen-grid.offy-h

    ω = reshape( state.Γ[:, 1], grid.nx-1, grid.ny-1 ) / h^2;
    heatmap(x, y, ω', c=:RdBu,
        size=(4*grid.nx, 4*grid.ny), colorbar=:false)
end


function plot_body(body)
    plot!(body.xb[:, 1], body.xb[:, 2], color=:black)
end
