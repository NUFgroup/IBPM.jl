"""
State types for IBPM
"""

abstract type State end

function circ2_st_vflx!( ψ::AbstractArray,
                         q::AbstractArray,
                         Γ::AbstractArray,
                         model::IBModel{UniformGrid, <:Body} )
    """
    Take vorticity and return velocity flux and streamfunction
    """
    # Solve for streamfcn  ψ = RCinv * Γ
    mul!(ψ, model.mats.RCinv, Γ)

    #--Get velocity flux from curl of stream function
    mul!(q, model.mats.C, ψ)
end


function circ2_st_vflx!( ψ::AbstractArray,
                         q::AbstractArray,
                         Γ::AbstractArray,
                         model::IBModel{MultiGrid, <:Body},
                         ngrids::Int )
         grid = model.grid

         for lev = ngrids:-1:1

             grid.stbc .*= 0.0   # Reset boundary conditions in pre-allocated memory
             #--Solve Poisson problem for ψ

             #BCs for Poisson problem (will be overwritten if mg > 1)
             if ( lev < ngrids )
                get_stfn_BCs!( grid.stbc, @view(ψ[:, lev+1]), model );

                # Solve for streamfcn
                # TODO: in place addition with @view macro
                ψ[:, lev] = Array( model.mats.RCinv * (Γ[:, lev] .+ sum(grid.stbc, dims=2)) )
             else  # don't need bcs for largest grid
             # TODO: Figure out how to use view for in-place
                ψ[:, lev] = model.mats.RCinv * Γ[:, lev]
                #@views mul!(ψ[:, lev], model.mats.RCinv, Γ[:, lev])
             end

             #--Get velocity on first grid from stream function
             @views curl!( q[:,lev], ψ[:,lev], grid.stbc, model );
         end
end

mutable struct IBState{T<:Grid} <: State
        q::Array{Float64, 2}
        q0::Array{Float64, 2}
        Γ::Array{Float64, 2}     # Circulation
        ψ::Array{Float64, 2}     # Streamfunction
        nonlin::Array{Array{Float64, 2}, 1}  # Memory of nonlinear terms
        fb::Array{Array{Float64, 1}, 1}          # Surface stresses
        F̃b::Array{Float64, 1}                    # Body forces * dt
        CD::Array{Float64, 1}    # Drag coefficient
        CL::Array{Float64, 1}    # Lift coefficient
        cfl::Float64
        slip::Float64
end

function init_state(prob::IBProblem)
        grid = prob.model.grid
        mg = (grid isa UniformGrid) ? 1 : grid.mg  # Number of grid levels

        nb, nf = get_body_info(prob.model.bodies)
        return IBState{typeof(grid)}(
        zeros(grid.nq, mg),   # Flux
        zeros(grid.nq, mg),   # Background flux
        zeros(grid.nΓ, mg),   # Circulation
        zeros(grid.nΓ, mg),   # Streamfunction
        [zeros(grid.nΓ, mg) for i=1:length(prob.scheme.β)],   # Nonlinear
        [zeros(nf[i]) for i=1:length(nf)],   # Surface stresses on each body
        zeros(sum(nf)),
        zeros(length(prob.model.bodies)),
        zeros(length(prob.model.bodies)),
        0.0,
        0.0
        )
end


function base_flux!(state::IBState{UniformGrid},
                    grid::UniformGrid,
                    Uinf::Float64,
                    α::Float64)
        """
        Initialize irrotational freestream flux
        """
    m = grid.nx;
    n = grid.ny;
    state.q0[ 1:(m-1)*n ] .= Uinf * grid.h * cos(α);  # x-flux
    state.q0[ (m-1)*n+1:end ] .= Uinf * grid.h * sin(α);  # y-flux
end


function base_flux!(state::IBState{MultiGrid},
                    grid::MultiGrid,
                    Uinf::Float64,
                    α::Float64)
        """
        Initialize irrotational freestream flux
        """
    m = grid.nx;
    n = grid.ny;
    for lev = 1 : grid.mg
        # Coarse grid spacing
        hc = grid.h * 2^( lev - 1 );

        # write fluid velocity flux in body-fixed frame
        state.q0[ 1:(m-1)*n, lev ] .= Uinf * hc * cos(α);      # x-flux
        state.q0[ (m-1)*n+1:end, lev ] .= Uinf * hc * sin(α);  # y-flux
    end
end
