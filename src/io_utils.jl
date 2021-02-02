using MAT

function save_model(filename, model::IBModel)
      matwrite(filename, Dict(
          "Re" => model.Re,
          "nx" => model.grid.nx,
          "ny" => model.grid.ny,
          "offx" => model.grid.offx,
          "offy" => model.grid.offy,
          "len" => model.grid.len,
          "mg" => model.grid.mg
      ); compress = true)
end

function save_state(filename, state::State)
      matwrite(filename, Dict(
          "q" => state.q,
          "circ" => state.Γ,
          "stfn" => state.ψ
      ); compress = true)
end


function load(filename)
    return matread(filename)
end
