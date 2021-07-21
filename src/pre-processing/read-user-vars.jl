function read_user_vars(Δt::Union{Missing,Float64},
    Δx::Union{Missing,Float64},
    freestream::NamedTuple,
    Re::Float64,
    T::Union{Float64,Missing})

    if ismissing(Δx) == true
        Δx = 2.0/Re #default to a grid Re of 2
    end

    #--freestream info
    #(This could probably be made less clunky...)
        xkey = :Ux in keys(freestream)
        ykey = :Uy in keys(freestream)
        anglekey = :inclination in keys(freestream)

        if xkey==false
            Ux = t->t^0.0
        else
            Ux = freestream.Ux
        end

        if ykey==false
            Uy=t->0.0*t^0.0
        else
            Uy = freestream.Uy
        end

        if anglekey==false
            inclination=t->0.0*t^0.0
        else
            inclination = freestream.inclination
        end

        tnames = (:Ux, :Uy, :inclination)
        tvals = (Ux, Uy, inclination)

        freestream = (;zip(tnames, tvals)...)
    #--

    #approximate Δt using Δx if not provided by user
    if ismissing(Δt) == true
        if ismissing(T)==false
            tvect = range(0.0,T,length=5000)
        else
            tvect = range(0.0,20.0,length=5000)
        end
        Umax = sqrt(maximum(freestream.Ux.(tvect))^2.0 .+
            maximum(freestream.Uy.(tvect))^2.0 )

        Δt = 0.1 * Δx/(5.0*Umax) #satisfy CFL 0.2 constraint, with
                                 #safety factor on max velocity

        if ismissing(T)==true
            T = 20.0*Δt
        end
    end

    return Δx, Δt, T, freestream
end
