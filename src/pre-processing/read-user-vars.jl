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
            if (typeof(freestream.Ux)==Float64 ||
                typeof(freestream.Ux)==Int64) == true #convert to function if constant
                Ux = t->freestream.Ux*t^0.0
            else
                Ux = freestream.Ux
            end
        end

        if ykey==false
            Uy=t->0.0*t^0.0
        else
            if (typeof(freestream.Uy)==Float64 ||
                typeof(freestream.Uy)==Int64) == true #convert to function if constant
                Uy = t->freestream.Uy*t^0.0
            else
                Uy = freestream.Uy
            end
        end

        if anglekey==false
            inclination=t->0.0*t^0.0
        else
            if (typeof(freestream.inclination)==Float64 ||
                typeof(freestream.inclination)==Int64) == true
                inclination = t->freestream.inclination*t^0.0
            else
                inclination = freestream.inclination
            end
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
