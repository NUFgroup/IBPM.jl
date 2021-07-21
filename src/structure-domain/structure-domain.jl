function make_body( bdy_in, Δx )

    #Currently has functionality for piping into the supported bodies within
    #sample-bodies.jl. Would be nice to add the ability for the user to be able
    #to prescribe the body as a function of arclength...

    #Ugly way of initializing array. Have to replace this...
    #Not sure that the code can handle a sim with multiple Motion types...
    if (:motion in keys(bdy_in[1]))==true
        if bdy_in[1].motion==:static
            bodies = RigidBody{Static}[]
        elseif bdy_in[1].motion==:movinggrid
            bodies = RigidBody{MovingGrid}[]
        elseif bdy_in[1].motion==:movingbody
            bodies = RigidBody{MotionFunction}[]
        end
    end


    for bdy_j in bdy_in
        #read user vals
        if (:motion in keys(bdy_j))==true
            if bdy_j.motion==:static
                motion = Static()
            elseif bdy_j.motion==:movinggrid
                #back out desired grid motion
                if (:motionfcn in keys(bdy_j))==false
                    motion = MovingGrid(t -> 1.0, t -> 0.0, t -> 0.0, t -> 0.0,
                        xc=0.0, yc=0.0 )
                else
                    mfcn = bdy_j.motionfcn

                    if (:U in keys(bdy_j.motionfcn))==false
                        Uf = (t -> 1.0)
                    else
                        Uf = mfcn.U
                    end

                    if (:V in keys(bdy_j.motionfcn))==false
                        Vf = (t -> 0.0)
                    else
                        Vf = mfcn.V
                    end

                    if (:θ in keys(bdy_j.motionfcn))==false
                        θf = (t -> 0.0)
                    else
                        θf = mfcn.θ
                    end

                    if (:θ̇ in keys(bdy_j.motionfcn))==false
                        θ̇f = (t -> 0.0)
                    else
                        θ̇f = mfcn.θ̇
                    end

                    if (:xc in keys(bdy_j.motionfcn))==false
                        xcf = 0.0
                    else
                        xcf = mfcn.xc
                    end

                    if (:yc in keys(bdy_j.motionfcn))==false
                        ycf = 0.0
                    else
                        ycf = mfcn.yc
                    end
                    @show Uf, Vf, θf, θ̇f
                    #put together:
                    motion = MovingGrid(Uf, Vf, θf, θ̇f, xc=xcf, yc=ycf)
                end
            elseif bdy_j.motion==:movingbody
                if (:motionfcn in keys(bdy_j))==false
                    motion = MotionFunction([t -> -t, t -> 0.0, t -> 0.0],
                        [t -> -1.0, t->0.0, t->0.0])
                else
                    mfcn = bdy_j.motionfcn

                    if (:xc in keys(bdy_j.motionfcn))==false
                        xcf = (t -> -t)
                    else
                        xcf = mfcn.xc
                    end

                    if (:yc in keys(bdy_j.motionfcn))==false
                        ycf = (t -> 0.0)
                    else
                        ycf = mfcn.yc
                    end

                    if (:θ in keys(bdy_j.motionfcn))==false
                        θcf = (t -> 0.0)
                    else
                        θcf = mfcn.θ
                    end

                    if (:uc in keys(bdy_j.motionfcn))==false
                        ucf = (t -> -1.0)
                    else
                        ucf = mfcn.uc
                    end

                    if (:vc in keys(bdy_j.motionfcn))==false
                        vcf = (t -> 0.0)
                    else
                        vcf = mfcn.vc
                    end

                    if (:θ̇ in keys(bdy_j.motionfcn))==false
                        θ̇cf = (t -> 0.0)
                    else
                        θ̇cf = mfcn.θ̇
                    end

                    motion = MotionFunction([xcf, ycf, θcf],
                        [ucf, vcf, θ̇cf])
                end
            end
        else
            motion = Static()
        end

        if (:center in keys(bdy_j))==true
            center = bdy_j.center
        else
            center = [0.0; 0.0]
        end

        if bdy_j.type == :cylinder
            body = make_cylinder( bdy_j.lengthscale, Δx,
                center[1], center[2]; motion=motion )
        end

        if bdy_j.type == :plate
            body = make_plate( bdy_j.lengthscale, bdy_j.spec.α, Δx,
                center[1], center[2]; motion=motion )
        end

        bodies = [bodies; body]
        # bodies=[body]
    end

    return bodies
end
