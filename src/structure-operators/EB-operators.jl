function get_structural_mats(bodies)
    #If there's a deforming body, there can only be one body (for now...)
    btypes = [typeof(bodies[i]) for i=1:length(bodies)]
    if sum(btypes .<: DeformingBody) > 0
        @show cp = 1
        @assert length(bodies) == 1
    end

    #get structural matrices (or fill with Missing if rigid body/ies)
    get_structural_mats(typeof(bodies[1]), bodies[1])
end

function get_structural_mats(::Type{RigidBody{W}} where W <: Motion, body)
    return missing, missing, missing
end

function get_structural_mats(::Type{DeformingBody{W}} where W <: Motion, body)

    nb = length(body.xb[:,1]) # number of body points
    nel = nb - 1 #number of finite elements
    m, kb, ke = body.m, body.kb, body.ke


    Mmat = zeros(Float64,3*nb, 3*nb)
    Kmat = zeros(Float64,3*nb, 3*nb)
    Qmat = zeros(Float64,3*nb, 3*nb)
    Rmat = zeros(Float64,3*nb, 3*nb)

    #We will build these matrices by element and assemble in a loop
    for i_el = 1:nel

        Δx = body.ds[i_el]

        el_ind = @. (i_el-1)*3 + ( 1:6 ) #indices corresponding with the 6 unknowns
                                      #associated w/ each element

        #--build local mass and stiffnes matrices
            #Consistent
            M_e = m[i_el]*Δx/420.0 * [ 140.0 0.0 0.0 70.0 0.0 0.0;
            0.0 156.0 22.0*Δx 0.0 54.0 -13.0*Δx;
            0.0 22.0*Δx 4.0*Δx^2.0 0.0 13.0*Δx -3.0*Δx^2.0;
            70.0 0.0 0.0 140.0 0.0 0.0;
            0.0 54.0 13.0*Δx 0.0 156.0 -22.0*Δx;
            0.0 -13.0*Δx -3.0*Δx^2.0 0.0 -22.0*Δx 4*Δx^2.0
            ]

            #in local coord sys. In the geometrically nonlinear sims, rotation
            #matrices will be used to express global displacements in a local
            #coordinate system, and dynamics will be linear about that
            K_e = 1.0/(Δx^3.0) *
            [ ke[i_el]*Δx^2.0 0.0 0.0 -ke[i_el]*Δx^2.0 0.0 0.0;
              0.0 kb[i_el]*12.0 kb[i_el]*6.0*Δx 0.0 -kb[i_el]*12.0 kb[i_el]*6*Δx;
              0.0 kb[i_el]*6.0*Δx kb[i_el]*4.0*Δx^2.0 0.0 -kb[i_el]*6.0*Δx kb[i_el]*2.0*Δx^2.0;
              -ke[i_el]*Δx^2.0 0.0 0.0 ke[i_el]*Δx^2.0 0.0 0.0;
              0.0 -kb[i_el]*12.0 -kb[i_el]*6.0*Δx 0.0 kb[i_el]*12.0 -kb[i_el]*6.0*Δx;
              0.0 kb[i_el]*6.0*Δx kb[i_el]*2.0*Δx^2.0 0.0 -kb[i_el]*6.0*Δx kb[i_el]*4.0*Δx^2.0
            ]
        #--

        #--Assemble into global matrices
            for i = 1:6 #Add contributions for each DOF in the element
                i_ind = el_ind[i]

                for j = 1:6
                    j_ind = el_ind[j]
                    Mmat[i_ind,j_ind] = Mmat[i_ind,j_ind] + M_e[i,j]
                    Kmat[i_ind,j_ind] = Kmat[i_ind,j_ind] + K_e[i,j]
                    Qmat[i_ind,j_ind] = Mmat[i_ind,j_ind] / m[i_el]
                end
            end
        #--

        #--Rotation matrix
        θ = atan2( body.xb[i_el,2]-body.xref[i_el,2],
            body.xb[i_el,1]-body.xref[i_el,1] )
        Rind = 3*(i_el-1)+1
        Rmat[Rind:Rind+1, Rind:Rind+1] = [cos(θ) -sin(θ); sin(θ) cos(θ)]
        #--

    end

    #--Rotation matrix for last body point
        θ = atan2( body.xb[nb,2]-body.xref[nb,2],
            body.xb[nb,1]-body.xref[nb,1] )
        Rind = 3*(nb-1)+1
        Rmat[Rind:Rind+1, Rind:Rind+1] = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    #--


    #--

    Kmat = Rmat' * Kmat * Rmat

    #account for BCs
    for BCinf in body.BCinfo
        for j in 1:length(BCinf.BCvec)
            Mmat[3*(BCinf.BCnode-1)+j,:] = 0.0
            Kmat[3*(BCinf.BCnode-1)+j,:] = 0.0
            Qmat[3*(BCinf.BCnode-1)+j,:] = 0.0
            Kmat[3*(BCinf.BCnode-1)+j,3*(BCinf.BCnode-1)+j] = 1.0

        end
    end

    # MiKmat = Mmat \ Kmat


    return Mmat, Kmat, Qmat
end
