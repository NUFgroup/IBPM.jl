mutable struct user_var{T}
	val :: T
	t :: Union{Vector{Float64},
		StepRangeLen{ Float64, Base.TwicePrecision{Float64},
		Base.TwicePrecision{Float64} }
		}
	fcn :: Any
end


function init_save_data( t, save_info, state )

	"""
	Initialize a data structure that will store the desired user data as a named
	Tuple, with names prescribed by user and vectors/arrays deduced by user info
	"""

	#Initial save variables as data structure

	#if data structure unspecified, the full state data struct will be saved
	#at the last time instance to facilitate a restart
	if ismissing(save_info) == true
			data = user_var( Any[]  ,
					t[end], (t,state)->state )

	else
		#Initialize data structure data as a vector of user_var structs
		# data = user_var[]
		data = user_var[]

		Δt = t[2]-t[1] #assumes constant time step size for now

		for j = 1:length( save_info.save_fcns )
			if (:save_times in keys(save_info))==true
				sv_ti = save_info.save_times[j]

				#save every time
				if (ismissing(sv_ti) == true) | (Float64(sv_ti) == 0.0)
						tsv = t

				#--save over subintervals

				#user prescribes single Float
				elseif length(sv_ti) == 1

						#round sv_ti to the nearest value that will be encountered by time
						#advancement given choice of Δt
						sv_ti_rd = round(sv_ti/Δt) * Δt
						tsv = t[1] : sv_ti_rd : t[end]

				#user prescribes array of values
				elseif length(sv_ti) > 1

						#for now, just picking the time instance during the sim
						#that is nearest to the desired save instance. would be cool
						#to eventually interpolate at the desired save instances
						t_sv = Float64[]
						for ttmp in sv_ti
							ind_min = argmin( t .- ttmp )
							push!( t_sv, t[ind_min] )
						end

						tsv = transpose(unique( tsv ))

				end
				#--

			else
				tsv = [t[end]]
			end

			#Initialize data structure here...
			#Utilize knowledge of data type if specified by user
			if (:save_types in keys(save_info))==true
				data = [data; user_var( save_info.save_types[j][]  ,
						tsv, save_info.save_fcns[j] )
						]
			else
				data = [data; user_var( Any[]  ,
						tsv, save_info.save_fcns[j] )
						]
			end

		end

	end

	return data


end

function save_data!(tnow, t, prob, state, data)

		for j = 1 : length(data)

			#Saving?
			tdiff = @. abs( data[j].t - tnow )
			ind_diff = @. (tdiff+1.0 ≈ 1.0) #artificially add 1 to set sig figs
			if sum(ind_diff) > 0
				push!(data[j].val, deepcopy(data[j].fcn( tnow, state )) )
			end

		  end

end #function
