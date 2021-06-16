mutable struct user_var{T}
	val :: T
	t :: Union{Vector{Float64},
		StepRangeLen{ Float64, Base.TwicePrecision{Float64},
		Base.TwicePrecision{Float64} }
		}
	fcn :: Any
end


# mutable struct user_data{T}
#
# 	uvar1 :: user_var{T}
# 	uvar2 :: user_var{T}
# 	uvar3 :: user_var{T}
# 	uvar4 :: user_var{T}
# 	uvar5 :: user_var{T}
# 	uvar6 :: user_var{T}
# 	uvar7 :: user_var{T}
# 	uvar8 :: user_var{T}
# 	uvar9 :: user_var{T}
# 	uvar10 :: user_var{T}
# 	uvar11 :: user_var{T}
# 	uvar12 :: user_var{T}
# 	uvar13 :: user_var{T}
# 	uvar14 :: user_var{T}
# 	uvar15 :: user_var{T}
#
# end


function init_save_data( t, save_info, state )

	"""
	Initialize a data structure that will store the desired user data as a named
	Tuple, with names prescribed by user and vectors/arrays deduced by user info
	"""

	#Initial save variables as data structure
	if ismissing(save_info) == true
			data = Float64[]

	else
		# #Initialize data structure data as a vector of user_var structs
		data = user_var[]

		Δt = t[2]-t[1] #assumes constant time step size for now

		for j = 1:length( save_info.save_times )
			sv_ti = save_info.save_times[j]

			#save every time
			if (isempty(sv_ti) == true) | (ismissing(sv_ti) == true) |
				(sv_ti == 0) | (sv_ti == 0.0)
					t_cnt = length(t)
					tsv = t

			#--save over subintervals

			#user prescribes single Float
			elseif length(sv_ti) == 1

					#round sv_ti to the nearest value that will be encountered by time
					#advancement given choice of Δt
					sv_ti_rd = round(sv_ti/Δt) * Δt
					tsv = t[1] : sv_ti_rd : t[end]
					t_cnt = length( tsv )

			#user prescribes array of values
			elseif length(sv_ti) > 1

					t_sv = Float64[]
					for ttmp in sv_ti
							ind_min = argmin( t .- ttmp )
							push!( t_sv, t[ind_min] )
					end

					tsv = transpose(unique( tsv ))
					t_cnt = length( tsv )

			end
			#--

			#Put info about data in here...

			@show [save_info.save_types[j][] ]

			push!(data,
				user_var(
					[ save_info.save_types[j][] ] ,
					tsv, save_info.save_fcns[j] )
					)
			# push!(data,
			# 	user_var(
			# 		zeros(size(
			# 		save_info.save_fcns[j](0.0,state))...,t_cnt
			# 		) , tsv, save_info.save_fcns[j] )
			# 		)

		end

	end

	return data


end

# function init_save_data( t, save_info, state )
#
# 	"""
# 	Initialize a data structure that will store the desired user data as a named
# 	Tuple, with names prescribed by user and vectors/arrays deduced by user info
# 	"""
#
# 	#Initial save variables as data structure
# 	if ismissing(save_info) == true
# 			data = Float64[]
#
# 	else
# 		#Initialize data structure data as an array of user_var structs
# 		data :: Array{user_var,2}
#
# 		Δt = t[2]-t[1] #assumes constant time step size for now
#
# 		#define the number of save instances for each
# 		tsv = [] #not specifying type; typically a mix of StepRange and Vector{Float}
# 		tcnt = Int64[]
#
# 		for j = 1:length( save_info.save_times )
# 			sv_ti = save_info.save_times[j]
# 			@show sv_ti
#
# 			#save every time
# 			if (isempty(sv_ti) == true) | (ismissing(sv_ti) == true) |
# 				(sv_ti == 0) | (sv_ti == 0.0)
# 					tcnt = [tcnt; length(t)]
# 					tsv = [tsv; [t]]
#
# 			#--save over subintervals
#
# 			#user prescribes single Float
# 			elseif length(sv_ti) == 1
#
# 					#round sv_ti to the nearest value that will be encountered by time
# 					#advancement given choice of Δt
# 					sv_ti_rd = round(sv_ti/Δt) * Δt
# 					tsv_v = t[1] : sv_ti_rd : t[end]
# 					tsv = [tsv; [tsv_v] ]
# 					t_cnt = [tcnt; length( tsv_v )]
#
# 			#user prescribes array of values
# 			elseif length(sv_ti) > 1
#
# 					t_sv = Float64[]
# 					for ttmp in sv_ti
# 							ind_min = argmin( t .- ttmp )
# 							push!( t_sv, t[ind_min] )
# 					end
#
# 					t_sv = transpose(unique( t_sv ))
# 					tsv = [tsv; [t_sv] ]
# 					t_cnt = [tcnt; length( t_sv )]
#
# 			end
# 			#--
#
# 			# #Put info about data in here...
# 			data[j] = user_var{}
#
# 		end
#
# 	end
#
# end

#
# function init_save_data( t, save_info, state )
#
# 	"""
# 	Initialize a data structure that will store the desired user data as a named
# 	Tuple, with names prescribed by user and vectors/arrays deduced by user info
# 	"""
#
# 	#Initial save variables as data structure
# 	if ismissing(save_info) == true
# 			data = Float64[]
#
# 	else
# 		Δt = t[2]-t[1] #assumes constant time step size for now
#
# 	  #define the number of save instances for each
# 		tsv = [] #not specifying type; will be a mix of StepRange and Vector{Float}
# 		tcnt = Int64[]
#
# 	  for sv_ti in save_info.save_times
#
# 			@show sv_ti
#
# 			#save every time
# 			if (isempty(sv_ti) == true) | (ismissing(sv_ti) == true) |
# 				(sv_ti == 0) | (sv_ti == 0.0)
# 					tcnt = [tcnt; length(t)]
# 					tsv = [tsv; [t]]
#
# 			#--save over subintervals
#
# 			#user prescribes single Float
# 			elseif length(sv_ti) == 1
#
# 					#round sv_ti to the nearest value that will be encountered by time
# 					#advancement given choice of Δt
# 					sv_ti_rd = round(sv_ti/Δt) * Δt
# 					tsv_v = t[1] : sv_ti_rd : t[end]
# 					tsv = [tsv; [tsv_v] ]
# 					t_cnt = [tcnt; length( tsv_v )]
#
# 			#user prescribes array of values
# 			elseif length(sv_ti) > 1
#
# 					t_sv = Float64[]
# 					for ttmp in sv_ti
# 							ind_min = argmin( t .- ttmp )
# 							push!( t_sv, t[ind_min] )
# 					end
#
# 					t_sv = transpose(unique( t_sv ))
# 					tsv = [tsv; [t_sv] ]
# 					t_cnt = [tcnt; length( t_sv )]
#
# 			end
# 			#--
#
# 	  end
#
# 		#Create a vector of NamedTuples for each quantity the user wants to save
# 		in_tup = []
# 		for j in 1 : length(save_info.save_fcns)
# 				svfc = save_info.save_fcns[j]
# 				@show zeros(size(svfc(0.0,state)))
# 				@show tsv[j]
# 				@show [ [zeros(size(svfc(0.0,state)))]; [tsv[j]] ]
#
# 				push!( in_tup, NamedTuple{(:val, :t)}(
# 				[ [zeros(size(svfc(0.0,state)))]; [tsv[j]] ]
# 				) )
# 		end
#
# 		data = NamedTuple{ Tuple(save_info.save_names) }( in_tup )
#     # data = NamedTuple{ Tuple(save_info.save_names) }( [
#     #   zeros(size(svfc(0.0,state))) for svfc in save_info.save_names;
# 		#
# 		# 	])
#
#   end
#
# 	@show data
#
#   return data
#
# end

function save_data!(tnow, t, prob, state, data)

		for j = 1 : length(data)

			#Saving?
			tdiff = @. abs( data[j].t - tnow )
			ind_diff = @. (tdiff ≈ 0.0)
			if sum(ind_diff) > 0
				@show typeof(data[j].fcn( tnow, state ))
				@show typeof(data[j].val)
				push!( data[j].val, data[j].fcn( tnow, state ) )
			end

		  end

end #function
