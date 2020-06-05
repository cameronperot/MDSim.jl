struct Parameters
	N                        ::Int
	d                        ::Int
	m                        ::Float64
	L                        ::Float64
	T                        ::Float64
	τ                        ::Float64
	ε_LJ                     ::Float64
	σ_LJ                     ::Float64
	k_B                      ::Float64
	seed                     ::Int
	n_iter                   ::Int
	ensemble_type            ::Symbol
	rescale_velocity_interval::Int
	rng                      ::AbstractRNG
	diag_indices             ::Vector{CartesianIndex{3}}

	function Parameters(params_dict)
		new(
			params_dict[:N],
			params_dict[:d],
			params_dict[:m],
			params_dict[:L],
			params_dict[:T],
			params_dict[:τ],
			params_dict[:ε_LJ],
			params_dict[:σ_LJ],
			params_dict[:k_B],
			params_dict[:seed],
			params_dict[:n_iter],
			params_dict[:ensemble_type],
			params_dict[:rescale_velocity_interval],
			MersenneTwister(params_dict[:seed]),
			[CartesianIndex((i, 1, i)) for i in 1:params_dict[:N]]
		)
	end
end


mutable struct State
	x                        ::Matrix{Float64}
	v                        ::Matrix{Float64}
	a                        ::Matrix{Float64}
	a_prev                   ::Matrix{Float64}
	r_nearest_image          ::Array{Float64, 3}
	r_nearest_image_magnitude::Array{Float64, 3}
	E_kin                    ::Float64
	E_pot                    ::Float64

	function State(P::Parameters, params_dict::Dict)
		x = get(params_dict, :x, initialize_random_positions(P))
		v = get(params_dict, :v, initialize_random_velocities(P))
		a                         = zeros((P.N, P.d))
		a_prev                    = zeros((P.N, P.d))
		r_nearest_image           = zeros((P.N, P.d, P.N))
		r_nearest_image_magnitude = zeros((P.N, 1, P.N))
		E_kin                     = 0
		E_pot                     = 0

		new(
			x,
			v,
			a,
			a_prev,
			r_nearest_image,
			r_nearest_image_magnitude,
			E_kin,
			E_pot,
		)
	end
end


mutable struct History
	x    ::Array{Float64, 3}
	v    ::Array{Float64, 3}
	E_kin::Vector{Float64}
	E_pot::Vector{Float64}

	function History(P::Parameters)
		new(
			zeros((P.N, P.d, P.n_iter + 1)),
			zeros((P.N, P.d, P.n_iter + 1)),
			zeros(P.n_iter + 1),
			zeros(P.n_iter + 1),
		)
	end
end


mutable struct MDSimulation
	params  ::Parameters
	state   ::State
	history ::History

	function MDSimulation(params_dict)
		params  = Parameters(params_dict)
		state   = State(params, params_dict)
		history = History(params)

		mdsim = new(params, state, history)
		compute_accelerations!(mdsim)
		compute_energies!(mdsim)
		update_history!(mdsim, 1)

		return mdsim
	end
end
