function initialize_random_positions(P::Parameters)
	if P.d == 2
		x_range = range(P.σ_LJ / 2, stop=P.L - P.σ_LJ / 2, length=ceil(Int, sqrt(P.N)))
		x = [
			[i for i in x_range, j in x_range],
			[j for i in x_range, j in x_range],
			]
		x = [x[j][i] for i in 1:length(x[1]), j in 1:length(x)][1:P.N, :]
	elseif P.d == 3
		x_range = range(P.σ_LJ / 2, stop=P.L - P.σ_LJ / 2, length=ceil(Int, cbrt(P.N)))
		x = [
			[i for i in x_range, j in x_range, k in x_range],
			[j for i in x_range, j in x_range, k in x_range],
			[k for i in x_range, j in x_range, k in x_range],
			]
		x = [x[j][i] for i in 1:length(x[1]), j in 1:length(x)][1:P.N, :]
	else
		error("d must be ∈ {2, 3}")
	end

	x += rand(P.rng, Normal(0, 0.01 * P.σ_LJ), size(x))
	return x
end


function initialize_random_velocities(P::Parameters)
	μ = 0
	σ = sqrt(P.k_B * P.T / P.m)
	v = rand(P.rng, Normal(μ, σ), (P.N, P.d))
	λ = sqrt(P.d * P.k_B * (P.N - 1) * P.T / (P.m * sum(v.^2)))
	v *= λ
	v .- mean(v)
	return v
end
