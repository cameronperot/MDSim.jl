norm(A::Array; dims::Int) = sqrt.(sum(x->x^2, A, dims=dims))


function compute_r_nearest_image!(mdsim::MDSimulation)
	P = mdsim.params
	x = mdsim.state.x .% P.L

	mdsim.state.r_nearest_image = reshape(x', (1, P.d, P.N)) .- reshape(x, (P.N, P.d, 1))
	mdsim.state.r_nearest_image[mdsim.state.r_nearest_image .> P.L / 2]  .-= P.L
	mdsim.state.r_nearest_image[mdsim.state.r_nearest_image .< -P.L / 2] .+= P.L

	mdsim.state.r_nearest_image_magnitude = norm(mdsim.state.r_nearest_image, dims=2)
	mdsim.state.r_nearest_image_magnitude[P.diag_indices] .= Inf
end


function compute_row_vector_magnitudes(r)
	if size(r, 2) == 2
		return @. sqrt(r[:, 1]^2 + r[:, 2]^2)
	elseif size(r, 2) == 3
		return @. sqrt(r[:, 1]^2 + r[:, 2]^2 + r[:, 3]^2)
	end
end


function compute_energies!(mdsim::MDSimulation)
	P = mdsim.params
	α = (P.σ_LJ ./ mdsim.state.r_nearest_image_magnitude).^6
	mdsim.state.E_pot = sum((4 * P.ε_LJ) .* (α.^2 .- α)) / 2
	mdsim.state.E_kin = 0.5P.m * sum(mdsim.state.v.^2)
end


function compute_accelerations!(mdsim::MDSimulation)
	P = mdsim.params
	mdsim.state.a_prev = mdsim.state.a
	compute_r_nearest_image!(mdsim)

	α = (P.σ_LJ ./ mdsim.state.r_nearest_image_magnitude).^6
	mdsim.state.a = sum(
		((48 * P.ε_LJ / P.m) * (α.^2 - 0.5α) ./ mdsim.state.r_nearest_image_magnitude.^2)
		.* mdsim.state.r_nearest_image,
		dims=1
	)[1, :, :]'
end
