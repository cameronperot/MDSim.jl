function rescale_velocity!(mdsim::MDSimulation)
	P = mdsim.params
	λ = sqrt(P.d * P.k_B * (P.N - 1) * P.T / (P.m * sum(mdsim.state.v.^2)))
	mdsim.state.v *= λ
end


function update_history!(mdsim::MDSimulation, n::Int)
	mdsim.history.x[:, :, n] = mdsim.state.x
	mdsim.history.v[:, :, n] = mdsim.state.v
	mdsim.history.E_kin[n]   = mdsim.state.E_kin
	mdsim.history.E_pot[n]   = mdsim.state.E_pot
end


function velocity_verlet!(mdsim::MDSimulation)
	P = mdsim.params
	for n in 1:P.n_iter
		@. mdsim.state.x += P.τ * mdsim.state.v + (P.τ^2 / 2) * mdsim.state.a
		compute_accelerations!(mdsim)
		@. mdsim.state.v += (P.τ / 2) * (mdsim.state.a_prev + mdsim.state.a)

		if n % P.rescale_velocity_interval == 0 && P.ensemble_type == :NVT
			rescale_velocity!(mdsim)
		end

		compute_energies!(mdsim)
		update_history!(mdsim, n + 1)
	end
end
