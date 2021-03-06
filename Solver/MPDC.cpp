#include "MPDC.h"
#include <Converger/Converger.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

MPDC::MPDC(const unsigned T)
	: Solver(T) {}

int MPDC::analyze() {
	auto& C = get_converger();
	auto& G = get_integrator();
	auto& W = G->get_domain().lock()->get_factory();

	suanpan_info("current analysis time: %.5f.\n", W->get_trial_time());

	const auto& max_iteration = C->get_max_iteration();

	// ninja anchor
	auto& ninja = get_ninja(W);

	// get column index for each nonzero dof
	// uvec load_ref_idx = find(load_ref);
	// for(auto I = 0; I < load_ref_idx.n_elem; ++I) load_ref_idx(I) -= I * load_ref.n_rows;

	auto& idx = W->get_reference_dof();

	mat disp_a;

	// iteration counter
	unsigned counter = 0;

	while(true) {
		// process modifiers
		if(G->process_modifier() != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		// assemble resistance
		G->assemble_resistance();
		// assemble stiffness
		G->assemble_matrix();
		// process loads
		if(G->process_load() != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		// process constraints
		if(G->process_constraint() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		// solve ninja
		if(G->solve(ninja, G->get_displacement_residual()) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		// solve reference displacement
		if(G->solve_trs(disp_a, W->get_reference_load()) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		if(const auto n_size = W->get_size(); 0 != W->get_mpc()) {
			mat right, kernel;
			auto& border = W->get_auxiliary_stiffness();
			if(G->solve_trs(right, border) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
			ninja -= right * solve(kernel = border.t() * right.head_rows(n_size), border.t() * ninja.head_rows(n_size) - W->get_auxiliary_load());
			disp_a -= right * solve(kernel, border.t() * disp_a.head_rows(n_size));
		}

		const vec incre_lambda = solve(mat(disp_a.rows(idx)), W->get_trial_settlement()(idx) - W->get_trial_displacement()(idx) - ninja.rows(idx));

		ninja += disp_a * incre_lambda;

		// avoid machine error accumulation
		G->erase_machine_error();
		// update internal variable
		G->update_internal(ninja);
		// update trial load factor
		W->update_trial_load_factor(W->get_trial_load_factor() + incre_lambda);
		// update trial displacement
		W->update_trial_displacement(W->get_trial_displacement() + ninja);
		// update for nodes and elements
		if(G->update_trial_status() != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		// exit if converged
		if(C->is_converged()) return SUANPAN_SUCCESS;
		// exit if maximum iteration is hit
		if(++counter > max_iteration) return SUANPAN_FAIL;
	}
}
