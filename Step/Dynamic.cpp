////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "Dynamic.h"
#include <Converger/AbsIncreDisp.h>
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Load/GroupNodalDisplacement.h>
#include <Solver/BFGS.h>
#include <Solver/MPDC.h>
#include <Solver/Newton.h>
#include <Solver/Ramm.h>
#include <Solver/Integrator/LeeNewmark.h>
#include <Solver/Integrator/LeeNewmarkFull.h>

Dynamic::Dynamic(const unsigned T, const double P)
	: Step(T, P) {}

int Dynamic::initialize() {
	if(sparse_mat) factory->set_storage_scheme(StorageScheme::SPARSE);
	else if(symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BANDSYMM);
	else if(!symm_mat && band_mat) factory->set_storage_scheme(StorageScheme::BAND);
	else if(symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::SYMMPACK);
	else if(!symm_mat && !band_mat) factory->set_storage_scheme(StorageScheme::FULL);

	factory->set_analysis_type(AnalysisType::DYNAMICS);
	factory->set_precision(precision);
	factory->set_tolerance(tolerance);
	factory->set_solver(system_solver);

	const auto& t_domain = database.lock();

	if(SUANPAN_SUCCESS != t_domain->restart()) return SUANPAN_FAIL;

	if(nullptr == tester) tester = make_shared<AbsIncreDisp>();
	if(nullptr == modifier) modifier = make_shared<Newmark>();

	// avoid arc length solver
	if(nullptr != solver) {
		const auto& t_solver = *solver;
		if(typeid(t_solver) == typeid(Ramm)) solver = nullptr;
	}

	// automatically enable displacement controlled solver
	if(nullptr == solver) {
		auto flag = false;
		for(const auto& I : t_domain->get_load_pool())
			if(const auto& J = *I; (typeid(J) == typeid(NodalDisplacement) || typeid(J) == typeid(GroupNodalDisplacement)) && I->get_start_step() == get_tag()) {
				flag = true;
				break;
			}
		flag ? solver = make_shared<MPDC>() : solver = make_shared<Newton>();
	}

	const auto& t_solver = *solver;
	const auto& t_modifier = *modifier;
	if(typeid(t_solver) == typeid(BFGS) && typeid(t_modifier) == typeid(LeeNewmark) && typeid(t_modifier) == typeid(LeeNewmarkFull)) {
		suanpan_error("currently BFGS solver is not supported by Lee damping model.\n");
		return SUANPAN_FAIL;
	}

	tester->set_domain(t_domain);
	modifier->set_domain(t_domain);
	solver->set_converger(tester);
	solver->set_integrator(modifier);

	if(SUANPAN_SUCCESS != modifier->initialize()) return SUANPAN_FAIL;

	modifier->update_parameter(get_ini_step_size());
	modifier->update_compatibility();

	return SUANPAN_SUCCESS;
}

int Dynamic::analyze() {
	auto& S = get_solver();
	auto& G = get_integrator();

	auto remain_time = get_time_period();
	auto step_time = get_ini_step_size();

	// record initial state
	// if(W->get_current_time() == 0.) G->record();

	unsigned num_increment = 0, num_converged_step = 0;

	while(true) {
		// check if the target time point is hit
		if(remain_time <= 1E-7) return SUANPAN_SUCCESS;
		// check if the maximum substep number is hit
		if(++num_increment > get_max_substep()) {
			suanpan_warning("analyze() reaches maximum substep number %u.\n", get_max_substep());
			return SUANPAN_FAIL;
		}
		// update incremental and trial time
		G->update_incre_time(step_time);
		// call solver
		const auto code = S->analyze();
		if(code == SUANPAN_SUCCESS) {
			// success step
			// commit converged iteration
			G->commit_status();
			// record response
			G->record();
			// eat current increment
			set_time_left(remain_time -= step_time);
			if(!is_fixed_step_size() && ++num_converged_step > 5 && G->allow_to_change_time_step()) {
				step_time = std::min(get_max_step_size(), step_time * time_step_amplification);
				num_converged_step = 0;
			}
			// check if time overflows
			if(step_time > remain_time) step_time = remain_time;
		} else if(code == SUANPAN_FAIL) {
			// failed step
			// reset to the start of current substep
			G->reset_status();
			// check if minimum step size is hit
			if(step_time <= get_min_step_size()) {
				suanpan_error("analyze() reaches minimum step size %.3E.\n", get_min_step_size());
				return SUANPAN_FAIL;
			}
			// check if fixed step size
			if(is_fixed_step_size() || !G->allow_to_change_time_step()) {
				suanpan_error("analyze() does not converge for given fixed step size %.3E.\n", step_time);
				return SUANPAN_FAIL;
			}
			// step size is allowed to decrease
			step_time *= .5;
			set_max_substep(num_increment + static_cast<unsigned>(remain_time / step_time) + 1);
			num_converged_step = 0;
		} else return SUANPAN_FAIL; // positive codes are from lapack subroutines
	}
}
