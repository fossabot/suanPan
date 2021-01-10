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

#include "MultiplierBC.h"
#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Domain/Node.h>
#include <Step/Frequency.h>

/**
 * \brief method to apply the BC to the system.
 * \param D `Domain`
 * \return 0
 */
int MultiplierBC::process(const shared_ptr<DomainBase>& D) {
	// for eignevalue problem, only use penalty method
	auto& t_step = *D->get_current_step();
	if(typeid(t_step) == typeid(Frequency)) return PenaltyBC::process(D);

	auto& t_matrix = D->get_factory()->get_stiffness();
	auto& t_load = get_trial_load(D->get_factory());
	auto& t_sushi = get_sushi(D->get_factory());

	for(const auto& I : nodes)
		if(D->find<Node>(I))
			if(auto& t_node = D->get<Node>(I); t_node->is_active()) {
				auto& t_dof = t_node->get_reordered_dof();
				for(const auto& J : dofs)
					if(J <= t_dof.n_elem)
						if(const auto& t_idx = t_dof(J - 1); D->insert_restrained_dof(t_idx)) {
							t_matrix->unify(t_idx);
							t_sushi(t_idx) = t_load(t_idx) = 0.;
						}
			}

	return SUANPAN_SUCCESS;
}
