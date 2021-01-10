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

#include "Damper01.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>

uvec Damper01::IS{0, 1};
uvec Damper01::JS{2, 3};

Damper01::Damper01(const unsigned T, uvec&& NT, const unsigned D)
	: MaterialElement1D(T, d_node, d_dof, std::forward<uvec>(NT), uvec{D}, false) {}

void Damper01::initialize(const shared_ptr<DomainBase>& D) {
	damper = D->get<Material>(material_tag(0))->get_copy();

	const mat coord = get_coordinate(d_dof).t();

	access::rw(direction_cosine) = normalise(coord.col(1) - coord.col(0));

	const auto t_disp = get_current_displacement();
	const auto t_vec = get_current_velocity();

	access::rw(damper->get_current_strain()) = vec{dot(direction_cosine, t_disp(JS) - t_disp(IS))};
	access::rw(damper->get_current_strain_rate()) = vec{dot(direction_cosine, t_vec(JS) - t_vec(IS))};

	initial_damping.set_size(d_size, d_size);
	initial_damping(IS, IS) = direction_cosine * damper->get_initial_damping() * direction_cosine.t();
	initial_damping(IS, JS) = -initial_damping(IS, IS);
	initial_damping(JS, JS) = initial_damping(IS, IS);
	initial_damping(JS, IS) = initial_damping(IS, JS);

	trial_damping = current_damping = initial_damping;
}

int Damper01::update_status() {
	const auto t_disp = get_trial_displacement();
	const auto t_vec = get_trial_velocity();

	if(damper->update_trial_status(dot(direction_cosine, t_disp(JS) - t_disp(IS)), dot(direction_cosine, t_vec(JS) - t_vec(IS))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	trial_resistance.set_size(d_size);
	trial_resistance(JS) = direction_cosine * damper->get_trial_stress();
	trial_resistance(IS) = -trial_resistance(JS);

	trial_damping(IS, IS) = direction_cosine * damper->get_trial_damping() * direction_cosine.t();
	trial_damping(IS, JS) = -trial_damping(IS, IS);
	trial_damping(JS, JS) = trial_damping(IS, IS);
	trial_damping(JS, IS) = trial_damping(IS, JS);

	return SUANPAN_SUCCESS;
}

int Damper01::commit_status() { return damper->commit_status(); }

int Damper01::clear_status() { return damper->clear_status(); }

int Damper01::reset_status() { return damper->reset_status(); }

vector<vec> Damper01::record(const OutputType P) { return damper->record(P); }

void Damper01::print() { suanpan_info("A viscous damper element using displacement and velocity as basic quantities.\n"); }
