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

#include "Damper02.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Viscosity/Maxwell.h>

uvec Damper02::IS{0, 1};
uvec Damper02::JS{2, 3};

Damper02::Damper02(const unsigned T, uvec&& NT, const unsigned DT, const unsigned ST, const bool UM, const unsigned PC, const double BT)
	: MaterialElement1D(T, d_node, d_dof, std::forward<uvec>(NT), {}, false)
	, device(make_unique<Maxwell>(0, DT, ST, UM, PC, BT)) {}

void Damper02::initialize(const shared_ptr<DomainBase>& D) {
	device->Material::initialize(D);
	device->initialize(D);

	const mat coord = get_coordinate(d_dof).t();

	access::rw(direction_cosine) = normalise(coord.col(1) - coord.col(0));

	const auto t_disp = get_current_displacement();
	const auto t_vec = get_current_velocity();

	access::rw(device->get_current_strain()) = vec{dot(direction_cosine, t_disp(JS) - t_disp(IS))};
	access::rw(device->get_current_strain_rate()) = vec{dot(direction_cosine, t_vec(JS) - t_vec(IS))};

	initial_damping.set_size(d_size, d_size);
	initial_damping(IS, IS) = direction_cosine * device->get_initial_damping() * direction_cosine.t();
	initial_damping(IS, JS) = -initial_damping(IS, IS);
	initial_damping(JS, JS) = initial_damping(IS, IS);
	initial_damping(JS, IS) = initial_damping(IS, JS);

	initial_stiffness.set_size(d_size, d_size);
	initial_stiffness(IS, IS) = direction_cosine * device->get_initial_stiffness() * direction_cosine.t();
	initial_stiffness(IS, JS) = -initial_stiffness(IS, IS);
	initial_stiffness(JS, JS) = initial_stiffness(IS, IS);
	initial_stiffness(JS, IS) = initial_stiffness(IS, JS);

	trial_damping = current_damping = initial_damping;
	trial_stiffness = current_stiffness = initial_stiffness;
}

int Damper02::update_status() {
	const auto t_disp = get_trial_displacement();
	const auto t_vec = get_trial_velocity();

	if(device->update_trial_status(dot(direction_cosine, t_disp(JS) - t_disp(IS)), dot(direction_cosine, t_vec(JS) - t_vec(IS))) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	trial_resistance.set_size(d_size);
	trial_resistance(JS) = direction_cosine * device->get_trial_stress();
	trial_resistance(IS) = -trial_resistance(JS);

	trial_damping.set_size(d_size, d_size);
	trial_damping(IS, IS) = direction_cosine * device->get_trial_damping() * direction_cosine.t();
	trial_damping(IS, JS) = -trial_damping(IS, IS);
	trial_damping(JS, JS) = trial_damping(IS, IS);
	trial_damping(JS, IS) = trial_damping(IS, JS);

	trial_stiffness.set_size(d_size, d_size);
	trial_stiffness(IS, IS) = direction_cosine * device->get_trial_stiffness() * direction_cosine.t();
	trial_stiffness(IS, JS) = -trial_stiffness(IS, IS);
	trial_stiffness(JS, JS) = trial_stiffness(IS, IS);
	trial_stiffness(JS, IS) = trial_stiffness(IS, JS);

	return SUANPAN_SUCCESS;
}

int Damper02::commit_status() { return device->commit_status(); }

int Damper02::clear_status() { return device->clear_status(); }

int Damper02::reset_status() { return device->reset_status(); }

vector<vec> Damper02::record(const OutputType P) { return device->record(P); }

void Damper02::print() { suanpan_info("A viscous damper element using displacement and velocity as basic quantities.\n"); }
