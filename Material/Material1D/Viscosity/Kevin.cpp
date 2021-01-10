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

#include "Kevin.h"
#include <Domain/DomainBase.h>
#include <Recorder/OutputType.h>

Kevin::Kevin(const unsigned T, const unsigned DT, const unsigned ST)
	: Material1D(T, 0.)
	, damper_tag(DT)
	, spring_tag(ST) {}

Kevin::Kevin(const Kevin& old_obj)
	: Material1D(old_obj)
	, damper_tag(old_obj.damper_tag)
	, spring_tag(old_obj.spring_tag)
	, damper(suanpan::make_copy(old_obj.damper))
	, spring(suanpan::make_copy(old_obj.spring)) {}

void Kevin::initialize(const shared_ptr<DomainBase>& D) {
	if(nullptr == D || !D->find_material(damper_tag) || !D->find_material(spring_tag)) {
		D->disable_material(get_tag());
		return;
	}

	damper = get_material(D, damper_tag)->get_copy();
	spring = get_material(D, spring_tag)->get_copy();

	damper->Material::initialize(D);
	damper->initialize(D);
	spring->Material::initialize(D);
	spring->initialize(D);

	trial_strain_rate = current_strain_rate = incre_strain_rate.zeros(1);

	trial_damping = current_damping = initial_damping = damper->get_initial_damping();
	trial_stiffness = current_stiffness = initial_stiffness = spring->get_initial_stiffness();
}

unique_ptr<Material> Kevin::get_copy() { return make_unique<Kevin>(*this); }

int Kevin::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
	incre_strain = (trial_strain = t_strain) - current_strain;
	incre_strain_rate = (trial_strain_rate = t_strain_rate) - current_strain_rate;

	if(fabs(incre_strain(0) + fabs(incre_strain_rate(0))) <= datum::eps) return SUANPAN_SUCCESS;

	spring->update_trial_status(trial_strain, trial_strain_rate);
	damper->update_trial_status(trial_strain, trial_strain_rate);

	trial_stiffness = spring->get_trial_stiffness();
	trial_damping = damper->get_trial_damping();

	return SUANPAN_SUCCESS;
}

int Kevin::clear_status() {
	trial_strain = current_strain.zeros();
	trial_stress = current_stress.zeros();
	trial_strain_rate = current_strain_rate.zeros();
	trial_damping = current_damping = initial_damping;
	trial_stiffness = current_stiffness = initial_stiffness;
	return spring->clear_status() + damper->clear_status();
}

int Kevin::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_strain_rate = trial_strain_rate;
	current_damping = trial_damping;
	current_stiffness = trial_stiffness;
	return spring->commit_status() + damper->commit_status();
}

int Kevin::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_strain_rate = current_strain_rate;
	trial_damping = current_damping;
	trial_stiffness = current_stiffness;
	return spring->reset_status() + damper->reset_status();
}

vector<vec> Kevin::record(const OutputType P) {
	vector<vec> data;

	if(OutputType::SD == P || OutputType::SS == P || OutputType::S == P) data.emplace_back(current_stress);
	else if(OutputType::ED == P) data.emplace_back(damper->get_current_strain());
	else if(OutputType::VD == P) data.emplace_back(damper->get_current_strain_rate());
	else if(OutputType::ES == P) data.emplace_back(spring->get_current_strain());
	else if(OutputType::VS == P) data.emplace_back(spring->get_current_strain_rate());
	else if(OutputType::E == P) data.emplace_back(current_strain);
	else if(OutputType::V == P) data.emplace_back(current_strain_rate);

	return data;
}

void Kevin::print() { suanpan_info("A Kevin material model.\n"); }
