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

#include "PolyElastic1D.h"

PolyElastic1D::PolyElastic1D(const unsigned T, vec&& P, const double R)
	: DataPolyElastic1D{std::forward<vec>(P)}
	, Material1D(T, R) {}

void PolyElastic1D::initialize(const shared_ptr<DomainBase>&) { trial_stiffness = current_stiffness = initial_stiffness = pool(0); }

unique_ptr<Material> PolyElastic1D::get_copy() { return make_unique<PolyElastic1D>(*this); }

int PolyElastic1D::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(fabs(incre_strain(0)) <= datum::eps) return SUANPAN_SUCCESS;

	auto pow_strain = 1., exponent = 1.;
	trial_stress = 0.;
	trial_stiffness = 0.;
	for(auto I : pool) {
		trial_stiffness += exponent * I * pow_strain;
		exponent += 2.;
		pow_strain *= trial_strain(0);
		trial_stress += I * pow_strain;
		pow_strain *= trial_strain(0);
	}

	return SUANPAN_SUCCESS;
}

int PolyElastic1D::clear_status() {
	current_strain.zeros();
	current_stress.zeros();
	current_stiffness = initial_stiffness;
	return reset_status();
}

int PolyElastic1D::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_stiffness = trial_stiffness;
	return SUANPAN_SUCCESS;
}

int PolyElastic1D::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_stiffness = current_stiffness;
	return SUANPAN_SUCCESS;
}

void PolyElastic1D::print() {
	suanpan_info("An elastic model based on polynomial.\n");
	suanpan_info("current Strain: %.3E\tcurrent Stress: %.3E\n", current_strain(0), current_stress(0));
}
