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

#include "Section3D.h"
#include <Material/Material.h>

Section3D::IntegrationPoint::IntegrationPoint(const double CY, const double CZ, const double W, unique_ptr<Material>&& M)
	: coor_y(CY)
	, coor_z(CZ)
	, weight(W)
	, s_material(std::forward<unique_ptr<Material>>(M)) {}

Section3D::IntegrationPoint::IntegrationPoint(const IntegrationPoint& old_obj)
	: coor_y(old_obj.coor_y)
	, coor_z(old_obj.coor_z)
	, weight(old_obj.weight)
	, s_material(old_obj.s_material->get_copy()) {}

Section3D::Section3D(const unsigned T, const unsigned MT, const double A, vec&& E)
	: Section(T, SectionType::D3, MT, A, std::forward<vec>(E)) {}

int Section3D::update_trial_status(const vec& t_deformation) {
	const vec incre_deformation = (trial_deformation = t_deformation) - current_deformation;

	if(norm(incre_deformation) <= datum::eps) return SUANPAN_SUCCESS;

	trial_stiffness.zeros();
	trial_resistance.zeros();

	for(const auto& I : int_pt) {
		const auto arm_y = eccentricity(0) - I.coor_y;
		const auto arm_z = I.coor_z - eccentricity(1);
		if(I.s_material->update_trial_status(trial_deformation(0) + trial_deformation(1) * arm_y + trial_deformation(2) * arm_z) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		const auto tmp_a = I.s_material->get_trial_stiffness().at(0) * I.weight;
		trial_stiffness(0, 0) += tmp_a;
		trial_stiffness(0, 1) += tmp_a * arm_y;
		trial_stiffness(0, 2) += tmp_a * arm_z;
		trial_stiffness(1, 1) += tmp_a * arm_y * arm_y;
		trial_stiffness(1, 2) += tmp_a * arm_y * arm_z;
		trial_stiffness(2, 2) += tmp_a * arm_z * arm_z;
		const auto tmp_b = I.s_material->get_trial_stress().at(0) * I.weight;
		trial_resistance(0) += tmp_b;
		trial_resistance(1) += tmp_b * arm_y;
		trial_resistance(2) += tmp_b * arm_z;
	}
	trial_stiffness(1, 0) = trial_stiffness(0, 1);
	trial_stiffness(2, 0) = trial_stiffness(0, 2);
	trial_stiffness(2, 1) = trial_stiffness(1, 2);

	return SUANPAN_SUCCESS;
}

int Section3D::clear_status() {
	current_deformation = trial_deformation.zeros();
	current_resistance = trial_resistance.zeros();
	current_stiffness = trial_stiffness = initial_stiffness;
	auto code = 0;
	for(const auto& I : int_pt) code += I.s_material->clear_status();
	return code;
}

int Section3D::commit_status() {
	current_deformation = trial_deformation;
	current_resistance = trial_resistance;
	current_stiffness = trial_stiffness;
	auto code = 0;
	for(const auto& I : int_pt) code += I.s_material->commit_status();
	return code;
}

int Section3D::reset_status() {
	trial_deformation = current_deformation;
	trial_resistance = current_resistance;
	trial_stiffness = current_stiffness;
	auto code = 0;
	for(const auto& I : int_pt) code += I.s_material->reset_status();
	return code;
}
