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

#include "Rotation2D.h"
#include <Domain/DomainBase.h>

void Rotation2D::form_transformation(mat&& R) {
	const auto &R1 = R(0, 0), &R2 = R(0, 1);
	const auto &R3 = R(1, 0), &R4 = R(1, 1);

	right.set_size(3, 3);
	left.set_size(3, 3);

	right(span(0, 1), span(0, 1)) = square(R);

	left(span(0, 1), span(0, 1)) = right(span(0, 1), span(0, 1)).t();

	right(0, 2) = 2. * R1 * R2;
	right(1, 2) = 2. * R3 * R4;

	right(2, 0) = R1 * R3;
	right(2, 1) = R2 * R4;
	right(2, 2) = R1 * R4 + R2 * R3;

	left(0, 2) = 2. * R1 * R3;
	left(1, 2) = 2. * R2 * R4;

	left(2, 0) = R1 * R2;
	left(2, 1) = R3 * R4;
	left(2, 2) = R1 * R4 + R2 * R3;
}

Rotation2D::Rotation2D(const unsigned T, const unsigned MT, const double A)
	: Material2D(T, PlaneType::N, 0.)
	, mat_tag(MT) {
	const auto C = cos(A), S = sin(A);

	form_transformation({{C, -S}, {S, C}});
}

Rotation2D::Rotation2D(const Rotation2D& old_obj)
	: Material2D(old_obj)
	, mat_tag(old_obj.mat_tag)
	, mat_obj(old_obj.mat_obj == nullptr ? nullptr : old_obj.mat_obj->get_copy())
	, left(old_obj.left)
	, right(old_obj.right) {}

void Rotation2D::initialize(const shared_ptr<DomainBase>& D) {
	if(nullptr == D || !D->find_material(mat_tag) || D->get_material(mat_tag)->get_material_type() != MaterialType::D2) {
		D->disable_material(get_tag());
		return;
	}

	mat_obj = D->get_material(mat_tag)->get_copy();

	mat_obj->Material::initialize(D);
	mat_obj->initialize(D);
	access::rw(density) = mat_obj->get_parameter(ParameterType::DENSITY);

	trial_stiffness = current_stiffness = initial_stiffness = left * mat_obj->get_initial_stiffness() * right;
}

double Rotation2D::get_parameter(const ParameterType P) const { return mat_obj->get_parameter(P); }

unique_ptr<Material> Rotation2D::get_copy() { return make_unique<Rotation2D>(*this); }

int Rotation2D::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

	if(mat_obj->update_trial_status(right * trial_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	trial_stress = left * mat_obj->get_trial_stress();
	trial_stiffness = left * mat_obj->get_trial_stiffness() * right;

	return SUANPAN_SUCCESS;
}

int Rotation2D::clear_status() {
	current_strain.zeros();
	trial_strain.zeros();
	current_stress.zeros();
	trial_stress.zeros();
	trial_stiffness = current_stiffness = initial_stiffness;
	return mat_obj->clear_status();
}

int Rotation2D::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_stiffness = trial_stiffness;
	return mat_obj->commit_status();
}

int Rotation2D::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_stiffness = current_stiffness;
	return mat_obj->reset_status();
}

vector<vec> Rotation2D::record(const OutputType P) { return mat_obj->record(P); }
