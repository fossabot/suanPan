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

#include "Rotation3D.h"
#include <Domain/DomainBase.h>

void Rotation3D::form_transformation(mat&& R) {
	const auto &R1 = R(0, 0), &R2 = R(0, 1), &R3 = R(0, 2);
	const auto &R4 = R(1, 0), &R5 = R(1, 1), &R6 = R(1, 2);
	const auto &R7 = R(2, 0), &R8 = R(2, 1), &R9 = R(2, 2);

	right.set_size(6, 6);
	left.set_size(6, 6);

	right(span(0, 2), span(0, 2)) = square(R);

	left(span(0, 2), span(0, 2)) = right(span(0, 2), span(0, 2)).t();

	right(0, 3) = 2. * R1 * R2;
	right(0, 4) = 2. * R2 * R3;
	right(0, 5) = 2. * R3 * R1;

	right(1, 3) = 2. * R4 * R5;
	right(1, 4) = 2. * R5 * R6;
	right(1, 5) = 2. * R6 * R4;

	right(2, 3) = 2. * R7 * R8;
	right(2, 4) = 2. * R8 * R9;
	right(2, 5) = 2. * R9 * R7;

	right(3, 0) = R1 * R4;
	right(3, 1) = R2 * R5;
	right(3, 2) = R3 * R6;
	right(3, 3) = R1 * R5 + R2 * R4;
	right(3, 4) = R2 * R6 + R3 * R5;
	right(3, 5) = R1 * R6 + R3 * R4;

	right(4, 0) = R4 * R7;
	right(4, 1) = R5 * R8;
	right(4, 2) = R6 * R9;
	right(4, 3) = R4 * R8 + R5 * R7;
	right(4, 4) = R5 * R9 + R6 * R8;
	right(4, 5) = R4 * R9 + R6 * R7;

	right(5, 0) = R1 * R7;
	right(5, 1) = R2 * R8;
	right(5, 2) = R3 * R9;
	right(5, 3) = R1 * R8 + R2 * R7;
	right(5, 4) = R2 * R9 + R3 * R8;
	right(5, 5) = R1 * R9 + R3 * R7;

	left(0, 3) = 2. * R1 * R4;
	left(0, 4) = 2. * R4 * R7;
	left(0, 5) = 2. * R7 * R1;

	left(1, 3) = 2. * R2 * R5;
	left(1, 4) = 2. * R5 * R8;
	left(1, 5) = 2. * R8 * R2;

	left(2, 3) = 2. * R3 * R6;
	left(2, 4) = 2. * R6 * R9;
	left(2, 5) = 2. * R9 * R3;

	left(3, 0) = R1 * R2;
	left(3, 1) = R4 * R5;
	left(3, 2) = R7 * R8;
	left(3, 3) = R1 * R5 + R2 * R4;
	left(3, 4) = R4 * R8 + R5 * R7;
	left(3, 5) = R1 * R8 + R2 * R7;

	left(4, 0) = R2 * R3;
	left(4, 1) = R5 * R6;
	left(4, 2) = R8 * R9;
	left(4, 3) = R2 * R6 + R3 * R5;
	left(4, 4) = R5 * R9 + R6 * R8;
	left(4, 5) = R2 * R9 + R3 * R8;

	left(5, 0) = R1 * R3;
	left(5, 1) = R4 * R6;
	left(5, 2) = R7 * R9;
	left(5, 3) = R1 * R6 + R3 * R4;
	left(5, 4) = R4 * R9 + R6 * R7;
	left(5, 5) = R1 * R9 + R3 * R7;
}

Rotation3D::Rotation3D(const unsigned T, const unsigned MT, const double A1, const double A2, const double A3)
	: Material3D(T, 0.)
	, mat_tag(MT) {
	const auto C1 = cos(A1), C2 = cos(A2), C3 = cos(A3);
	const auto S1 = sin(A1), S2 = sin(A2), S3 = sin(A3);

	const auto R1 = C1 * C2 * C3 - S1 * S3;
	const auto R2 = -C1 * C2 * S3 - S1 * C3;
	const auto R3 = C1 * S2;
	const auto R4 = S1 * C2 * C3 + C1 * S3;
	const auto R5 = -S1 * C2 * S3 + C1 * C3;
	const auto R6 = S1 * S2;
	const auto R7 = -S2 * C3;
	const auto R8 = S2 * S3;
	const auto R9 = C2;

	form_transformation({{R1, R2, R3}, {R4, R5, R6}, {R7, R8, R9}});
}

Rotation3D::Rotation3D(const unsigned T, const unsigned MT, mat&& R)
	: Material3D(T, 0.)
	, mat_tag(MT) { form_transformation(std::forward<mat>(R)); }

Rotation3D::Rotation3D(const Rotation3D& old_obj)
	: Material3D(old_obj)
	, mat_tag(old_obj.mat_tag)
	, mat_obj(nullptr == old_obj.mat_obj ? nullptr : old_obj.mat_obj->get_copy())
	, left(old_obj.left)
	, right(old_obj.right) {}

void Rotation3D::initialize(const shared_ptr<DomainBase>& D) {
	if(nullptr == D || !D->find_material(mat_tag) || D->get_material(mat_tag)->get_material_type() != MaterialType::D3) {
		D->disable_material(get_tag());
		return;
	}

	mat_obj = D->get_material(mat_tag)->get_copy();

	mat_obj->Material::initialize(D);
	mat_obj->initialize(D);
	access::rw(density) = mat_obj->get_parameter(ParameterType::DENSITY);

	trial_stiffness = current_stiffness = initial_stiffness = left * mat_obj->get_initial_stiffness() * right;
}

unique_ptr<Material> Rotation3D::get_copy() { return make_unique<Rotation3D>(*this); }

int Rotation3D::update_trial_status(const vec& t_strain) {
	incre_strain = (trial_strain = t_strain) - current_strain;

	if(norm(incre_strain) <= datum::eps) return SUANPAN_SUCCESS;

	if(mat_obj->update_trial_status(right * trial_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

	trial_stress = left * mat_obj->get_trial_stress();
	trial_stiffness = left * mat_obj->get_trial_stiffness() * right;

	return SUANPAN_SUCCESS;
}

int Rotation3D::clear_status() {
	current_strain.zeros();
	trial_strain.zeros();
	current_stress.zeros();
	trial_stress.zeros();
	trial_stiffness = current_stiffness = initial_stiffness;
	return mat_obj->clear_status();
}

int Rotation3D::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_stiffness = trial_stiffness;
	return mat_obj->commit_status();
}

int Rotation3D::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_stiffness = current_stiffness;
	return mat_obj->reset_status();
}

vector<vec> Rotation3D::record(const OutputType P) { return mat_obj->record(P); }
