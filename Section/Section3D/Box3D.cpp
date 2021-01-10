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

#include "Box3D.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>

Box3D::Box3D(const unsigned T, const double B, const double H, const double TH, const unsigned M, const unsigned S, const double E1, const double E2)
	: Section3D(T, M, 2. * TH * (B + H - 4. * TH), vec{E1, E2})
	, width(B)
	, height(H)
	, thickness(TH)
	, int_pt_num(S > 20 ? 20 : S) {}

Box3D::Box3D(const unsigned T, vec&& D, const unsigned M, const unsigned S, vec&& EC)
	: Section3D(T, M, 2. * D(2) * (D(0) + D(1) - 4. * D(2)), std::forward<vec>(EC))
	, width(D(0))
	, height(D(1))
	, thickness(D(2))
	, int_pt_num(S > 20 ? 20 : S) {}

void Box3D::initialize(const shared_ptr<DomainBase>& D) {
	auto& material_proto = D->get_material(material_tag);

	linear_density = area * material_proto->get_parameter(ParameterType::DENSITY);

	const IntegrationPlan plan(1, int_pt_num, IntegrationType::GAUSS);

	const auto net_flange = width - 2. * thickness;
	const auto net_web = height - 2. * thickness;
	const auto web_area = net_web * thickness;
	const auto flange_area = net_flange * thickness;
	const auto web_middle = .5 * (width - thickness);
	const auto flange_middle = .5 * (height - thickness);

	int_pt.clear();
	int_pt.reserve(4 * static_cast<size_t>(int_pt_num));
	for(unsigned I = 0; I < int_pt_num; ++I) {
		int_pt.emplace_back(.5 * plan(I, 0) * net_web, web_middle, .5 * plan(I, 1) * web_area, material_proto->get_copy());
		int_pt.emplace_back(.5 * plan(I, 0) * net_web, -web_middle, .5 * plan(I, 1) * web_area, material_proto->get_copy());
		int_pt.emplace_back(flange_middle, .5 * plan(I, 0) * net_flange, .5 * plan(I, 1) * flange_area, material_proto->get_copy());
		int_pt.emplace_back(-flange_middle, .5 * plan(I, 0) * net_flange, .5 * plan(I, 1) * flange_area, material_proto->get_copy());
	}

	initial_stiffness.zeros(3, 3);
	for(const auto& I : int_pt) {
		const auto tmp_a = I.s_material->get_initial_stiffness().at(0) * I.weight;
		const auto arm_y = eccentricity(0) - I.coor_y;
		const auto arm_z = I.coor_z - eccentricity(1);
		initial_stiffness(0, 0) += tmp_a;
		initial_stiffness(0, 1) += tmp_a * arm_y;
		initial_stiffness(0, 2) += tmp_a * arm_z;
		initial_stiffness(1, 1) += tmp_a * arm_y * arm_y;
		initial_stiffness(1, 2) += tmp_a * arm_y * arm_z;
		initial_stiffness(2, 2) += tmp_a * arm_z * arm_z;
	}

	initial_stiffness(1, 0) = initial_stiffness(0, 1);
	initial_stiffness(2, 0) = initial_stiffness(0, 2);
	initial_stiffness(2, 1) = initial_stiffness(1, 2);

	trial_stiffness = current_stiffness = initial_stiffness;
}

unique_ptr<Section> Box3D::get_copy() { return make_unique<Box3D>(*this); }

void Box3D::print() { suanpan_info("A Box3D Section.\n"); }
