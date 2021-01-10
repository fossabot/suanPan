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

#include "ISection2D.h"
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Domain/DomainBase.h>

ISection2D::ISection2D(const unsigned T, const double TFW, const double TFT, const double BFW, const double BFT, const double WH, const double WT, const unsigned MT, const unsigned IP, const double EC)
	: Section2D(T, MT, TFW * TFT + BFW * BFT + WH * WT, EC)
	, top_flange_width(TFW)
	, top_flange_thickness(TFT)
	, bottom_flange_width(BFW)
	, bottom_flange_thickness(BFT)
	, web_height(WH)
	, web_thickness(WT)
	, int_pt_num(IP > 20 ? 20 : IP) {}

ISection2D::ISection2D(const unsigned T, vec&& D, const unsigned MT, const unsigned IP, const double EC)
	: Section2D(T, MT, D(0) * D(1) + D(2) * D(3) + D(4) * D(5), EC)
	, top_flange_width(D(0))
	, top_flange_thickness(D(1))
	, bottom_flange_width(D(2))
	, bottom_flange_thickness(D(3))
	, web_height(D(4))
	, web_thickness(D(5))
	, int_pt_num(IP > 20 ? 20 : IP) {}

void ISection2D::initialize(const shared_ptr<DomainBase>& D) {
	const auto& mat_proto = D->get_material(material_tag)->get_copy();

	linear_density = mat_proto->get_parameter(ParameterType::DENSITY) * area;

	const auto height = top_flange_thickness + bottom_flange_thickness + web_height;
	const auto web_area = web_height * web_thickness;
	const auto b_flange_area = bottom_flange_width * bottom_flange_thickness;
	const auto t_flange_area = top_flange_width * top_flange_thickness;
	auto total_q = .5 * bottom_flange_thickness * b_flange_area;
	total_q += (.5 * web_height + bottom_flange_thickness) * web_area;
	total_q += (.5 * top_flange_thickness + bottom_flange_thickness + web_height) * t_flange_area;

	access::rw(eccentricity) += total_q / area;

	const IntegrationPlan plan_flange(1, 2, IntegrationType::GAUSS);
	const IntegrationPlan plan_web(1, int_pt_num, IntegrationType::GAUSS);

	int_pt.clear();
	int_pt.reserve(int_pt_num + 2 * plan_flange.n_rows);
	for(unsigned I = 0; I < int_pt_num; ++I) int_pt.emplace_back(bottom_flange_thickness + web_height * (.5 + .5 * plan_web(I, 0)), .5 * plan_web(I, 1) * web_area, mat_proto->get_copy());
	if(b_flange_area != 0.) for(unsigned I = 0; I < plan_flange.n_rows; ++I) int_pt.emplace_back((.5 + .5 * plan_flange(I, 0)) * bottom_flange_thickness, .5 * plan_flange(I, 1) * b_flange_area, mat_proto->get_copy());
	if(t_flange_area != 0.) for(unsigned I = 0; I < plan_flange.n_rows; ++I) int_pt.emplace_back(height - (.5 + .5 * plan_flange(I, 0)) * top_flange_thickness, .5 * plan_flange(I, 1) * t_flange_area, mat_proto->get_copy());

	initial_stiffness.zeros(2, 2);
	for(const auto& I : int_pt) {
		auto tmp_a = I.s_material->get_initial_stiffness().at(0) * I.weight;
		const auto tmp_b = eccentricity(0) - I.coor;
		initial_stiffness(0, 0) += tmp_a;
		initial_stiffness(0, 1) += tmp_a *= tmp_b;
		initial_stiffness(1, 1) += tmp_a *= tmp_b;
	}
	initial_stiffness(1, 0) = initial_stiffness(0, 1);

	trial_stiffness = current_stiffness = initial_stiffness;
}

unique_ptr<Section> ISection2D::get_copy() { return make_unique<ISection2D>(*this); }

void ISection2D::print() {
	suanpan_info("A I-shape section with following inetgeration points.\n");
	auto J = 1;
	for(const auto& I : int_pt) {
		suanpan_info("IP %u: %.4E.\n", J++, I.coor);
		I.s_material->print();
	}
}
