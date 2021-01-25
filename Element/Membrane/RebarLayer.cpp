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

#include "RebarLayer.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>

RebarLayer::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
	: coor(std::forward<vec>(C))
	, weight(W)
	, m_material(std::forward<unique_ptr<Material>>(M))
	, strain_mat(m_size, fill::zeros) {}

RebarLayer::RebarLayer(const unsigned T, uvec&& N, vec&& C, const double TH, const char DR)
	: MaterialElement1D(T, m_node, m_dof, std::forward<uvec>(N), conv_to<uvec>::from(C(span(6, 8))), false)
	, thickness(TH)
	, direction(suanpan::to_upper(DR))
	, length_a(C(0))
	, length_b(C(1))
	, length_c(C(2))
	, rho_a(C(3))
	, rho_b(C(4))
	, rho_c(C(5)) {}

void RebarLayer::initialize(const shared_ptr<DomainBase>& D) {
	const auto ele_coor = get_coordinate(2);

	const auto LX1 = ele_coor(1, 1) - ele_coor(0, 1);
	const auto LX2 = ele_coor(2, 1) - ele_coor(1, 1);
	const auto LX3 = ele_coor(3, 1) - ele_coor(2, 1);
	const auto LX4 = ele_coor(0, 1) - ele_coor(3, 1);
	const auto LY1 = ele_coor(0, 0) - ele_coor(1, 0);
	const auto LY2 = ele_coor(1, 0) - ele_coor(2, 0);
	const auto LY3 = ele_coor(2, 0) - ele_coor(3, 0);
	const auto LY4 = ele_coor(3, 0) - ele_coor(0, 0);

	const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

	int_pt.clear();
	int_pt.reserve(3llu * plan.n_rows);

	if(rho_a > 0. && length_a > 0.) {
		auto& material_a = D->get<Material>(material_tag(0));
		for(unsigned I = 0; I < plan.n_rows; ++I) int_pt.emplace_back(vec{length_a - 1. + length_a * plan(I, 0), plan(I, 1)}, plan(I, 2) * rho_a * length_a, material_a->get_copy());
	}
	if(rho_b > 0. && length_b > 0.) {
		auto& material_b = D->get<Material>(material_tag(1));
		for(unsigned I = 0; I < plan.n_rows; ++I) int_pt.emplace_back(vec{length_a - length_c + length_b * plan(I, 0), plan(I, 1)}, plan(I, 2) * rho_b * length_b, material_b->get_copy());
	}
	if(rho_c > 0. && length_c > 0.) {
		auto& material_c = D->get<Material>(material_tag(2));
		for(unsigned I = 0; I < plan.n_rows; ++I) int_pt.emplace_back(vec{1. - length_c + length_c * plan(I, 0), plan(I, 1)}, plan(I, 2) * rho_c * length_c, material_c->get_copy());
	}

	if(direction == 'X') for(auto&& I : int_pt) I.coor = flipud(I.coor);

	mat pnt(2, 8);

	for(auto&& I : int_pt) {
		const auto pn = shape::quad(I.coor, 1);
		const mat jacob = pn * ele_coor;
		const mat pn_pxy = solve(jacob, pn);
		I.weight *= det(jacob);

		const auto &X = I.coor(0), &Y = I.coor(1);

		const auto X2 = 2. * X, Y2 = 2. * Y;

		const auto AA = X + 1., BB = X - 1., CC = Y + 1., DD = Y - 1.;

		pnt(0, 0) = DD * (LX4 * CC - LX1 * X2);
		pnt(0, 1) = DD * (LX2 * CC + LX1 * X2);
		pnt(0, 2) = CC * (LX3 * X2 - LX2 * DD);
		pnt(0, 3) = -CC * (LX3 * X2 + LX4 * DD);
		pnt(0, 4) = DD * (LY4 * CC - LY1 * X2);
		pnt(0, 5) = DD * (LY2 * CC + LY1 * X2);
		pnt(0, 6) = CC * (LY3 * X2 - LY2 * DD);
		pnt(0, 7) = -CC * (LY3 * X2 + LY4 * DD);
		pnt(1, 0) = BB * (LX4 * Y2 - LX1 * AA);
		pnt(1, 1) = AA * (LX1 * BB + LX2 * Y2);
		pnt(1, 2) = AA * (LX3 * BB - LX2 * Y2);
		pnt(1, 3) = -BB * (LX3 * AA + LX4 * Y2);
		pnt(1, 4) = BB * (LY4 * Y2 - LY1 * AA);
		pnt(1, 5) = AA * (LY1 * BB + LY2 * Y2);
		pnt(1, 6) = AA * (LY3 * BB - LY2 * Y2);
		pnt(1, 7) = -BB * (LY3 * AA + LY4 * Y2);

		const mat pnt_pxy = solve(jacob, pnt / 16.);

		if(direction == 'Y')
			for(unsigned J = 0, K = 1, L = 2, M = 4; J < m_node; ++J, K += m_dof, L += m_dof, ++M) {
				I.strain_mat(K) = pn_pxy(1, J);
				I.strain_mat(L) = pnt_pxy(1, M);
			}
		else
			for(unsigned J = 0, K = 0, L = 2; J < m_node; ++J, K += m_dof, L += m_dof) {
				I.strain_mat(K) = pn_pxy(0, J);
				I.strain_mat(L) = pnt_pxy(0, J);
			}
	}

	initial_stiffness.zeros(m_size, m_size);
	for(const auto& I : int_pt) initial_stiffness += I.strain_mat * I.m_material->get_trial_stiffness() * I.strain_mat.t() * I.weight * thickness;
	trial_stiffness = current_stiffness = initial_stiffness;
}

int RebarLayer::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness.zeros(m_size, m_size);
	trial_resistance.zeros(m_size);
	for(const auto& I : int_pt) {
		if(SUANPAN_SUCCESS != I.m_material->update_trial_status(dot(I.strain_mat, t_disp))) return SUANPAN_FAIL;
		trial_stiffness += I.weight * thickness * I.strain_mat * I.m_material->get_trial_stiffness() * I.strain_mat.t();
		trial_resistance += I.weight * thickness * I.strain_mat * I.m_material->get_trial_stress();
	}

	return SUANPAN_SUCCESS;
}

int RebarLayer::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	return code;
}

int RebarLayer::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	return code;
}

int RebarLayer::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	return code;
}

vector<vec> RebarLayer::record(const OutputType T) {
	vector<vec> output;
	for(const auto& I : int_pt) for(const auto& J : I.m_material->record(T)) output.emplace_back(J);
	return output;
}

void RebarLayer::print() {
	if(!is_initialized()) return;

	suanpan_info("Material model response:\n");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("Integration Point %lu:\n", I + 1);
		int_pt[I].m_material->print();
	}
}
