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

#include "RCP4.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>

RCP4::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M, mat&& PNPXY)
	: coor(std::forward<vec>(C))
	, weight(W)
	, m_material(std::forward<unique_ptr<Material>>(M))
	, pn_pxy(std::forward<mat>(PNPXY))
	, B(3, m_size, fill::zeros) {}

void RCP4::stack_stiffness(mat& K, const mat& D, const mat& N, const double F) {
	const auto D11 = F * D(0, 0);
	const auto D12 = F * D(0, 1);
	const auto D13 = F * D(0, 2);
	const auto D21 = F * D(1, 0);
	const auto D22 = F * D(1, 1);
	const auto D23 = F * D(1, 2);
	const auto D31 = F * D(2, 0);
	const auto D32 = F * D(2, 1);
	const auto D33 = F * D(2, 2);

	const auto& NX1 = N(0, 0);
	const auto& NY1 = N(1, 1);
	const auto& NX2 = N(0, 2);
	const auto& NY2 = N(1, 3);
	const auto& NX3 = N(0, 4);
	const auto& NY3 = N(1, 5);
	const auto& NX4 = N(0, 6);
	const auto& NY4 = N(1, 7);

	const auto D11NX1 = D11 * NX1;
	const auto D11NX2 = D11 * NX2;
	const auto D11NX3 = D11 * NX3;
	const auto D11NX4 = D11 * NX4;

	const auto D12NX1 = D12 * NX1;
	const auto D12NX2 = D12 * NX2;
	const auto D12NX3 = D12 * NX3;
	const auto D12NX4 = D12 * NX4;

	const auto D13NX1 = D13 * NX1;
	const auto D13NX2 = D13 * NX2;
	const auto D13NX3 = D13 * NX3;
	const auto D13NX4 = D13 * NX4;

	const auto D21NY1 = D21 * NY1;
	const auto D21NY2 = D21 * NY2;
	const auto D21NY3 = D21 * NY3;
	const auto D21NY4 = D21 * NY4;

	const auto D22NY1 = D22 * NY1;
	const auto D22NY2 = D22 * NY2;
	const auto D22NY3 = D22 * NY3;
	const auto D22NY4 = D22 * NY4;

	const auto D23NY1 = D23 * NY1;
	const auto D23NY2 = D23 * NY2;
	const auto D23NY3 = D23 * NY3;
	const auto D23NY4 = D23 * NY4;

	const auto D31NX1 = D31 * NX1;
	const auto D31NX2 = D31 * NX2;
	const auto D31NX3 = D31 * NX3;
	const auto D31NX4 = D31 * NX4;
	const auto D31NY1 = D31 * NY1;
	const auto D31NY2 = D31 * NY2;
	const auto D31NY3 = D31 * NY3;
	const auto D31NY4 = D31 * NY4;

	const auto D32NX1 = D32 * NX1;
	const auto D32NX2 = D32 * NX2;
	const auto D32NX3 = D32 * NX3;
	const auto D32NX4 = D32 * NX4;
	const auto D32NY1 = D32 * NY1;
	const auto D32NY2 = D32 * NY2;
	const auto D32NY3 = D32 * NY3;
	const auto D32NY4 = D32 * NY4;

	const auto D33NX1 = D33 * NX1;
	const auto D33NX2 = D33 * NX2;
	const auto D33NX3 = D33 * NX3;
	const auto D33NX4 = D33 * NX4;
	const auto D33NY1 = D33 * NY1;
	const auto D33NY2 = D33 * NY2;
	const auto D33NY3 = D33 * NY3;
	const auto D33NY4 = D33 * NY4;

	const auto D11NX1D31NY1 = D11NX1 + D31NY1;
	const auto D13NX1D33NY1 = D13NX1 + D33NY1;
	const auto D12NX1D32NY1 = D12NX1 + D32NY1;
	const auto D31NX1D21NY1 = D31NX1 + D21NY1;
	const auto D33NX1D23NY1 = D33NX1 + D23NY1;
	const auto D32NX1D22NY1 = D32NX1 + D22NY1;
	const auto D11NX2D31NY2 = D11NX2 + D31NY2;
	const auto D13NX2D33NY2 = D13NX2 + D33NY2;
	const auto D12NX2D32NY2 = D12NX2 + D32NY2;
	const auto D31NX2D21NY2 = D31NX2 + D21NY2;
	const auto D33NX2D23NY2 = D33NX2 + D23NY2;
	const auto D32NX2D22NY2 = D32NX2 + D22NY2;
	const auto D11NX3D31NY3 = D11NX3 + D31NY3;
	const auto D13NX3D33NY3 = D13NX3 + D33NY3;
	const auto D12NX3D32NY3 = D12NX3 + D32NY3;
	const auto D31NX3D21NY3 = D31NX3 + D21NY3;
	const auto D33NX3D23NY3 = D33NX3 + D23NY3;
	const auto D32NX3D22NY3 = D32NX3 + D22NY3;
	const auto D11NX4D31NY4 = D11NX4 + D31NY4;
	const auto D13NX4D33NY4 = D13NX4 + D33NY4;
	const auto D12NX4D32NY4 = D12NX4 + D32NY4;
	const auto D31NX4D21NY4 = D31NX4 + D21NY4;
	const auto D33NX4D23NY4 = D33NX4 + D23NY4;
	const auto D32NX4D22NY4 = D32NX4 + D22NY4;

	K(0, 0) += NX1 * D11NX1D31NY1 + NY1 * D13NX1D33NY1;
	K(0, 1) += NX1 * D13NX1D33NY1 + NY1 * D12NX1D32NY1;
	K(0, 2) += NX2 * D11NX1D31NY1 + NY2 * D13NX1D33NY1;
	K(0, 3) += NX2 * D13NX1D33NY1 + NY2 * D12NX1D32NY1;
	K(0, 4) += NX3 * D11NX1D31NY1 + NY3 * D13NX1D33NY1;
	K(0, 5) += NX3 * D13NX1D33NY1 + NY3 * D12NX1D32NY1;
	K(0, 6) += NX4 * D11NX1D31NY1 + NY4 * D13NX1D33NY1;
	K(0, 7) += NX4 * D13NX1D33NY1 + NY4 * D12NX1D32NY1;
	K(1, 0) += NX1 * D31NX1D21NY1 + NY1 * D33NX1D23NY1;
	K(1, 1) += NX1 * D33NX1D23NY1 + NY1 * D32NX1D22NY1;
	K(1, 2) += NX2 * D31NX1D21NY1 + NY2 * D33NX1D23NY1;
	K(1, 3) += NX2 * D33NX1D23NY1 + NY2 * D32NX1D22NY1;
	K(1, 4) += NX3 * D31NX1D21NY1 + NY3 * D33NX1D23NY1;
	K(1, 5) += NX3 * D33NX1D23NY1 + NY3 * D32NX1D22NY1;
	K(1, 6) += NX4 * D31NX1D21NY1 + NY4 * D33NX1D23NY1;
	K(1, 7) += NX4 * D33NX1D23NY1 + NY4 * D32NX1D22NY1;
	K(2, 0) += NX1 * D11NX2D31NY2 + NY1 * D13NX2D33NY2;
	K(2, 1) += NX1 * D13NX2D33NY2 + NY1 * D12NX2D32NY2;
	K(2, 2) += NX2 * D11NX2D31NY2 + NY2 * D13NX2D33NY2;
	K(2, 3) += NX2 * D13NX2D33NY2 + NY2 * D12NX2D32NY2;
	K(2, 4) += NX3 * D11NX2D31NY2 + NY3 * D13NX2D33NY2;
	K(2, 5) += NX3 * D13NX2D33NY2 + NY3 * D12NX2D32NY2;
	K(2, 6) += NX4 * D11NX2D31NY2 + NY4 * D13NX2D33NY2;
	K(2, 7) += NX4 * D13NX2D33NY2 + NY4 * D12NX2D32NY2;
	K(3, 0) += NX1 * D31NX2D21NY2 + NY1 * D33NX2D23NY2;
	K(3, 1) += NX1 * D33NX2D23NY2 + NY1 * D32NX2D22NY2;
	K(3, 2) += NX2 * D31NX2D21NY2 + NY2 * D33NX2D23NY2;
	K(3, 3) += NX2 * D33NX2D23NY2 + NY2 * D32NX2D22NY2;
	K(3, 4) += NX3 * D31NX2D21NY2 + NY3 * D33NX2D23NY2;
	K(3, 5) += NX3 * D33NX2D23NY2 + NY3 * D32NX2D22NY2;
	K(3, 6) += NX4 * D31NX2D21NY2 + NY4 * D33NX2D23NY2;
	K(3, 7) += NX4 * D33NX2D23NY2 + NY4 * D32NX2D22NY2;
	K(4, 0) += NX1 * D11NX3D31NY3 + NY1 * D13NX3D33NY3;
	K(4, 1) += NX1 * D13NX3D33NY3 + NY1 * D12NX3D32NY3;
	K(4, 2) += NX2 * D11NX3D31NY3 + NY2 * D13NX3D33NY3;
	K(4, 3) += NX2 * D13NX3D33NY3 + NY2 * D12NX3D32NY3;
	K(4, 4) += NX3 * D11NX3D31NY3 + NY3 * D13NX3D33NY3;
	K(4, 5) += NX3 * D13NX3D33NY3 + NY3 * D12NX3D32NY3;
	K(4, 6) += NX4 * D11NX3D31NY3 + NY4 * D13NX3D33NY3;
	K(4, 7) += NX4 * D13NX3D33NY3 + NY4 * D12NX3D32NY3;
	K(5, 0) += NX1 * D31NX3D21NY3 + NY1 * D33NX3D23NY3;
	K(5, 1) += NX1 * D33NX3D23NY3 + NY1 * D32NX3D22NY3;
	K(5, 2) += NX2 * D31NX3D21NY3 + NY2 * D33NX3D23NY3;
	K(5, 3) += NX2 * D33NX3D23NY3 + NY2 * D32NX3D22NY3;
	K(5, 4) += NX3 * D31NX3D21NY3 + NY3 * D33NX3D23NY3;
	K(5, 5) += NX3 * D33NX3D23NY3 + NY3 * D32NX3D22NY3;
	K(5, 6) += NX4 * D31NX3D21NY3 + NY4 * D33NX3D23NY3;
	K(5, 7) += NX4 * D33NX3D23NY3 + NY4 * D32NX3D22NY3;
	K(6, 0) += NX1 * D11NX4D31NY4 + NY1 * D13NX4D33NY4;
	K(6, 1) += NX1 * D13NX4D33NY4 + NY1 * D12NX4D32NY4;
	K(6, 2) += NX2 * D11NX4D31NY4 + NY2 * D13NX4D33NY4;
	K(6, 3) += NX2 * D13NX4D33NY4 + NY2 * D12NX4D32NY4;
	K(6, 4) += NX3 * D11NX4D31NY4 + NY3 * D13NX4D33NY4;
	K(6, 5) += NX3 * D13NX4D33NY4 + NY3 * D12NX4D32NY4;
	K(6, 6) += NX4 * D11NX4D31NY4 + NY4 * D13NX4D33NY4;
	K(6, 7) += NX4 * D13NX4D33NY4 + NY4 * D12NX4D32NY4;
	K(7, 0) += NX1 * D31NX4D21NY4 + NY1 * D33NX4D23NY4;
	K(7, 1) += NX1 * D33NX4D23NY4 + NY1 * D32NX4D22NY4;
	K(7, 2) += NX2 * D31NX4D21NY4 + NY2 * D33NX4D23NY4;
	K(7, 3) += NX2 * D33NX4D23NY4 + NY2 * D32NX4D22NY4;
	K(7, 4) += NX3 * D31NX4D21NY4 + NY3 * D33NX4D23NY4;
	K(7, 5) += NX3 * D33NX4D23NY4 + NY3 * D32NX4D22NY4;
	K(7, 6) += NX4 * D31NX4D21NY4 + NY4 * D33NX4D23NY4;
	K(7, 7) += NX4 * D33NX4D23NY4 + NY4 * D32NX4D22NY4;
}

RCP4::RCP4(const unsigned T, uvec&& N, const unsigned M, const double TH, const double RX, const unsigned MX, const double RY, const unsigned MY, const bool R, const bool F)
	: MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, F)
	, thickness(TH)
	, reduced_scheme(R)
	, rho_x(RX)
	, rho_y(RY)
	, mat_x(MX)
	, mat_y(MY) {}

void RCP4::initialize(const shared_ptr<DomainBase>& D) {
	auto& material_proto = D->get<Material>(material_tag(0));

	if(static_cast<double>(PlaneType::E) == material_proto->get_parameter(ParameterType::PLANETYPE)) suanpan::hacker(thickness) = 1.;

	const auto ele_coor = get_coordinate(2);

	access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

	if(reduced_scheme) {
		hourglassing.zeros(m_size, m_size);
		const auto area = .5 * ((ele_coor(2, 0) - ele_coor(0, 0)) * (ele_coor(3, 1) - ele_coor(1, 1)) + (ele_coor(1, 0) - ele_coor(3, 0)) * (ele_coor(2, 1) - ele_coor(0, 1)));
		vec b1(4), b2(4);
		b1(0) = ele_coor(1, 1) - ele_coor(3, 1);
		b1(1) = ele_coor(2, 1) - ele_coor(0, 1);
		b1(2) = ele_coor(3, 1) - ele_coor(1, 1);
		b1(3) = ele_coor(0, 1) - ele_coor(2, 1);
		b2(0) = ele_coor(3, 0) - ele_coor(1, 0);
		b2(1) = ele_coor(0, 0) - ele_coor(2, 0);
		b2(2) = ele_coor(1, 0) - ele_coor(3, 0);
		b2(3) = ele_coor(2, 0) - ele_coor(0, 0);
		const vec h{std::initializer_list<double>{1., -1., 1., -1.}};
		vec gamma = 2. * area * h - dot(h, ele_coor.col(0)) * b1 - dot(h, ele_coor.col(1)) * b2;
		mat t_hourglassing = gamma * gamma.t();
		for(uword I = 0, K = 0, M = 1; I < m_node; ++I, K += m_dof, M += m_dof) for(uword J = 0, L = 0, N = 1; J < m_node; ++J, L += m_dof, N += m_dof) hourglassing(M, N) = hourglassing(K, L) = t_hourglassing(I, J);
	}

	auto ini_stiffness = material_proto->get_initial_stiffness();

	if(reinforced_x) {
		if(D->find_material(mat_x)) ini_stiffness(0, 0) += D->get_material(mat_x)->get_initial_stiffness().at(0) * rho_x;
		else {
			suanpan_error("RCP4: cannot find the reinforcement material.\n");
			D->disable_element(get_tag());
			return;
		}
	}
	if(reinforced_y) {
		if(D->find_material(mat_y)) ini_stiffness(1, 1) += D->get_material(mat_y)->get_initial_stiffness().at(0) * rho_y;
		else {
			suanpan_error("RCP4: cannot find the reinforcement material.\n");
			D->disable_element(get_tag());
			return;
		}
	}

	const IntegrationPlan plan(2, reduced_scheme ? 1 : 2, IntegrationType::GAUSS);

	initial_stiffness.zeros(m_size, m_size);

	int_pt.clear();
	int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1)};
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), plan(I, 2) * det(jacob), material_proto->get_copy(), solve(jacob, pn));

		auto& c_pt = int_pt.back();

		if(reinforced_x) c_pt.rebar_x = D->get_material(mat_x)->get_copy();
		if(reinforced_y) c_pt.rebar_y = D->get_material(mat_y)->get_copy();

		for(auto J = 0u, K = 0u, L = 1u; J < m_node; ++J, K += m_dof, L += m_dof) {
			c_pt.B(0, K) = c_pt.B(2, L) = c_pt.pn_pxy(0, J);
			c_pt.B(2, K) = c_pt.B(1, L) = c_pt.pn_pxy(1, J);
		}
		initial_stiffness += c_pt.weight * thickness * c_pt.B.t() * ini_stiffness * c_pt.B;
	}
	trial_stiffness = current_stiffness = initial_stiffness;

	if(nlgeom)
		for(auto&& I : int_pt) {
			I.BG.zeros(m_dof * m_dof, m_size);
			for(uword J = 0; J < m_node; ++J) {
				I.BG(0, m_dof * J) = I.BG(2, m_dof * J + 1) = I.pn_pxy(0, J);
				I.BG(1, m_dof * J) = I.BG(3, m_dof * J + 1) = I.pn_pxy(1, J);
			}
		}

	if(const auto t_density = material_proto->get_parameter(); t_density > 0.) {
		initial_mass.zeros(m_size, m_size);
		for(const auto& I : int_pt) {
			const auto n_int = shape::quad(I.coor, 0);
			const auto t_factor = t_density * I.weight * thickness;
			for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
		}
		for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
			initial_mass(K, K) = initial_mass(I, I);
			for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
		}
		ConstantMass(this);
	}
}

int RCP4::update_status() {
	auto code = 0;

	const auto t_disp = get_trial_displacement();
	const mat ele_disp = reshape(t_disp, m_dof, m_node);

	vec t_strain(3);

	trial_stiffness.zeros(m_size, m_size);
	trial_resistance.zeros(m_size);

	if(nlgeom) {
		trial_geometry.zeros(m_size, m_size);

		mat sigma(4, 4, fill::zeros);
		mat BN(3, m_size, fill::zeros);
		for(auto& I : int_pt) {
			mat gradient = I.pn_pxy * ele_disp.t();
			gradient(0, 0) += 1.;
			gradient(1, 1) += 1.;
			mat t_mat = .5 * gradient * gradient.t();
			t_strain(0) = t_mat(0, 0) - .5;
			t_strain(1) = t_mat(1, 1) - .5;
			t_strain(2) = t_mat(0, 1) + t_mat(1, 0);
			for(auto J = 0u, K = 0u, L = 1u; J < m_node; ++J, K += m_dof, L += m_dof) {
				BN(0, K) = I.pn_pxy(0, J) * gradient(0, 0);
				BN(1, K) = I.pn_pxy(1, J) * gradient(1, 0);
				BN(2, K) = I.pn_pxy(0, J) * gradient(1, 0) + I.pn_pxy(1, J) * gradient(0, 0);
				BN(0, L) = I.pn_pxy(0, J) * gradient(0, 1);
				BN(1, L) = I.pn_pxy(1, J) * gradient(1, 1);
				BN(2, L) = I.pn_pxy(0, J) * gradient(1, 1) + I.pn_pxy(1, J) * gradient(0, 1);
			}

			code += I.m_material->update_trial_status(t_strain);

			const auto t_factor = I.weight * thickness;

			auto t_stress = I.m_material->get_trial_stress();
			auto t_stiff = I.m_material->get_trial_stiffness();

			if(reinforced_x) {
				I.rebar_x->update_trial_status(t_strain(0));
				t_stress(0) += I.rebar_x->get_trial_stress().at(0) * rho_x;
				t_stiff(0, 0) += I.rebar_x->get_trial_stiffness().at(0) * rho_x;
			}
			if(reinforced_y) {
				I.rebar_y->update_trial_status(t_strain(1));
				t_stress(1) += I.rebar_y->get_trial_stress().at(0) * rho_y;
				t_stiff(1, 1) += I.rebar_y->get_trial_stiffness().at(0) * rho_y;
			}

			sigma(0, 0) = sigma(2, 2) = t_stress(0);
			sigma(1, 1) = sigma(3, 3) = t_stress(1);
			sigma(0, 1) = sigma(1, 0) = sigma(2, 3) = sigma(3, 2) = t_stress(2);
			trial_geometry += t_factor * I.BG.t() * sigma * I.BG;
			trial_stiffness += t_factor * BN.t() * t_stiff * BN;
			trial_resistance += t_factor * BN.t() * t_stress;
		}
	} else
		for(const auto& I : int_pt) {
			t_strain.zeros();
			for(auto J = 0u, K = 0u, L = 1u; J < m_node; ++J, K += m_dof, L += m_dof) {
				t_strain(0) += t_disp(K) * I.pn_pxy(0, J);
				t_strain(1) += t_disp(L) * I.pn_pxy(1, J);
				t_strain(2) += t_disp(K) * I.pn_pxy(1, J) + t_disp(L) * I.pn_pxy(0, J);
			}

			code += I.m_material->update_trial_status(t_strain);

			const auto t_factor = I.weight * thickness;

			auto t_stress = I.m_material->get_trial_stress();
			auto t_stiff = I.m_material->get_trial_stiffness();

			if(reinforced_x) {
				I.rebar_x->update_trial_status(t_strain(0));
				t_stress(0) += I.rebar_x->get_trial_stress().at(0) * rho_x;
				t_stiff(0, 0) += I.rebar_x->get_trial_stiffness().at(0) * rho_x;
			}
			if(reinforced_y) {
				I.rebar_y->update_trial_status(t_strain(1));
				t_stress(1) += I.rebar_y->get_trial_stress().at(0) * rho_y;
				t_stiff(1, 1) += I.rebar_y->get_trial_stiffness().at(0) * rho_y;
			}

			stack_stiffness(trial_stiffness, t_stiff, I.B, t_factor);
			trial_resistance += t_factor * I.B.t() * t_stress;
		}

	if(reduced_scheme) {
		trial_stiffness += hourglassing;
		trial_resistance += hourglassing * t_disp;
	}

	return code;
}

int RCP4::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->commit_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->commit_status();
	return code;
}

int RCP4::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->clear_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->clear_status();
	return code;
}

int RCP4::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->reset_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->reset_status();
	return code;
}

vector<vec> RCP4::record(const OutputType P) {
	vector<vec> output;
	output.reserve(int_pt.size());

	for(const auto& I : int_pt) for(const auto& J : I.m_material->record(P)) output.emplace_back(J);

	return output;
}

void RCP4::print() {
	suanpan_info("Element %u is a four-node membrane element (RCP4)%s.\n", get_tag(), nlgeom ? " with nonlinear geomotry (TL formulation)" : "");
	suanpan_info("The nodes connected are:\n");
	node_encoding.t().print();
	if(!is_initialized()) return;
	suanpan_info("\nMaterial model response:");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("\nIntegration Point %lu:\t", I + 1);
		int_pt[I].coor.t().print();
		int_pt[I].m_material->print();
	}
}
