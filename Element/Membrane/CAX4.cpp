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

#include "CAX4.h"
#include <Domain/DomainBase.h>
#include <Material/Material2D/Material2D.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/utility.h>

CAX4::IntegrationPoint::IntegrationPoint(vec&& C, const double W, unique_ptr<Material>&& M)
	: coor(std::forward<vec>(C))
	, weight(W)
	, m_material(std::forward<unique_ptr<Material>>(M))
	, strain_mat(4, m_size, fill::zeros) {}

vec CAX4::isoparametric_mapping(const vec& in) {
	vec out(4);

	out(0) = .25 * (in(0) + in(1) + in(2) + in(3));
	out(1) = .25 * (in(1) + in(2) - in(3) - in(0));
	out(2) = .25 * (in(2) + in(3) - in(0) - in(1));
	out(3) = .25 * (in(0) + in(2) - in(1) - in(3));

	return out;
}

CAX4::CAX4(const unsigned T, uvec&& N, const unsigned M, const bool F)
	: MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(N), uvec{M}, F) {}

void CAX4::initialize(const shared_ptr<DomainBase>& D) {
	auto& material_proto = D->get<Material>(material_tag(0));

	if(!suanpan::approx_equal(static_cast<double>(PlaneType::A), material_proto->get_parameter(ParameterType::PLANETYPE))) {
		D->disable_element(get_tag());
		return;
	}

	const auto ele_coor = get_coordinate(2);

	auto& ini_stiffness = material_proto->get_initial_stiffness();

	const IntegrationPlan plan(2, 2, IntegrationType::GAUSS);

	initial_stiffness.zeros(m_size, m_size);

	int_pt.clear(), int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		vec t_vec{plan(I, 0), plan(I, 1)};
		const auto n = shape::quad(t_vec, 0);
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		const mat pn_pxy = solve(jacob, pn);
		const auto gx = dot(vec{1., t_vec(0), t_vec(1), t_vec(0) * t_vec(1)}, isoparametric_mapping(ele_coor.col(0)));
		int_pt.emplace_back(std::move(t_vec), 2. * datum::pi * gx * plan(I, 2) * det(jacob), material_proto->get_copy());

		auto& c_pt = int_pt.back();

		for(unsigned J = 0, K = 0, L = 1; J < m_node; ++J, K += m_dof, L += m_dof) {
			c_pt.strain_mat(0, K) = c_pt.strain_mat(3, L) = pn_pxy(0, J);
			c_pt.strain_mat(3, K) = c_pt.strain_mat(1, L) = pn_pxy(1, J);
			c_pt.strain_mat(2, K) = n(J) / gx;
		}
		initial_stiffness += c_pt.weight * c_pt.strain_mat.t() * ini_stiffness * c_pt.strain_mat;
	}
	trial_stiffness = current_stiffness = initial_stiffness;

	// trial_mass = current_mass = initial_mass.zeros(m_size, m_size);
}

int CAX4::update_status() {
	const auto t_disp = get_trial_displacement();

	trial_stiffness.zeros(m_size, m_size);
	trial_resistance.zeros(m_size);

	for(const auto& I : int_pt) {
		vec t_strain(4, fill::zeros);
		for(unsigned J = 0, K = 1; J < m_size; J += m_dof, K += m_dof) {
			t_strain(0) += t_disp(J) * I.strain_mat(0, J);
			t_strain(1) += t_disp(K) * I.strain_mat(3, J);
			t_strain(2) += t_disp(J) * I.strain_mat(2, J);
			t_strain(3) += t_disp(J) * I.strain_mat(3, J) + t_disp(K) * I.strain_mat(0, J);
		}

		if(I.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;

		trial_stiffness += I.weight * I.strain_mat.t() * I.m_material->get_trial_stiffness() * I.strain_mat;
		trial_resistance += I.weight * I.strain_mat.t() * I.m_material->get_trial_stress();
	}

	return SUANPAN_SUCCESS;
}

int CAX4::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	return code;
}

int CAX4::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	return code;
}

int CAX4::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	return code;
}

vector<vec> CAX4::record(const OutputType P) {
	vector<vec> output;
	for(const auto& I : int_pt) for(const auto& J : I.m_material->record(P)) output.emplace_back(J);
	return output;
}

void CAX4::print() {
	suanpan_info("Element %u is a four-node axisymmteric element (CAX4)%s.\n", get_tag(), nlgeom ? " with nonlinear geomotry (TL formulation)" : "");
	suanpan_info("The nodes connected are:\n");
	node_encoding.t().print();
	if(!is_initialized()) return;
	suanpan_info("Material model response:\n");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("Integration Point %lu:\t", I + 1);
		int_pt[I].coor.t().print();
		int_pt[I].m_material->print();
	}
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void CAX4::Setup() {
	vtk_cell = vtkSmartPointer<vtkQuad>::New();
	const auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < m_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void CAX4::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
	mat t_disp(6, m_node, fill::zeros);

	if(OutputType::A == type) t_disp.rows(0, 1) = reshape(get_current_acceleration(), m_dof, m_node);
	else if(OutputType::V == type) t_disp.rows(0, 1) = reshape(get_current_velocity(), m_dof, m_node);
	else if(OutputType::U == type) t_disp.rows(0, 1) = reshape(get_current_displacement(), m_dof, m_node);

	for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

void CAX4::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * reshape(get_current_displacement(), m_dof, m_node).t();
	for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
