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

#include "SRGCMQ.h"
#include <Domain/DomainBase.h>
#include <Material/Material1D/Material1D.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>

const mat SRGCMQ::iso_mapping;

SRGCMQ::IntegrationPoint::IntegrationPoint(vec&& C, const double F, unique_ptr<Material>&& MAT)
	: coor(std::forward<vec>(C))
	, factor(F)
	, m_material(std::forward<unique_ptr<Material>>(MAT)) {}

mat SRGCMQ::form_transformation(const mat& jacobian) {
	mat trans_mat(3, 3);

	trans_mat(0, 0) = jacobian(0, 0) * jacobian(0, 0);
	trans_mat(1, 0) = jacobian(0, 1) * jacobian(0, 1);
	trans_mat(2, 0) = jacobian(0, 0) * jacobian(0, 1);

	trans_mat(0, 1) = jacobian(1, 0) * jacobian(1, 0);
	trans_mat(1, 1) = jacobian(1, 1) * jacobian(1, 1);
	trans_mat(2, 1) = jacobian(1, 0) * jacobian(1, 1);

	trans_mat(0, 2) = 2. * jacobian(0, 0) * jacobian(1, 0);
	trans_mat(1, 2) = 2. * jacobian(1, 0) * jacobian(1, 1);
	trans_mat(2, 2) = jacobian(0, 0) * jacobian(1, 1) + jacobian(0, 1) * jacobian(1, 0);

	return trans_mat / accu(square(trans_mat));
}

mat SRGCMQ::form_drilling_mass(const vec& coor, const vec& lxy) {
	mat poly_mass(2, m_size, fill::zeros);

	auto &X = coor(0), &Y = coor(1);

	auto &LX1 = lxy(0), &LX2 = lxy(1), &LX3 = lxy(2), &LX4 = lxy(3);
	auto &LY1 = lxy(4), &LY2 = lxy(5), &LY3 = lxy(6), &LY4 = lxy(7);

	const auto XX = X * X - 1., YY = Y * Y - 1., YP = 1. + Y, YM = 1. - Y, XP = 1. + X, XM = 1. - X;

	poly_mass(0, 2) = LX1 * XX * YM - LX4 * YY * XM;
	poly_mass(0, 5) = LX2 * YY * XP - LX1 * XX * YM;
	poly_mass(0, 8) = LX3 * XX * YP - LX2 * YY * XP;
	poly_mass(0, 11) = LX4 * YY * XM - LX3 * XX * YP;
	poly_mass(1, 2) = LY1 * XX * YM - LY4 * YY * XM;
	poly_mass(1, 5) = LY2 * YY * XP - LY1 * XX * YM;
	poly_mass(1, 8) = LY3 * XX * YP - LY2 * YY * XP;
	poly_mass(1, 11) = LY4 * YY * XM - LY3 * XX * YP;

	poly_mass /= 16.;

	return poly_mass;
}

mat SRGCMQ::form_drilling_displacement(const vec& coor, const vec& lxy) {
	mat poly_drilling(2, 8);

	auto &X = coor(0), &Y = coor(1);

	auto &LX1 = lxy(0), &LX2 = lxy(1), &LX3 = lxy(2), &LX4 = lxy(3);
	auto &LY1 = lxy(4), &LY2 = lxy(5), &LY3 = lxy(6), &LY4 = lxy(7);

	const auto X2 = 2. * X, Y2 = 2. * Y, XP = X + 1., XM = X - 1., YP = Y + 1., YM = Y - 1.;

	poly_drilling(0, 0) = YM * (LX4 * YP - LX1 * X2);
	poly_drilling(0, 1) = YM * (LX2 * YP + LX1 * X2);
	poly_drilling(0, 2) = YP * (LX3 * X2 - LX2 * YM);
	poly_drilling(0, 3) = -YP * (LX3 * X2 + LX4 * YM);
	poly_drilling(0, 4) = YM * (LY4 * YP - LY1 * X2);
	poly_drilling(0, 5) = YM * (LY2 * YP + LY1 * X2);
	poly_drilling(0, 6) = YP * (LY3 * X2 - LY2 * YM);
	poly_drilling(0, 7) = -YP * (LY3 * X2 + LY4 * YM);
	poly_drilling(1, 0) = XM * (LX4 * Y2 - LX1 * XP);
	poly_drilling(1, 1) = XP * (LX1 * XM + LX2 * Y2);
	poly_drilling(1, 2) = XP * (LX3 * XM - LX2 * Y2);
	poly_drilling(1, 3) = -XM * (LX3 * XP + LX4 * Y2);
	poly_drilling(1, 4) = XM * (LY4 * Y2 - LY1 * XP);
	poly_drilling(1, 5) = XP * (LY1 * XM + LY2 * Y2);
	poly_drilling(1, 6) = XP * (LY3 * XM - LY2 * Y2);
	poly_drilling(1, 7) = -XM * (LY3 * XP + LY4 * Y2);

	poly_drilling /= 16.;

	return poly_drilling;
}

mat SRGCMQ::form_displacement(const mat& pn_pxy, const mat& pnt_pxy) {
	mat poly_disp(3, m_size, fill::zeros);

	for(unsigned J = 0, K = 0, L = 1, M = 2; J < m_node; ++J, K += m_dof, L += m_dof, M += m_dof) {
		poly_disp(0, K) = poly_disp(2, L) = pn_pxy(0, J);
		poly_disp(2, K) = poly_disp(1, L) = pn_pxy(1, J);
		poly_disp(0, M) = pnt_pxy(0, J);
		poly_disp(1, M) = pnt_pxy(1, J + 4llu);
		poly_disp(2, M) = pnt_pxy(0, J + 4llu) + pnt_pxy(1, J);
	}

	return poly_disp;
}

vec SRGCMQ::form_stress_mode(const double X, const double Y) { return vec{0., X, Y, X * Y}; }

SRGCMQ::SRGCMQ(const unsigned T, uvec&& NN, const unsigned MAT, const double TH, const char IP, const double RHOX, const unsigned MATX, const double RHOY, const unsigned MATY)
	: MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(NN), uvec{MAT}, false)
	, thickness(TH)
	, int_scheme(IP)
	, rho_x(RHOX)
	, rho_y(RHOY)
	, mat_x(MATX)
	, mat_y(MATY) {
	if(iso_mapping.is_empty()) {
		mat t_mapping(4, 4);
		t_mapping.fill(.25);
		t_mapping(1, 0) = t_mapping(1, 3) = t_mapping(2, 0) = t_mapping(2, 1) = t_mapping(3, 1) = t_mapping(3, 3) = -.25;
		access::rw(iso_mapping) = t_mapping;
	}
}

void SRGCMQ::initialize(const shared_ptr<DomainBase>& D) {
	const auto ele_coor = get_coordinate(2);

	access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

	const mat trans_mat = trans(iso_mapping * ele_coor);

	const auto jacob_trans = form_transformation(shape::quad(vec{0., 0.}, 1) * ele_coor);

	vec diff_coor(8);
	diff_coor(0) = ele_coor(1, 1) - ele_coor(0, 1);
	diff_coor(1) = ele_coor(2, 1) - ele_coor(1, 1);
	diff_coor(2) = ele_coor(3, 1) - ele_coor(2, 1);
	diff_coor(3) = ele_coor(0, 1) - ele_coor(3, 1);
	diff_coor(4) = ele_coor(0, 0) - ele_coor(1, 0);
	diff_coor(5) = ele_coor(1, 0) - ele_coor(2, 0);
	diff_coor(6) = ele_coor(2, 0) - ele_coor(3, 0);
	diff_coor(7) = ele_coor(3, 0) - ele_coor(0, 0);

	auto& mat_proto = D->get<Material>(material_tag(0));

	if(static_cast<double>(PlaneType::E) == mat_proto->get_parameter(ParameterType::PLANETYPE)) suanpan::hacker(thickness) = 1.;

	auto mat_stiff = mat_proto->get_initial_stiffness();

	if(reinforced_x) {
		if(D->find_material(mat_x) && std::dynamic_pointer_cast<Material1D>(D->get_material(mat_x)) != nullptr) mat_stiff(0, 0) += D->get_material(mat_x)->get_initial_stiffness().at(0) * rho_x;
		else {
			suanpan_error("SRGCMQ: cannot find the reinforcement material.\n");
			D->disable_element(get_tag());
			return;
		}
	}
	if(reinforced_y) {
		if(D->find_material(mat_y) && std::dynamic_pointer_cast<Material1D>(D->get_material(mat_y)) != nullptr) mat_stiff(1, 1) += D->get_material(mat_y)->get_initial_stiffness().at(0) * rho_y;
		else {
			suanpan_error("SRGCMQ: cannot find the reinforcement material.\n");
			D->disable_element(get_tag());
			return;
		}
	}

	const IntegrationPlan plan(2, int_scheme == 'I' ? 2 : 3, int_scheme == 'I' ? IntegrationType::IRONS : int_scheme == 'L' ? IntegrationType::LOBATTO : IntegrationType::GAUSS);

	mat N(11, 12, fill::zeros), H(11, 11, fill::zeros), HT(11, 11, fill::zeros);

	int_pt.clear();
	int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		const auto &X = plan(I, 0), &Y = plan(I, 1);

		vec t_vec{X, Y};
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), det(jacob) * plan(I, 2) * thickness, mat_proto->get_copy());

		auto& c_pt = int_pt.back();

		if(reinforced_x) c_pt.rebar_x = D->get_material(mat_x)->get_copy();
		if(reinforced_y) c_pt.rebar_y = D->get_material(mat_y)->get_copy();

		c_pt.poly_strain = shape::stress11(trans_mat * form_stress_mode(X, Y));
		const auto& poly_stress = c_pt.poly_strain;

		N += c_pt.factor * poly_stress.t() * form_displacement(solve(jacob, pn), solve(jacob, form_drilling_displacement(c_pt.coor, diff_coor)));
		H += c_pt.factor * poly_stress.t() * c_pt.poly_strain;
		HT += c_pt.factor * c_pt.poly_strain.t() * mat_stiff * c_pt.poly_strain;
	}

	const mat NT = solve(H, N);

	trial_stiffness = current_stiffness = initial_stiffness = NT.t() * HT * NT;

	for(auto& I : int_pt) I.poly_strain *= NT;

	if(const auto t_density = mat_proto->get_parameter(); t_density > 0.) {
		initial_mass.zeros(m_size, m_size);
		for(const auto& I : int_pt) {
			const auto n_int = shape::quad(I.coor, 0);
			const auto t_factor = t_density * I.factor;
			for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, M = L; K < m_node; ++K, M += m_dof) initial_mass(L, M) += t_factor * n_int(J) * n_int(K);
		}
		for(auto I = 0u, K = 1u; I < m_size; I += m_dof, K += m_dof) {
			initial_mass(K, K) = initial_mass(I, I);
			for(auto J = I + m_dof, L = K + m_dof; J < m_size; J += m_dof, L += m_dof) initial_mass(J, I) = initial_mass(K, L) = initial_mass(L, K) = initial_mass(I, J);
		}
		for(const auto& I : int_pt) {
			const auto n_int = form_drilling_mass(I.coor, diff_coor);
			initial_mass += n_int.t() * n_int * t_density * I.factor;
		}
		ConstantMass(this);
	}

	body_force.zeros(m_size, m_dof);
	for(const auto& I : int_pt) {
		const auto n_int = shape::quad(I.coor, 0);
		const auto n_int_drill = form_drilling_mass(I.coor, diff_coor);
		for(auto J = 0u, L = 0u, M = 2u; J < m_node; ++J, L += m_dof, M += m_dof) {
			for(auto K = 0llu; K < 2; ++K) body_force(L + K, K) += I.factor * n_int(J);
			body_force(M, 2) += I.factor * n_int_drill(J);
		}
	}
}

int SRGCMQ::update_status() {
	const auto trial_disp = get_trial_displacement();

	trial_resistance.zeros(m_size);
	trial_stiffness.zeros(m_size, m_size);
	for(const auto& t_pt : int_pt) {
		const vec t_strain = t_pt.poly_strain * trial_disp;
		if(t_pt.m_material->update_trial_status(t_strain) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		if(!reinforced_x && !reinforced_y) {
			trial_resistance += t_pt.poly_strain.t() * t_pt.factor * t_pt.m_material->get_trial_stress();
			trial_stiffness += t_pt.poly_strain.t() * t_pt.factor * t_pt.m_material->get_trial_stiffness() * t_pt.poly_strain;
		} else {
			auto t_stress = t_pt.m_material->get_trial_stress();
			auto t_stiffness = t_pt.m_material->get_trial_stiffness();
			if(reinforced_x) {
				t_pt.rebar_x->update_trial_status(t_strain(0));
				t_stress(0) += t_pt.rebar_x->get_trial_stress().at(0) * rho_x;
				t_stiffness(0, 0) += t_pt.rebar_x->get_trial_stiffness().at(0) * rho_x;
			}
			if(reinforced_y) {
				t_pt.rebar_y->update_trial_status(t_strain(1));
				t_stress(1) += t_pt.rebar_y->get_trial_stress().at(0) * rho_y;
				t_stiffness(1, 1) += t_pt.rebar_y->get_trial_stiffness().at(0) * rho_y;
			}
			trial_resistance += t_pt.poly_strain.t() * t_pt.factor * t_stress;
			trial_stiffness += t_pt.poly_strain.t() * t_pt.factor * t_stiffness * t_pt.poly_strain;
		}
	}

	return SUANPAN_SUCCESS;
}

int SRGCMQ::commit_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->commit_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->commit_status();
	return code;
}

int SRGCMQ::clear_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->clear_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->clear_status();
	return code;
}

int SRGCMQ::reset_status() {
	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	if(reinforced_x) for(const auto& I : int_pt) code += I.rebar_x->reset_status();
	if(reinforced_y) for(const auto& I : int_pt) code += I.rebar_y->reset_status();
	return code;
}

vector<vec> SRGCMQ::record(const OutputType T) {
	vector<vec> data;

	if(T == OutputType::REBARE) {
		if(reinforced_x) for(const auto& I : int_pt) for(const auto& J : I.rebar_x->record(OutputType::E)) data.emplace_back(J);
		if(reinforced_y) for(const auto& I : int_pt) for(const auto& J : I.rebar_y->record(OutputType::E)) data.emplace_back(J);
	} else if(T == OutputType::REBARS) {
		if(reinforced_x) for(const auto& I : int_pt) for(const auto& J : I.rebar_x->record(OutputType::S)) data.emplace_back(J);
		if(reinforced_y) for(const auto& I : int_pt) for(const auto& J : I.rebar_y->record(OutputType::S)) data.emplace_back(J);
	} else for(const auto& I : int_pt) for(const auto& J : I.m_material->record(T)) data.emplace_back(J);

	return data;
}

void SRGCMQ::print() {
	suanpan_info("SRGCMQ mixed quad element %u connects nodes:\n", get_tag());
	node_encoding.t().print();
	if(!is_initialized()) return;
	suanpan_info("\nMaterial model response:");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("\nIntegration Point %lu:\t", I + 1);
		int_pt[I].coor.t().print();
		int_pt[I].m_material->print();
	}
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void SRGCMQ::Setup() {
	vtk_cell = vtkSmartPointer<vtkQuad>::New();
	const auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < m_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void SRGCMQ::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
	mat t_disp(6, m_node, fill::zeros);

	if(OutputType::A == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_acceleration(), m_dof, m_node);
	else if(OutputType::V == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_velocity(), m_dof, m_node);
	else if(OutputType::U == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), m_dof, m_node);

	for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

mat SRGCMQ::GetData(const OutputType P) {
	if(P == OutputType::S) {
		mat B(3 * int_pt.size(), 1);
		mat A(3 * int_pt.size(), 11);
		for(uword I = 0, J = 0; I < int_pt.size(); ++I, J = J + 3) {
			B.rows(J, J + 2llu) = int_pt[I].m_material->get_current_stress();
			A.rows(J, J + 2llu) = shape::stress11(int_pt[I].coor);
		}
		mat X;
		if(!solve(X, A, B)) return {};

		mat t_stress(6, m_node, fill::zeros);

		t_stress(uvec{0, 1, 3}, uvec{0}) = shape::stress11(-1., -1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{1}) = shape::stress11(1., -1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{2}) = shape::stress11(1., 1.) * X;
		t_stress(uvec{0, 1, 3}, uvec{3}) = shape::stress11(-1., 1.) * X;

		return t_stress;
	}

	mat A(int_pt.size(), 9);
	mat B(int_pt.size(), 6, fill::zeros);

	for(size_t I = 0; I < int_pt.size(); ++I) {
		if(const auto C = int_pt[I].m_material->record(P); !C.empty()) B(I, 0, size(C[0])) = C[0];
		A.row(I) = interpolation::quadratic(int_pt[I].coor);
	}

	mat data(m_node, 9);

	data.row(0) = interpolation::quadratic(-1., -1.);
	data.row(1) = interpolation::quadratic(1., -1.);
	data.row(2) = interpolation::quadratic(1., 1.);
	data.row(3) = interpolation::quadratic(-1., 1.);

	return (data * solve(A, B)).t();
}

void SRGCMQ::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * mat(reshape(get_current_displacement(), m_dof, m_node)).rows(0, 1).t();
	for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
