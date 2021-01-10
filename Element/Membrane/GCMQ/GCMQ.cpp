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

#include "GCMQ.h"
#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Material/Material2D/Material2D.h>
#include <Recorder/OutputType.h>
#include <Toolbox/IntegrationPlan.h>
#include <Toolbox/shapeFunction.h>
#include <Toolbox/tensorToolbox.h>
#include <Toolbox/utility.h>

const mat GCMQ::mapping;

/**
 * \brief create converter for resultant forces
 * \param E edge label
 * \param T element thickness
 * \param NP vector of node pointers
 * \param IP integration plan
 * \param TRANS transformation matrix from parent to global
 */
GCMQ::ResultantConverter::ResultantConverter(const unsigned E, const double T, const vector<weak_ptr<Node>>& NP, const IntegrationPlan& IP, const mat& TRANS)
	: direction_cosine(2, 3) {
	const auto &X1 = IP(0, 0), &X2 = IP(1, 0);

	vec node_i, node_j, int_pt_a, int_pt_b;

	if(E == 1) {
		node_i = NP[0].lock()->get_coordinate();
		node_j = NP[1].lock()->get_coordinate();
		int_pt_a = TRANS * form_stress_mode(X1, -1.);
		int_pt_b = TRANS * form_stress_mode(X2, -1.);
	} else if(E == 2) {
		node_i = NP[1].lock()->get_coordinate();
		node_j = NP[2].lock()->get_coordinate();
		int_pt_a = TRANS * form_stress_mode(1., X1);
		int_pt_b = TRANS * form_stress_mode(1., X2);
	} else if(E == 3) {
		node_i = NP[2].lock()->get_coordinate();
		node_j = NP[3].lock()->get_coordinate();
		int_pt_a = TRANS * form_stress_mode(X1, 1.);
		int_pt_b = TRANS * form_stress_mode(X2, 1.);
	} else if(E == 4) {
		node_i = NP[3].lock()->get_coordinate();
		node_j = NP[0].lock()->get_coordinate();
		int_pt_a = TRANS * form_stress_mode(-1., X1);
		int_pt_b = TRANS * form_stress_mode(-1., X2);
	} else throw logic_error("need valid edge tag");

	const auto incre_y = node_j(1) - node_i(1);
	const auto incre_x = node_j(0) - node_i(0);

	const auto edge_length = sqrt(incre_x * incre_x + incre_y * incre_y);

	const auto angle = 2. * atan2(incre_y, incre_x) - datum::pi;

	const auto sin_angle = sin(angle);
	const auto cos_angle = cos(angle);

	// transformation to local reference frame
	direction_cosine(0, 0) = .5 + .5 * cos_angle;
	direction_cosine(0, 1) = .5 - .5 * cos_angle;
	direction_cosine(0, 2) = sin_angle;
	direction_cosine(1, 0) = -(direction_cosine(1, 1) = .5 * sin_angle);
	direction_cosine(1, 2) = cos_angle;

	auto weight = .5 * edge_length * T;

	const mat part_a = shape::stress11(int_pt_a) * IP(0, 1) * weight;
	const mat part_b = shape::stress11(int_pt_b) * IP(1, 1) * weight;
	converter_a = part_a + part_b;
	converter_b = .5 * edge_length * (part_a * X1 + part_b * X2);
}

double GCMQ::ResultantConverter::F(const vec& alpha) const { return dot(direction_cosine.row(0), converter_a * alpha); }

double GCMQ::ResultantConverter::V(const vec& alpha) const { return dot(direction_cosine.row(1), converter_a * alpha); }

double GCMQ::ResultantConverter::M(const vec& alpha) const { return dot(direction_cosine.row(0), converter_b * alpha); }

GCMQ::IntegrationPoint::IntegrationPoint(vec&& C, const double F, unique_ptr<Material>&& MAT)
	: coor(std::forward<vec>(C))
	, factor(F)
	, m_material(std::forward<unique_ptr<Material>>(MAT)) {}

mat GCMQ::form_transformation(const mat& jacobian) {
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

mat GCMQ::form_drilling_mass(const vec& coor, const vec& lxy) {
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

mat GCMQ::form_drilling_displacement(const vec& coor, const vec& lxy) {
	mat poly_drilling(2, 8);

	auto &X = coor(0), &Y = coor(1);

	auto &LX1 = lxy(0), &LX2 = lxy(1), &LX3 = lxy(2), &LX4 = lxy(3);
	auto &LY1 = lxy(4), &LY2 = lxy(5), &LY3 = lxy(6), &LY4 = lxy(7);

	const auto X2 = 2. * X, Y2 = 2. * Y, XP = X + 1., XM = X - 1., YP = Y + 1., YM = Y - 1.;

	poly_drilling(0, 0) = +YM * (LX4 * YP - LX1 * X2);
	poly_drilling(0, 1) = +YM * (LX2 * YP + LX1 * X2);
	poly_drilling(0, 2) = -YP * (LX2 * YM - LX3 * X2);
	poly_drilling(0, 3) = -YP * (LX4 * YM + LX3 * X2);
	poly_drilling(0, 4) = +YM * (LY4 * YP - LY1 * X2);
	poly_drilling(0, 5) = +YM * (LY2 * YP + LY1 * X2);
	poly_drilling(0, 6) = -YP * (LY2 * YM - LY3 * X2);
	poly_drilling(0, 7) = -YP * (LY4 * YM + LY3 * X2);
	poly_drilling(1, 0) = -XM * (LX1 * XP - LX4 * Y2);
	poly_drilling(1, 1) = +XP * (LX1 * XM + LX2 * Y2);
	poly_drilling(1, 2) = +XP * (LX3 * XM - LX2 * Y2);
	poly_drilling(1, 3) = -XM * (LX3 * XP + LX4 * Y2);
	poly_drilling(1, 4) = -XM * (LY1 * XP - LY4 * Y2);
	poly_drilling(1, 5) = +XP * (LY1 * XM + LY2 * Y2);
	poly_drilling(1, 6) = +XP * (LY3 * XM - LY2 * Y2);
	poly_drilling(1, 7) = -XM * (LY3 * XP + LY4 * Y2);

	poly_drilling /= 16.;

	return poly_drilling;
}

mat GCMQ::form_displacement(const mat& pn_pxy, const mat& pnt_pxy) {
	mat poly_disp(3, m_size, fill::zeros);

	for(unsigned J = 0, K = 0, L = 1, M = 2, N = 4; J < m_node; ++J, K += m_dof, L += m_dof, M += m_dof, ++N) {
		poly_disp(0, K) = poly_disp(2, L) = pn_pxy(0, J);
		poly_disp(2, K) = poly_disp(1, L) = pn_pxy(1, J);
		poly_disp(0, M) = pnt_pxy(0, J);
		poly_disp(1, M) = pnt_pxy(1, N);
		poly_disp(2, M) = pnt_pxy(0, N) + pnt_pxy(1, J);
	}

	return poly_disp;
}

mat GCMQ::form_enhanced_strain(const vec& coor, const int num_enhanced_mode) {
	mat poly_enhanced_strain(3, num_enhanced_mode, fill::zeros);

	auto &X = coor(0), &Y = coor(1);

	if(1 == num_enhanced_mode) {
		poly_enhanced_strain(0, 0) = 3. * X * X - 1.;
		poly_enhanced_strain(1, 0) = 3. * Y * Y - 1.;
	} else if(2 == num_enhanced_mode) {
		poly_enhanced_strain(2, 1) = poly_enhanced_strain(0, 0) = 3. * X * X - 1.;
		poly_enhanced_strain(2, 0) = poly_enhanced_strain(1, 1) = 3. * Y * Y - 1.;
	} else throw invalid_argument("not supported");

	return poly_enhanced_strain;
}

vec GCMQ::form_stress_mode(const double X, const double Y) { return vec{0., X, Y, X * Y}; }

GCMQ::GCMQ(const unsigned T, uvec&& NN, const unsigned MAT, const double TH, const char IP)
	: MaterialElement2D(T, m_node, m_dof, std::forward<uvec>(NN), uvec{MAT}, false)
	, thickness(TH)
	, int_scheme(IP) {
	if(mapping.is_empty()) {
		mat t_mapping(4, 4);
		t_mapping.fill(.25);
		t_mapping(1, 0) = t_mapping(1, 3) = t_mapping(2, 0) = t_mapping(2, 1) = t_mapping(3, 1) = t_mapping(3, 3) = -.25;
		access::rw(mapping) = t_mapping;
	}
}

void GCMQ::initialize(const shared_ptr<DomainBase>& D) {
	const auto material_proto = std::dynamic_pointer_cast<Material2D>(D->get<Material>(material_tag(0)));

	access::rw(mat_stiffness) = material_proto->get_initial_stiffness();

	if(PlaneType::E == material_proto->plane_type) suanpan::hacker(thickness) = 1.;

	const auto ele_coor = get_coordinate(2);

	access::rw(characteristic_length) = sqrt(area::shoelace(ele_coor));

	vec diff_coor(8);
	diff_coor(0) = ele_coor(1, 1) - ele_coor(0, 1);
	diff_coor(1) = ele_coor(2, 1) - ele_coor(1, 1);
	diff_coor(2) = ele_coor(3, 1) - ele_coor(2, 1);
	diff_coor(3) = ele_coor(0, 1) - ele_coor(3, 1);
	diff_coor(4) = ele_coor(0, 0) - ele_coor(1, 0);
	diff_coor(5) = ele_coor(1, 0) - ele_coor(2, 0);
	diff_coor(6) = ele_coor(2, 0) - ele_coor(3, 0);
	diff_coor(7) = ele_coor(3, 0) - ele_coor(0, 0);

	access::rw(iso_mapping) = trans(mapping * ele_coor);

	const auto jacob_trans = form_transformation(shape::quad(vec{0., 0.}, 1) * ele_coor);

	const IntegrationPlan edge_plan(1, 2, IntegrationType::GAUSS);
	edge.clear();
	edge.reserve(4);
	for(unsigned I = 1; I <= 4; ++I) edge.emplace_back(I, thickness, node_ptr, edge_plan, iso_mapping);

	const IntegrationPlan plan(2, int_scheme == 'I' ? 2 : 3, int_scheme == 'I' ? IntegrationType::IRONS : int_scheme == 'L' ? IntegrationType::LOBATTO : IntegrationType::GAUSS);

	mat H(11, 11, fill::zeros), HTT(11, 11, fill::zeros);

	N.zeros(11, 12);
	M.zeros(11, enhanced_mode);

	int_pt.clear();
	int_pt.reserve(plan.n_rows);
	for(unsigned I = 0; I < plan.n_rows; ++I) {
		const auto &X = plan(I, 0), &Y = plan(I, 1);

		vec t_vec{X, Y};
		const auto pn = shape::quad(t_vec, 1);
		const mat jacob = pn * ele_coor;
		int_pt.emplace_back(std::move(t_vec), det(jacob) * plan(I, 2) * thickness, material_proto->get_copy());

		auto& c_pt = int_pt.back();

		const vec coord = iso_mapping * form_stress_mode(X, Y);

		c_pt.poly_stress = shape::stress11(coord);
		c_pt.poly_strain = solve(mat_stiffness, c_pt.poly_stress);

		M += c_pt.factor * c_pt.poly_stress.t() * jacob_trans * form_enhanced_strain(c_pt.coor, enhanced_mode);
		N += c_pt.factor * c_pt.poly_stress.t() * form_displacement(solve(jacob, pn), solve(jacob, form_drilling_displacement(c_pt.coor, diff_coor)));
		H += c_pt.factor * c_pt.poly_stress.t() * c_pt.poly_strain;
		HTT += c_pt.factor * c_pt.poly_strain.t() * mat_stiffness * c_pt.poly_strain;
	}

	HT = trans(H);

	if(!solve(NT, H, N) || !solve(MT, H, M)) {
		suanpan_error("GCMQ: fail to initialize, disable and return.\n");
		D->disable_element(get_tag());
		return;
	}

	const mat T = HTT * MT, W = NT.t() * T;

	trial_stiffness = current_stiffness = initial_stiffness = NT.t() * HTT * NT - W * (trial_viwt = current_viwt = initial_viwt = solve(mat(MT.t() * T), W.t()));

	pre_disp.zeros(m_size);

	trial_vif = current_vif.zeros(enhanced_mode);
	trial_zeta = current_zeta.zeros(enhanced_mode);
	trial_beta = current_beta.zeros(11);
	trial_alpha = current_alpha.zeros(11);
	trial_q = current_q.zeros(11);

	if(const auto t_density = material_proto->get_parameter(); t_density > 0.) {
		initial_mass.zeros(m_size, m_size);
		for(const auto& I : int_pt) {
			const auto n_int = shape::quad(I.coor, 0);
			const auto t_factor = t_density * I.factor;
			for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) for(auto K = J, P = L; K < m_node; ++K, P += m_dof) initial_mass(L, P) += t_factor * n_int(J) * n_int(K);
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
		const mat n_int = I.factor * shape::quad(I.coor, 0);
		const mat n_int_drill = I.factor * form_drilling_mass(I.coor, diff_coor);
		for(auto J = 0u, L = 0u; J < m_node; ++J, L += m_dof) {
			for(auto K = 0llu; K < 2; ++K) body_force(L + K, K) += n_int(J);
			body_force(L + 2llu, 2) += n_int_drill(J);
		}
	}
}

int GCMQ::update_status() {
	vec incre_disp = -pre_disp;
	incre_disp += pre_disp = get_incre_displacement();

	const vec incre_zeta = -trial_viwt * incre_disp - trial_vif;

	trial_zeta += incre_zeta;
	trial_beta += NT * incre_disp + MT * incre_zeta;

	vec local_stress(11, fill::zeros);
	mat local_stiffness(11, 11, fill::zeros);
	for(const auto& t_pt : int_pt) {
		if(t_pt.m_material->update_trial_status(t_pt.poly_strain * trial_beta) != SUANPAN_SUCCESS) return SUANPAN_FAIL;
		local_stress += t_pt.factor * t_pt.poly_strain.t() * t_pt.m_material->get_trial_stress();
		local_stiffness += t_pt.factor * t_pt.poly_strain.t() * t_pt.m_material->get_trial_stiffness() * t_pt.poly_strain;
	}

	const mat T = NT.t() * local_stiffness, V = MT.t() * local_stiffness * MT, W = T * MT;

	trial_alpha = solve(HT, local_stress);

	if(!solve(trial_viwt, V, W.t())) return SUANPAN_FAIL;
	if(!solve(trial_vif, V, M.t() * trial_alpha)) return SUANPAN_FAIL;

	trial_resistance = N.t() * trial_alpha - W * trial_vif;
	trial_stiffness = T * (NT - MT * trial_viwt);

	return SUANPAN_SUCCESS;
}

int GCMQ::commit_status() {
	current_zeta = trial_zeta;
	current_beta = trial_beta;
	current_alpha = trial_alpha;
	current_q = trial_q;
	current_vif = trial_vif;
	current_viwt = trial_viwt;

	pre_disp.zeros();

	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->commit_status();
	return code;
}

int GCMQ::clear_status() {
	trial_zeta = current_zeta.zeros();
	trial_beta = current_beta.zeros();
	trial_alpha = current_alpha.zeros();
	trial_q = current_q.zeros();
	trial_vif = current_vif.zeros();

	pre_disp.zeros();

	current_viwt = trial_viwt = initial_viwt;

	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->clear_status();
	return code;
}

int GCMQ::reset_status() {
	trial_zeta = current_zeta;
	trial_beta = current_beta;
	trial_alpha = current_alpha;
	trial_q = current_q;
	trial_vif = current_vif;
	trial_viwt = current_viwt;

	pre_disp.zeros();

	auto code = 0;
	for(const auto& I : int_pt) code += I.m_material->reset_status();
	return code;
}

vector<vec> GCMQ::record(const OutputType T) {
	vector<vec> data;

	if(T == OutputType::S) for(const auto& I : int_pt) data.emplace_back(I.poly_stress * current_alpha);
	else if(T == OutputType::S11) for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(0) * current_alpha);
	else if(T == OutputType::S22) for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(1) * current_alpha);
	else if(T == OutputType::S12) for(const auto& I : int_pt) data.emplace_back(I.poly_stress.row(2) * current_alpha);
	else if(T == OutputType::SP) for(const auto& I : int_pt) data.emplace_back(transform::stress::principal(I.poly_stress * current_alpha));
	else if(T == OutputType::SP1) for(const auto& I : int_pt) data.emplace_back(vec{transform::stress::principal(I.poly_stress * current_alpha).at(0)});
	else if(T == OutputType::SP2) for(const auto& I : int_pt) data.emplace_back(vec{transform::stress::principal(I.poly_stress * current_alpha).at(1)});
	else if(T == OutputType::MISES)
		for(const auto& I : int_pt) {
			const vec t_stress = I.poly_stress * current_alpha;
			data.emplace_back(vec{sqrt(t_stress(0) * t_stress(0) - t_stress(0) * t_stress(1) + t_stress(1) * t_stress(1) + 3. * t_stress(2) * t_stress(2))});
		}
	else if(T == OutputType::SINT) data.emplace_back(current_alpha);
	else if(T == OutputType::E) for(const auto& I : int_pt) data.emplace_back(I.poly_strain * current_beta);
	else if(T == OutputType::E11) for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(0) * current_beta);
	else if(T == OutputType::E22) for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(1) * current_beta);
	else if(T == OutputType::E12) for(const auto& I : int_pt) data.emplace_back(I.poly_strain.row(2) * current_beta);
	else if(T == OutputType::EP) for(const auto& I : int_pt) data.emplace_back(transform::strain::principal(I.poly_strain * current_beta));
	else if(T == OutputType::EP1) for(const auto& I : int_pt) data.emplace_back(vec{transform::strain::principal(I.poly_strain * current_beta).at(0)});
	else if(T == OutputType::EP2) for(const auto& I : int_pt) data.emplace_back(vec{transform::strain::principal(I.poly_strain * current_beta).at(1)});
	else if(T == OutputType::EINT) data.emplace_back(current_beta);
	else if(T == OutputType::PE) for(const auto& I : int_pt) data.emplace_back(I.poly_strain * current_beta - solve(mat_stiffness, I.poly_stress * current_alpha));
	else if(T == OutputType::PEP) for(const auto& I : int_pt) data.emplace_back(transform::strain::principal(I.poly_strain * current_beta - solve(mat_stiffness, I.poly_stress * current_alpha)));
	else if(T == OutputType::RESULTANT) for(const auto& I : edge) data.emplace_back(vec{I.F(current_alpha), I.V(current_alpha), I.M(current_alpha)});
	else if(T == OutputType::AXIAL) data.emplace_back(vec{edge[0].F(current_alpha), edge[1].F(current_alpha), edge[2].F(current_alpha), edge[3].F(current_alpha)});
	else if(T == OutputType::SHEAR) data.emplace_back(vec{edge[0].V(current_alpha), edge[1].V(current_alpha), edge[2].V(current_alpha), edge[3].V(current_alpha)});
	else if(T == OutputType::MOMENT) data.emplace_back(vec{edge[0].M(current_alpha), edge[1].M(current_alpha), edge[2].M(current_alpha), edge[3].M(current_alpha)});
	else if(T == OutputType::K) data.emplace_back(vectorise(current_stiffness));
	else if(T == OutputType::M) data.emplace_back(vectorise(current_mass));
	else for(const auto& I : int_pt) for(const auto& J : I.m_material->record(T)) data.emplace_back(J);

	return data;
}

void GCMQ::print() {
	suanpan_info("GCMQ mixed quad element %u connects nodes:\n", get_tag());
	node_encoding.t().print();
	if(!is_initialized()) return;
	suanpan_info("Material model response:\n");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("Integration Point %lu:\t", I + 1);
		int_pt[I].coor.t().print();
		int_pt[I].m_material->print();
	}
	suanpan_info("Element model response:\n");
	for(size_t I = 0; I < int_pt.size(); ++I) {
		suanpan_info("Integration Point %lu:\t", I + 1);
		int_pt[I].coor.t().print();
		suanpan_info("strain: ");
		(int_pt[I].poly_strain * current_beta).t().print();
		suanpan_info("stress: ");
		(int_pt[I].poly_stress * current_alpha).t().print();
	}
}

#ifdef SUANPAN_VTK
#include <vtkQuad.h>

void GCMQ::Setup() {
	vtk_cell = vtkSmartPointer<vtkQuad>::New();
	const auto ele_coor = get_coordinate(2);
	for(unsigned I = 0; I < m_node; ++I) {
		vtk_cell->GetPointIds()->SetId(I, node_encoding(I));
		vtk_cell->GetPoints()->SetPoint(I, ele_coor(I, 0), ele_coor(I, 1), 0.);
	}
}

void GCMQ::GetData(vtkSmartPointer<vtkDoubleArray>& arrays, const OutputType type) {
	mat t_disp(6, m_node, fill::zeros);

	if(OutputType::A == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_acceleration(), m_dof, m_node);
	else if(OutputType::V == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_velocity(), m_dof, m_node);
	else if(OutputType::U == type) t_disp.rows(uvec{0, 1, 5}) = reshape(get_current_displacement(), m_dof, m_node);

	for(unsigned I = 0; I < m_node; ++I) arrays->SetTuple(node_encoding(I), t_disp.colptr(I));
}

mat GCMQ::GetData(const OutputType P) {
	if(OutputType::S == P) {
		mat t_stress(6, m_node, fill::zeros);
		t_stress(uvec{0, 1, 3}, uvec{0}) = shape::stress11(iso_mapping * form_stress_mode(-1., -1.)) * current_alpha;
		t_stress(uvec{0, 1, 3}, uvec{1}) = shape::stress11(iso_mapping * form_stress_mode(1., -1.)) * current_alpha;
		t_stress(uvec{0, 1, 3}, uvec{2}) = shape::stress11(iso_mapping * form_stress_mode(1., 1.)) * current_alpha;
		t_stress(uvec{0, 1, 3}, uvec{3}) = shape::stress11(iso_mapping * form_stress_mode(-1., 1.)) * current_alpha;
		return t_stress;
	}
	if(OutputType::E == P) {
		mat t_strain(6, m_node, fill::zeros);
		t_strain(uvec{0, 1, 3}, uvec{0}) = solve(mat_stiffness, shape::stress11(iso_mapping * form_stress_mode(-1., -1.)) * current_beta);
		t_strain(uvec{0, 1, 3}, uvec{1}) = solve(mat_stiffness, shape::stress11(iso_mapping * form_stress_mode(1., -1.)) * current_beta);
		t_strain(uvec{0, 1, 3}, uvec{2}) = solve(mat_stiffness, shape::stress11(iso_mapping * form_stress_mode(1., 1.)) * current_beta);
		t_strain(uvec{0, 1, 3}, uvec{3}) = solve(mat_stiffness, shape::stress11(iso_mapping * form_stress_mode(-1., 1.)) * current_beta);
		return t_strain;
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

void GCMQ::SetDeformation(vtkSmartPointer<vtkPoints>& nodes, const double amplifier) {
	const mat ele_disp = get_coordinate(2) + amplifier * mat(reshape(get_current_displacement(), m_dof, m_node)).rows(0, 1).t();
	for(unsigned I = 0; I < m_node; ++I) nodes->SetPoint(node_encoding(I), ele_disp(I, 0), ele_disp(I, 1), 0.);
}

#endif
