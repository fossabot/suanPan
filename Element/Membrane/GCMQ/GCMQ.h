/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
/**
 * @class GCMQ
 * @brief A GCMQ class.
 * 
 * Reference:
 *   1. A new drilling quadrilateral membrane element with high coarse-mesh accuracy using a modified Hu-Washizu principle
 *      https://doi.org/10.1002/nme.6066
 *
 * @author tlc
 * @date 07/03/2019
 * @version 0.6.0
 * @file GCMQ.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef GCMQ_H
#define GCMQ_H

#include <Element/MaterialElement.h>

class IntegrationPlan;

class GCMQ final : public MaterialElement2D {
	class ResultantConverter final {
		mat converter_a, converter_b;
		mat direction_cosine;
	public:
		ResultantConverter(unsigned, double, const vector<weak_ptr<Node>>&, const IntegrationPlan&, const mat&);
		[[nodiscard]] double F(const vec&) const;
		[[nodiscard]] double V(const vec&) const;
		[[nodiscard]] double M(const vec&) const;
	};

	struct IntegrationPoint final {
		vec coor;
		const double factor;
		unique_ptr<Material> m_material;
		mat poly_stress, poly_strain;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
	};

	static constexpr unsigned m_node = 4, m_dof = 3, m_size = m_dof * m_node;

	static constexpr int enhanced_mode = 2;

	static const mat mapping;

	const double thickness;

	const char int_scheme;

	vector<IntegrationPoint> int_pt;
	vector<ResultantConverter> edge;

	const mat mat_stiffness, iso_mapping;

	mat HT, NT, MT, N, M; // constant matrices

	mat initial_viwt, trial_viwt, current_viwt;
	vec trial_vif, current_vif;
	vec trial_zeta, current_zeta;   // enhanced strain
	vec trial_beta, current_beta;   // strain
	vec trial_alpha, current_alpha; // stress
	vec trial_q, current_q;

	vec pre_disp;

	static mat form_transformation(const mat&);
	static mat form_drilling_mass(const vec&, const vec&);
	static mat form_drilling_displacement(const vec&, const vec&);
	static mat form_displacement(const mat&, const mat&);
	static mat form_enhanced_strain(const vec&, int);
	static vec form_stress_mode(double, double);
public:
	GCMQ(unsigned,    // element tag
	     uvec&&,      // node tag
	     unsigned,    // material tag
	     double = 1., // thickness
	     char = 'I'   // integration type
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int commit_status() override;
	int clear_status() override;
	int reset_status() override;

	vector<vec> record(OutputType) override;

	void print() override;

#ifdef SUANPAN_VTK
	void Setup() override;
	void GetData(vtkSmartPointer<vtkDoubleArray>&, OutputType) override;
	mat GetData(OutputType) override;
	void SetDeformation(vtkSmartPointer<vtkPoints>&, double) override;
#endif
};

#endif

//! @}
