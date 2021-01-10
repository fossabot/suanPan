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
 * @class SRGCMQ
 * @brief A SRGCMQ class.
 * 
 * Reference:
 *   1. A new drilling quadrilateral membrane element with high coarse-mesh accuracy using a modified Hu-Washizu principle
 *      https://doi.org/10.1002/nme.6066
 *
 * @author tlc
 * @date 07/03/2019
 * @version 0.2.0
 * @file SRGCMQ.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef SRGCMQ_H
#define SRGCMQ_H

#include <Element/MaterialElement.h>

class IntegrationPlan;

class SRGCMQ final : public MaterialElement2D {
	struct IntegrationPoint final {
		vec coor;
		const double factor;
		unique_ptr<Material> m_material;
		unique_ptr<Material> rebar_x = nullptr, rebar_y = nullptr;
		mat poly_strain;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
	};

	static constexpr unsigned m_node = 4, m_dof = 3, m_size = m_dof * m_node;

	static const mat iso_mapping;

	const double thickness;

	const char int_scheme;

	const double rho_x, rho_y;
	const unsigned mat_x, mat_y;

	const bool reinforced_x = rho_x != 0. && mat_x != 0;
	const bool reinforced_y = rho_y != 0. && mat_y != 0;

	vector<IntegrationPoint> int_pt;

	static mat form_transformation(const mat&);
	static mat form_drilling_mass(const vec&, const vec&);
	static mat form_drilling_displacement(const vec&, const vec&);
	static mat form_displacement(const mat&, const mat&);
	static vec form_stress_mode(double, double);
public:
	SRGCMQ(unsigned,     // element tag
	       uvec&&,       // node tag
	       unsigned,     // material tag
	       double = 1.,  // thickness
	       char = 'I',   // integration type
	       double = 0.,  // reinforcement ratio along x-axis
	       unsigned = 0, // material tag of x-reinforcement
	       double = 0.,  // reinforcement ratio along y-axis
	       unsigned = 0  // material tag of y-reinforcement
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
