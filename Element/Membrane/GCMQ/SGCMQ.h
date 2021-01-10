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
 * @class SGCMQ
 * @brief A SGCMQ class.
 * 
 * Reference:
 *   1. A new drilling quadrilateral membrane element with high coarse-mesh accuracy using a modified Hu-Washizu principle
 *      https://doi.org/10.1002/nme.6066
 *
 * @author tlc
 * @date 07/03/2019
 * @version 0.6.0
 * @file SGCMQ.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef SGCMQ_H
#define SGCMQ_H

#include <Element/MaterialElement.h>

class IntegrationPlan;

class SGCMQ final : public MaterialElement2D {
	struct IntegrationPoint final {
		vec coor;
		const double factor;
		unique_ptr<Material> m_material;
		mat poly_strain;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
	};

	static constexpr unsigned m_node = 4, m_dof = 3, m_size = m_dof * m_node;

	static const mat mapping;

	const double thickness;

	const char int_scheme;

	vector<IntegrationPoint> int_pt;

	static mat form_drilling_mass(const vec&, const vec&);
	static mat form_drilling_displacement(const vec&, const vec&);
	static mat form_displacement(const mat&, const mat&);
	static vec form_stress_mode(double, double);
public:
	SGCMQ(unsigned,    // element tag
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
