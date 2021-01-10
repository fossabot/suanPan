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
 * @class RCP4
 * @brief The RCP4 class handles RCPS4, RCPE4, RCPS4R and RCPE4R elements. It is
 * a four node constant strain membrane element with optional reduced
 * integration for both plane stress and plane strain problems and optional
 * switch for TL nonlinear geometry formulation.
 * @author tlc
 * @date 27/09/2017
 * @version 0.1.3
 * @file RCP4.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef RCP4_H
#define RCP4_H

#include <Element/MaterialElement.h>

class RCP4 final : public MaterialElement2D {
	struct IntegrationPoint final {
		vec coor;
		double weight;
		unique_ptr<Material> m_material;
		unique_ptr<Material> rebar_x, rebar_y;
		mat pn_pxy;
		mat B, BG;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&, mat&&);
	};

	static constexpr unsigned m_node = 4, m_dof = 2, m_size = m_dof * m_node;

	const double thickness;

	const bool reduced_scheme;

	const double rho_x, rho_y;
	const unsigned mat_x, mat_y;

	const bool reinforced_x = rho_x != 0. && mat_x != 0;
	const bool reinforced_y = rho_y != 0. && mat_y != 0;

	vector<IntegrationPoint> int_pt;

	mat hourglassing;

	static void stack_stiffness(mat&, const mat&, const mat&, double);
public:
	RCP4(unsigned,     // tag
	     uvec&&,       // node tag
	     unsigned,     // material tag
	     double = 1.,  // thickness
	     double = 0.,  // reinforcement ratio x
	     unsigned = 0, // material tag x
	     double = 0.,  // reinforcement ratio y
	     unsigned = 0, // material tag y
	     bool = false, // reduced integration
	     bool = false  // nonlinear geometry switch
	);

	void initialize(const shared_ptr<DomainBase>&) override;

	int update_status() override;

	int commit_status() override;
	int clear_status() override;
	int reset_status() override;

	vector<vec> record(OutputType) override;

	void print() override;
};

#endif

//! @}
