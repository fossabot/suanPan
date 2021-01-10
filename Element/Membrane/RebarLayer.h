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
 * @class RebarLayer
 * @brief The RebarLayer class.
 *
 * A three--segment RebarLayer element based on the GQ12 element.
 *
 * @author tlc
 * @date 26/01/2018
 * @version 0.1.2
 * @file RebarLayer.h
 * @addtogroup Membrane
 * @ingroup Element
 * @{
 */

#ifndef REBARLAYER_H
#define REBARLAYER_H

#include <Element/MaterialElement.h>

class RebarLayer final : public MaterialElement1D {
	struct IntegrationPoint final {
		vec coor;
		double weight;
		unique_ptr<Material> m_material;
		vec strain_mat;
		IntegrationPoint(vec&&, double, unique_ptr<Material>&&);
	};

	static constexpr unsigned m_node = 4, m_dof = 3, m_size = m_dof * m_node;

	const double thickness;

	const char direction;

	const double length_a, length_b, length_c;
	const double rho_a, rho_b, rho_c;

	vector<IntegrationPoint> int_pt;
public:
	RebarLayer(unsigned, uvec&&, vec&&, double = 1., char = 'Y');

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
