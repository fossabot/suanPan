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
 * @class Rebar2D
 * @brief A Rebar2D material class.
 * @author tlc
 * @date 03/10/2017
 * @version 0.1.1
 * @file Rebar2D.h
 * @addtogroup Material-2D
 * @{
 */

#ifndef REBAR2D_H
#define REBAR2D_H

#include <Material/Material2D/Material2D.h>

class Rebar2D final : public Material2D {
	const unsigned tag_major, tag_minor;

	const double ratio_major, ratio_minor;

	const double inclination;

	const mat trans_mat;

	unique_ptr<Material> rebar_major, rebar_minor;
public:
	Rebar2D(unsigned,   // tag
	        unsigned,   // material tag along major axis
	        unsigned,   // material tag along minor axis
	        double,     // reinforcement ratio along major axis
	        double,     // reinforcement ratio along minor axis
	        double = 0. // inclination
	);
	Rebar2D(const Rebar2D&);

	void initialize(const shared_ptr<DomainBase>&) override;

	unique_ptr<Material> get_copy() override;
	int update_trial_status(const vec&) override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	void print() override;
};

#endif

//! @}
