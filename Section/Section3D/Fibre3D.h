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
 * @class Fibre3D
 * @brief A Fibre3D class.
 * @author tlc
 * @date 10/06/2018
 * @version 0.1.0
 * @file Fibre3D.h
 * @addtogroup Section-2D
 * @ingroup Section
 * @{
 */

#ifndef Fibre3D_H
#define Fibre3D_H

#include <Section/Section3D/Section3D.h>

class Fibre3D final : public Section3D {
	podarray<unsigned> fibre_tag;

	vector<unique_ptr<Section>> fibre;
public:
	explicit Fibre3D(unsigned, podarray<unsigned>&&);
	Fibre3D(const Fibre3D&);
	Fibre3D(Fibre3D&&) noexcept = delete;
	Fibre3D& operator=(const Fibre3D&) = delete;
	Fibre3D& operator=(Fibre3D&&) noexcept = delete;
	~Fibre3D() override = default;

	void initialize(const shared_ptr<DomainBase>&) override;

	unique_ptr<Section> get_copy() override;

	int update_trial_status(const vec&) override;

	int clear_status() override;
	int commit_status() override;
	int reset_status() override;

	void print() override;
};

#endif

//! @}
