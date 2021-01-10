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
 * @class Cosine
 * @brief An Amplitude class that can generate Amplitude pattern.
 *
 * @author tlc
 * @date 26/10/2017
 * @version 0.1.0
 * @file Cosine.h
 * @addtogroup Amplitude
 * @{
 */

#ifndef COSINE_H
#define COSINE_H

#include <Load/Amplitude/Amplitude.h>

class Cosine final : public Amplitude {
	const double period;
	const std::vector<double> amp;
public:
	explicit Cosine(unsigned, double, std::vector<double>&&, unsigned);

	double get_amplitude(double) override;

	void print() override;
};

#endif

//! @}
