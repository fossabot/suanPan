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
 * @fn sort_color
 * @brief A four color sorting algorithm.
 * @author tlc
 * @date 22/02/2020
 * @version 0.1.2
 * @file sort_color.h
 * @addtogroup Utility
 * @{
 */

#ifndef COLOR_H
#define COLOR_H

#include <suanPan.h>

using std::vector;
using std::multimap;

unsigned sort_color(vector<vector<unsigned>>&, unsigned);

#endif

//! @}
