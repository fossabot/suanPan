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

#ifndef DEBUG_H
#define DEBUG_H

#include <functional>

void suanpan_info(const char*, ...);
void suanpan_debug(const char*, ...);
void suanpan_extra_debug(const char*, ...);
void suanpan_warning(const char*, ...);
void suanpan_error(const char*, ...);
void suanpan_fatal(const char*, ...);
void suanpan_print(const char*, ...);

void suanpan_debug(const std::function<void()>&);

bool check_debugger();

void swap_constant();

#endif
