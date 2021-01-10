////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2017-2021 Theodore Chang
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
////////////////////////////////////////////////////////////////////////////////

#include "VisualisationRecorder.h"
#include <Domain/DomainBase.h>

VisualisationRecorder::VisualisationRecorder(const unsigned T, const OutputType L, const unsigned I, const unsigned W, const double S)
	: Recorder(T, {}, L, I, false, false)
	, width(W) {
#ifdef SUANPAN_VTK
	config.save_file = true;
	config.type = get_variable_type();
	config.scale = S;

	if(const auto t_name = vtk_get_name(config.type); 'U' == t_name[0] || 'V' == t_name[0] || 'A' == t_name[0]) funtion_handler = &vtk_plot_node_quantity;
	else funtion_handler = &vtk_plot_element_quantity;
#endif
}

void VisualisationRecorder::record(const shared_ptr<DomainBase>& D) {
#ifdef SUANPAN_VTK
	if(1 != interval && counter++ != interval) return;

	counter = 1;

	ostringstream file_name;
	file_name << to_char(get_variable_type()) << '-' << std::setw(width) << std::setfill('0') << ++total_counter << ".vtk";

	config.file_name = file_name.str();

	(*funtion_handler)(D, config);
#endif
}

void VisualisationRecorder::save() {}

void VisualisationRecorder::print() { suanpan_info("A Visualisation Recorder.\n"); }
