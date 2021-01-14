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

#include <suanPan.h>
#include <Material/ExternalMaterialData.h>

void allocate_material(ExternalMaterialData* data, int* info) {
	if(nullptr == data) {
		*info = -1;
		return;
	}

	data->size = 1;

	data->c_strain = 0;
	data->t_strain = 1;
	data->c_stress = 2;
	data->t_stress = 3;
	data->i_stiffness = 4;
	data->c_stiffness = 5;
	data->t_stiffness = 6;

	data->pool = new double[7];

	if(nullptr == data) *info = -1;

	for(auto I = 0; I < 7; ++I) data->pool[I] = 0.;

	data->pool[data->c_stiffness] = data->pool[data->t_stiffness] = data->pool[data->i_stiffness] = data->constant[0];

	if(data->constant_size > 1) data->density = data->constant[1];
}

void deallocate_material(ExternalMaterialData* data, int* info) {
	if(nullptr == data) {
		*info = -1;
		return;
	}

	delete[] data->pool;

	*info = 0;
}

void static_update(ExternalMaterialData* data, int* info) {
	data->pool[data->t_stress] = data->pool[data->t_stiffness] * data->pool[data->t_strain];

	*info = 0;
}

void dynamic_update(ExternalMaterialData* data, int* info) { static_update(data, info); }

void commit(ExternalMaterialData* data, int* info) {
	data->pool[data->c_strain] = data->pool[data->t_strain];
	data->pool[data->c_stress] = data->pool[data->t_stress];

	*info = 0;
}

void reset(ExternalMaterialData* data, int* info) {
	data->pool[data->t_strain] = data->pool[data->c_strain];
	data->pool[data->t_stress] = data->pool[data->c_stress];

	*info = 0;
}

void clear(ExternalMaterialData* data, int* info) {
	data->pool[data->c_strain] = data->pool[data->t_strain] = 0.;
	data->pool[data->c_stress] = data->pool[data->t_stress] = 0.;

	*info = 0;
}

void validate(ExternalMaterialData* data, int* info) { *info = 0 == data->constant_size ? -1 : 0; }

SUANPAN_EXPORT void elasticexternal_handler(ExternalMaterialData* data, int* info) {
	if(0 == *info) allocate_material(data, info);
	else if(1 == *info) deallocate_material(data, info);
	else if(2 == *info) static_update(data, info);
	else if(3 == *info) dynamic_update(data, info);
	else if(4 == *info) commit(data, info);
	else if(5 == *info) reset(data, info);
	else if(6 == *info) clear(data, info);
	else if(7 == *info) validate(data, info);
	else {
		suanpan_error("unknown flag recieved.\n");
		*info = -1;
	}
}
