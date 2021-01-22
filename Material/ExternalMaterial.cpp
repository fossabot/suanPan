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

#include "ExternalMaterial.h"

MaterialType ExternalMaterial::get_type(const ExternalMaterialData& D) {
	if(D.size == 6) return MaterialType::D3;
	if(D.size == 3) return MaterialType::D2;
	if(D.size == 1) return MaterialType::D1;
	return MaterialType::D0;
}

ExternalMaterial::ExternalMaterial(const unsigned T, vector<double>&& P, void* H)
	: Material(T, MaterialType::D0, 0.)
	, constant(std::forward<vector<double>>(P))
	// ReSharper disable once CppFunctionalStyleCast
	, cooker(Interface(H)) {
	data.constant = constant.data();
	data.constant_size = static_cast<unsigned>(constant.size());

	auto info = 0;

	cooker(&data, &info);

	access::rw(density) = data.density;
	access::rw(material_type) = get_type(data);
}

ExternalMaterial::ExternalMaterial(const ExternalMaterial& old_obj)
	: Material(old_obj)
	, constant(old_obj.constant)
	, cooker(old_obj.cooker) {
	data.constant = constant.data();
	data.constant_size = static_cast<unsigned>(constant.size());

	auto info = 0;

	cooker(&data, &info);

	// need to reinitialize to setup containers
	initialize(nullptr);
}

ExternalMaterial::~ExternalMaterial() {
	auto info = 1;

	cooker(&data, &info);
}

bool ExternalMaterial::validate() {
	auto info = 7;

	cooker(&data, &info);

	return 0 == info;
}

void ExternalMaterial::initialize(const shared_ptr<DomainBase>&) {
	PureWrapper(this);

	// ! very tricky implementation
	// ! we need to make all variables wrappers of some external memory
	// ! armadillo has no existing method to do so
	// ! can only make a wrapper and steal everything from it to ensure compatibility

	if(-1 != data.c_strain) {
		vec t_holder(&data.pool[data.c_strain], data.size, false);
		current_strain.steal_mem(t_holder);
	} // current status
	if(-1 != data.c_strain_rate) {
		vec t_holder(&data.pool[data.c_strain_rate], data.size, false);
		current_strain_rate.steal_mem(t_holder);
	} // current status
	if(-1 != data.c_stress) {
		vec t_holder(&data.pool[data.c_stress], data.size, false);
		current_stress.steal_mem(t_holder);
	} // current status 

	if(-1 != data.t_strain) {
		vec t_holder(&data.pool[data.t_strain], data.size, false);
		trial_strain.steal_mem(t_holder);
	} // trial status
	if(-1 != data.t_strain_rate) {
		vec t_holder(&data.pool[data.t_strain_rate], data.size, false);
		trial_strain_rate.steal_mem(t_holder);
	} // trial status
	if(-1 != data.t_stress) {
		vec t_holder(&data.pool[data.t_stress], data.size, false);
		trial_stress.steal_mem(t_holder);
	} // trial status 

	if(-1 != data.i_stiffness) {
		mat t_holder(&data.pool[data.i_stiffness], data.size, data.size, false);
		initial_stiffness.steal_mem(t_holder);
	} // stiffness matrix
	if(-1 != data.c_stiffness) {
		mat t_holder(&data.pool[data.c_stiffness], data.size, data.size, false);
		current_stiffness.steal_mem(t_holder);
	} // stiffness matrix
	if(-1 != data.t_stiffness) {
		mat t_holder(&data.pool[data.t_stiffness], data.size, data.size, false);
		trial_stiffness.steal_mem(t_holder);
	} // stiffness matrix

	if(-1 != data.i_damping) {
		mat t_holder(&data.pool[data.i_damping], data.size, data.size, false);
		initial_damping.steal_mem(t_holder);
	} // damping matrix
	if(-1 != data.c_damping) {
		mat t_holder(&data.pool[data.c_damping], data.size, data.size, false);
		current_damping.steal_mem(t_holder);
	} // damping matrix
	if(-1 != data.t_damping) {
		mat t_holder(&data.pool[data.t_damping], data.size, data.size, false);
		trial_damping.steal_mem(t_holder);
	} // damping matrix
}

void ExternalMaterial::initialize_history(unsigned) {}

void ExternalMaterial::set_initial_history(const vec&) {}

unique_ptr<Material> ExternalMaterial::get_copy() { return make_unique<ExternalMaterial>(*this); }

int ExternalMaterial::update_trial_status(const vec& t_strain) {
	if(-1 == data.t_strain) return SUANPAN_FAIL;

	trial_strain = t_strain;

	auto info = 2;

	cooker(&data, &info);

	return info;
}

int ExternalMaterial::update_trial_status(const vec& t_strain, const vec& t_strain_rate) {
	if(-1 == data.t_strain || -1 == data.t_strain_rate) return SUANPAN_FAIL;

	trial_strain = t_strain;
	trial_strain_rate = t_strain_rate;

	auto info = 3;

	cooker(&data, &info);

	return info;
}

int ExternalMaterial::commit_status() {
	auto info = 4;

	cooker(&data, &info);

	return info;
}

int ExternalMaterial::reset_status() {
	auto info = 5;

	cooker(&data, &info);

	return info;
}

int ExternalMaterial::clear_status() {
	auto info = 6;

	cooker(&data, &info);

	return info;
}

vector<vec> ExternalMaterial::record(OutputType) { return {}; }
