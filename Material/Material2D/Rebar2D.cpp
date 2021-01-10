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

#include "Rebar2D.h"
#include <Domain/DomainBase.h>
#include <Toolbox/tensorToolbox.h>

Rebar2D::Rebar2D(const unsigned T, const unsigned XT, const unsigned YT, const double RX, const double RY, const double A)
	: Material2D(T, PlaneType::S, 0.)
	, tag_major(XT)
	, tag_minor(YT)
	, ratio_major(RX)
	, ratio_minor(RY)
	, inclination(A)
	, trans_mat(transform::strain::trans(inclination)) {}

Rebar2D::Rebar2D(const Rebar2D& old_obj)
	: Material2D(old_obj)
	, tag_major(old_obj.tag_major)
	, tag_minor(old_obj.tag_minor)
	, ratio_major(old_obj.ratio_major)
	, ratio_minor(old_obj.ratio_minor)
	, inclination(old_obj.inclination)
	, trans_mat(old_obj.trans_mat)
	, rebar_major(old_obj.rebar_major == nullptr ? nullptr : old_obj.rebar_major->get_copy())
	, rebar_minor(old_obj.rebar_minor == nullptr ? nullptr : old_obj.rebar_minor->get_copy()) {}

void Rebar2D::initialize(const shared_ptr<DomainBase>& D) {
	if(nullptr == D || !D->find_material(tag_major) || !D->find_material(tag_minor)) {
		D->disable_material(get_tag());
		return;
	}

	rebar_major = get_material(D, tag_major)->get_copy();
	rebar_minor = get_material(D, tag_minor)->get_copy();

	access::rw(density) = ratio_major * rebar_major->get_parameter() + ratio_minor * rebar_minor->get_parameter();

	initial_stiffness.zeros(3, 3);

	initial_stiffness(0, 0) = ratio_major * rebar_major->get_initial_stiffness().at(0);
	initial_stiffness(1, 1) = ratio_minor * rebar_minor->get_initial_stiffness().at(0);

	trial_stiffness = current_stiffness = initial_stiffness = trans_mat.t() * diagmat(initial_stiffness) * trans_mat;
}

unique_ptr<Material> Rebar2D::get_copy() { return make_unique<Rebar2D>(*this); }

int Rebar2D::update_trial_status(const vec& t_strain) {
	trial_strain = t_strain;

	const vec main_strain = trans_mat * trial_strain;

	// update status
	rebar_major->update_trial_status(main_strain(0));
	rebar_minor->update_trial_status(main_strain(1));

	vec main_stress(3);

	// collect main stress components
	main_stress(0) = ratio_major * rebar_major->get_trial_stress().at(0);
	main_stress(1) = ratio_minor * rebar_minor->get_trial_stress().at(0);
	main_stress(2) = 0.;

	// collect principal stiffness components
	trial_stiffness(0, 0) = ratio_major * rebar_major->get_trial_stiffness().at(0);
	trial_stiffness(1, 1) = ratio_minor * rebar_minor->get_trial_stiffness().at(0);
	trial_stiffness(2, 2) = 0.;

	// transform back to nominal direction
	trial_stress = trans_mat.t() * main_stress;
	trial_stiffness = trans_mat.t() * diagmat(trial_stiffness) * trans_mat;

	return SUANPAN_SUCCESS;
}

int Rebar2D::clear_status() {
	current_strain.zeros();
	trial_strain.zeros();
	current_stress.zeros();
	trial_stress.zeros();
	trial_stiffness = current_stiffness = initial_stiffness;
	return rebar_major->clear_status() + rebar_minor->clear_status();
}

int Rebar2D::commit_status() {
	current_strain = trial_strain;
	current_stress = trial_stress;
	current_stiffness = trial_stiffness;
	return rebar_major->commit_status() + rebar_minor->commit_status();
}

int Rebar2D::reset_status() {
	trial_strain = current_strain;
	trial_stress = current_stress;
	trial_stiffness = current_stiffness;
	return rebar_major->reset_status() + rebar_minor->reset_status();
}

void Rebar2D::print() {
	suanpan_info("A rebar layer with major/minor reinforcement ratio of %.3E and %.3E.\n", ratio_major, ratio_minor);
	suanpan_info("Major: ");
	rebar_major->print();
	suanpan_info("Minor: ");
	rebar_minor->print();
}
