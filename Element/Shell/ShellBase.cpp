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

#include "ShellBase.h"
#include <Domain/Node.h>

const uvec ShellBase::m_dof{0, 1, 5};
const uvec ShellBase::p_dof{2, 3, 4};

vec ShellBase::reshuffle(const vec& membrane_resistance, const vec& plate_resistance) {
	suanpan_debug([&]() { if(membrane_resistance.n_elem != plate_resistance.n_elem) throw invalid_argument("size conflicts"); });

	const auto t_size = 2 * plate_resistance.n_elem;

	vec total_resistance(t_size);

	for(unsigned I = 0, J = 0; I < t_size; I += 6, J += 3) {
		const span K(J, J + 2llu);
		total_resistance(I + m_dof) = membrane_resistance(K);
		total_resistance(I + p_dof) = plate_resistance(K);
	}

	return total_resistance;
}

mat ShellBase::reshuffle(const mat& membrane_stiffness, const mat& plate_stiffness) {
	suanpan_debug([&]() {
		if(membrane_stiffness.n_cols != membrane_stiffness.n_rows) throw invalid_argument("size conflicts");
		if(membrane_stiffness.n_cols != plate_stiffness.n_cols) throw invalid_argument("size conflicts");
		if(plate_stiffness.n_cols != plate_stiffness.n_rows) throw invalid_argument("size conflicts");
	});

	const auto t_size = 2 * plate_stiffness.n_cols;

	mat total_stiffness(t_size, t_size, fill::zeros);

	for(unsigned I = 0, K = 0; I < t_size; I += 6, K += 3) {
		const span M(K, K + 2llu);
		for(unsigned J = 0, L = 0; J < t_size; J += 6, L += 3) {
			const span N(L, L + 2llu);
			total_stiffness(I + m_dof, J + m_dof) = membrane_stiffness(M, N);
			total_stiffness(I + p_dof, J + p_dof) = plate_stiffness(M, N);
		}
	}

	return total_stiffness;
}

void ShellBase::direction_cosine() {
	const mat coor = get_coordinate(3).t();

	trans_mat.set_size(3, 3);

	trans_mat.col(0) = normalise(coor.col(1) - coor.col(0));
	trans_mat.col(2) = normalise(cross(trans_mat.col(0), coor.col(2) - coor.col(0)));
	trans_mat.col(1) = cross(trans_mat.col(2), trans_mat.col(0));
}

mat ShellBase::get_local_coordinate() const {
	mat l_coordinate = get_coordinate(3).t();
	l_coordinate = trans_mat.t() * (l_coordinate - repmat(l_coordinate.col(0), 1, l_coordinate.n_cols));
	if(norm(l_coordinate.row(2)) > 1E-6) suanpan_warning("non-planar shell geomtry detected.\n");
	return l_coordinate.rows(0, 1).t();
}

vec& ShellBase::transform_from_local_to_global(vec& resistance) const {
	for(unsigned I = 0; I < resistance.n_elem; I += 3) {
		const span t_span(I, I + 2llu);
		resistance(t_span) = trans_mat * resistance(t_span);
	}

	return resistance;
}

vec& ShellBase::transform_from_global_to_local(vec& displacement) const {
	for(unsigned I = 0; I < displacement.n_elem; I += 3) {
		const span t_span(I, I + 2llu);
		displacement(t_span) = trans_mat.t() * displacement(t_span);
	}

	return displacement;
}

mat& ShellBase::transform_from_local_to_global(mat& stiffness) const {
	suanpan_debug([&]() { if(stiffness.n_cols != stiffness.n_rows) throw invalid_argument("size conflicts"); });

	for(unsigned I = 0; I < stiffness.n_cols; I += 3) {
		const span i_span(I, I + 2llu);
		for(unsigned J = 0; J < stiffness.n_cols; J += 3) {
			const span j_span(J, J + 2llu);
			stiffness(i_span, j_span) = trans_mat * stiffness(i_span, j_span) * trans_mat.t();
		}
	}

	return stiffness;
}

vec ShellBase::transform_from_local_to_global(vec&& resistance) const { return std::move(transform_from_local_to_global(resistance)); }

vec ShellBase::transform_from_global_to_local(vec&& displacement) const { return std::move(transform_from_global_to_local(displacement)); }

mat ShellBase::transform_from_local_to_global(mat&& stiffness) const { return std::move(transform_from_local_to_global(stiffness)); }
