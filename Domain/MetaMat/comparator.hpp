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

#ifndef COMPARATOR_HPP
#define COMPARATOR_HPP

template<typename index_t> class sparse_comparator {
protected:
	const index_t* row_idx;
	const index_t* col_idx;
public:
	sparse_comparator(const index_t* const in_row_idx, const index_t* const in_col_idx)
		: row_idx(in_row_idx)
		, col_idx(in_col_idx) {}

	virtual ~sparse_comparator() = default;
	virtual bool operator()(index_t, index_t) const = 0;
};

template<typename index_t> class csr_comparator final : public sparse_comparator<index_t> {
public:
	using sparse_comparator<index_t>::sparse_comparator;
	using sparse_comparator<index_t>::row_idx;
	using sparse_comparator<index_t>::col_idx;

	bool operator()(const index_t idx_a, const index_t idx_b) const override {
		if(row_idx[idx_a] < row_idx[idx_b]) return true;
		if(row_idx[idx_a] > row_idx[idx_b]) return false;
		if(col_idx[idx_a] < col_idx[idx_b]) return true;
		return false;
	}
};

template<typename index_t> class csc_comparator final : public sparse_comparator<index_t> {
public:
	using sparse_comparator<index_t>::sparse_comparator;
	using sparse_comparator<index_t>::row_idx;
	using sparse_comparator<index_t>::col_idx;

	bool operator()(const index_t idx_a, const index_t idx_b) const override {
		if(col_idx[idx_a] < col_idx[idx_b]) return true;
		if(col_idx[idx_a] > col_idx[idx_b]) return false;
		if(row_idx[idx_a] < row_idx[idx_b]) return true;
		return false;
	}
};

template<typename index_t> class abs_comparator final : public sparse_comparator<index_t> {
public:
	using sparse_comparator<index_t>::row_idx;
	using sparse_comparator<index_t>::col_idx;

	explicit abs_comparator(const index_t* const location)
		: sparse_comparator<index_t>(location, nullptr) {}

	bool operator()(const index_t idx_a, const index_t idx_b) const override { return row_idx[idx_a] < row_idx[idx_b]; }
};

#endif
