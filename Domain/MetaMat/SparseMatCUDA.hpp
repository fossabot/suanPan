/*******************************************************************************
 * Copyright (C) 2017-2019 Theodore Chang
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
 * @class SparseMatCUDA
 * @brief A SparseMatCUDA class that holds matrices.
 *
 * @author tlc
 * @date 14/01/2021
 * @version 0.1.0
 * @file SparseMatCUDA.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMATCUDA_HPP
#define SPARSEMATCUDA_HPP

#ifdef SUANPAN_CUDA

#include <cusolverSp.h>
#include <cusparse.h>

template<typename T> class SparseMatCUDA final : public SparseMat<T> {
	cusolverSpHandle_t handle;
	cusparseMatDescr_t descrA;
public:
	using SparseMat<T>::tolerance;
	using SparseMat<T>::triplet_mat;

	SparseMatCUDA();
	SparseMatCUDA(uword, uword, uword = 0);

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> SparseMatCUDA<T>::SparseMatCUDA()
	: SparseMat<T>() {
	cusolverSpCreate(&handle);
	cusparseCreateMatDescr(&descrA);
}

template<typename T> SparseMatCUDA<T>::SparseMatCUDA(const uword in_row, const uword in_col, const uword in_elem)
	: SparseMat<T>(in_row, in_col, in_elem) {
	cusolverSpCreate(&handle);
	cusparseCreateMatDescr(&descrA);
}

template<typename T> unique_ptr<MetaMat<T>> SparseMatCUDA<T>::make_copy() { return make_unique<SparseMatCUDA<T>>(*this); }

template<typename T> int SparseMatCUDA<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
	csr_form<T, int> csr_mat(triplet_mat);

	out_mat.set_size(size(in_mat));

	int singularity;

	const cusolverStatus_t info = cusolverSpDcsrlsvluHost(handle, csr_mat.n_rows, csr_mat.c_size, descrA, csr_mat.val_idx, csr_mat.row_ptr, csr_mat.col_idx, in_mat.memptr(), tolerance, 3, out_mat.memptr(), &singularity);

	return CUSOLVER_STATUS_SUCCESS == info ? SUANPAN_SUCCESS : SUANPAN_FAIL;
}

#endif

#endif

//! @}
