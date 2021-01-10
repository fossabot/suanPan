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
 * @class FullMatMAGMA
 * @brief A FullMatMAGMA class that holds matrices.
 *
 * @author tlc
 * @date 06/09/2017
 * @version 0.1.0
 * @file FullMatMAGMA.hpp
 * @addtogroup MetaMat
 * @{
 */

// ReSharper disable CppCStyleCast
#ifndef FULLMATMAGMA_HPP
#define FULLMATMAGMA_HPP

#ifdef SUANPAN_MAGMA

template<typename T> class FullMatMAGMA final : public FullMat<T> {
public:
	using FullMat<T>::IPIV;
	using FullMat<T>::n_rows;
	using FullMat<T>::FullMat;

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;

	int solve_trs(Mat<T>&, const Mat<T>&) override;
};

template<typename T> unique_ptr<MetaMat<T>> FullMatMAGMA<T>::make_copy() { return make_unique<FullMatMAGMA<T>>(*this); }

template<typename T> int FullMatMAGMA<T>::solve(Mat<T>& X, const Mat<T>& B) {
	X.set_size(B.n_rows, B.n_cols);

	auto INFO = 0;

	const auto N = static_cast<magma_int_t>(n_rows);
	const auto NRHS = static_cast<magma_int_t>(B.n_cols);
	const auto LDA = N;
	const auto LDB = static_cast<magma_int_t>(B.n_rows);
	IPIV.zeros(N);

	magma_queue_t QUEUE;
	magma_queue_create(0, &QUEUE);
	const auto LDDA = magma_roundup(N, 64);

	if(std::is_same<T, float>::value) {
		using E = float;
		magmaFloat_ptr GPUA, GPUB;
		magma_smalloc(&GPUA, LDDA * N);
		magma_smalloc(&GPUB, LDDA * NRHS);
		magma_ssetmatrix(N, N, (E*)this->memptr(), LDA, GPUA, LDDA, QUEUE);
		magma_ssetmatrix(N, NRHS, (E*)B.memptr(), LDB, GPUB, LDDA, QUEUE);
		magma_sgesv_gpu(N, NRHS, GPUA, LDDA, IPIV.memptr(), GPUB, LDDA, &INFO);
		magma_sgetmatrix(N, NRHS, GPUB, LDDA, (E*)X.memptr(), LDB, QUEUE);
		magma_free(GPUA);
		magma_free(GPUB);
	} else if(std::is_same<T, double>::value) {
		using E = double;
		magmaDouble_ptr GPUA, GPUB;
		magma_dmalloc(&GPUA, LDDA * N);
		magma_dmalloc(&GPUB, LDDA * NRHS);
		magma_dsetmatrix(N, N, (E*)this->memptr(), LDA, GPUA, LDDA, QUEUE);
		magma_dsetmatrix(N, NRHS, (E*)B.memptr(), LDB, GPUB, LDDA, QUEUE);
		magma_dgesv_gpu(N, NRHS, GPUA, LDDA, IPIV.memptr(), GPUB, LDDA, &INFO);
		magma_dgetmatrix(N, NRHS, GPUB, LDDA, (E*)X.memptr(), LDB, QUEUE);
		magma_free(GPUA);
		magma_free(GPUB);
	}

	if(INFO != 0) suanpan_error("solve() receives error code %u from the base driver, the matrix is probably singular.\n", INFO);

	return INFO;
}

template<typename T> int FullMatMAGMA<T>::solve_trs(Mat<T>& X, const Mat<T>& B) { return solve(X, B); }

#endif

#endif

//! @}
