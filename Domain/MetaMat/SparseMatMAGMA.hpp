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
 * @class SparseMatMAGMA
 * @brief A SparseMatMAGMA class that holds matrices.
 *
 * @author tlc
 * @date 14/08/2020
 * @version 0.1.0
 * @file SparseMatMAGMA.hpp
 * @addtogroup MetaMat
 * @{
 */

#ifndef SPARSEMATMAGMA_HPP
#define SPARSEMATMAGMA_HPP

#ifdef SUANPAN_MAGMA

template<typename T> class SparseMatMAGMA final : public SparseMat<T> {
	magma_queue_t queue{};

	magma_dopts doption{};
	magma_d_matrix A{}, dA{}, b{}, db{}, x{}, dx{};

	magma_sopts soption{};
	magma_s_matrix sA{}, sdA{}, sb{}, sdb{}, sx{}, sdx{};

	void init_magma();
public:
	using SparseMat<T>::precision;
	using SparseMat<T>::tolerance;
	using SparseMat<T>::triplet_mat;

	SparseMatMAGMA();
	SparseMatMAGMA(uword, uword, uword = 0);
	~SparseMatMAGMA();

	unique_ptr<MetaMat<T>> make_copy() override;

	int solve(Mat<T>&, const Mat<T>&) override;
};

template<typename T> void SparseMatMAGMA<T>::init_magma() {
	magma_queue_create(0, &queue);

	if(std::is_same<T, float>::value) {
		soption.solver_par.solver = Magma_PIDRMERGE;
		soption.solver_par.restart = 8;
		soption.solver_par.maxiter = 500;
		soption.solver_par.rtol = 1e-7f;
		soption.solver_par.atol = 1e-7f;
		soption.solver_par.maxiter = 200;
		soption.precond_par.solver = Magma_ILU;
		soption.precond_par.levels = 0;
		soption.precond_par.trisolver = Magma_CUSOLVE;

		magma_ssolverinfo_init(&soption.solver_par, &soption.precond_par, queue);
	} else {
		doption.solver_par.solver = Magma_PIDRMERGE;
		doption.solver_par.restart = 8;
		doption.solver_par.maxiter = 500;
		doption.solver_par.rtol = 1e-14;
		doption.solver_par.atol = 1e-14;
		doption.solver_par.maxiter = 200;
		doption.precond_par.solver = Magma_ILU;
		doption.precond_par.levels = 0;
		doption.precond_par.trisolver = Magma_CUSOLVE;

		magma_dsolverinfo_init(&doption.solver_par, &doption.precond_par, queue);
	}
}

template<typename T> SparseMatMAGMA<T>::SparseMatMAGMA() { init_magma(); }

template<typename T> SparseMatMAGMA<T>::SparseMatMAGMA(const uword in_row, const uword in_col, const uword in_elem)
	: SparseMat(in_row, in_col, in_elem) { init_magma(); }

template<typename T> SparseMatMAGMA<T>::~SparseMatMAGMA() {
	if(std::is_same<T, float>::value) {
		magma_smfree(&sdx, queue);
		magma_smfree(&sdb, queue);
		magma_smfree(&sdA, queue);
		magma_smfree(&sx, queue);
		magma_smfree(&sb, queue);
		magma_smfree(&sA, queue);
	} else {
		magma_dmfree(&dx, queue);
		magma_dmfree(&db, queue);
		magma_dmfree(&dA, queue);
		magma_dmfree(&x, queue);
		magma_dmfree(&b, queue);
		magma_dmfree(&A, queue);
	}

	magma_queue_destroy(queue);
}

template<typename T> unique_ptr<MetaMat<T>> SparseMatMAGMA<T>::make_copy() { return make_unique<SparseMatMAGMA<T>>(*this); }

template<typename T> int SparseMatMAGMA<T>::solve(Mat<T>& out_mat, const Mat<T>& in_mat) {
	magma_int_t l_n_rows, l_n_cols;

	csr_form<T, uword> csr_mat(triplet_mat);

	podarray<magma_index_t> s_row_ptr(csr_mat.n_rows + 1);
	podarray<magma_index_t> s_col_idx(csr_mat.c_size);

	tbb::parallel_for(static_cast<int>(0), static_cast<int>(csr_mat.n_rows + 1), [&](const int I) { s_row_ptr[I] = static_cast<magma_index_t>(csr_mat.row_ptr[I]); });
	tbb::parallel_for(static_cast<int>(0), static_cast<int>(csr_mat.c_size), [&](const int I) { s_col_idx[I] = static_cast<magma_index_t>(csr_mat.col_idx[I]); });

	if(Precision::DOUBLE == precision) {
		// double precision

		out_mat = in_mat;

		magma_dcsrset(static_cast<magma_int_t>(csr_mat.n_rows), static_cast<magma_int_t>(csr_mat.n_cols), s_row_ptr.memptr(), s_col_idx.memptr(), csr_mat.val_idx, &A, queue);

		magma_dmtransfer(A, &dA, Magma_CPU, Magma_DEV, queue);

		for(uword I = 0; I < out_mat.n_cols; ++I) {
			magma_dvset(static_cast<magma_int_t>(out_mat.n_rows), 1, out_mat.colptr(I), &b, queue);

			magma_dmtransfer(b, &db, Magma_CPU, Magma_DEV, queue);

			magma_d_precondsetup(dA, db, &doption.solver_par, &doption.precond_par, queue);

			magma_dvinit(&dx, Magma_DEV, A.num_rows, 1, 0., queue);

			const auto info = magma_d_solver(dA, db, &dx, &doption, queue);

			if(0 != info) return SUANPAN_FAIL;

			magma_dmtransfer(dx, &x, Magma_DEV, Magma_CPU, queue);
			magma_dvcopy(x, &l_n_rows, &l_n_cols, out_mat.colptr(I), queue);
		}

		suanpan_debug([&]() { magma_dsolverinfo(&doption.solver_par, &doption.precond_par, queue); });
	} else {
		// single precision (mixed precision)

		podarray<float> s_val_idx(csr_mat.c_size);

		tbb::parallel_for(static_cast<int>(0), static_cast<int>(csr_mat.c_size), [&](const int I) { s_val_idx[I] = static_cast<float>(csr_mat.val_idx[I]); });

		out_mat.zeros(arma::size(in_mat));

		magma_scsrset(static_cast<magma_int_t>(csr_mat.n_rows), static_cast<magma_int_t>(csr_mat.n_cols), s_row_ptr.memptr(), s_col_idx.memptr(), s_val_idx.memptr(), &sA, queue);

		magma_smtransfer(sA, &sdA, Magma_CPU, Magma_DEV, queue);

		for(uword I = 0; I < out_mat.n_cols; ++I) {
			mat out_col(out_mat.colptr(I), in_mat.n_rows, 1, false, true);

			vec full_residual = in_mat.col(I);
			fvec residual(full_residual.n_elem);
			tbb::parallel_for(static_cast<uword>(0), full_residual.n_elem, [&](const uword J) { residual(J) = static_cast<float>(full_residual(J)); });
			fvec incre = arma::zeros<fvec>(residual.n_rows);

			auto multiplier = 1.;

			auto counter = 0;
			while(++counter < 10) {
				magma_svset(static_cast<magma_int_t>(residual.n_rows), 1, residual.memptr(), &sb, queue);

				magma_smtransfer(sb, &sdb, Magma_CPU, Magma_DEV, queue);

				magma_s_precondsetup(sdA, sdb, &soption.solver_par, &soption.precond_par, queue);

				magma_svinit(&sdx, Magma_DEV, sA.num_rows, 1, 0., queue);

				const auto info = magma_s_solver(sdA, sdb, &sdx, &soption, queue);

				if(0 != info) return SUANPAN_FAIL;

				magma_smtransfer(sdx, &sx, Magma_DEV, Magma_CPU, queue);
				magma_svcopy(sx, &l_n_rows, &l_n_cols, incre.memptr(), queue);

				suanpan_debug([&]() { magma_ssolverinfo(&soption.solver_par, &soption.precond_par, queue); });

				tbb::parallel_for(static_cast<uword>(0), static_cast<uword>(incre.n_rows), [&](const uword J) { out_col(J) += multiplier * incre(J); });

				full_residual = in_mat.col(I) - csr_mat * out_col;

				multiplier = norm(full_residual);

				suanpan_debug("mixed precision algorithm multiplier: %.5E\n", multiplier);

				if(multiplier < tolerance) break;

				tbb::parallel_for(static_cast<uword>(0), static_cast<uword>(residual.n_rows), [&](const uword J) { residual(J) = static_cast<float>(full_residual(J) / multiplier); });
			}
		}
	}

	return SUANPAN_SUCCESS;
}

#endif

#endif

//! @}
