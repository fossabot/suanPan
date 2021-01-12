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
 * @class Quaternion
 * @brief An Quaternion class.
 *
 * @author tlc
 * @date 05/09/2020
 * @version 0.1.0
 * @file Quaternion.hpp
 * @addtogroup Utility
 * @{
 */

#ifndef QUATERNION_H
#define QUATERNION_H

#include <suanPan.h>
#include <Toolbox/tensorToolbox.h>

template<typename T> class Quaternion {
	T re;
	Col<T> im;
public:
	Quaternion();
	Quaternion(T, T, T, T);
	Quaternion(T, const Col<T>&);
	Quaternion(T, Col<T>&&);

	const T& real() const;
	const Col<T>& imag() const;

	T norm() const;
	Quaternion& normalise();

	Quaternion inv() const;
	Quaternion conj() const;

	Quaternion operator+(const Quaternion&) const;
	Quaternion& operator+=(const Quaternion&);
	Quaternion operator-(const Quaternion&) const;
	Quaternion& operator-=(const Quaternion&);
	Quaternion operator*(const Quaternion&) const;
	Quaternion& operator*=(const Quaternion&);
	Quaternion operator/(const Quaternion&) const;
	Quaternion& operator/=(const Quaternion&);

	Mat<T> operator*(const Mat<T>&) const;
};

template<typename T> Quaternion<T>::Quaternion()
	: re(T(0))
	, im(arma::zeros<Col<T>>(3)) {}

template<typename T> Quaternion<T>::Quaternion(const T R, const T I, const T J, const T K)
	: re(R)
	, im({I, J, K}) {}

template<typename T> Quaternion<T>::Quaternion(const T R, const Col<T>& I)
	: re(R)
	, im(I) {}

template<typename T> Quaternion<T>::Quaternion(const T R, Col<T>&& I)
	: re(R)
	, im(std::forward<Col<T>>(I)) {}

template<typename T> const T& Quaternion<T>::real() const { return re; }

template<typename T> const Col<T>& Quaternion<T>::imag() const { return im; }

template<typename T> T Quaternion<T>::norm() const { return re * re + arma::dot(im, im); }

template<typename T> Quaternion<T>& Quaternion<T>::normalise() {
	const auto magnitude = std::sqrt(norm());
	re /= magnitude;
	im /= magnitude;
	return *this;
}

template<typename T> Quaternion<T> Quaternion<T>::inv() const {
	const auto L = norm();

	return Quaternion<T>(re / L, -im / L);
}

template<typename T> Quaternion<T> Quaternion<T>::conj() const { return Quaternion<T>(re, -im); }

template<typename T> Quaternion<T> Quaternion<T>::operator+(const Quaternion& B) const {
	Quaternion<T> A = *this;

	return A += B;
}

template<typename T> Quaternion<T>& Quaternion<T>::operator+=(const Quaternion& B) {
	re += B.re;
	im += B.im;

	return *this;
}

template<typename T> Quaternion<T> Quaternion<T>::operator-(const Quaternion& B) const {
	Quaternion<T> A = *this;

	return A -= B;
}

template<typename T> Quaternion<T>& Quaternion<T>::operator-=(const Quaternion& B) {
	re -= B.re;
	im -= B.im;

	return *this;
}

template<typename T> Quaternion<T> Quaternion<T>::operator*(const Quaternion& B) const {
	Quaternion<T> A;

	A.re = re * B.re - arma::dot(im, B.im);
	A.im = re * B.im + B.re * im + arma::cross(im, B.im);

	return A;
}

template<typename T> Quaternion<T>& Quaternion<T>::operator*=(const Quaternion& B) { return *this = *this * B; }

template<typename T> Quaternion<T> Quaternion<T>::operator/(const Quaternion& B) const { return *this * B.inv(); }

template<typename T> Quaternion<T>& Quaternion<T>::operator/=(const Quaternion& B) { return *this = *this * B.inv(); }

template<typename T> Mat<T> Quaternion<T>::operator*(const Mat<T>& I) const {
	const Mat<T> R = 2. * re * transform::skew_symm(im) + 2. * im * im.t() + (re * re - arma::dot(im, im)) * eye(3, 3);

	return R * I;
}

#endif

//! @}
