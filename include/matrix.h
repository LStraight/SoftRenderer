#pragma once

#include <math.h>
#include <assert.h>
#include <iostream>
#include <ostream>
#include <sstream>

#include <string>

template<size_t R, size_t C, typename T>
struct Matrix {
	T m_mat[R][C];

	Matrix() {}
	Matrix(const Matrix<R,C,T>& m) {
		for (size_t i = 0; i < R; ++i) {
			for (size_t j = 0; j < C; ++j) {
				m_mat[i][j] = m.m_mat[i][j];
			}
		}
	}
	Matrix(const std::initializer_list<Vector<C, T>>& u) {
		auto it = u.begin();
		for (size_t i = 0; i < R; ++i) SetRow(i, *it++);
	}

	const T* operator[](size_t row)const {
		assert(row < R);
		return m_mat[row];
	}
	T* operator[](size_t row) {
		assert(row < R);
		return m_mat[row];
	}

	Vector<C, T> Row(size_t row) const {
		assert(row < R);
		Vector<C, T> a;
		for (size_t i; i < C; ++i)
			a[i] = m_mat[row][i];
		return a;
	}

	Vector<R, T> Col(size_t col) const {
		assert(col < C);
		Vector<R, T> a;
		for (size_t i; i < R; ++i)
			a[i] = m_mat[i][col];
		return a;
	}

	void SetRow(size_t row, const Vector<C, T>& a) {
		assert(row < R);
		for (size_t i = 0; i < C; ++i)
			m_mat[row][i] = a[i];
	}

	void SetCol(size_t col, const Vector<R, T>& a) {
		assert(col < C);
		for (size_t i = 0; i < R; ++i)
			m_mat[i][col] = a[i];
	}

	Matrix<R - 1, C - 1, T> GetMinor(size_t row, size_t col) const {
		Matrix<R - 1, C - 1, T> res;
		for (size_t i = 0; i < R - 1; ++i) {
			for (size_t j = 0; j < C - 1; ++j) {
				res.m_mat[i][j] = m_mat[i < row ? i : i + 1][j < col ? j : j + 1];
			}
		}
		return res;
	}

	Matrix<C, R, T> Transpose() const {
		Matrix<C, R, T> res;
		for (size_t i = 0; i < R; ++i) {
			for (size_t j = 0; j < C; ++j) {
				res.m_mat[i][j] = m_mat[j][i];
			}
		}
		return res;
	}

	static Matrix<R, C, T> GetZero() {
		Matrix<R, C, T> res;
		for (size_t i = 0; i < R; ++i) {
			for (size_t j = 0; j < C; ++j) {
				res.m_mat[i][j] = 0;
			}
		}
		return res;
	}

	static Matrix<R, C, T> GetIndentity() {
		Matrix<R, C, T> res;
		for (size_t i = 0; i < R; ++i) {
			for (size_t j = 0; j < C; ++j) {
				res.m_mat[i][j] = (i == j) ? 1 : 0;
			}
		}
		return res;
	}
};


typedef Matrix<4, 4, float> Mat44f;
typedef Matrix<3, 3, float> Mat33f;
typedef Matrix<4, 3, float> Mat43f;
typedef Matrix<3, 4, float> Mat34f;