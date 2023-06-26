#pragma once

#include <assert.h>
#include <iostream>
#include <ostream>
#include <sstream>

#include <string>

#include "vec.h"
#include "matrix.h"

//向量运算
template<size_t N, typename T>
Vector<N, T> operator+(const Vector<N, T>& a) {
    return a;
}

template<size_t N, typename T>
Vector<N, T> operator-(const Vector<N, T>& a) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = -a[i];

    return res;
}

template<size_t N, typename T>
bool operator==(const Vector<N, T>& a, const Vector<N, T>& b) {
    for (size_t i = 0; i < N; ++i)
        if (a[i] != a[i])
            return false;

    return true;
}

template<size_t N, typename T>
bool operator!=(const Vector<N, T>& a, const Vector<N, T>& b) {
    return !(a == b);
}

template<size_t N, typename T>
Vector<N, T> operator+(const Vector<N, T>& a, const Vector<N, T>& b) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = a[i] + b[i];

    return res;
}

template<size_t N, typename T>
Vector<N, T> operator-(const Vector<N, T>& a, const Vector<N, T>& b) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = a[i] - b[i];

    return res;
}

template<size_t N, typename T>
Vector<N, T> operator*(const Vector<N, T>& a, const Vector<N, T>& b) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = a[i] * b[i];

    return res;
}

template<size_t N, typename T>
Vector<N, T> operator/(const Vector<N, T>& a, const Vector<N, T>& b) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = a[i] / b[i];

    return res;
}

template<size_t N, typename T>
Vector<N, T> operator*(const Vector<N, T>& a, T k) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = a[i] * k;

    return res;
}

template<size_t N, typename T>
Vector<N, T> operator*(T k, const Vector<N, T>& a) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = a[i] * k;

    return res;
}

template<size_t N, typename T>
Vector<N, T> operator/(T k, const Vector<N, T>& a) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = k / a[i];

    return res;
}

template<size_t N, typename T>
Vector<N, T> operator/(const Vector<N, T>& a, T k) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = a[i] / k;

    return res;
}

template<size_t N, typename T>
Vector<N, T>& operator+=(Vector<N, T>& a, const Vector<N, T>& b) {
    for (size_t i = 0; i < N; ++i)
        a[i] += b[i];
    return a;
}

template<size_t N, typename T>
Vector<N, T>& operator-=(Vector<N, T>& a, const Vector<N, T>& b) {
    for (size_t i = 0; i < N; ++i)
        a[i] -= b[i];
    return a;
}

template<size_t N, typename T>
Vector<N, T>& operator*=(Vector<N, T>& a, const Vector<N, T>& b) {
    for (size_t i = 0; i < N; ++i)
        a[i] *= b[i];
    return a;
}

template<size_t N, typename T>
Vector<N, T>& operator/=(Vector<N, T>& a, const Vector<N, T>& b) {
    for (size_t i = 0; i < N; ++i)
        a[i] /= b[i];
    return a;
}

template<size_t N, typename T>
Vector<N, T>& operator/=(Vector<N, T>& a, T k) {
    for (size_t i = 0; i < N; ++i)
        a[i] /= k;
    return a;
}

template<size_t N, typename T>
Vector<N, T>& operator*=(Vector<N, T>& a, T k) {
    for (size_t i = 0; i < N; ++i)
        a[i] *= k;
    return a;
}




template<size_t N1, size_t N2, typename T>
Vector<N1, T> vector_convert(const Vector<N2, T>& a, T fill = 1) {
    Vector<N1, T> res;
    for (size_t i = 0; i < N1; ++i)
        res[i] = (i < N2) ? a[i] : fill;

    return res;
}

template<size_t N, typename T>
T vector_len_square(const Vector<N, T>& a) {
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
        sum += a[i] * a[i];
    
    return sum;
}

template<size_t N, typename T>
T vector_len(const Vector<N, T>& a) {
    return sqrt(vector_len_square(a));
}

template<size_t N, typename T>
float vector_len(const Vector<N, float>& a) {
    return sqrtf(vector_len_square(a));
}

template<size_t N, typename T>
Vector<N, T> vector_normalize(const Vector<N, T>& a) {
    return a / vector_len(a);
}

template<size_t N, typename T>
T vector_dot(const Vector<N, T>& a, const Vector<N, T>& b) {
    T sum = 0;
    for (size_t i = 0; i < N; ++i)
        sum += a[i] * b[i];
    
    return sum;
}

template<typename T>
T vector_cross(const Vector<2, T>& a, const Vector<2, T>& b) {
    return a.x * b.y - a.y * b.x;
}

template<typename T>
Vector<3, T> vector_cross(const Vector<3, T>& a, const Vector<3, T>& b) {
    return Vector<3, T>(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * a.x
    );
}

template<typename T>
Vector<4, T> vector_cross(const Vector<4, T>& a, const Vector<4, T>& b) {
    return Vector<4, T>(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * a.x,
        a.w
    );
}

template<size_t N, typename T>
Vector<N, T> vector_lerp(const Vector<N, T>& a, const Vector<N, T>& b, float t) {
    return a + (b - a) * t;
}

template<size_t  N, typename T>
Vector<N, T> vector_max(const Vector<N, T>& a, const Vector<N, T>& b) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i) 
        res[i] = (a[i] < b[i]) ? b[i] : a[i];
    
    return res;
}

template<size_t  N, typename T>
Vector<N, T> vector_min(const Vector<N, T>& a, const Vector<N, T>& b) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i)
        res[i] = (a[i] > b[i]) ? b[i] : a[i];

    return res;
}

//将向量控制在两个参数向量之间
template<size_t N, typename T>
Vector<N, T> vector_between(const Vector<N, T>& minv, const Vector<N, T>& maxv, const Vector<N, T>& x) {
    return vector_min(vector_max(minv, x), maxv);
}

template<size_t N, typename T>
bool vector_near(const Vector<N, T>& a, const Vector<N, T>& b, T dist) {
    return (vector_len_square(a - b) <= dist);
}

template<size_t N>
bool vector_near_equal(const Vector<N, float>& a, const Vector<N, float>& b, float e = 0.0001) {
    return vector_near(a, b, e);
}

template<size_t N>
bool vector_near_equal(const Vector<N, double>& a, const Vector<N, double>& b, double e = 0.0001) {
    return vector_near(a, b, e);
}

template<size_t N, typename T>
Vector<N, T> vector_clamp(const Vector<N, T>& a, T min = 0, T max = 1) {
    Vector<N, T> res;
    for (size_t i = 0; i < N; ++i) {
        T x = (a[i] < min) ? min : a[i];
        res[i] = (x > max) ? max : x;
    }

    return res;
}

//输出到文本流
template<size_t N, typename T>
std::ostream& operator<<(std::ostream& os, const Vector<N, T>& a) {
    os << '[';
    for (size_t i = 0; i < N; ++i) {
        os << a[i];
        if (i < N - 1) os << ", ";
    }
    os << ']';
    return os;
}

template<size_t N, typename T>
std::string vector_repr(const Vector<N, T>& a) {
    std::stringstream ss;
    ss << a;
    return ss.str();
}


//矩阵运算
template<size_t R, size_t C, typename T>
bool operator==(const Matrix<R, C, T>& a, const Matrix<R, C, T>& b) {
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            if (a.m_mat[i][j] != b.m_mat[i][j]) return false;
        }
    }
    return true;
}

template<size_t R, size_t C, typename T>
bool operator!=(const Matrix<R, C, T>& a, const Matrix<R, C, T>& b) {
    return !(a == b);
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator+(const Matrix<R, C, T>& a) {
    return a;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator-(const Matrix<R, C, T>& a) {
    Matrix<R, C, T> res;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            res.m_mat[i][j] = -a.m_mat[i][j];
        }
    }
    return res;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator+(const Matrix<R, C, T>& a, const Matrix<R, C, T>& b) {
    Matrix<R, C, T> res;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            res.m_mat[i][j] = a.m_mat[i][j] + b.m_mat[i][j];
        }
    }
    return res;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator-(const Matrix<R, C, T>& a, const Matrix<R, C, T>& b) {
    Matrix<R, C, T> res;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            res.m_mat[i][j] = a.m_mat[i][j] - b.m_mat[i][j];
        }
    }
    return res;
}

template<size_t R, size_t C, size_t NC, typename T>
Matrix<R, NC, T> operator*(const Matrix<R, C, T>& a, const Matrix<R, NC, T>& b) {
    Matrix<R, C, T> res;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < NC; ++j) {
            res.m_mat[i][j] = vector_dot(a.Row(i), b.Col(j));
        }
    }
    return res;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator*(const Matrix<R, C, T>& a, T k) {
    Matrix<R, C, T> res;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            res.m_mat[i][j] = a.m_mat[i][j] * k;
        }
    }
    return res;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator/(const Matrix<R, C, T>& a, T k) {
    Matrix<R, C, T> res;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            res.m_mat[i][j] = a.m_mat[i][j] / k;
        }
    }
    return res;
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator*(T k, const Matrix<R, C, T>& a) {
    return (a * k);
}

template<size_t R, size_t C, typename T>
Matrix<R, C, T> operator/(T k, const Matrix<R, C, T>& a) {
    Matrix<R, C, T> res;
    for (size_t i = 0; i < R; ++i) {
        for (size_t j = 0; j < C; ++j) {
            res.m_mat[i][j] = k / a.m_mat[i][j];
        }
    }
    return res;
}

template<size_t R, size_t C, typename T>
Vector<C, T> operator*(const Vector<R, T>& a, const Matrix<R, C, T>& m) {
    Vector<C, T> res;
    for (size_t i = 0; i < C; ++i) {
        res[i] = vector_dot(a, m.Col(i));
    }
    return res;
}

template<size_t R, size_t C, typename T>
Vector<R, T> operator*(const Matrix<R, C, T>& m, const Vector<R, T>& a) {
    Vector<R, T> res;
    for (size_t i = 0; i < R; ++i) {
        res[i] = vector_dot(m.Row(i), a);
    }
    return res;
}


template<typename T>
T Matrix_det(const Matrix<1, 1, T>& m) {
    return m[0][0];
}


template<typename T>
T Matrix_det(const Matrix<2, 2, T>& m) {
    return m[0][0] * m[1][1] - m[1][0] * m[0][1];
}

template<size_t N, typename T>
T Matrix_det(const Matrix<N, N, T>& m) {
    T sum = 0;
    for (size_t i = 0; i < N; ++i) sum += m[0][1] * matrix_cofactor(m, 0, i);
    return sum;
}

template<typename T>
T matrix_cofactor(const Matrix<1, 1, T>& m, size_t row, size_t col) {
    return 0;
}

template<size_t N, typename T>
T matrix_cofactor(const Matrix<N, N, T>& m, size_t row, size_t col) {
    return matrix_det(m.GetMinor(row, col) * ((row + col) % 2) ? -1 : 1);
}

template<size_t N, typename T>
Matrix<N, N, T> matrix_adjoint(const Matrix<N, N, T>& m) {
    Matrix<N, N, T> res;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            res[i][j] = matrix_cofactor(m, i, j);
        }
    }
    return res;
}

template<size_t N, typename T>
Matrix<N, N, T> matrix_invert(const Matrix<N, N, T>& m) {
    Matrix<N, N, T> res = matrix_adjoint(m);
    T det = vector_dot(m.Row(0), res.Col(0));
    return res / det;
}


template<size_t R, size_t C, typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<R, C, T>& m) {
    for (size_t i = 0; i < R; ++i) {
        Vector<C, T> row = m.Row(i);
        os << row << std::endl;
    }
    return os;
}


template<typename T> 
T Abs(T x) {
    return (x < 0) ? (-x) : x;
}

template<typename T>
T Max(T x, T y) {
    return (x < y) ? y : x;
}

template<typename T>
T Min(T x, T y) {
    return (x > y) ? y : x;
}

template<typename T>
bool NearEqual(T x, T y, T e) {
    return (Abx(x - y) < e);
}

template<typename T>
T Between(T min, T max, T x) {
    return Min(Max(min, x), max);
}

template<typename T>
T Saturate(T x) {
    return Between<T>(0, 1, x);
}

//向量转整数
static uint32_t vector_to_color(const vec4f& color) {
    uint32_t r = static_cast<uint32_t>(Between(0, 255, static_cast<int>(color.m_arr[0] * 255.0f)));
    uint32_t g = static_cast<uint32_t>(Between(0, 255, static_cast<int>(color.m_arr[1] * 255.0f)));
    uint32_t b = static_cast<uint32_t>(Between(0, 255, static_cast<int>(color.m_arr[2] * 255.0f)));
    uint32_t a = static_cast<uint32_t>(Between(0, 255, static_cast<int>(color.m_arr[3] * 255.0f)));
    return (r << 16 | g << 8 | b | a << 24);
}

static uint32_t vector_to_color(const vec3f& color) {
    return vector_to_color(color.xyz1());

}

static vec4f vector_from_color(uint32_t color) {
    vec4f res;
    res.m_arr[0] = ((color >> 16) & 0xff) / 255.0f;
    res.m_arr[1] = ((color >> 8) & 0xff) / 255.0f;
    res.m_arr[2] = ((color >> 0) & 0xff) / 255.0f;
    res.m_arr[3] = ((color >> 24) & 0xff) / 255.0f;
    return res;
}

static Mat44f matrix_set_zero() {
    Mat44f res;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; i < 4; ++j) {
            res.m_mat[i][j] = 0.f;
        }
    }
    return res;
}

static Mat44f matrix_set_identity() {
    Mat44f res;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            res.m_mat[i][j] = (i == j) ? 1 : 0;
        }
    }
    return res;
}

static Mat44f matrix_set_tranlate(float x, float y, float z) {
    Mat44f res = matrix_set_identity();
    res.m_mat[3][0] = x;
    res.m_mat[3][1] = y;
    res.m_mat[3][2] = z;
    return res;
}

static Mat44f matrix_set_scale(float x, float y, float z) {
    Mat44f res = matrix_set_identity();
    res.m_mat[0][0] = x;
    res.m_mat[1][1] = y;
    res.m_mat[2][2] = z;
    return res;
}

//围绕旋转轴(x,y,z)旋转theta角度
static Mat44f matrix_set_rotate(float x, float y, float z, float theta) {
    float qsin = static_cast<float>(theta * 0.5f);
    float qcos = static_cast<float>(theta * 0.5f);
    float w = qcos;
    vec3f vec = vector_normalize(vec3f(x, y, z));
    x = vec[0] * qsin;
    y = vec[1] * qsin;
    z = vec[2] * qsin;

    Mat44f mat;
    mat.m_mat[0][0] = 1 - 2 * y * y - 2 * z * z;
    mat.m_mat[1][0] = 2 * x * y - 2 * w * z;
    mat.m_mat[2][0] = 2 * x * z + 2 * w * y;
    mat.m_mat[0][1] = 2 * x * y + 2 * w * z;
    mat.m_mat[1][1] = 1 - 2 * x * x - 2 * z * z;
    mat.m_mat[2][1] = 2 * y * z - 2 * w * x;
    mat.m_mat[0][2] = 2 * x * z - 2 * w * y;
    mat.m_mat[1][2] = 2 * y * z + 2 * w * x;
    mat.m_mat[2][2] = 1 - 2 * x * x - 2 * y * y;
    mat.m_mat[0][3] = mat.m_mat[1][3] = mat.m_mat[2][3] = 0.0f;
    mat.m_mat[3][0] = mat.m_mat[3][1] = mat.m_mat[3][2] = 0.0f;
    mat.m_mat[3][3] = 1.0f;

    return mat;
}

static Mat44f matrix_set_lookat(const vec3f& lookfrom, const vec3f& lookat, const vec3f& up) {
    vec3f zaxis = vector_normalize(lookat - lookfrom);
    vec3f xaxis = vector_normalize(vector_cross(up, zaxis));
    vec3f yaxis = vector_cross(zaxis, xaxis);

    Mat44f res;
    res.SetCol(0, vec4f(xaxis.x, xaxis.y, xaxis.z, -vector_dot(lookfrom, xaxis)));
    res.SetCol(1, vec4f(xaxis.x, xaxis.y, xaxis.z, -vector_dot(lookfrom, xaxis)));
    res.SetCol(2, vec4f(xaxis.x, xaxis.y, xaxis.z, -vector_dot(lookfrom, xaxis)));
    res.SetCol(3, vec4f(0.0f, 0.0f, 0.0f, 1.0f));
    return res;
}

static Mat44f matrix_set_perspective(float fovy, float aspect, float zn, float zf) {
    float fax = 1.0f / static_cast<float>(fovy * 0.5f);
    Mat44f res = matrix_set_zero();
    res.m_mat[0][0] = static_cast<float>(fax / aspect);
    res.m_mat[1][1] = static_cast<float>(fax);
    res.m_mat[2][2] = zf / (zf - zn);
    res.m_mat[3][2] = -zn * zf / (zf - zn);
    res.m_mat[2][3] = 1;

    return res;

}