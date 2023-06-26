#pragma once

#include <math.h>
#include <assert.h>


//ÏòÁ¿Ä£°å
template<size_t N, typename T>
struct Vector {
    T m_arr[N];
    Vector() {
        for (size_t i = 0; i < N; ++i)
            m_arr[i] = T();
    }
    Vector(const T* ptr) {
        for (size_t i = 0; i < N; ++i)
            m_arr[i] = ptr[i];
    }
    Vector(const Vector<N, T>& v) {
        for (size_t i = 0; i < N; ++i)
            m_arr[i] = v.m_arr[i];
    }
    Vector(const std::initializer_list<T>& v) {
        auto it = v.begin();
        for (size_t i = 0; i < N; ++i)
            m_arr[i] = *it++;
    }

    const T& operator[](size_t i)const {
        assert(i < N);
        return m_arr[i];
    }
    T& operator[](size_t i) {
        assert(i < N);
        return m_arr[i];
    }

    void load(const T* ptr) {
        for (size_t i = 0; i < N; ++i) {
            m_arr[i] = ptr[i];
        }
    }
    void save(T* ptr) {
        for (size_t i = 0; i < N; ++i)
            ptr[i] = m_arr[i];
    }
};

template<typename T>
struct Vector<2, T> {
    union {
        struct { T x, y; };
        struct { T u, v; };
        T m_arr[2];
    };

    inline Vector(): x(T()), y(T()) {}
    Vector(T _x, T _y) :x(_x), y(_y) {}
    Vector(const Vector<2, T>& u) :x(u.x), y(u.y) {}
    Vector(const T* ptr) :x(ptr[0]), y(ptr[1]) {}

    const T& operator[](size_t i) const {
        assert(i < 2);
        return m_arr[i];
    }
    T& operator[](size_t i) {
        assert(i < 2);
        return m_arr[i];
    }

    void load(const T* ptr) {
        m_arr[0] = ptr[0];
        m_arr[1] = ptr[1];
    }
    void save(T* ptr) {
        ptr[0] = m_arr[0];
        ptr[1] = m_arr[1];
    }

    Vector<2, T> xy() const { return *this; }
    Vector<3, T> xy1() const { return Vector<3, T>(x, y, 1); }
    Vector<4, T> xy11() const { return Vector<4, T>(x, y, 1, 1); }
};

template<typename T>
struct Vector<3, T> {
    union {
        struct { T x, y, z; };
        struct { T r, g, b; };
        T m_arr[3];
    };

    Vector() :x(T()), y(T()), z(T()) {}
    Vector(T _x, T _y, T _z) :x(_x), y(_y), z(_z) {}
    Vector(const Vector<3, T>& u) :x(u.x), y(u.y), z(u.z) {}
    Vector(const T* ptr) :x(ptr[0]), y(ptr[1]), z(ptr[2]) {}

    const T& operator[](size_t i) const {
        assert(i < 3);
        return m_arr[i];
    }
    T& operator[](size_t i) {
        assert(i < 3);
        return m_arr[i];
    }

    void load(const T* ptr) {
        m_arr[0] = ptr[0];
        m_arr[1] = ptr[1];
        m_arr[2] = ptr[2];
    }
    void save(T* ptr) {
        ptr[0] = m_arr[0];
        ptr[1] = m_arr[1];
        ptr[2] = m_arr[2];
    }

    Vector<2, T> xy() const { return Vector<2, T>(x, y); }
    Vector<3, T> xyz() const { return *this; }
    Vector<4, T> xyz1() const { return Vector<4, T>(x, y, z, 1); }
};

template<typename T>
struct Vector<4, T> {
    union {
        struct { T x, y, z, w; };
        struct { T r, g, b, a; };
        T m_arr[4];
    };

    Vector() :x(T()), y(T()), z(T()), w(T()) {}
    Vector(T _x, T _y, T _z, T _w) :x(_x), y(_y), z(_z), w(_w) {}
    Vector(const Vector<4, T>& u) :x(u.x), y(u.y), z(u.z), w(u.w) {}
    Vector(const T* ptr) :x(ptr[0]), y(ptr[1]), z(ptr[2]), w(ptr[3]) {}

    const T& operator[](size_t i) const {
        assert(i < 4);
        return m_arr[i];
    }
    T& operator[](size_t i) {
        assert(i < 4);
        return m_arr[i];
    }

    void load(const T* ptr) {
        for (size_t i = 0; i < 4; ++i)
            m_arr[i] = ptr[i];
    }
    void save(T* ptr) {
        for (size_t i = 0; i < 4; ++i)
            ptr[i] = m_arr[i];
    }

    Vector<2, T> xy() const { return Vector<2, T>(x, y); }
    Vector<3, T> xyz() const { return Vector<3, T>(x, y, z); }
    Vector<4, T> xyzw() const { return *this; }
};


typedef Vector<2, float>    vec2f;
typedef Vector<2, double>   vec2d;
typedef Vector<2, int>      vec2i;
typedef Vector<3, float>    vec3f;
typedef Vector<3, double>   vec3d;
typedef Vector<3, int>      vec3i;
typedef Vector<4, float>    vec4f;
typedef Vector<4, double>   vec4d;
typedef Vector<4, int>      vec4i;





