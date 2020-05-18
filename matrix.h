#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>

//...................................CHECK_COMPOSITE................................................

template<unsigned I, unsigned M>
class check {
public:
    static const unsigned n = std::max((int)(check<I - 1, M>::n), (int)(M % I == 0));
};

template<unsigned M>
class check<1, M> {
public:
    static const unsigned n = 0;
};

template<unsigned M>
class check<0, M> {
public:
    static const unsigned n = 0;
};

//..............................ASSERT_N_IS_EQUAL_TO_M.........................................................

template<unsigned N, unsigned M>
void ASSERT_A_IS_EQUAL_TO_B() {
    int a[(N == M) ? 1 : -1];
    a[0] = 3;
    if (a[0] == 5) {
        return;
    }
    else {
        return;
    }
}

//.............................ASSERT_COMPOSITE....................................................................

template<unsigned M>
void ASSERT_COMPOSITE() {
    int a[(0 == check<int(sqrt(M)), M>::n) ? 1 : -1];
    a[0] = 3;
    if (a[0] == 5) {
        return;
    }
    else {
        return;
    }
};

//..................................FINITE............................................................
template<unsigned M>
class Finite {
public:
    Finite();//+
    Finite(const int& other);//+
    Finite(const Finite<M>& other);//+
    ~Finite();//+

    Finite& operator =(const Finite& other);//+

    bool operator == (const Finite& other) const;
    bool operator != (const Finite& other) const;
    Finite operator +(const Finite& other) const;//+
    Finite operator -(const Finite& other) const;//+
    Finite operator *(const Finite& other) const;//+
    Finite operator /(const Finite& other) const;//+
    Finite& operator ++();
    Finite operator ++(int);
    Finite& operator +=(const Finite& other);//+
    Finite& operator -=(const Finite& other);//+
    Finite& operator *=(const Finite& other);//+
    Finite& operator /=(const Finite& other);//+
    operator int() const;

    friend std::istream& operator >>(std::istream& in, Finite& finite) {
        in >> finite.num;
        finite.norm();
        return in;
    }
    friend std::ostream& operator <<(std::ostream& out, const Finite& finite) {
        out << finite.num;
        return out;
    }

    Finite power(const unsigned& pow) const;//+
    Finite inverse() const;//+
private:
    void norm();//+
    int num = 0;
};

template<unsigned M>
Finite<M> abs(Finite<M> F) {
    return F;
}

template<unsigned M>
void Finite<M>::norm() {
    num %= (int)M;
    num += (int)M;
    num %= (int)M;
}

template<unsigned M>
Finite<M>::operator int() const {
    return num;
}

template<unsigned M>
Finite<M>::Finite() {
    num = 0;
}

template<unsigned M>
Finite<M>::Finite(const int& other) {
    num = other;
    norm();
}

template<unsigned M>
Finite<M>::Finite(const Finite<M>& other) {
    num = other.num;
    norm();
}

template<unsigned M>
Finite<M>::~Finite() {}

template<unsigned M>
Finite<M>& Finite<M>::operator =(const Finite<M>& other) {
    num = other.num;
    return *this;
}

template<unsigned M>
bool Finite<M>::operator ==(const Finite<M>& other) const {
    return num == other.num;
}

template<unsigned M>
bool Finite<M>::operator !=(const Finite<M>& other) const {
    return num != other.num;
}

template<unsigned M>
Finite<M>& Finite<M>::operator +=(const Finite<M>& other) {
    num += other.num;
    norm();
    return *this;
}

template<unsigned M>
Finite<M>& Finite<M>::operator -=(const Finite<M>& other) {
    num += M - other.num;
    norm();
    return *this;
}

template<unsigned M>
Finite<M>& Finite<M>::operator *=(const Finite<M>& other) {
    num *= other.num;
    norm();
    return *this;
}

template<unsigned M>
Finite<M>& Finite<M>::operator /=(const Finite<M>& other) {
    *this *= other.inverse();
    return *this;
}

template<unsigned M>
Finite<M>& Finite<M>::operator ++() {
    return *this += 1;
}

template<unsigned M>
Finite<M> Finite<M>::operator ++(int) {
    Finite<M> newFinite = *this;
    *this += 1;
    return newFinite;
}

template<unsigned M>
Finite<M> Finite<M>::operator +(const Finite<M>& other) const {
    Finite<M> newFinite = *this;
    return newFinite += other;
}

template<unsigned M>
Finite<M> Finite<M>::operator -(const Finite<M>& other) const {
    Finite<M> newFinite = *this;
    return newFinite -= other;
}

template<unsigned M>
Finite<M> Finite<M>::operator *(const Finite<M>& other) const {
    Finite<M> newFinite = *this;
    return newFinite *= other;
}

template<unsigned M>
Finite<M> Finite<M>::operator /(const Finite<M>& other) const {
    Finite<M> newFinite = *this;
    return newFinite /= other;
}

template<unsigned M>
Finite<M> Finite<M>::power(const unsigned& pow) const {
    if (pow == 0) {
        Finite<M> res = 1;
        return res;
    }
    if (pow == 1) {
        Finite<M> res = *this;
        return res;
    }
    if (pow % 2 == 0) {
        return (power(pow / 2) * power(pow / 2));
    }
    else {
        return (power(pow - 1) * (*this));
    }
}

template<unsigned M>
Finite<M> Finite<M>::inverse() const {
    ASSERT_COMPOSITE<M>();
    return power(M - 2);
}

//....................................FIND_DEGREE.............................................................

namespace FindDegree {
    template<unsigned N, unsigned BASE>
    class findLog {
    public:
        static const unsigned log = findLog<N / BASE, BASE>::log + 1;
    };

    template<unsigned BASE>
    class findLog<0, BASE> {
    public:
        static const unsigned log = -1;
    };

    template<unsigned N, unsigned BASE>
    class lowerNum {
    public:
        static const unsigned num = lowerNum<N - 1, BASE>::num * BASE;
    };

    template<unsigned BASE>
    class lowerNum <0, BASE> {
    public:
        static const unsigned num = 1;
    };

    template<unsigned N, unsigned BASE>
    class findDeg {
    public:
        static const unsigned logN = findLog<N, BASE>::log;
        static const unsigned lowerAns = lowerNum<logN, BASE>::num;
        static const unsigned ans = (lowerAns == N) ? N : lowerAns * BASE;
    };

    template<unsigned BASE>
    class findDeg<1, BASE> {
    public:
        static const unsigned ans = BASE;
    };
}

//...........................................ARRAY_FUNCTIONS..................................................

template<unsigned N, unsigned M, typename Field = Rational>
Field** generateArray() {
    Field** arr = new Field * [N];
    for (unsigned i = 0; i < N; ++i) {
        arr[i] = new Field[M];
        for (unsigned j = 0; j < M; ++j) {
            arr[i][j] = 0;
        }
    }
    return arr;
}

template<unsigned N, unsigned M, typename Field = Rational>
void fillArray(Field**& arr, Field** otherArr) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            arr[i][j] = otherArr[i][j];
        }
    }
}

template<unsigned N, unsigned M, typename Field = Rational>
void deleteArray(Field**& a) {
    for (unsigned i = 0; i < N; ++i) {
        delete[] a[i];
    }
    delete[] a;
}

//.........................RATIONAL.....................................

Rational abs(const Rational& num) {
    return (num >= 0) ? num : -num;
}

//.........................................BASEMATRIX.............................................................

template<unsigned N, unsigned M, typename Field = Rational> //finite<K>, K != simple
class BaseMatrix {
protected:
    Field** a;
    unsigned rowsCount();
    unsigned colCount();
public:
    BaseMatrix();
    BaseMatrix(const int& other);
    BaseMatrix(const BaseMatrix<N, M, Field>& other);
    ~BaseMatrix();
    BaseMatrix& operator = (const BaseMatrix<N, M, Field>& other);

    bool operator == (const BaseMatrix<N, M, Field>& other) const;
    bool operator != (const BaseMatrix<N, M, Field>& other) const;

    BaseMatrix& operator += (const BaseMatrix<N, M, Field>& other);
    BaseMatrix& operator -= (const BaseMatrix<N, M, Field>& other);

    BaseMatrix operator + (const BaseMatrix<N, M, Field>& other) const;
    BaseMatrix operator - (const BaseMatrix<N, M, Field>& other) const;

    BaseMatrix& operator *= (const Field& other);
    BaseMatrix operator * (const Field& other) const;

    template<unsigned K>
    BaseMatrix<N, K, Field> operator * (const BaseMatrix<M, K, Field>& other) const;//if can

    BaseMatrix<N, M, Field>& operator *= (const BaseMatrix<M, M, Field>& other);//if can

    Field det(); //if can

    BaseMatrix<M, N, Field> transposed() const;

    unsigned rank() const;

    BaseMatrix inverted() const;//if can
    BaseMatrix& invert();//if can

    Field trace() const;//if can

    std::vector<Field> getRow(const unsigned& cnt) const;
    std::vector<Field> getColumn(const unsigned& cnt) const;

    Field*& operator [](const unsigned& cnt);
    Field* operator [](const unsigned& cnt) const;

    BaseMatrix<N / 2, N / 2, Field>** quater() const;
};

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field> operator * (const Field& F, const BaseMatrix<N, M, Field>& Matr) {
    return Matr * F;
}

template<unsigned N, unsigned M, typename Field>
unsigned BaseMatrix<N, M, Field>::rowsCount() {
    return N;
}

template<unsigned N, unsigned M, typename Field>
unsigned BaseMatrix<N, M, Field>::colCount() {
    return M;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>::BaseMatrix() {
    a = generateArray<N, M, Field>();
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>::BaseMatrix(const int& other) {
    a = generateArray<N, M, Field>();
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            a[i][j] = other;
        }
    }
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>::BaseMatrix(const BaseMatrix& other) {
    a = generateArray<N, M, Field>();
    fillArray<N, M, Field>(a, other.a);
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>::~BaseMatrix() {
    deleteArray<N, M, Field>(a);
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>& BaseMatrix<N, M, Field>::operator = (const BaseMatrix& other) {
    fillArray<N, M, Field>(a, other.a);
    return *this;
}

template<unsigned N, unsigned M, typename Field>
bool BaseMatrix<N, M, Field>::operator == (const BaseMatrix<N, M, Field>& other) const {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            if (a[i][j] != other[i][j]) return 0;
        }
    }
    return 1;
}

template<unsigned N, unsigned M, typename Field>
bool BaseMatrix<N, M, Field>::operator != (const BaseMatrix<N, M, Field>& other) const {
    return !(*this == other);
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>& BaseMatrix<N, M, Field>::operator += (const BaseMatrix& other) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            a[i][j] += other[i][j];
        }
    }
    return *this;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>& BaseMatrix<N, M, Field>::operator -= (const BaseMatrix& other) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            a[i][j] -= other[i][j];
        }
    }
    return *this;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field> BaseMatrix<N, M, Field>::operator + (const BaseMatrix& other) const {
    BaseMatrix New = *this;
    New += other;
    return New;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field> BaseMatrix<N, M, Field>::operator - (const BaseMatrix& other) const {
    BaseMatrix New = *this;
    New -= other;
    return New;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>& BaseMatrix<N, M, Field>::operator *= (const Field& other) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            a[i][j] *= other;
        }
    }
    return *this;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field> BaseMatrix<N, M, Field>::operator * (const Field& other) const {
    BaseMatrix New = *this;
    New *= other;
    return New;
}

template<unsigned N, unsigned M, typename Field>
Field BaseMatrix<N, M, Field>::trace() const {
    ASSERT_A_IS_EQUAL_TO_B<N, M>();
    Field Trace = 0;
    for (unsigned i = 0; i < N; ++i) {
        Trace += a[i][i];
    }
    return Trace;
}

template<unsigned N, unsigned M, typename Field>
std::vector<Field> BaseMatrix<N, M, Field>::getRow(const unsigned& cnt) const {
    std::vector<Field> New;
    for (unsigned i = 0; i < M; ++i) {
        New.push_back(a[cnt][i]);
    }
    return New;
}

template<unsigned N, unsigned M, typename Field>
std::vector<Field> BaseMatrix<N, M, Field>::getColumn(const unsigned& cnt) const {
    std::vector<Field> New;
    for (unsigned i = 0; i < N; ++i) {
        New.push_back(a[i][cnt]);
    }
    return New;
}

template<unsigned N, unsigned M, typename Field>
Field*& BaseMatrix<N, M, Field>::operator [](const unsigned& cnt) {
    return a[cnt];
}

template<unsigned N, unsigned M, typename Field>
Field* BaseMatrix<N, M, Field>::operator [](const unsigned& cnt) const {
    return a[cnt];
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<M, N, Field> BaseMatrix<N, M, Field>::transposed() const {
    BaseMatrix<M, N, Field> New;
    for (unsigned j = 0; j < M; ++j) {
        for (unsigned i = 0; i < N; ++i) {
            New[j][i] = a[i][j];
        }
    }
    return New;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N / 2, N / 2, Field>** BaseMatrix<N, M, Field>::quater() const {
    BaseMatrix<N / 2, N / 2, Field>** Res = generateArray<2, 2, BaseMatrix<N / 2, N / 2, Field>>();
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < N; ++j) {
            if (i < N / 2 && j < N / 2) {
                Res[0][0][i][j] = a[i][j];
            }
            if (i >= N / 2 && j < N / 2) {
                Res[1][0][i - N / 2][j] = a[i][j];
            }
            if (i >= N / 2 && j >= N / 2) {
                Res[1][1][i - N / 2][j - N / 2] = a[i][j];
            }
            if (i < N / 2 && j >= N / 2) {
                Res[0][1][i][j - N / 2] = a[i][j];
            }
        }
    }
    return Res;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, N, Field> baseMultiply(const BaseMatrix<N, N, Field>& A, const BaseMatrix<N, N, Field>& B) {
    BaseMatrix<N, N, Field> Res;
    Res[0][0] = (A[0][0] * B[0][0]) + (A[0][1] * B[1][0]);
    Res[0][1] = (A[0][0] * B[0][1]) + (A[0][1] * B[1][1]);
    Res[1][0] = (A[1][0] * B[0][0]) + (A[1][1] * B[1][0]);
    Res[1][1] = (A[1][0] * B[0][1]) + (A[1][1] * B[1][1]);
    return Res;
}

template<unsigned N, typename Field>
BaseMatrix<N, N, Field> getMatrixFromQuater(BaseMatrix<N / 2, N / 2, Field>** A) {
    BaseMatrix<N, N, Field> Res;
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < N; ++j) {
            if (i < N / 2 && j < N / 2) {
                Res[i][j] = A[0][0][i][j];
            }
            if (i >= N / 2 && j < N / 2) {
                Res[i][j] = A[1][0][i - N / 2][j];
            }
            if (i >= N / 2 && j >= N / 2) {
                Res[i][j] = A[1][1][i - N / 2][j - N / 2];
            }
            if (i < N / 2 && j >= N / 2) {
                Res[i][j] = A[0][1][i][j - N / 2];
            }
        }
    }
    return Res;
}

template<unsigned N, typename Field>
BaseMatrix<N, N, Field> multiply(const BaseMatrix<N, N, Field>& matrixA, const BaseMatrix<N, N, Field>& matrixB) {
    if (N == 2) {
        return baseMultiply<N, N, Field>(matrixA, matrixB);
    }
    BaseMatrix<N, N, Field> Res;
    BaseMatrix<N / 2, N / 2, Field>** A = matrixA.quater();
    BaseMatrix<N / 2, N / 2, Field>** B = matrixB.quater();
    BaseMatrix<N / 2, N / 2, Field>** C = generateArray<2, 2, BaseMatrix<N / 2, N / 2, Field>>();
    BaseMatrix<N / 2, N / 2, Field>* P = new BaseMatrix<N / 2, N / 2, Field>[7];
    P[0] = multiply<N / 2, Field>((A[0][0] + A[1][1]), (B[0][0] + B[1][1]));
    P[1] = multiply<N / 2, Field>((A[1][0] + A[1][1]), B[0][0]);
    P[2] = multiply<N / 2, Field>(A[0][0], (B[0][1] - B[1][1]));
    P[3] = multiply<N / 2, Field>(A[1][1], (B[1][0] - B[0][0]));
    P[4] = multiply<N / 2, Field>((A[0][0] + A[0][1]), B[1][1]);
    P[5] = multiply<N / 2, Field>((A[1][0] - A[0][0]), (B[0][0] + B[0][1]));
    P[6] = multiply<N / 2, Field>((A[0][1] - A[1][1]), (B[1][0] + B[1][1]));
    C[0][0] = P[0] + P[3] - P[4] + P[6];
    C[0][1] = P[2] + P[4];
    C[1][0] = P[1] + P[3];
    C[1][1] = P[0] - P[1] + P[2] + P[5];
    Res = getMatrixFromQuater<N, Field>(C);
    deleteArray<2, 2, BaseMatrix<N / 2, N / 2, Field>>(A);
    deleteArray<2, 2, BaseMatrix<N / 2, N / 2, Field>>(B);
    deleteArray<2, 2, BaseMatrix<N / 2, N / 2, Field>>(C);
    delete[] P;
    return Res;
}

template<unsigned N, unsigned M, typename Field>
template<unsigned K>
BaseMatrix<N, K, Field> BaseMatrix<N, M, Field>::operator * (const BaseMatrix<M, K, Field>& other) const {
    using FindDegree::findDeg;
    const unsigned newSize = findDeg<std::max(N, std::max(M, K)), 2>::ans;
    BaseMatrix<newSize, newSize, Field> A, B, C;
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            A[i][j] = a[i][j];
        }
    }
    for (unsigned i = 0; i < M; ++i) {
        for (unsigned j = 0; j < K; ++j) {
            B[i][j] = other[i][j];
        }
    }
    C = multiply(A, B);
    BaseMatrix<N, K, Field> Res;
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < K; ++j) {
            Res[i][j] = C[i][j];
        }
    }
    return Res;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>& BaseMatrix<N, M, Field>::operator *= (const BaseMatrix<M, M, Field>& other) {
    *this = *this * other;
    return *this;
}


template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field>& BaseMatrix<N, M, Field>::invert() {
    ASSERT_A_IS_EQUAL_TO_B<N, M>();
    if (this->rank() < N) {
        std::cout << "can't be inverted";
        return *this;
    }
    BaseMatrix<N, N + N, Field> AE;
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < N + N; ++j) {
            if (j < N) AE[i][j] = a[i][j];
            else if (i == j - N) AE[i][j] = 1;
        }
    }
    int* used = new int[N];
    int* place = new int[N];
    for (unsigned i = 0; i < N; ++i) {
        used[i] = -1;
    }
    for (unsigned row = 0; row < N; ++row) {
        unsigned curRow = row;
        for (unsigned i = 0; i < N; i++) {
            if (used[i] == -1 && (abs(AE[i][row]) > abs(AE[curRow][row]))) {
                curRow = i;
            }
        }
        used[curRow] = row;
        place[row] = curRow;
        for (unsigned i = 0; i < N; ++i) {
            if (i != curRow) {
                Field c = AE[i][row] / AE[curRow][row];
                for (unsigned j = row; j < N + N; ++j) {
                    AE[i][j] -= AE[curRow][j] * c;
                }
            }
        }
        for (unsigned i = row + 1; i < N + N; i++) {
            AE[curRow][i] /= AE[curRow][row];
        }
        AE[curRow][row] = 1;
    }
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < N; ++j) {
            a[i][j] = AE[place[i]][j + N];
        }
    }
    delete[] used;
    delete[] place;
    return *this;
}

template<unsigned N, unsigned M, typename Field>
BaseMatrix<N, M, Field> BaseMatrix<N, M, Field>::inverted() const {
    BaseMatrix New = *this;
    New.invert();
    return New;
}

template<unsigned N, unsigned M, typename Field>
Field BaseMatrix<N, M, Field>::det() {
    ASSERT_A_IS_EQUAL_TO_B<N, M>();
    if (this->rank() < N) {
        return 0;
    }
    BaseMatrix<N, M, Field> Res = *this;
    Field determinant = 1;
    for (unsigned col = 0, row = 0; row < N && col < M; ++col) {
        unsigned curRow = row;
        for (unsigned i = row; i < N; i++) {
            if (abs(Res[i][col]) > abs(Res[curRow][col])) {
                curRow = i;
            }
        }
        if (row != curRow) {
            determinant *= -1;
        }
        for (unsigned j = col; j < M; ++j) {
            Field tmp = Res[curRow][j];
            Res[curRow][j] = Res[row][j];
            Res[row][j] = tmp;
        }
        determinant *= Res[row][col];
        for (unsigned i = col + 1; i < N; ++i) {
            Res[row][i] /= Res[row][col];
        }
        Res[row][col] = 1;
        for (unsigned i = 0; i < N; ++i) {
            if (i != row) {
                Field c = Res[i][col] / Res[row][col];
                for (unsigned j = col; j < M; ++j) {
                    Res[i][j] -= Res[row][j] * c;
                }
            }
        }
        ++row;
    }
    return determinant;
}

template<unsigned N, unsigned M, typename Field>
std::istream& operator >> (std::istream& in, BaseMatrix<N, M, Field>& object) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            in >> object[i][j];
        }
    }
    return in;
}

template<unsigned N, unsigned M, typename Field>
std::ostream& operator << (std::ostream& out, const BaseMatrix<N, M, Field>& object) {
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            out << object[i][j] << " ";
        }
        out << std::endl;
    }
    return out;
}

template<unsigned N, unsigned M, typename Field>
unsigned BaseMatrix<N, M, Field>::rank() const {
    const double EPS = 1e-9;
    BaseMatrix<N, M, Field> Res = *this;
    unsigned rk = 0;
    for (unsigned col = 0, row = 0; row < N && col < M; ++col) {
        unsigned curRow = row;
        for (unsigned i = row; i < N; i++) {
            if (abs(Res[i][col]) > abs(Res[curRow][col])) {
                curRow = i;
            }
        }
        if (abs(Res[curRow][col]) < EPS) {
            continue;
        }
        rk++;
        for (unsigned j = col; j < M; ++j) {
            Field tmp = Res[curRow][j];
            Res[curRow][j] = Res[row][j];
            Res[row][j] = tmp;
        }
        for (unsigned i = col + 1; i < M; ++i) {
            Res[row][i] /= Res[row][col];
        }
        Res[row][col] = 1;
        for (unsigned i = 0; i < N; ++i) {
            if (i != row) {
                Field c = Res[i][col] / Res[row][col];
                for (unsigned j = col; j < M; ++j) {
                    Res[i][j] -= Res[row][j] * c;
                }
            }
        }
        ++row;
    }
    return rk;
}

//............................SQUARE_MATRIX....................................................................

template<unsigned N, typename Field = Rational>
class SquareMatrix : public BaseMatrix<N, N, Field> {
public:
    SquareMatrix();
    ~SquareMatrix();
    SquareMatrix(const std::initializer_list<std::vector<Field>>& list);
    SquareMatrix(const BaseMatrix<N, N, Field>& other);
    SquareMatrix& operator = (const SquareMatrix<N, Field>& other);
};

template<unsigned N, typename Field>
SquareMatrix<N, Field>::SquareMatrix() {
    this->a = generateArray<N, N, Field>();
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < N; ++j) {
            this->a[i][j] = (i == j);
        }
    }
}

template<unsigned N, typename Field>
SquareMatrix<N, Field>::~SquareMatrix() {}

template<unsigned N, typename Field>
SquareMatrix<N, Field>::SquareMatrix(const std::initializer_list<std::vector<Field>>& list) {
    this->a = generateArray<N, N, Field>();
    unsigned i = 0;
    for (auto& elem : list) {
        for (unsigned j = 0; j < N; ++j) {
            this->a[i][j] = elem[j];
        }
        i++;
    }
}

template<unsigned N, typename Field>
SquareMatrix<N, Field>::SquareMatrix(const BaseMatrix<N, N, Field>& other) {
    this->a = generateArray<N, N, Field>();
    fillArray<N, N, Field>(this->a, other.a);
}

template<unsigned N, typename Field>
SquareMatrix<N, Field>& SquareMatrix<N, Field>::operator = (const SquareMatrix& other) {
    fillArray<N, N, Field>(this->a, other.a);
    return *this;
}

//...................................MATRIX.............................................................

template<unsigned N, unsigned M, typename Field>
class Matrix : public BaseMatrix<N, M, Field> {
public:
    Matrix();
    ~Matrix();
    Matrix(const std::initializer_list<std::vector<Field>>& list);
    Matrix(const BaseMatrix<N, M, Field>& other);
    Matrix& operator = (const Matrix<N, M, Field>& other);
};

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>::Matrix() {
    this->a = generateArray<N, M, Field>();
};

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>::Matrix(const std::initializer_list<std::vector<Field>>& list) {
    this->a = generateArray<N, M, Field>();
    unsigned i = 0;
    for (auto elem : list) {
        for (unsigned j = 0; j < N; ++j) {
            this->a[i][j] = elem[j];
        }
        i++;
    }
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>::Matrix(const BaseMatrix<N, M, Field>& other) {
    this->a = generateArray<N, M, Field>();
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            this->a[i][j] = other[i][j];
        }
    }
}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>::~Matrix() {}

template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator = (const Matrix<N, M, Field>& other) {
    fillArray<N, M, Field>(this->a, other.a);
    return *this;
}

//........................................FINITE_MATRIX..........................................................

template<unsigned N, unsigned M, unsigned P>
class Matrix<N, M, Finite<P>> : public BaseMatrix<N, M, Finite<P>> {
public:
    Matrix();
    ~Matrix();
    Matrix(const std::initializer_list<std::vector<Finite<P>>>& list);
    Matrix(const BaseMatrix<N, M, Finite<P>>& other);
    Matrix& operator = (const Matrix<N, M, Finite<P>>& other);
};

template<unsigned N, unsigned M, unsigned P>
Matrix<N, M, Finite<P>>::Matrix() {
    ASSERT_COMPOSITE<P>();
    this->a = generateArray<N, M, Finite<P>>();
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            this->a[i][j] = (i == j);
        }
    }
};

template<unsigned N, unsigned M, unsigned P>
Matrix<N, M, Finite<P>>::Matrix(const std::initializer_list<std::vector<Finite<P>>>& list) {
    this->a = generateArray<N, M, Finite<P>>();
    unsigned i = 0;
    for (auto elem : list) {
        for (unsigned j = 0; j < M; ++j) {
            this->a[i][j] = elem[j];
        }
        i++;
    }
}
template<unsigned N, unsigned M, unsigned P>
Matrix<N, M, Finite<P>>::Matrix(const BaseMatrix<N, M, Finite<P>>& other) {
    this->a = generateArray<N, M, Finite<P>>();
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            this->a[i][j] = other[i][j];
        }
    }
}

template<unsigned N, unsigned M, unsigned P>
Matrix<N, M, Finite<P>>::~Matrix() {}

template<unsigned N, unsigned M, unsigned P>
Matrix<N, M, Finite<P>>& Matrix<N, M, Finite<P>>::operator = (const Matrix<N, M, Finite<P>>& other) {
    fillArray<N, M, Finite<P>>(this->a, other.a);
    return *this;
}

//....................MAIN.......................................................................................


