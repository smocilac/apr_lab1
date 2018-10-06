#ifndef MATRICA_H
#define MATRICA_H

#include <cstdlib>
#include <utility>
#include <memory>
#include <vector>
#include <valarray>
#include <ostream>
#include <istream>
#include <sstream>
#include <fstream>


template <class T = double>
class Matrica
{
public:
    Matrica(const std::pair<size_t, size_t>&);
    Matrica(size_t x, size_t y);
    Matrica(const Matrica&);
    Matrica(const std::string&);

    ~Matrica();

    int LU_decompose();
    static Matrica forward_substitution(const Matrica<T> &A, const Matrica<T> &b);
    static Matrica backward_substitution(const Matrica<T> &A, const Matrica<T> &b);

    Matrica& operator=(const Matrica&);

    T& operator()(size_t x, size_t y);
    T operator()(size_t x, size_t y) const;

    Matrica& operator+=(T);
    Matrica operator+(T);
    Matrica& operator+=(const Matrica&);
    Matrica operator+(const Matrica&);

    Matrica& operator-=(T);
    Matrica operator-(T);
    Matrica& operator-=(const Matrica&);
    Matrica operator-(const Matrica&);

    Matrica& operator*=(T);
    Matrica operator*(T);
    Matrica operator*(const Matrica&);

    bool operator==(const Matrica&);

    Matrica operator~();


    template <typename U>
    friend std::ostream& operator<< (std::ostream& os, const Matrica<U>& mat);

    template <typename U>
    friend std::istream& operator>> (std::istream& is, Matrica<U>& mat);

    std::valarray<T> &array() {return values;}
    std::valarray<T> array() const {return values;}

    std::pair<size_t, size_t> get_dimensions() const{return std::pair<size_t, size_t>(m_x, m_y);}

    size_t row_size() const { return m_y;}
    size_t col_size() const { return m_x;}



private:

    std::valarray<T> values;

    size_t m_x;
    size_t m_y;
};


template<class T>
Matrica<T>::Matrica(const std::pair<size_t, size_t> &dimension) : m_x(dimension.first), m_y(dimension.second), values(std::valarray<T>(dimension.first * dimension.second))
{

}


template<class T>
Matrica<T>::Matrica(size_t x, size_t y) : m_x(x), m_y(y), values(std::valarray<T>(x * y))
{

}


template<class T>
Matrica<T>::Matrica(const Matrica<T> &mat) : m_x(mat.col_size()), m_y(mat.row_size()), values(std::valarray<T>(mat.array()))
{

}


template<class T>
Matrica<T>::Matrica(const std::string &file_name)
{
    std::vector<T> vals;
    std::ifstream input_file;
    input_file.open (file_name);
    std::string line;
    size_t row = 0;
    while (std::getline(input_file, line))
    {
        std::istringstream iss(line);
        T val;
        while(iss >> val){
            vals.push_back(val);
        }

        if (!row)
            row = vals.size();
    }

    values = std::valarray<T>(vals.data(), vals.size());

    m_x = row;
    m_y = values.size() / row;
}


template<class T>
Matrica<T>::~Matrica()
{

}


template <class T>
int Matrica<T>::LU_decompose(){
    for (size_t i = 0; i < m_y - 1; i++){
        for (size_t j = i + 1; j < m_y; j++){
            (*this)(i, j) /= (*this)(i, i);
            for (size_t k = i + 1; k < m_x; k++){
                (*this)(k, j) -= (*this)(i,j) * (*this)(k,i);
            }
        }
    }
}

template <class T>
Matrica<T> Matrica<T>::forward_substitution(const Matrica<T> &A, const Matrica<T> &b){
    Matrica<T> y(b);

    for (size_t i = 0; i < y.row_size() - 1; i++){
        for (size_t j = i + 1; j < y.row_size() ; j++){
            y(0, j) -= A(i, j) * y(0, i);
        }
    }
    return y;
}


template <class T>
Matrica<T> Matrica<T>::backward_substitution(const Matrica<T> &A, const Matrica<T> &y){
    Matrica<T> x(y);

    for (int i = x.row_size() - 1; i >= 0; i--){

        for (int j = i + 1; j < x.row_size(); j++){
            x(0, i) -= A(j, i) * x(0, j);
        }
        x(0, i) /= A(i, i);

    }
    return x;
}


template<class T>
Matrica<T> &Matrica<T>::operator=(const Matrica<T> &mat)
{
    values = std::valarray<T>(mat.array());

    return *this;
}


template<class T>
T &Matrica<T>::operator()(size_t x, size_t y)
{
    return values[y * m_x + x];
}


template<class T>
T Matrica<T>::operator()(size_t x, size_t y) const
{
    return values[y * m_x + x];
}


template<class T>
Matrica<T> &Matrica<T>::operator+=(T t)
{
    values += t;

    return *this;
}


template<class T>
Matrica<T> Matrica<T>::operator+(T t)
{
    Matrica<T> mat(*this);
    mat.array() += t;
    return mat;
}


template<class T>
Matrica<T> &Matrica<T>::operator+=(const Matrica<T> &mat)
{
    values += mat.array();

    return *this;
}


template<class T>
Matrica<T> Matrica<T>::operator+(const Matrica<T> &mat)
{
    Matrica<T> mat_out(*this);

    mat_out.values += mat.values;

    return mat_out;
}


template<class T>
Matrica<T> &Matrica<T>::operator-=(T t)
{
    values -= t;

    return *this;
}


template<class T>
Matrica<T> Matrica<T>::operator-(T t)
{
    Matrica<T> mat(*this);

    mat -= t;

    return mat;
}


template<class T>
Matrica<T> &Matrica<T>::operator-=(const Matrica<T> &mat)
{
    values -= mat.array();

    return *this;
}


template<class T>
Matrica<T> Matrica<T>::operator-(const Matrica<T> &mat)
{
    Matrica<T> mat_out(*this);

    mat_out.array() -= mat.array();

    return mat_out;
}


template<class T>
Matrica<T> &Matrica<T>::operator*=(T t)
{
    values *= t;

    return *this;
}


template<class T>
Matrica<T> Matrica<T>::operator*(T t)
{
    Matrica<T> mat(*this);

    mat.array() *= t;

    return mat;
}



//    a b c      q w       (a*q+b*z+c*v)  (a*w+b*x+c*b)
//    d e f      z x   =   (d*q+e*z+f*v)  (d*w+e*x+f*b)
//               v b
template<class T>
Matrica<T> Matrica<T>::operator*(const Matrica<T> &mat)
{
    Matrica<T> mat_out(mat.m_x, m_y);

    for (int k = 0; k < m_y; k ++){
        for (size_t i = 0; i < mat.m_x; i++){
            T val = 0;
            for (size_t j = 0; j < mat.m_y; j++){
                val += (*this)(j, k) * mat(i, j);
            }
            mat_out(i, k) = val;
        }
    }

    return mat_out;
}

template <class T>
Matrica<T> Matrica<T>::operator~(){
    Matrica<T> mat(m_y, m_x);
    for (size_t i = 0; i < m_y; i++){
        for (size_t j = 0; j < m_x; j++){
            mat(i, j) = (*this)(j, i);
        }
    }

    return mat;
}


template<class T>
bool Matrica<T>::operator==(const Matrica<T> &mat)
{
    return values == mat.array();
}


template <class T>
std::ostream& operator<<(std::ostream& os, const Matrica<T>& mat){
    for (size_t i = 0; i < mat.m_y; i++){
        for (size_t j = 0; j < mat.m_x; j++){
            os << mat.values[i * mat.m_x + j];
            if (j < mat.m_x - 1)
                os << " ";
        }
        os << std::endl;
    }

    return os;
}


template <class T>
std::istream& operator>> (std::istream& is, Matrica<T>& mat){
    for (size_t i = 0; i < mat.m_y; i++){
        for (size_t j = 0; j < mat.m_x; j++){
            is >> mat(j, i);
        }
    }

    return is;
}


#endif // MATRICA_H
