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
    Matrica(const std::valarray<T> &va, size_t x, size_t y) ;
    Matrica(const std::pair<size_t, size_t>&);
    Matrica(size_t x, size_t y);
    Matrica(const Matrica&);
    Matrica(const std::string&);

    ~Matrica();

    void init_permutations();
    void replace_rows(size_t, size_t);
    void replace_columns(size_t, size_t);

    int LU_decompose();
    int LUP_decompose();
    size_t find_max_below_elem(size_t, size_t);
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


    Matrica& operator/=(T);
    Matrica operator/(T);
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

    static Matrica<T> get_permutations(const Matrica<T> &matrica);

    std::pair<size_t, size_t> get_dimensions() const{return std::pair<size_t, size_t>(m_x, m_y);}

    size_t row_size() const { return m_y;}
    size_t col_size() const { return m_x;}

    static const T epsilon;

private:

    std::valarray<T> values;
    std::valarray<T> permutations;

    size_t m_x;
    size_t m_y;
};


/**
 * Proizvoljno veliki epsilon. Defaultna vrijednost 1e-10
 */
template<class T>
const T Matrica<T>::epsilon = 0.0000000001;


/**
 *  Funkcija dohvaca permutacijsku matricu.
 */
template<class T>
Matrica<T> Matrica<T>::get_permutations(const Matrica<T> &matrica){
    return Matrica<T>(matrica.permutations, matrica.m_x, matrica.m_y);
}


/**
 * Konstruktor, prima inicijalno polje te x i y dimenziju matrice.
 * Korisnik se mora osigurati da je ulazno polje (barem) velicine x * y
 */
template<class T>
Matrica<T>::Matrica(const std::valarray<T> &va, size_t x, size_t y) : m_x(x), m_y(y), values(va), permutations(x*y)
{
    init_permutations();
}


/**
 * Konstruktor, prima x i y dimenziju matrice.
 * Matrica se inicijalizira nulama.
 */
template<class T>
Matrica<T>::Matrica(const std::pair<size_t, size_t> &dimension) : permutations(std::valarray<T>(dimension.first * dimension.second)), m_x(dimension.first), m_y(dimension.second), values(std::valarray<T>(dimension.first * dimension.second))
{
    init_permutations();
}


/**
 * Konstruktor, prima x i y dimenziju matrice.
 * Matrica se inicijalizira nulama.
 */
template<class T>
Matrica<T>::Matrica(size_t x, size_t y) : m_x(x), m_y(y), values(std::valarray<T>(x * y)),permutations(std::valarray<T>(x * y))
{
    init_permutations();
}


/**
 * Copy constructor, kao argument prima matricu cije vrijednosti i permutacijsku tablicu kopira.
 * Napomena: kopiraju se samo vrijednosti, ne i permutacijska matrica!
 */
template<class T>
Matrica<T>::Matrica(const Matrica<T> &mat) : m_x(mat.col_size()), m_y(mat.row_size()), values(std::valarray<T>(mat.array())), permutations(std::valarray<T>(mat.array().size()))
{
    init_permutations();
}


/**
 * Konstruktor, argument je path do filea koji treba otvoriti.
 * Vrijednosti unutar filea moraju biti formata:
 * x11 x12 ... x1n \n
 *       ....
 * xm1 xm2 ... xmn \n
 *
 */
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
    permutations = std::valarray<T>(vals.size());
    m_x = row;
    m_y = values.size() / row;

    init_permutations();
}


template<class T>
Matrica<T>::~Matrica()
{

}


/**
 * Funkcija koja inicijalizira permutacijsku matricu, odnosno na dijagonalama postavlja jedinice.
 * Sve vrijednost van glavne dijagonale su nula.
 */
template<class T>
void Matrica<T>::init_permutations(){
    for(size_t i = 0, j = 0; i < m_x && j < m_y; i++, j++)
        permutations[j*m_x + i] = 1;
}


/**
 * Zamjenjuje redak r1 sa redkom r2 u matrici. Azurira permutacijsku matricu sukladno tome.
 */
template<class T>
void Matrica<T>::replace_rows(size_t r1, size_t r2){
    if (r1 >= m_y || r2 >= m_y || r1 == r2)
        return;

    for (size_t i = 0; i < m_x; i++){
        // replace permutation matrix rows
        T tmp = permutations[r1*m_x + i];
        permutations[r1*m_x + i] = permutations[r2*m_x + i];
        permutations[r2*m_x + i] = tmp;
        // replace values matrix rows
        tmp = values[r1*m_x + i];
        values[r1*m_x + i] = values[r2*m_x + i];
        values[r2*m_x + i] = tmp;
    }
}


/**
 * Zamjenjuje stupac c1 sa stupcom c2 u matrici. Azurira permutacijsku matricu sukladno tome.
 */
template<class T>
void Matrica<T>::replace_columns(size_t c1, size_t c2){
    if (c1 >= m_x || c2 >= m_x)
        return;

    for (size_t i = 0; i < m_x; i++){
        // replace permutation matrix rows
        T tmp = permutations[i*m_x + c1];
        permutations[i*m_x + c1] = permutations[i*m_x + c2];
        permutations[i*m_x + c2] = tmp;
        // replace values matrix rows
        tmp = values[i*m_x + c1];
        values[i*m_x + c1] = values[i*m_x + c2];
        values[i*m_x + c2] = tmp;
    }
}


/**
 * Funkcija radi LU dekompoziciju matrice.
 */
template <class T>
int Matrica<T>::LU_decompose(){
    int is_error = 0;

    for (size_t i = 0; i < m_y - 1; i++){
        for (size_t j = i + 1; j < m_y; j++){
            // divide with zero flag
            if (std::abs((*this)(i, i)) <= Matrica<T>::epsilon)
                is_error = 1;

            (*this)(i, j) /= (*this)(i, i);
            for (size_t k = i + 1; k < m_x; k++){
                (*this)(k, j) -= (*this)(i,j) * (*this)(k,i);
            }
        }
    }

    return is_error;
}


/**
 * Funkcija radi LUP dekompoziciju matrice.
 */
template <class T>
int Matrica<T>::LUP_decompose(){
    int is_error = 0;

    for (size_t i = 0; i < m_y - 1; i++){

        size_t max_row = find_max_below_elem(i, i);
        replace_rows(i, max_row);

        for (size_t j = i + 1; j < m_y; j++){

            // divide with zero flag
            if (std::abs((*this)(i, i)) <= Matrica<T>::epsilon)
                is_error = 1;

            (*this)(i, j) /= (*this)(i, i);
            for (size_t k = i + 1; k < m_x; k++){
                (*this)(k, j) -= (*this)(i,j) * (*this)(k,i);
            }
        }
    }

    return is_error;
}


/**
 * Funkcija pronalazi najveci element u stupcu x ispod redka y.
 */
template <class T>
size_t Matrica<T>::find_max_below_elem(size_t x, size_t y){

    if (x >= m_x || y >= m_y)
        return m_y;
    size_t row_index = y;

    T _max = std::abs((*this)(x, y));

    for (size_t i = x, j = y + 1; j < m_y & i < m_x; j++){
        if (std::abs((*this)(x, j)) > _max){
            _max = std::abs((*this)(x, j));
            row_index = j;
        }
    }
    return row_index;
}


/**
 * Unaprijedna supstitucija (koristi se dio matrice A ispod glavne dijagonale.)
 */
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


/**
 * Supstitucija unatrag (koristi se dio matrice A iznad glavne dijagonale i glavna dijagonala)
 */
template <class T>
Matrica<T> Matrica<T>::backward_substitution(const Matrica<T> &A, const Matrica<T> &y){
    Matrica<T> x(y);

    for (int i = x.row_size() - 1; i >= 0; i--){

        for (int j = i + 1; j < x.row_size(); j++){
            x(0, i) -= A(j, i) * x(0, j);
        }

        if (std::abs(A(i, i)) <= Matrica<T>::epsilon)
            return Matrica<>(0, 0);

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


/**
 * indeksira elemente matrice u 2d prostoru.
 */
template<class T>
T &Matrica<T>::operator()(size_t x, size_t y)
{
    return values[y * m_x + x];
}


/**
 * indeksira elemente matrice u 2d prostoru (const)
 */
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
Matrica<T> &Matrica<T>::operator/=(T t)
{
    values /= t;

    return *this;
}


template<class T>
Matrica<T> Matrica<T>::operator/(T t)
{
    Matrica<T> mat(*this);

    mat.array() /= t;

    return mat;
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


/**
 * Mnozenje matrica (kartezijev produkt)
 *
 *    a b c      q w       (a*q+b*z+c*v)  (a*w+b*x+c*b)
 *    d e f      z x   =   (d*q+e*z+f*v)  (d*w+e*x+f*b)
 *               v b
 */
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
    // valarray ima interno definiran komparator nad svim elementima.
    return (values == mat.array()).min();

    // Ako bude potrebe moze se usporedba svesti na ovisnost o nekom nasem epsilonu
    /*auto mat_array_ptr = mat.array();

    if (mat_array_ptr.size() != values.size())
        return false;

    for ( int i = 0 ; i < values.size(); i++){
        if (std::abs(values[i] - mat_array_ptr) > Matrica<T>::epsilon)
            return false;
    }
    return true;*/
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
