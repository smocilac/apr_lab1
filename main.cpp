#include <iostream>
#include <valarray>
#include <fstream>

#include "matrica.hpp"

using namespace std;

/* testa funkcija, demonstrira samo prvi zadatak iz vjezbe */
void zadatak1();

/* testna funkcija, usporeduje LU i LUP algoritam */
void rijesi_sustav(std::string file_A, std::string file_b);

/* teste funckije za LU i LUP */
void test_using_LU(Matrica<> A, Matrica<> b);
void test_using_LUP(Matrica<> A, Matrica<> b);


int main()
{
    zadatak1();

    //rijesi_sustav("../matrice/mat1_zad2.txt", "../matrice/mat2_zad2.txt");
    //rijesi_sustav("../matrice/mat1_zad3.txt", "../matrice/mat2_zad3.txt");
    //rijesi_sustav("../matrice/mat1_zad4.txt", "../matrice/mat2_zad4.txt"); // TODO: LUP daje tocne vrijednost bez decimala.Zasto?
    //rijesi_sustav("../matrice/mat1_zad5.txt", "../matrice/mat2_zad5.txt"); // TODO: zasto nije NUlA?
    //rijesi_sustav("../matrice/mat1_zad6.txt", "../matrice/mat2_zad6.txt"); // TODO: oboje radi, moze imati problema. Zasto?

    return 0;
}


void zadatak1(){
    Matrica<double> mat1("../matrice/mat1_zad2.txt");

    Matrica<double> mat = mat1 * 3;
    mat /= 3;

    std::cout << mat << mat1;
    std::cout << (mat == mat1) << std::endl;
}


void rijesi_sustav(std::string file_A, std::string file_b){
    Matrica<double> A(file_A);
    Matrica<double> b(file_b);

    cout << "--------------------LU-----------------------" << endl;
    test_using_LU(A, b);
    cout << endl <<  "--------------------LUP-----------------------" << endl;
    test_using_LUP(A, b);
}


void test_using_LU(Matrica<> A, Matrica<> b){
    cout << "A = " << std::endl << A;
    cout << endl << "b = " << std::endl << b;

    if (A.LU_decompose()){
        std::cerr << "WARNING: could not decompose using LU algorithm, solution may not be valid." << std::endl;
    }

    cout << endl << "After LU: " << std::endl << A;

    Matrica<> y = Matrica<>::forward_substitution(A, b);
    cout << endl << "Forward substitution vector: " << std::endl << y;

    cout << endl << "Backward substitution vector: " << std::endl;

    Matrica<> x = Matrica<double>::backward_substitution(A, y);
    if (x.array().size() == 0)
        std::cerr << "WARNING: Backward supstitution :: divide by zero !" << std::endl;
    else
        std::cout << endl << Matrica<double>::backward_substitution(A, y);
}


void test_using_LUP(Matrica<> A, Matrica<> b){
    cout << "A = " << std::endl << A;
    cout << endl << "b = " << std::endl << b;

    if (A.LUP_decompose()){
        std::cerr << "WARNING: could not decompose using LUP algorithm, solution may not be valid." << std::endl;
    }

    cout << endl << "After LUP: " << std::endl << A;
    cout << endl << "Permutation matrix: " << std::endl << Matrica<>::get_permutations(A);
    b = Matrica<>::get_permutations(A) * b;
    Matrica<> y = Matrica<>::forward_substitution(A, b);
    cout << endl << "Forward substitution vector: " << std::endl << y;

    Matrica<> x = Matrica<double>::backward_substitution(A, y);

    if (x.array().size() == 0)
        std::cerr << "WARNING: Backward supstitution :: divide by zero !" << std::endl;
    else
        cout << endl << "Backward substitution vector: " << std::endl << x;
}

