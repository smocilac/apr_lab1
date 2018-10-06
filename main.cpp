#include <iostream>
#include <valarray>
#include <fstream>

#include "matrica.hpp"

using namespace std;

int main()
{

    Matrica<> mat1("mat1.txt");
    Matrica<> b("b.txt");
    mat1.LU_decompose();
    cout << Matrica<double>::backward_substitution(mat1, Matrica<>::forward_substitution(mat1, b));
    //cout << mat1;

//    Matrica<> m(4, 4);
//    Matrica<> m2(4, 4);


//    m2 += 3;
//    m += 2;
//    m2 *= 3;

//    cout << m2 ;

//    ifstream myfile;
//    myfile.open ("example.txt");
//    myfile >> m;

//    Matrica<> m5("example.txt");
//    cout << m5;

//    m3 = ~m5;
//    std::cout << endl;

//    cout << m3;

//    ofstream myOutFile;
//    myOutFile.open("stored_matrix.txt");
//    myOutFile << m5;

    return 0;
}
