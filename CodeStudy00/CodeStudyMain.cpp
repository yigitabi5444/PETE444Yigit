#include <iostream>
#include <fstream>

#include "Matrix.hpp"

using namespace std;

using namespace LinearAlgebra;

int main()
{
   cout << "Hello world" << endl;

   Matrix A(10, 10);

   cout << "I'm a random matrix:" << endl;
   A.Randomize(-2.0, 3.0, 0.9);

   A(0, 0) = -1;
   A(0, 3) = -1;
   A[1][1] = -1;

   cout << A;
}
