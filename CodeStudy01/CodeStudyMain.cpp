#include <iostream>
#include <fstream>
#include <time.h>

#include "Matrix.hpp"

using namespace std;

using namespace LinearAlgebra;

void subfunc(ofstream& o);

int main()
{
   ofstream out("Matrix.out");
   out.precision(2);

   subfunc(out);

   out.close();
   cout << "done!" << endl;
}

void subfunc(ofstream& o)
{
   Vector V({ 1, 1, 1 });

   Matrix Y({ Vector({ 1, 0, 0}),
              Vector({ 0, 1, 0}),
              Vector({ 0, 0, 1}) });

   Matrix Z({ Vector({ 0,  0,  1}),
              Vector({ 2, -1,  3}),
              Vector({ 1,  1,  4}) });

   cout << V << endl;
   cout << endl;
   cout << Y << endl;
   cout << endl;
   cout << Z << endl;
   cout << V * Z << endl;
   cout << endl;
   cout << Z * V << endl;
   cout << endl;
   cout << Z * Y << endl;
   cout << endl;
   cout << det(Z) << endl;
   cout << endl;

   V.Resize(2);
   Y.Resize(2, 2);
   Z.Resize(2, 2);

   Z(0, 0) = 1;
   Z(0, 1) = 1;
   Z(1, 0) = 3;
   Z(1, 1) = -1;
   V(0) = 2;
   V(1) = 2;

   cout << Z << endl;
   cout << V << endl;
   cout << endl;
   Vector U = V / Z;
   cout << U << endl;
   cout << endl;

   int N;
   cout << "Enter Matrix Size for Computation = ";
   cin >> N;

   o << "Matrix Size for Computation = " << N << endl << endl;

   Matrix A(N, N), C(N, N);
   Vector B(N), X(N);

   A.MakeZero(); B.MakeZero();
   // Set the matrix randomly as band matrix with 5 entries in a row
   A.RandomizeBand(-1.0, 5);
   B.Randomize(-1.0, 0.6);
   A.setprecision(5);
   B.setprecision(5);

   A(2, 6) = 1.5;

   if ( N <= 10 ) {
      o << A << endl;
      o << B << endl;
   }

   Matrix Acopy(A);
   Vector Bcopy(B);

   clock_t tStart, tFinish;

   cout << "Find solution of A*X = B using Gauss Elimination" << endl;
   o << "Find solution of A*X = B using Gauss Elimination" << endl;
   tStart = clock();
   X = B / A;
   tFinish = clock();
   X.setprecision(5);
   o << "The solution took " << tFinish-tStart << " cpu ticks spent, or about "
             << (tFinish-tStart)/float(CLOCKS_PER_SEC) << " seconds." << endl;
   o << "Infinity Norm of A*X - B = " << scientific << InfNorm(Acopy*X-Bcopy) << endl;
   o << X << endl;
   o << endl;

   cout << endl;

   cout << "Set a band matrix in a loop" << endl;
   N = 6;
   A.Resize(N, N);
   X.Resize(N);

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
         if (i == j)
            A(i, j) = -2;
         else if (j == i+1)
            A(i, j) = 1;
         else if (j == i-1)
            A(i, j) = 1;
      }
   }

   X = -1;

   cout << A;
   cout << "The right hand side vector:" << endl;
   cout << X << endl;
   cout << "The solution vector:" << endl;
   cout << X/A << endl;

   std::ifstream infile;
   infile.open("Matrix.dat");

   cout << endl;
   cout << "Reading the matrix from a file" << endl;
   infile >> N;

   A.Resize(N, N);

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
         infile >> A(i, j);
      }
   }

   cout << A << endl;

   cout << endl;
}
