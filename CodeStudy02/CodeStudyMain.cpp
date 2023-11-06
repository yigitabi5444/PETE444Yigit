#include <iostream>
#include <math.h>

#include "Matrix.hpp"

using namespace std;
using namespace LinearAlgebra;
// An arbitrary function
double func(double x) {
   double y = x * sin(x);
   return y;
}

// Derivative of the function above
double dfunc(double x) {
   double y = sin(x) + x * cos(x);
   return y;
}

int main()
{
   double a = -10;
   double b = 10;

   int n = 100;

   Vector X(n+1);
   Vector f(n+1);
   Vector dfa(n+1); // Analytical derivative
   Vector dfnf(n+1); // Numerical derivative by forward difference
   Vector dfnc(n+1); // Numerical derivative by central difference

   double dx = (b - a) / n; // uniform incremental length

   for (int i = 0; i <= n; i++) {
      double x = a + dx * i;
      X(i) = x;
      f(i) = func(x);   // Calculate and store function value
      dfa(i) = dfunc(x); // Calculate and store analytical derivative of function
   }

   // Calculate derivative using forward difference
   for (int i = 0; i < n; i++) {
      dfnf(i) = (f(i + 1) - f(i)) / dx;
      cout << X(i) << " " << dfnf(i) << " " << dfa(i) << endl;
   }

   cout << endl;
   cout << endl;

   // Calculate derivative using central difference
   dfnc(0) = (f(1) - f(0)) / dx;
   for (int i = 1; i < n; i++) {
      dfnc(i) = (f(i + 1) - f(i-1)) / (2.0*dx);
      cout << X(i) << " " << dfnc(i) << " " << dfa(i) << endl;
   }
}
