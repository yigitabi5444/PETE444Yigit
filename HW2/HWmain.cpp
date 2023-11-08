#include <iostream>
#include <math.h>
#include <fstream>

#include "Matrix.hpp"

using namespace std;
using namespace LinearAlgebra;

// Function to calculate f(x)
double func(double x) {
   double y = x * acos(sin(x));
   return y;
}

// Derivative of the function above
double dfunc(double x) {
   double y = acos(sin(x)) - x * sin(x) / sqrt(1 - pow(sin(x), 2));
   return y;
}

// Second derivative of the function above
double d2func(double x) {
   double y = -2 * cos(x) / sqrt(1 - pow(sin(x), 2)) - x * cos(x) / pow(1 - pow(sin(x), 2), 1.5) - pow(sin(x), 2) / pow(1 - pow(sin(x), 2), 1.5);
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
   Vector dfnf2(n+1); // Numerical derivative by central difference with 2nd order accuracy
   Vector dfnf4(n+1); // Numerical derivative by central difference with 4th order accuracy

   double dx = (b - a) / n; // uniform incremental length

   // Calculate function values and analytical derivative
   for (int i = 0; i <= n; i++) {
      double x = a + dx * i;
      X(i) = x;
      f(i) = func(x);   // Calculate and store function value
      dfa(i) = dfunc(x); // Calculate and store analytical derivative of function
   }

   // Calculate derivative using central difference with 2nd order accuracy
   dfnf2(0) = (f(1) - f(0)) / dx;
   dfnf2(n) = (f(n) - f(n-1)) / dx;
   for (int i = 1; i < n; i++) {
      dfnf2(i) = (f(i + 1) - f(i-1)) / (2.0*dx);
   }

   // Calculate derivative using central difference with 4th order accuracy
   dfnf4(0) = (-25*f(0) + 48*f(1) - 36*f(2) + 16*f(3) - 3*f(4)) / (12*dx);
   dfnf4(1) = (-3*f(0) - 10*f(1) + 18*f(2) - 6*f(3) + f(4)) / (12*dx);
   dfnf4(n-1) = (3*f(n) + 10*f(n-1) - 18*f(n-2) + 6*f(n-3) - f(n-4)) / (12*dx);
   dfnf4(n) = (25*f(n) - 48*f(n-1) + 36*f(n-2) - 16*f(n-3) + 3*f(n-4)) / (12*dx);
   for (int i = 2; i < n-1; i++) {
      dfnf4(i) = (-f(i+2) + 8*f(i+1) - 8*f(i-1) + f(i-2)) / (12*dx);
   }

   // Calculate second derivative using central difference with 2nd order accuracy
   Vector d2fnf2(n+1);
   d2fnf2(0) = (dfnf2(1) - dfnf2(0)) / dx;
   d2fnf2(n) = (dfnf2(n) - dfnf2(n-1)) / dx;
   for (int i = 1; i < n; i++) {
      d2fnf2(i) = (dfnf2(i + 1) - 2*dfnf2(i) + dfnf2(i-1)) / pow(dx, 2);
   }

   // Calculate second derivative using central difference with 4th order accuracy
   Vector d2fnf4(n+1);
   d2fnf4(0) = (dfnf4(2) - 2*dfnf4(1) + dfnf4(0)) / pow(dx, 2);
   d2fnf4(1) = (dfnf4(3) - 2*dfnf4(2) + dfnf4(1)) / pow(dx, 2);
   d2fnf4(n-1) = (dfnf4(n-3) - 2*dfnf4(n-2) + dfnf4(n-1)) / pow(dx, 2);
   d2fnf4(n) = (dfnf4(n-2) - 2*dfnf4(n-1) + dfnf4(n)) / pow(dx, 2);
   for (int i = 2; i < n-1; i++) {
      d2fnf4(i) = (-dfnf4(i+2) + 16*dfnf4(i+1) - 30*dfnf4(i) + 16*dfnf4(i-1) - dfnf4(i-2)) / (12*pow(dx, 2));
   }

   // Output results to file
   ofstream outfile;
   outfile.open("results.txt");
   outfile << "x\tf(x)\tdf(x)\tdf(x) (2nd order)\tdf(x) (4th order)\td2f(x) (2nd order)\td2f(x) (4th order)\n";
   for (int i = 0; i <= n; i++) {
      outfile << X(i) << "\t" << f(i) << "\t" << dfa(i) << "\t" << dfnf2(i) << "\t" << dfnf4(i) << "\t" << d2fnf2(i) << "\t" << d2fnf4(i) << "\n";
   }
   outfile.close();

   return 0;
}
