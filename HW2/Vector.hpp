#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <iostream>
#include <vector>

#include "gendefs.h"

/*
Copyright (c) 1996-2022, Ismail Durgut

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimers.
* Redistributions in binary form must reproduce the above copyright notice in
the documentation and/or other materials provided with the distribution.
* The name of Ismail Durgut may not be used to endorse or promote
products derived from this Software without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

namespace LinearAlgebra {

   class Matrix;

   enum VecError
   {
      NOERROR, NOMEM, IDERR, NEVS
   }; // VecError

   class Vector
   {
      friend class Matrix;
   private:
      long N;
      std::vector<double> Array;
      friend std::ostream& operator<<(std::ostream& s, const Vector& V);
      friend std::istream& operator>>(std::istream& s, Vector& V);
      unsigned widthPrint, precisionPrint;
      bool bColVec;

   public:
      Vector();
      Vector(long n);
      Vector(long n, const double t);
      Vector(const Vector& Source);
      Vector(const std::vector<double>& Source);

      ~Vector();

      void Setup(long n);
      void Setup(long n, const double t);
      void Setup(const Vector& Source);
      void Resize(long n);
      void Assign(long n, const double t=0.0);

      void SetColumnVector() { bColVec = true; }
      void SetRowVector() { bColVec = false; }

      inline std::vector<double>::pointer Data() { return Array.data(); }

      long Size() const { return N; }

      double& operator[](long j);
      const double& operator[](long j) const;

      double& operator()(long j);
      const double& operator()(long j) const;

      Vector& operator=(const Vector& Source);
      Vector& operator=(const std::vector<double>& Source);
      Vector& operator=(const double& t);
      Vector& operator+=(const Vector& V);
      Vector& operator-=(const Vector& V);
      Vector& operator*=(const double& t);

      Vector& Randomize(const double& range, const double sparsity=1.0, int seed=0);
      Vector& Randomize(const double& range, long j1, long j2);
      Vector& Scale(Vector& minvec, Vector& maxvec);
      Vector& Upscale(Vector& minvec, Vector& maxvec);
      Vector DownSample(long n);
      double Distance(const Vector& V);
      double Difference(const Vector& V);
      void Swap(long j1, long j2) { double x = THIS[j1]; THIS[j1] = THIS[j2]; THIS[j2] = x; }

      inline unsigned setwidth(int iW) { unsigned pW = widthPrint; widthPrint = iW; return pW; }
      inline unsigned setprecision(int iP) { unsigned pP = precisionPrint; precisionPrint = iP; return pP; }

      void MakeZero();

      void Error(VecError vErr, const char* Msg = "") const;
   };

   Vector operator-(const Vector& V);
   Vector operator+(const Vector& V);
   Vector operator-(const Vector& V, const Vector& W);
   Vector operator+(const Vector& V, const Vector& W);
   Vector operator*(const double& t, const Vector& V);
   Vector operator*(const Vector& V, const double& t);
   Vector operator/(const Vector& V, const double& t);
   double operator*(const Vector& V, const Vector& W);
   Vector operator+(const Vector& V, const double& S);
   Vector operator+(const double& S, const Vector& V);
   Vector operator&(const Vector& V, const Vector& W);
   Vector operator&(const double& S, const Vector& V);
   Vector operator&(const Vector& V, const double& S);
   double EuclidNorm(const Vector& V);
   double CityNorm(const Vector& V);
   double MaxNorm(const Vector& V);
   double InfNorm(const Vector& V);
   double Sum(const Vector& V);
   long MaxMod(const Vector& V);

}  // namespace LinearAlgebra

#endif
