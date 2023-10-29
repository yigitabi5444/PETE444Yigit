#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <iostream>

#include "gendefs.h"

#include "Vector.hpp"

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

   class Matrix
   {
   private:
      long NumRow, NumCol;       // number of rows and columns
      std::vector<Vector> Rows;

      unsigned widthPrint, precisionPrint;
      bool bPrintIndex, bUsePivoting;
      bool bPrintNonZerosOnly;

      void RowSwap(long i1, long i2);

   public:
      Matrix();
      Matrix(long nrow, long ncol);
      Matrix(Matrix& Source);
      Matrix(const Matrix& Source);
      Matrix(const std::vector<Vector>& Source);

      ~Matrix();

      void Setup(long nrow, long ncol);
      void Setup(const Matrix& Source);

      void Resize(long nrow, long ncol);

      double& operator()(long i, long j);

      Vector& operator[](long i);
      const Vector& operator[](long i) const;

      friend Matrix operator-(Matrix& A);
      friend Matrix operator+(Matrix& A);
      friend Matrix operator+(Matrix& A, Matrix& B);
      friend Matrix operator-(Matrix& A, Matrix& B);
      friend Vector operator*(Vector& V, Matrix& A);
      friend Vector operator*(Matrix& A, Vector& V);
      friend Matrix operator*(Matrix& A, Matrix& B);
      friend Matrix operator*(double& t, Matrix& A);
      friend Matrix operator*(Matrix& A, double& t);
      friend Vector operator/(Vector& B, Matrix& A);
      //friend Matrix operator!(Matrix& A);
      friend Matrix operator~(Matrix& A);

      Matrix& operator=(Matrix& Source);
      Matrix& operator=(const Matrix& Source);
      Matrix& operator+=(Matrix& A);
      Matrix& operator-=(Matrix& A);
      Matrix& operator*=(Matrix& A);
      Matrix& operator*=(double& t);

      Matrix& Randomize(double rangeStr, double rangeEnd, double sparsity=1.0, int seed=0);
      Matrix& RandomizeBand(double range, long bandwidth, int seed=0);

      Vector SolveFor(Vector B);

      Matrix& MakeIdentity();
      void MakeZero();

      void Fill(Vector& V);

      void PrintNonzeros(std::ostream& s);

      inline void UsePivoting() { bUsePivoting = true; }
      inline void NoPivoting() { bUsePivoting = false; }
      inline void SetPrintNonZerosOnly() { bPrintNonZerosOnly = true; }

      inline unsigned setwidth(int iW) { unsigned pW = widthPrint; widthPrint = iW; return pW; }
      inline unsigned setprecision(int iP) { unsigned pP = precisionPrint; precisionPrint = iP; return pP; }

      friend double det(Matrix A);

      friend double Sum(Matrix& A);

      friend double InfNorm(const Matrix& A);
      friend double RowNorm(const Matrix& A);

      void Error(const char* Msg="");

      friend std::ostream& operator<<(std::ostream& s, Matrix& A);
      friend std::istream& operator>>(std::istream& s, Matrix& A);
   };

}  // namespace LinearAlgebra

#endif
