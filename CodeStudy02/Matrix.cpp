#include <limits>
#include <iomanip>
#include <string>

#include "Matrix.hpp"

namespace LinearAlgebra {

   /* * Matrix Member Functions and Operators * */

   void Matrix::Error(const char* Msg)
   {
      std::cerr << "Matrix Error: " << Msg << std::endl;
      exit(1);
   }

   double& Matrix::operator()(long i, long j)
   {
      return THIS[i][j];
   }

   Vector& Matrix::operator[](long i)
   {
      return Rows[i];
   }

   const Vector& Matrix::operator[](long i) const
   {
      return Rows[i];
   }

   Matrix::Matrix()
   {
      Setup(0, 0);
   }

   Matrix::Matrix(long nrow, long ncol)
   {
      Setup(nrow, ncol);
   }

   Matrix::Matrix(const Matrix& Source)
   {
      if (&Source == this) return;
      Setup(Source);
   }

   Matrix::Matrix(Matrix& Source)
   {
      if (&Source == this) return;
      Setup(Source);
   }

   Matrix::Matrix(const std::vector<Vector>& Source)
   {
      Setup((long)Source.size(), Source[0].N);

      for (long i = 0; i < NumRow; i++) {
         if (Source[i].N != NumCol)
            Error("Number of entries in all rows must be the same.");
         THIS[i].Setup(Source[i]);
      }
   }

   Matrix::~Matrix()
   {
      // no need!
   }

   void Matrix::Setup(long nrow, long ncol)
   {
      NumRow = nrow;
      NumCol = ncol;

      Rows.assign(NumRow, Vector(NumCol));

      widthPrint = 11;
      precisionPrint = 2;
      bPrintIndex = false;
      bUsePivoting = true;
      bPrintNonZerosOnly = false;
   }

   void Matrix::Setup(const Matrix& Source)
   {
      NumRow = Source.NumRow;
      NumCol = Source.NumCol;

      Rows.assign(NumRow, Vector(NumCol));

      for (long i = 0; i<NumRow; i++)
         THIS[i] = Source[i];

      widthPrint = Source.widthPrint;
      precisionPrint = Source.precisionPrint;
      bPrintIndex = Source.bPrintIndex;
      bUsePivoting = Source.bUsePivoting;
      bPrintNonZerosOnly = Source.bPrintNonZerosOnly;
   }

   void Matrix::Resize(long nrow, long ncol)
   {
      Setup(nrow, ncol);
   }

   void Matrix::Fill(Vector& Source)
   {
      long i, j, k=0;
      for (i=0; i<NumRow; i++) {
         for (j=0; j<NumCol; j++) {
            if (k == Source.Size())
               continue;
            THIS[i][j] = Source[k];
            k++;
         }
      }
   }

   Matrix operator~(Matrix& Source)
   {
      long NumRow = Source.NumRow;
      long NumCol = Source.NumCol;
      Matrix T(NumCol, NumRow);
      for ( long i=0; i<NumCol; i++ )
         for ( long j=0; j<NumRow; j++ )
            T[i][j] = Source[j][i];
      T.bUsePivoting = Source.bUsePivoting;
      T.widthPrint = Source.widthPrint;
      T.precisionPrint = Source.precisionPrint;
      T.bPrintIndex = Source.bPrintIndex;
      T.bPrintNonZerosOnly = Source.bPrintNonZerosOnly;
      return T;
   }

   Matrix& Matrix::operator=(Matrix& Source)
   {
      if (&Source == this) return THIS;
      Setup(Source);
      return THIS;
   }

   Matrix& Matrix::operator=(const Matrix& Source)
   {
      if ( &Source == this ) return THIS;
      Setup(Source);
      return THIS;
   }

   Matrix& Matrix::operator+=(Matrix& A)
   {
      return THIS = THIS + A;
   }

   Matrix& Matrix::operator-=(Matrix& A)
   {
      return THIS = THIS - A;
   }

   Matrix& Matrix::operator*=(Matrix& A)
   {
      return THIS = THIS * A;
   }

   Matrix& Matrix::operator*=(double& t)
   {
      return THIS = THIS * t;
   }

   void Matrix::PrintNonzeros(std::ostream& s)
   {
      for ( long i=0; i<NumRow; i++) {
         for ( long j=0; j<NumCol; j++) {
            if ( THIS[i][j] == 0.0 )
               s << ".";
            else
               s << "x";
          }
          s << std::endl;
      }
      s << std::endl;
   }

   std::ostream& operator<<(std::ostream& s, Matrix& A)
   {
      long i, j;
      long NumRow = A.NumRow;
      long NumCol = A.NumCol;

      if (A.bPrintIndex)
         s << NumRow << " " << NumCol << std::endl;

      if ( A.bPrintIndex ) {
         for (j=0; j<NumCol; j++)
            s << std::setw(A.widthPrint) << j;
         s << std::endl;
      }

      double Ainf = InfNorm(A);

      std::streamsize old_precision = s.precision(A.precisionPrint);

      for (i=0; i<NumRow; i++) {
         for (j=0; j<NumCol; j++) {
            double a = A[i][j];
            if ( (a != 0.0) && (abs(a) / Ainf < 1e-9) )
               a = 0.0;
            if (A.bPrintNonZerosOnly && a == 0.0)
               s << std::string(A.widthPrint, ' ');
            else
            s << std::setw(A.widthPrint) << std::fixed << std::setprecision(A.precisionPrint) << a+0.0;
         }
         s << std::endl;
      }
      s.precision(old_precision);
      return s;
   }

   std::istream& operator>>(std::istream& s, Matrix& A)
   {
      long NumRow, NumCol;
      s >> NumRow >> NumCol;
      A(NumRow, NumCol);
      for (long i=0; i<NumRow; i++) {
         for (long j=0; j<NumCol; j++) {
            double x;
            s >> x;
            if ( x != 0.0 )
               A[i][j] = x;
         }
      }
      return s;
   }

   Matrix operator-(Matrix& A)
   {
      Matrix B(A);
      for ( long i=0; i<A.NumRow; i++ )
         B[i] = -A[i];
      return B;
   }

   Matrix operator+(Matrix& A)
   {
      return A;
   }

   Matrix operator+(Matrix& A, Matrix& B)
   {
      long NumRow = A.NumRow;
      long NumCol = A.NumCol;
      if ( (NumRow != B.NumRow) || (NumCol != B.NumCol) )
         A.Error("must have equal dimensions");
      Matrix Sum(NumRow, NumCol);
      for ( long i=0; i<NumRow; i++ )
         Sum[i] = A[i] + B[i];
      return Sum;
   }

   Matrix operator-(Matrix& A, Matrix& B)
   {
      long NumRow = A.NumRow;
      long NumCol = A.NumCol;
      if ( (NumRow != B.NumRow) || (NumCol != B.NumCol) )
         A.Error("must have equal dimensions");
      Matrix Dif(NumRow, NumCol);
      for ( long i=0; i<NumRow; i++ )
         Dif[i] = A[i] - B[i];
      return Dif;
   }

   Matrix operator*(Matrix& A, double& t)
   {
      Matrix P(A);
      for ( long i=0; i<A.NumRow; i++ )
         P[i] = A[i] * t;
      return P;
   }

   Matrix operator*( double& t, Matrix& A )
   {
      return A * t;
   }

   Vector operator*(Vector& V, Matrix& A)
   {
      if ( V.Size() != A.NumRow )
         A.Error("# rows of the matrix must equal to the size of vector for vector-matrix multiplication");
      long N = V.Size();
      long NumCol = A.NumCol;
      Vector P(N);
      for (long j=0; j<NumCol; j++) {
         double s = 0.0;
         for (long i=0; i<N; i++)
            s += (V[i] * A[i][j]);
         P[j] = s;
      }
      return P;
   }

   Vector operator*(Matrix& A, Vector& V)
   {
      if (A.NumCol != V.Size())
         A.Error("# columns of the matrix must equal to the size of vector for matrix-vector multiplication");
      long N = V.Size();
      long NumRow = A.NumRow;
      Vector P(N);
      for (long i = 0; i < NumRow; i++) {
         double s = 0.0;
         for (long j = 0; j < N; j++)
            s += (A[i][j] * V[j]);
         P[i] = s;
      }
      P.SetColumnVector();
      return P;
   }

   Matrix operator*(Matrix& A, Matrix& B)
   {
      if ( A.NumCol != B.NumRow )
         A.Error("# rows of second matrix must equal # columns of first for multiplication");
      long NumRow = A.NumRow;
      long NumCol = B.NumCol;
      long N = A.NumCol;
      Matrix P(NumRow, NumCol);
      for (long i = 0; i < NumRow; i++) {
         for (long j = 0; j < NumCol; j++) {
            double s = 0.0;
            for (long k = 0; k < N; k++)
               s += (A[i][k] * B[k][j]);
            P[i][j] = s;
         }
      }
      return P;
   }

   Matrix& Matrix::Randomize(double rangeStr, double rangeEnd, double sparsity/*=1.0*/, int seed/*=0*/)
   {
      if (seed != 0)
         srand(seed);
      else
         srand((unsigned int)time(NULL));

      MakeZero();

      double range = rangeEnd - rangeStr;

      for ( long j=0; j<NumCol; j++ )
      {
         for ( long i=0; i<NumRow; i++ )
         {
            double rSpars = rand()/(RAND_MAX+1.0);
            if ( range!=0.0 && rSpars<sparsity )
            {
               int rnd = rand();
               double rX = rnd / (RAND_MAX + 1.0) * abs(range);
               rX = rangeStr + rX;
               THIS(i, j) = rX;
            }
         }
      }
      return THIS;
   }

   Matrix& Matrix::RandomizeBand(double range, long bandwidth, int seed/*=0*/)
   {
      if (seed != 0)
         srand(seed);
      else
         srand((unsigned int)time(NULL));

      MakeZero();

      for ( long j=0; j<NumCol; j++ )
      {
         long i1 = j - bandwidth/2; i1 = std::max((long)0, i1);
         long i2 = j + bandwidth/2; i2 = std::min((long)NumCol-1, i2);
         for ( long i=i1; i<=i2; i++ )
         {
            if ( range!=0.0 )
            {
               double rX = rand()/(RAND_MAX+1.0)*abs(range);
               if ( range < 0.0 ) rX = abs(range) - rX*2.0;
               if ( rX ) THIS(i, j) = rX;
            }
         }
      }
      return THIS;
   }

   Matrix& Matrix::MakeIdentity()
   {
      MakeZero();
      for (long i=0; i<NumRow; i++)
         THIS[i][i] = 1.0;
      return THIS;
   }

   void Matrix::RowSwap(long i1, long i2)
   {
      for ( long j=0; j<NumCol; j++ ) {
         // swap entries
         double x = THIS[i1][j];
         THIS[i1][j] = THIS[i2][j];
         THIS[i2][j] = x;
      }
   }

   Vector Matrix::SolveFor(Vector B)
   {
      Matrix A(*this);

      long n = A.NumRow;
      if ( n != A.NumCol )
         A.Error("must be square matrix for Gauss Solver");
      if (n != B.Size())
         A.Error("Right hand-side vector must have the same number of rows of coefficient matrix");

      Vector X(n);

      double Akk, Aik, Mik;
      long i, j, k, l;

      // Gaussian elimination with partial pivoting (optional)
      for ( k=0; k<n-1; k++ ) {
         if ( A.bUsePivoting ) {
            Akk = 0.0; l = 0;
            for ( i=k; i<n; i++ ) {
               Aik = A[i][k];
               if ( abs(Aik) > abs(Akk) ) {
                  Akk = Aik;
                  l = i;
               }
            }
            if ( Akk == 0.0 )
               A.Error("cannot solve due to singularity");

            if ( l > k ) {
               A.RowSwap(l, k);
               B.Swap(l, k);
            }
         }
         else {
            Akk = A[k][k];
         }

         if ( abs(Akk) < std::numeric_limits<double>::lowest())
            A.Error("cannot solve due to near-zero diagonal entry");

         for ( i=k+1; i<n; i++ ) {
            Mik = -A[i][k] / Akk;
            A[i][k] = Mik;
            for ( j=k+1; j<n; j++ ) {
               A[i][j] += Mik * A[k][j];
            }
         }
      }

      // Row operations on B
      //for ( k=0; k<n-1; k++ )          // bu dongu orijinal algoritma; sparse matrix icin 
      //   for ( i=k+1; i<n; i++ )    // kolon operasyonlarini non-zero entry'ler icin yapan asagidaki loop daha hizli calisir
      //   {
      //      B[i] -= A[i][k] * B[k]; Nmult++; Nsum++; 
      //   }
      for ( k=1; k<n; k++ ) {
         for ( i=0; i<k; i++ ) {
            B[k] += A[k][i] * B[i];
         }
      }

      // Back substitution
      double S;
      for ( i=n-1; i>=0; i-- ) {
         S = B[i];
         for ( j=i+1; j<n; j++ ) {
            S -= A[i][j] * X[j];
         }
         double Aii = A[i][i];
         if (abs(Aii) < std::numeric_limits<double>::lowest())
            A.Error("solution does not exist!");
         X[i] = S / A[i][i];
      }

      X.SetColumnVector();

      return X;
   }

   Vector operator/(Vector& B, Matrix& A)
   {
      long n = A.NumRow;
      if ( n != A.NumCol )
         A.Error("must be square matrix for Gauss Solver");
      if ( n != B.Size() )
         A.Error("Right hand-side vector must have the same number of rows of the coefficient matrix");

      Vector X(n);

      X = A.SolveFor(B);

      return X;
   }

   //Matrix operator!(Matrix& A)
   //{
   //   long m = A.NumRow;
   //   Matrix I(m, m);
   //   I.MakeIdentity();
   //   Matrix X = I/A;
   //   return X;
   //}

   double det(Matrix A)
   {
      long n = A.NumRow;
      if ( n != A.NumCol )
         A.Error("must be square matrix for Determinant Calculation");

      double Mik, detval = 1.0;
      double Akk, Aik;
      long i, j, k, l;

      for ( k=0; k<n; k++ ) {
         if ( A.bUsePivoting ) {
            Akk = 0.0; l = 0;
            for ( i=k; i<n; i++ ) {
               Aik = A[i][k];
               if ( abs(Aik) > abs(Akk) ) {
                  Akk = Aik;
                  l = i;
               }
            }
            if ( Akk == 0.0 ) 
               A.Error("cannot invert due to singularity");
            if ( l > k ) {
               A.RowSwap(l, k);
               detval = -detval; 
            }
         }
         else {
            Akk = A[k][k];
         }

         if ( abs(Akk) < std::numeric_limits<double>::lowest())
            A.Error("cannot solve due to near-zero or zero diagonal entry");

         detval *= Akk;
         for ( i=k+1; i<n; i++ ) {
            Mik = A[i][k] / Akk;
            for ( j=k+1; j<n; j++ ) {
               A[i][j] -= Mik * A[k][j];
            }
            A[i][k] = Mik;
         }
      }

      if ( detval == 0.0 ) 
         A.Error("cannot invert due to singularity");

      return detval;
   }

   void Matrix::MakeZero()
   {
      for ( long i=0; i<NumRow; i++ )
         THIS[i].MakeZero();
    }

   double Sum(Matrix& A)
   {
      long NumRow = A.NumRow;
      double S = Sum(A[0]);
      for (long i=1; i<NumRow; i++)
         S += Sum(A[i]);
      return S;
   }

   double RowNorm(const Matrix& A)
   {
      long NumRow = A.NumRow;
      double norm = CityNorm(A[0]);
      for (long i=1; i<NumRow; i++)
         norm = std::max(norm, CityNorm(A[i]));
      return norm;
   }

   double ColNorm(Matrix& A)
   {
      return RowNorm(~A);
   }

   double InfNorm(const Matrix& A)
   {
      long NumRow = A.NumRow;
      double norm = InfNorm(A[0]);
      for (long i=1; i<NumRow; i++)
         norm = std::max(norm, InfNorm(A[i]));
      return norm;
   }

}  // namespace LinearAlgebra

