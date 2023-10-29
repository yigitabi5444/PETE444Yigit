#include <algorithm>
#include <iomanip>

#include "Vector.hpp"

namespace LinearAlgebra {

   /* * Vector Member Functions and Operators * */

   void Vector::Error(VecError vErr, const char* Msg) const
   {
      std::cerr << "Vector Error = " << vErr << " : " <<   Msg << std::endl;
      exit(vErr);
   }

   double& Vector::operator[](long j)
   {
      #if defined(_DEBUG)
      if (j<0 || j >= N)
         Error(IDERR, "Index out of range");
      #endif
      return Array[j];
   }

   const double& Vector::operator[](long j) const
   {
      #if defined(_DEBUG)
      if ( j<0 || j>=N )
         Error(IDERR, "Index out of range");
      #endif
      return Array[j];
   }

   double& Vector::operator()(long j)
   {
      return THIS[j];
   }

   const double& Vector::operator()(long j) const
   {
      return THIS[j];
   }

   Vector::Vector()
   {
      Setup(0);
   }

   Vector::Vector(long n)
   {
      Setup(n);
   }

   Vector::Vector(long n, const double t)
   {
      Setup(n, t);
   }

   Vector::Vector(const Vector& Source)
   {
      Setup(Source);
   }

   Vector::Vector(const std::vector<double>& Source)
   {
      N = (long)Source.size();
      widthPrint = 11;
      precisionPrint = 2;
      Array.assign(Source.begin(), Source.end());
      bColVec = false;
   }

   Vector::~Vector()
   {
      // no need!
   }

   void Vector::Setup(long n)
   {
      if (n > MAXENTRY)
         Error(NOMEM, "Attempts to allocate memory for Vector more than absolute maximum!");
      N = n;
      widthPrint = 11;
      precisionPrint = 2;
      Array.assign(n, 0.0);
      bColVec = false;
   }

   void Vector::Setup(long n, const double t)
   {
      if (n > MAXENTRY)
         Error(NOMEM, "Attempts to allocate memory for Vector more than absolute maximum!");
      N = n;
      widthPrint = 11;
      precisionPrint = 2;
      Array.assign(n, t);
      bColVec = false;
   }

   void Vector::Setup(const Vector& Source)
   {
      if (&Source == this) return;
      N = Source.N;
      widthPrint = Source.widthPrint;
      precisionPrint = Source.precisionPrint;
      Array.assign(Source.Array.begin(), Source.Array.end());
      bColVec = Source.bColVec;
   }

   void Vector::Resize(long n)
   {
      if (n > MAXENTRY)
         Error(NOMEM, "Attempts to allocate memory for Vector more than absolute maximum!");
      N = n;
      Array.assign(n, 0.0);
   }

   void Vector::Assign(long n, const double t)
   {
      Setup(n, t);
   }

   // Operators

   Vector& Vector::operator=(const double& val)
   {
      Array.assign(N, val);
      bColVec = false;
      return THIS;
   }

   Vector& Vector::operator=(const std::vector<double>& Source)
   {
      N = (long)Source.size();
      widthPrint = 11;
      precisionPrint = 2;
      Array.assign(Source.begin(), Source.end());
      bColVec = false;
      return THIS;
   }

   Vector& Vector::operator+=(const Vector& V)
   {
      return THIS = THIS + V;
   }

   Vector& Vector::operator-=(const Vector& V)
   {
      return THIS = THIS - V;
   }

   Vector& Vector::operator*=(const double& t)
   {
      return THIS = THIS*t;
   }

   Vector& Vector::operator=(const Vector& Source)
   {
      Setup(Source);
      return THIS;
   }

   Vector operator-(const Vector& V)
   {
      Vector W(V.Size());
      for (long i=0; i<V.Size(); i++)
         W[i] = -V[i];
      return W;
   }

   Vector operator+(const Vector& V)
   {
      return V;
   }

   Vector operator+(const Vector& V, const Vector& W)
   {
      long M, N, L;
      M = V.Size();
      N = W.Size();
      L = std::min(M, N);
      Vector Sum(std::max(M, N));
      long i=0;
      for ( ; i<L; i++)
         Sum[i] = V[i] + W[i];
      for ( ; i<M; i++)
         Sum[i] = V[i];
      for ( ; i<N; i++)
         Sum[i] = W[i];
      return Sum;
   }

   Vector operator+(const double& S, const Vector& V)
   {
      long n = V.Size();
      Vector Sum(n);
      for ( long i=0; i<n; i++ )
         Sum[i] = S + V[i];
      return Sum;
   }

   Vector operator+(const Vector& V, const double& S)
   {
      long n = V.Size();
      Vector Sum(n);
      Sum = S + V;
      return Sum;
   }

   Vector operator-(const Vector& V, const Vector& W)
   {
      long M, N, L;
      M = V.Size();
      N = W.Size();
      L = std::min(M, N);
      Vector Sum(std::max(M, N));
      long i=0;
      for ( ; i<L; i++)
         Sum[i] = V[i] - W[i];
      for ( ; i<M; i++)
         Sum[i] = V[i];
      for ( ; i<N; i++)
         Sum[i] = -W[i];
      return Sum;
   }

   Vector operator*(const double& t, const Vector& V)
   {
      Vector P(V.Size());
      for (long i=0; i<V.Size(); i++)
         P[i] = t*V[i];
      return P;
   }

   Vector operator*(const Vector& V, const double& t)
   {
      return t*V;
   }

   Vector operator/(const Vector& V, const double& t)
   {
      Vector P(V.Size());
      for (long i=0; i<V.Size(); i++)
         P[i] = V[i]/t;
      return P;
   }

   double operator*(const Vector& V, const Vector& W)
   {
      long M, N;
      M = V.Size();
      N = W.Size();
      if ( M != N ) V.Error(NEVS, "Not equal vector sizes");
      double ScaPro = 0.0;
      for (long i=0; i<N; i++)
         ScaPro += V[i]*W[i];
      return ScaPro;
   }

   Vector operator&(const Vector& V, const Vector& W)
   {
      long i, M, N;
      M = V.Size();
      N = W.Size();
      Vector S(M+N);
      for ( i=0; i<M; i++ )
         S[i] = V[i];
      for ( i=0; i<N; i++ )
         S[i+M] = W[i];
      return S;
   }

   Vector operator&(const double& S, const Vector& V)
   {
      long i, N;
      N = V.Size();
      Vector W(N+1);
      W[0] = S;
      for ( i=0; i<N; i++ )
         W[i+1] = V[i];
      return W;
   }

   Vector operator&(const Vector& V, const double& S)
   {
      long i, N;
      N = V.Size();
      Vector W(N+1);
      for ( i=0; i<N; i++ )
         W[i] = V[i];
      W[N+1] = S;
      return W;
   }

   double Sum(const Vector& V)
   {
      double S = 0.0;
      for (long i=0; i<V.Size(); i++)
         S += V[i];
      return S;
   }

   double MaxNorm(const Vector& V)
   {
      double norm = 0.0;
      for (long i=0; i<V.Size(); i++)
         if (abs(V[i]) > norm ) norm = abs(V[i]);
      return norm;
   }

   double InfNorm(const Vector& V)
   {
      return MaxNorm(V);
   }

   double CityNorm(const Vector& V)
   {
      double norm = 0.0;
      for (long i=0; i<V.Size(); i++)
         norm += abs(V[i]);
      return norm;
   }

   double EuclidNorm(const Vector& V)
   {
      return sqrt(V*V);
   }

   std::ostream& operator<<(std::ostream& s, const Vector& V)
   {
      std::streamsize old_precision = s.precision(V.precisionPrint);
      long n = V.Size();
      for (long i=0; i<n; i++) {
         //s <<  setwidth(V.widthPrint) << V[i];
         s << std::setw(V.widthPrint) << std::fixed << std::setprecision(V.precisionPrint) << V[i]+0.0;
         if (V.bColVec)
            s << std::endl;
      }
      s.precision(old_precision);
      return s;
   }

   std::istream& operator>>(std::istream& s, Vector& V)
   {
      long n = V.Size();
      V.Resize(n);
      for (long i=0; i<n; i++)
         s >> V[i];
      return s;
   }

   Vector& Vector::Randomize(const double& range, const double sparsity/*=1.0*/, int seed/*=0*/)
   {
      if (seed != 0)
         srand(seed);
      else
         srand((unsigned int)time(NULL));

      double rX, rSpars;
      for (long i=0; i<N; i++) {
         Array[i] = 0.;
         rSpars = rand()/(RAND_MAX+1.0);
         if ( range!=0 && rSpars<sparsity ) {
            rX = rand()/(RAND_MAX+1.0)*abs(range);
            if ( range < 0.0 ) rX = abs(range) - rX*2.0;
            if ( rX  ) Array[i] = rX;
         }
      }
      return THIS;
   }

   Vector& Vector::Randomize(const double& range, long j1, long j2)
   {
      for (long i=j1; i<=j2; i++) {
         if ( range != 0 ) {
            double rnd = rand()/(RAND_MAX+1.0)*abs(range);
            if ( range < 0.0 ) rnd = abs(range) - rnd*2.0;
            if ( rnd != 0 ) THIS[i] = rnd;
         }
      }
      return THIS;
   }

   Vector& Vector::Scale(Vector& minvec, Vector& maxvec)
   {
      for (long i=0; i<N; i++) {
         if ( THIS[i] < minvec[i] )
            THIS[i] = 0.;
         else if ( THIS[i] > maxvec[i] )
            THIS[i] = 1.;
         else
            THIS[i] = (THIS[i]-minvec[i])/(maxvec[i]-minvec[i]);
      }
      return THIS;
   }

   Vector& Vector::Upscale(Vector& minvec, Vector& maxvec)
   {
      for (long i=0; i<N; i++) {
         if ( THIS[i] < 0.0 )
            THIS[i] = minvec[i];
         else if ( THIS[i] > 1.0 )
            THIS[i] = maxvec[i];
         else
            THIS[i] = minvec[i] + THIS[i]*(maxvec[i]-minvec[i]);
      }
      return THIS;
   }

   double Vector::Distance(const Vector& V)
   {
      double sum = 0.0, d;
      for (long i=0; i<N; i++) {
         d = THIS[i] - V[i];
         sum += (d*d);
      }
      return sqrt(sum);
   }

   double Vector::Difference(const Vector& V)
   {
      double sum = 0.0, d;
      for (long i=0; i<N; i++) {
         d = THIS[i] - V[i];
         sum += (d*d);
      }
      return sum;
   }

   Vector Vector::DownSample(long n)
   {
      #if defined(_DEBUG)
      if ( n<1 || n>N/2 )
         Error(IDERR, "Invalid downsampling frequency");
      #endif
      Vector X(N/n+(N%n==0?0:1));
      for ( long i=0,j=0; i<N; i+=n,j++)
         X[j] = THIS[i];
      X[X.Size()] = THIS[N];
      return X;
   }

   void Vector::MakeZero()
   {
      Array.assign(N, 0.0);
   }

   long MaxMod(const Vector& V)
   {
      double maxval = 0.0; long maxmod = 0;
      for (long i=0; i<V.Size(); i++)
         if (abs(V[i]) > maxval )
            maxmod = i;
      return maxmod;
   }

}  // namespace LinearAlgebra
