#ifndef _RMCP_H_
#define _RMCP_H_

#ifdef real
	#undef real
#endif

#include <mkl/mkl.h>


/*
s : real, single precision = float
d : real, double precision = double
c : complex, single precision = MKL_Complex8
z : complex, double precision = MKL_Complex16

blas
ge general matrix
gb general band matrix
sy symmetric matrix
sp symmetric matrix (packed storage)
sb symmetric band matrix
he Hermitian matrix
hp Hermitian matrix (packed storage)
hb Hermitian band matrix
tr triangular matrix
tp triangular matrix (packed storage)
tb triangular band matrix.

full storage
packed storage : symmetic, hermitian, triangular
band storage : band

lapack
ge general
gb general band
gt general tridiagonal
po symmetric or Hermitian positive-definite
pp symmetric or Hermitian positive-definite (packed storage)
pb symmetric or Hermitian positive-definite band
pt symmetric or Hermitian positive-definite tridiagonal
sy symmetric indefinite
sp symmetric indefinite (packed storage)
he Hermitian indefinite
hp Hermitian indefinite (packed storage)
tr triangular
tp triangular (packed storage)
tb triangular band
bd bidiagonal matrix
hs upper Hessenberg matrix
or (real) orthogonal matrix
op (real) orthogonal matrix (packed storage)
un (complex) unitary matrix
up (complex) unitary matrix (packed storage)
sb (real) symmetric band matrix
st (real) symmetric tridiagonal matrix
hb (complex) Hermitian band matrix

symmetric > Use Upper Triangle

quaternion : general quaternion
uquaternion : unit quaternion

*/

namespace rmcp {

	//////////////// structs ///////////////
	// underbar('_') means that that is just used internally.
	// Don't use the structs which is starting with '_' character.
	

	struct _sfullvector {
//		protected:
			int size_;
			float *element_;
	};

	struct _dfullvector {
//		protected:
			int size_;
			double *element_;
	};

	struct _sfullmatrix {
//		protected:
			int row_;
			int column_;
			int size_;
			float *element_;
	};

	struct _dfullmatrix {
//		protected:
			int row_;
			int column_;
			int size_;
			double *element_;
	};

	struct _dpackedmatrix {
//		protected:
			int row_;	/* column is equal to row */
			int size_;
			double *element_;
	};

	struct _dbandmatrix {
//		protected:
			int row_;
			int column_;
			int kl_;
			int ku_;
			double *element_;
	};

	struct _dtridmatrix {
//		protected:
			int row_;
	};


	//////////////// functions ////////////////
//	inline long min(long a, long b){ return ((a) <= (b) ? (a) : (b)); }
//	inline long max(long a, long b){ return ((a) >= (b) ? (a) : (b)); }

	/////////////// class prototype //////////	
	
	/* general vector, not sparse */ 
	class svector; 	/* the same as sgevector */
   	class dvector;  /* the same as dgevector */
   	class cvector;	/* the same as cgevector */
   	class zvector;  /* the same as zgevector */

	/* sparse vector */
	class sspvector;
   	class dspvector;
   	class cspvector;
   	class zspvector;

	/* general matrix, m by n matrix */
	class smatrix;	/* the same as sgematrix */
   	class dmatrix;	/* the same as dgematrix */
   	class cmatrix;	/* the same as cgematrix */
	class zmatrix;	/* the same as zgematrix */
	
	/* general square matrix */	
	class sgqmatrix;
	class dgqmatrix;	

	class dssmatrix;	/* skew-symmetric matrix */
	class dormatrix;	/* orthogonal matrix */
	class donmatrix;	/* orthonormal matrix */
	class dtrmatrix;	/* triangular matrix */
	class dsymatrix;	/* symmetric indefinite matrix */
	class dpomatrix;	/* symmetric positive-definite matrix */

	/* cubic spline */
	class scuspline;
	class dcuspline;

	/* optimization */	
	class ssqobjectfn;
	class dsqobjectfn;

	class ssqoptimize;
	class dsqoptimize;

	/* for robotics */

	/* 3-dimensional vector */
	class svec3;
	class dvec3;	
	/* 3 by 3, matrix */
	class smat3;
	class dmat3;
	/* homogeneous transformation matrix */
	class shtrans;
	class dhtrans;	
	
	// for Lie group
	class dSO3;		/* SO(3) */
	class dso3; 	/* so(3) */
	class dSE3; 	/* SE(3) */
	class dse3;		/* se(3) and se*(3) */

	// for Lie group dynamics, auxilary
	class djacobian;
	class dinertia;
	class dainertia;

	// quaternion
	class dquaternion;
	class duquaternion;

	////////////// type definition ///////////
	typedef svector sgevector;
	typedef dvector dgevector;
 	typedef cvector cgevector;
	typedef zvector zgevector;

	typedef smatrix sgematrix;
	typedef dmatrix dgematrix;
	typedef cmatrix cgematrix;
	typedef zmatrix zgematrix;

	///////////// const variable ///////////
	const double PI = 3.141592653589793238462643383279502884197;
	const double PI_2 = 1.570796326794896619231321691639751442099;
	const double PI_4 = 0.7853981633974483096156608458198757210495;
	const double EPS = 2.2204E-16;
	const double INF = 1.0E100;
	const double SMALL = 1.0e-9;
	const double TINY = 1.0E-12;
	const double SVD_EPS = 1.0e-6;
	const double DGEFA_EPS = 1.0E-6;
	const double LCP_EPS = 1.0E-6;

	const double MATRIX_NORM_EPS = 1.0E-20;	// matrix 계열 class
	const double VECTOR_NORM_EPS = 1.0E-20;	// vector 계열 class
	const double EXP_EPS = 1.0E-12;	// lie group class
	const double LOG_EPS = 1.0E-12;	// lie group class
	const double INV_EULER_EPS = 1.0E-16;	// lie group class

	const double SQOPTIMIZE_INF = 1.0E20;	// sqoptimize class

}

#endif /* _RMCP_H_ */


