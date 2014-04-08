#ifndef _RMCP_DVECTOR_H_
#define _RMCP_DVECTOR_H_

#include <rmcp/rmcp.h>
#include <iostream>

namespace rmcp {
	
	class dvector : public rmcp::_dfullvector {

		friend class rmcp::duquaternion;
		friend class rmcp::dmatrix;
		friend class rmcp::dso3;
		friend class rmcp::dSO3;
		friend class rmcp::dse3;
		friend class rmcp::dSE3;
		friend class rmcp::dinertia;
	
		public:
			dvector();	/* 3-dimensional vector with zero values */
			dvector(const int size);
			dvector(const int size, const double element);
			dvector(const int size, const double element[]);
			dvector(const double element0, const double element1, const double element2);	/* for 3-dimensional vector */
			dvector(const dvector &v);
//			dvector(const rmcp::dmatrix &m);
			virtual ~dvector();

			double &operator [] (const int i);
			const double &operator [] (int i) const;
			double &operator () (const int i);
			const double &operator () (int i) const;
			const dvector &operator + (void) const;
			dvector operator - (void) const;
			dvector &operator = (const dvector &v);
			dvector &operator = (const double element[]);
			dvector &operator += (const dvector &v);
			dvector &operator -= (const dvector &v);
			dvector &operator *= (double s);
			dvector &operator /= (double s);
			dvector operator * (double s) const;
			dvector operator / (double s) const;
			dvector operator + (const dvector &v) const;
			dvector operator - (const dvector &v) const;
			double operator * (const dvector &v) const;		// inner product
			dvector operator * (const dmatrix &m) const;
			dmatrix operator ^ (const dvector &v) const;  // matrix = col vec x row vec	
			dvector operator & (const dvector &v) const; // cross product, only 3 by 3 vector
			bool operator == (const dvector &v) const;

			int size(void) const { return size_; }
			double *element(void) const { return element_; }

			void set(const double element);
			void set(const double element[]);
			void set(const double element0, const double element1, const double element2); /* for 3-dimensional vector */
			void get(double element[]);
			void resize(const int size);
			void resize(const int size, const double element);
			void resize(const int size, double *element);

			bool zero(void) const;
			double norm(void) const;
			double square_sum(void) const;
			double normalize(void);
			double amax(void) const; // max of absolute values of elements
			double amin(void) const; // min of absolute values of elements
			int amax_index(void) const;
			int amin_index(void) const;
			void diagonal(rmcp::dmatrix &m);	// diagonal vector
			void sort(int dir = 0);
			void sort(unsigned int key[], int dir = 0);
			// dir = 0 : increasing order
			// dir = 1 : decreasing order
			// key[] : the size of the key is vector.size_

			double mean(void); // mean value
			double var(int n = 0); // variance, n = 0 : divided by size_-1, n = 1: divided by size_
			double std(int n = 0); // standard deviation, n = 0 : divided by size_-1, n = 1: divided by size_
			friend std::ostream &operator << (std::ostream &os, const dvector &v);
			friend dvector operator * (double s, const dvector &v);
			friend double inner(const dvector &x, const dvector &y);
			friend dvector cross(const dvector &x, const dvector &y);
			friend dvector rand(const int size);
			friend double cov(dvector &v1, dvector &v2, int n);
			friend double quadratic(const dvector &x, const dmatrix &A, const dvector &y);

			/* Lie Group Class(dso3, dSO3, dse3, dSE3) */
			dvector &operator %= (const dSE3 &T);		
			dvector &operator *= (const dSE3 &T);

			void rotate(const rmcp::dSE3 &T);
			void inv_rotate(const rmcp::dSE3 &T);
			void rotate(const rmcp::dSE3 &T, const dvector &v);
			void inv_rotate(const rmcp::dSE3 &T, const dvector &v);

/* deprecated */
//			friend dvector inv_Ad(const rmcp::dSE3 &T, const dvector &v);
//			friend rmcp::dse3 Ad(const dvector &p, const rmcp::dse3 &s);
//			friend dvector linear_Ad(const dvector &p, const rmcp::dse3 &s);	
//			friend dvector ad(const dvector &s1, const rmcp::dse3 &s2);
//			friend rmcp::dse3 dAd(const dvector &p, const rmcp::dse3 &t);

	};

	double inner(const dvector &x, const dvector &y);
	dvector cross(const dvector &x, const dvector &y);
	dvector rand(const int size);
	double quadratic(const dvector &x, const dmatrix &A, const dvector &y);

	double norm(const dvector &v);
	dvector zeros(const int size);
	dvector ones(const int size);
	double cov(dvector &v1, dvector &v2, int n = 0); // the sizes are equal.

	/* Lie Group Class */
	dvector rotate(const rmcp::dSE3 &T, const dvector &v);
	dvector inv_rotate(const rmcp::dSE3 &T, const dvector &v);

} // namespace rmcp

#endif // _RMCP_DVECTOR_H_
