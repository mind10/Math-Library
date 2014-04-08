#ifndef _RMCP_DSE3_H_
#define _RMCP_DSE3_H_

#include "rmcp/dvector.h"
#include "rmcp/dgqmatrix.h"

#include <iostream>

namespace rmcp {

	class dse3 : public rmcp::dvector {

		friend class dvector;
		friend class dSE3;
		friend class dinertia;
		friend class dainertia;
		friend class djacobian;

		public:
			dse3(void);
			dse3(double w0, double w1, double w2, double v0, double v1, double v2);
			dse3(const dse3 &S);
			dse3(const rmcp::dvector &w, const rmcp::dvector &v);
			dse3(const double element[]);
			
			const dse3 &operator + (void) const;
			dse3 operator - (void) const;
			dse3 &operator = (const dse3 &S);
			dse3 &operator += (const dse3 &s);
			dse3 &operator -= (const dse3 &s);
			dse3 &operator *= (double s);
			dse3 &operator /= (double s);
			dse3 operator + (const dse3 &S) const;
			dse3 operator - (const dse3 &S) const;
			dse3 operator * (double s) const;
			dse3 operator / (double s) const;
			double operator * (const dse3 &s);
			dgqmatrix operator * (dSE3 &T) const;
	
			void set(double w0, double w1, double w2, double v0, double v1, double v2);
			void set_angular(dvector &w);
			void set_linear(dvector &v);
			dvector get_angular(void);
			dvector get_linear(void);

			void log(const dSE3 &T);
			void Ad(const rmcp::dSE3 &T, const dse3 &s);
			void inv_Ad(const rmcp::dSE3 &T, const dse3 &s);
			void ad(const dse3 &s1, const dse3 &s2);
			void dAd(const rmcp::dSE3 &T, const dse3 &t);
			void inv_dAd(const dSE3 &T, const dse3 &t);
			void dad(const dse3 &s, const dse3 &t);
			void rotate(const dSE3 &T, const dse3 &S);
			void inv_rotate(const dSE3 &T, const dse3 &S);

			rmcp::dmatrix export_dmatrix(void);
//			void import_dmatrix(rmcp::dmatrix &M);

//			void multiply_T1sT2(const rmcp::dSE3 &T1, const dse3 &s, const dSE3 &T2); 
			// = T1 * s * T2,             T1 == dSE3, T2 == dSE3, s == dse3

			friend dse3 operator * (double d, const dse3 &s);

//			friend double square_sum(const dse3 &s);

			friend dainertia kronecker_product(const dse3 &x, const dse3 &y);
			friend dainertia kronecker_product(const dse3 &x, const dse3 &y, double d);
	};

	dainertia kronecker_product(const dse3 &x, const dse3 &y);
	dainertia kronecker_product(const dse3 &x, const dse3 &y, double d);
	
	// column order 4 X 4 homogeneous transformation matrix
	class dSE3 : public rmcp::dgqmatrix {

		friend class dinertia;
		friend class dainertia;
		friend class dse3;
		friend class dso3;
		friend class dSO3;
		
		public:
			dSE3(void);
			dSE3(const rmcp::dSE3 &T);
			dSE3(double R0, double R1, double R2,
				   	double R3, double R4, double R5,
				   	double R6, double R7, double R8,
				   	double p0, double p1, double p2);
			dSE3(const rmcp::dSO3 &R, const rmcp::dvector &p);
			dSE3(const rmcp::dvector &v0,
				   	const rmcp::dvector &v1,
				   	const rmcp::dvector &v2, 
					const rmcp::dvector &p);
			dSE3(const rmcp::dSO3 &R);
			dSE3(const rmcp::dvector &p);
			dSE3(const double s);
			dSE3(const double element[]);

			dSE3 &operator = (const dSE3& T);
			dSE3 operator * (const dSE3& T) const;
			dSE3 operator / (const dSE3& T) const;
			dSE3 operator % (const dSE3& T) const;
			rmcp::dvector operator % (const rmcp::dvector& p) const;
			rmcp::dvector operator * (const rmcp::dvector& p) const;
			dSE3 &operator *= (const dSE3& T);
			dSE3 &operator /= (const dSE3& T);
			dSE3 &operator %= (const dSE3& T);
			bool operator == (const dSE3& T) const;
			dgqmatrix operator * (dse3 &A) const;

			rmcp::dvector position(void) const;
			dSO3 rotation(void) const;

			dSE3 &set_rotation(const dSO3 &R);
			dSE3 &set_position(const rmcp::dvector &p);
			void set_identity(void);

			dSE3 &translate(const rmcp::dvector &p);
			dSE3 &rotate(const dSO3 &R);
			virtual void inv(void);
			void exp(const dse3 &S);
			void exp(const dse3 &S, double theta);
			void rot_x(double t);
			void rot_y(double t); 
			void rot_z(double t);

			rmcp::dmatrix Ad(void);
			rmcp::dmatrix inv_Ad(void);
	};

	dSE3 exp(const dse3 &s);
	dSE3 exp(const dse3 &s, double theta);
	dse3 log(const dSE3 &T);

	dse3 Ad(const dSE3 &T, const dse3 &s);
	dse3 inv_Ad(const dSE3 &T, const dse3 &s);
	dse3 ad(const dse3 &s1, const dse3 &s2);
	dse3 dAd(const dSE3 &T, const dse3 &t);
	dse3 inv_dAd(const dSE3 &T, const dse3 &t);
	dse3 dad(const dse3 &s, const dse3 &t);

	dSE3 inv(const dSE3 &T);

	dse3 rotate(const dSE3 &T, const dse3 &S);
	dse3 inv_rotate(const dSE3 &T, const dse3 &S);

	/* deprecated */
	//			friend dse3 linearize(const dSE3 &T);
	//			friend rmcp::dvector linear_Ad(const rmcp::dvector &p, const dse3 &s);

	//			friend rmcp::dvector ad(const rmcp::dvector &s1, const dse3 &s2);
	//			friend dse3 dAd(const rmcp::dvector &p, const dse3 &t);
	//			friend dse3 Ad(const rmcp::dvector &p, const dse3 &s);
	//			friend dse3 logR(const dSE3 &T);


} // namespace rmcp 

#endif // _RMCP_DSE3_H_

