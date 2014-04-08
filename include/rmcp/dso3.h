#ifndef _RMCP_DSO3_H_
#define _RMCP_DSO3_H_

#include "rmcp/dssmatrix.h"
#include "rmcp/donmatrix.h"

#include <iostream>

namespace rmcp {
/*	
 *	(w1, w2, w3)
 *
 *	  0  -w3   w2
 *	 w3    0  -w1
 *	-w2   w1    0
 *
 */
	
	class dso3 : public rmcp::dssmatrix {

		friend class rmcp::dvector;
		friend class rmcp::dse3;
		friend class rmcp::dSO3;
		friend class rmcp::dSE3;

		public:
			dso3(void);
			dso3(double d);
			dso3(const double v[]);	/* size of v[] = 3 */
			dso3(double v0, double v1, double v2);
			dso3(const rmcp::dvector &v); /* size of dvector = 3 */
			dso3(const dso3 &w);

			const dso3& operator + (void) const;
			dso3 operator - (void) const;
			dso3 &operator = (const dso3 &w);
			dso3 &operator = (const rmcp::dvector &v);
			dso3 &operator = (const double v[]); /* size of v[] = 3 */
			dso3 &operator += (const dso3 &w);
			dso3 &operator -= (const dso3 &w);
			dso3 &operator *= (double s);
			dso3 &operator /= (double s);
			dso3 operator + (const dso3 &w) const;
			dso3 operator - (const dso3 &w) const;
			dso3 operator * (double s) const;
			dso3 operator / (double s) const;
			bool operator == (const dso3 &w) const;

			void get(rmcp::dvector &w);		/* get [w], w : 3-dimensional vector */
			void set(double v0, double v1, double v2);
			virtual void set(const double v[]);
			virtual double norm(void) const;	/* vector norm, not matrix norm */
			
			void log(const rmcp::dSO3& R);

			friend rmcp::dso3 operator * (double s, const dso3 &v);
	};

	class dSO3 : public rmcp::donmatrix {

		friend class dso3;
		friend class dSE3;
		friend class duquaternion;

		public:
			enum EULER {
				ZYX = 0,
				ZYZ
			};

			dSO3(void);
			dSO3(const dSO3 &R);
			dSO3(const double element[]);
			dSO3(double R0, double R1, double R2, double R3, double R4,
					double R5, double R6, double R7, double R8);

			dSO3 &operator = (const dSO3 &R);
			dSO3 operator ~ (void) const;
			dSO3 &operator *= (const dSO3 &R);
			dSO3 operator * (const dSO3 &R) const;
			dSO3 operator % (const dSO3 &R) const;
			rmcp::dvector operator * (const rmcp::dvector &p) const;

			virtual void inv(void);
			void exp(const dso3 &w);
			/* theta = radian, w = normalized. */
			void exp(const dso3 &w, double theta);
			void euler_zyx(double z, double y, double x); // radian
			int inv_euler_zyx(double &z, double &y, double &x); // if singular, return -1. else, return 0. 
			void euler_zyz(double z1, double y, double z2);
			int inv_euler_zyz(double &z1, double &y, double &z2); // if singular, return -1. else, return 0.
			
			dSO3::EULER inv_euler(double &a, double &b, double &c);
			/* if return = ZYX, a = z, b = y, c = x, if return ZYZ, a = z1, b = y, c = z2 */

			void rollpitchyaw(double roll, double pitch, double yaw);
			void inv_rollpitchyaw(double &roll, double &pitch, double &yaw); /* alpha = roll, beta = pitch, gamma = yaw */

			void rot_x(double t);
			void rot_y(double t); 
			void rot_z(double t);
	};	

	rmcp::dSO3 exp(const rmcp::dso3 &w);
	/* theta = radian, w = normalized. */
	rmcp::dSO3 exp(const rmcp::dso3 &w, double theta); 
	rmcp::dSO3 inv(const rmcp::dSO3 &R);	
	rmcp::dso3 log(const rmcp::dSO3 &R);
	dSO3 rot_x(double t);
	dSO3 rot_y(double t); 
	dSO3 rot_z(double t);
	rmcp::dSO3 euler_zyx(double z, double y, double x);
	int inv_euler_zyx(rmcp::dSO3 &R, double &z, double &y, double &x);
	rmcp::dSO3 euler_zyz(double z1, double y, double z2);
	int inv_euler_zyz(rmcp::dSO3 &R, double &z1, double &y, double &z2);		

} // namespace rmcp 

#endif /* _RMCP_DSO3_H_ */

