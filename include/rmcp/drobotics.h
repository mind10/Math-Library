#ifndef _RMCP_DROBOTICS_H_
#define _RMCP_DROBOTICS_H_

#include "rmcp/dvector.h"
#include "rmcp/dgqmatrix.h"

namespace rmcp {

	enum ROT_TYPE {
		ROT_RPY = 0x3300, // Roll-Pitch-Yaw
		ROT_ZYZ, 		// Euler ZYZ
		ROT_ZXZ			// Euler ZXZ
	};

	/* three dimensional vector class */
	class dvec3 : public rmcp::dvector {

		friend class dmat3;
		friend class dhmat;

		public:
			dvec3();
			dvec3(const double element);
			dvec3(const double x, const double y, const double z);
			dvec3(const double element[]);
			dvec3(const dvec3& v);

			const dvec3 &operator + (void) const;
			dvec3 operator - (void) const;
			dvec3 &operator = (const dvec3 &v);
			dvec3 &operator = (const double element[]);
			dvec3 &operator += (const dvec3 &v);
			dvec3 &operator -= (const dvec3 &v);
			dvec3 &operator *= (double s);
			dvec3 &operator /= (double s);
			dvec3 operator + (const dvec3 &v) const;
			dvec3 operator - (const dvec3 &v) const;
			dvec3 operator * (double s) const;
			dvec3 operator / (double s) const;
			double operator * (const dvec3 &v) const;	// inner product
			dvec3 operator * (const dmat3 &v) const; // vector = vector * matrix
			dmat3 operator ^ (const dvec3 &v2) const; // matrix = column vector x row vector
			dvec3 operator & (const dvec3 &v) const; // cross product

			friend dvec3 operator * (double s, const dvec3 &v);
			friend dvec3 cross(const dvec3 &x, const dvec3 &y); // cross product (x x y)
			friend dmat3 cross(const dvec3 &v,  const dmat3 &m);  // cross product : v x m (each column of m)
			friend dmat3 cross(const dvec3 &v);      // get skew-symmetric matrix from vector v
			friend dvec3 inv_cross(const dmat3 &m);  // get vector from skew-symmetric matrix
		
			friend dmat3 rotate(const dvec3 &ang, rmcp::ROT_TYPE type);
			friend dmat3 rotate(const double ang, dvec3 &axis);	// angle-axis
			friend dvec3 quaternion_rotate(const dvec3 &A, const dvec3 &P, double theta);  // axis A에 대한 point P의 theta만큼의 회전
			// Angle form Matrix
			friend dvec3 inv_rotate(const dmat3 &R, rmcp::ROT_TYPE type);
			friend void inv_rotate(double &ang, dvec3 &axis, const dmat3 &R);	// angle-axis

	};

	dvec3 cross(const dvec3 &x, const dvec3 &y); // cross product (x x y)
	dmat3 cross(const dvec3 &v,  const dmat3 &m);  // cross product : v x m (each column of m)
	dmat3 cross(const dvec3 &v);      // get skew-symmetric matrix from vector v
	dvec3 inv_cross(const dmat3 &m);  // get vector from skew-symmetric matrix
	dmat3 rotate(const dvec3 &ang, rmcp::ROT_TYPE type);
	dmat3 rotate(const double ang, dvec3 &axis);	// angle-axis
	dvec3 quaternion_rotate(const dvec3 &A, const dvec3 &P, double theta);	// axis A에 대한 point P의 theta만큼의 회전
	dvec3 inv_rotate(const dmat3 &R, rmcp::ROT_TYPE type);
	void inv_rotate(double &ang, dvec3 &axis, const dmat3 &R);	// angle-axis

	/* three dimensional square matrix */
	class dmat3 : public rmcp::dgqmatrix {

		friend class dvec3;
		friend class dhmat;

		public:
			dmat3();
			dmat3(const double element);
			dmat3(const double x, const double y, const double z);
			dmat3(const double element[]);
			dmat3(const dmat3 &m);

			const dmat3 &operator + (void) const;
			dmat3 operator - (void) const;
			dmat3 operator ~ (void) const;	// transpose
			dmat3 &operator = (const dmat3 &m);
			dmat3 &operator = (const double element[]);
			dmat3 &operator += (const dmat3 &m);
			dmat3 &operator -= (const dmat3 &m);
			dmat3 &operator *= (double s);
			dmat3 &operator *= (const dmat3 &m);
			dmat3 &operator /= (double s);
			dmat3 operator + (const dmat3 &m) const;
			dmat3 operator - (const dmat3 &m) const;
			dmat3 operator * (const double s) const;
			dmat3 operator / (const double s) const;
			dmat3 operator * (const dmat3 &m) const;
			dvec3 operator * (const dvec3 &v) const;
			
			// Decomposition & Inverse
			//void    svd(dmat3 &u, dmat3 &s, dmat3 &v) const;
			virtual double det(void);
			virtual void transpose(void);
			virtual void inv(void);
			
			friend dmat3 operator * (double s, const dmat3 &m);
			friend dmat3 inv(const dmat3 &m);
		//	friend double det(dmat3 &M);
		
			// Matrix from Angle
			friend dmat3 rotate(const dvec3 &ang, rmcp::ROT_TYPE type);
			friend dmat3 rotate(const double ang1, const double ang2, const double ang3, rmcp::ROT_TYPE type);
			friend dmat3 rotate(const double ang, dvec3 &axis);	// angle-axis

			// Angle form Matrix
			friend dvec3 inv_rotate(const dmat3 &R, rmcp::ROT_TYPE type);
			friend void inv_rotate(double &ang1, double &ang2, double &ang3, const dmat3 &R, rmcp::ROT_TYPE type);
			friend void inv_rotate(double &ang, dvec3 &axis, const dmat3 &R);	// angle-axis

			friend dmat3 cross(const dvec3 &v,  const dmat3 &m);  // cross product : v x m (each column of m)
			friend dmat3 cross(const dvec3 &v);      // get skew-symmetric matrix from vector v
			friend dvec3 inv_cross(const dmat3 &m);  // get vector from skew-symmetric matrix

//			friend dvec3 quaternion_rotate(const dvec3 &A, const dvec3 &P, double theta);  // axis A에 대한 point P의 theta만큼의 회전
			// p.144 example 5.2
			/*
			// Transformation Matrix : Angle Dot (v) => Angular Velocity (w)
			dmat3 &dmRVelTrans(const dmat3 &R, int type);
			dmat3 &dmRVelTrans(const double &ang1, const double &ang2, int type);

			// Angle Dot (v) => Angular Velocity (w)
			dvec3 &dmRAngVel(const dvec3& v, const dmat3 &R, int type);
			dvec3 &dmRAngVel(const dvec3& v, const double ang1, const double ang2, int type);

			// Angular Velocity (w) => Angle Dot (v)
			dvec3 &dmRAngDot(const dvec3& w, const dmat3 &R, int type);
			dvec3 &dmRAngDot(const dvec3& w, const double ang1, const double ang2, int type);
			*/
			
	};

	dmat3 inv(const dmat3 &m);
		//	double det(dmat3 &M);
	dmat3 rotate(const dvec3 &ang, rmcp::ROT_TYPE type);
	dmat3 rotate(const double ang1, const double ang2, const double ang3, rmcp::ROT_TYPE type);
	dmat3 rotate(const double ang, dvec3 &axis);	// angle-axis
	dvec3 inv_rotate(const dmat3 &R, rmcp::ROT_TYPE type);
	void inv_rotate(double &ang1, double &ang2, double &ang3, const dmat3 &R, rmcp::ROT_TYPE type);
	void inv_rotate(double &ang, dvec3 &axis, const dmat3 &R);	// angle-axis
	dmat3 cross(const dvec3 &v,  const dmat3 &m);  // cross product : v x m (each column of m)
	dmat3 cross(const dvec3 &v);      // get skew-symmetric matrix from vector v
	dvec3 inv_cross(const dmat3 &m);  // get vector from skew-symmetric matrix
//	dvec3 quaternion_rotate(const dvec3 &A, const dvec3 &P, double theta);	// axis A에 대한 point P의 theta만큼의 회전

	/* homogeneous transformation matrix */
	class dhmat : public rmcp::dgqmatrix {
		
		public:
			dhmat();
			dhmat(const double element[]); // size of array = 12
			dhmat(const rmcp::dmat3 &R, const rmcp::dvec3 &v);
			dhmat(const dhmat &H);

			dhmat &operator = (const dhmat &m);
			dhmat &operator = (const double element[]); // size of array = 12
			dhmat &operator *= (const dhmat &m);
			dhmat &operator /= (const dhmat &m);
			dhmat operator * (const dhmat &m) const;
			dhmat operator / (const dhmat &m) const;
			dvec3 operator * (const dvec3 &v) const;
	
			virtual void set(const double element[]);
			void set(const dmat3 &R, const dvec3 &v);
			void set(const dhmat &H);
			void set(const dmat3 &R);
			void set(const dvec3 &v);

			virtual void inv(void);
			
			friend dhmat inv(const dhmat &m);
			
	};

	dhmat inv(const dhmat &m);

} /* namespace rmcp */

#endif /* _RMCP_DROBOTICS_H_ */

