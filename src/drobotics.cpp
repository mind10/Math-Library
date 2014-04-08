
#include "rmcp/drobotics.h"
#include <cmath>
#include <cassert>

using namespace rmcp;

static int i, j;

dvec3::dvec3() : rmcp::dvector(3, 0.0)
{
}

dvec3::dvec3(const double element) : rmcp::dvector(3, element)
{
}

dvec3::dvec3(const double x, const double y, const double z) : rmcp::dvector(3)
{
	element_[0] = x;
	element_[1] = y;
	element_[2] = z;
}

dvec3::dvec3(const double element[]) : rmcp::dvector(3, element)
{ 
}

dvec3::dvec3(const dvec3& v) : rmcp::dvector(3, v.element_)
{
}

const dvec3 &
dvec3::operator + (void) const
{
   	return *this;
}

dvec3 
dvec3::operator - (void) const	
{
	dvec3 vr;
	
	for (i = 0; i < 3; i++) {
		vr.element_[i] = - element_[i];
	}

   	return vr;
}

dvec3 &
dvec3::operator = (const dvec3 &v)
{
	cblas_dcopy(3, v.element_, 1, element_, 1);
   	return *this;
}

dvec3 &
dvec3::operator = (const double element[])
{
	cblas_dcopy(3, element, 1, element_, 1);
   	return *this;
}

dvec3 &
dvec3::operator += (const dvec3& v)
{
	for (i = 0; i < 3; i++) {
		element_[i] += v.element_[i];
	}

	return *this;
}

dvec3 &
dvec3::operator -= (const dvec3& v)
{
	for (i = 0; i < 3; i++) {
		element_[i] -= v.element_[i];
	}

   	return *this;
}

dvec3 & 
dvec3::operator *= (double s)
{
	cblas_dscal(3, s, element_, 1);
   	return *this;
}

dvec3 & 
dvec3::operator /= (double s)
{
	assert ( s != 0.0 &&
			"const dvec3 &dvec3::operator /= (double)");

	cblas_dscal(3, 1.0 / s, element_, 1);
   	return *this;
}

dvec3
dvec3::operator * (double s) const
{
	dvec3 vr(*this);
	cblas_dscal(3, s, vr.element_, 1);
   	return vr;
}

dvec3
dvec3::operator / (double s) const
{
	assert ( s != 0.0 &&
			"dvec3 dvec3::operator / (double)");

	dvec3 vr(*this);
	cblas_dscal(3, 1.0 / s, vr.element_, 1);
   	return vr;
}

dvec3
dvec3::operator + (const dvec3& v) const
{
	dvec3 vr(*this);

	for (i = 0; i < 3; i++) {
		vr.element_[i] += v.element_[i];
	}

   	return vr;
}

dvec3
dvec3::operator - (const dvec3& v) const
{
	dvec3 vr(*this);
	for (i = 0; i < 3; i++) {
		vr.element_[i] -= v.element_[i];
	}

   	return vr;
}

double
dvec3::operator * (const dvec3& v) const
{
	return cblas_ddot(3, element_, 1, v.element_, 1);
}

dvec3
dvec3::operator * (const dmat3 &m) const
{
	dvec3 vr;
	cblas_dgemv(CblasColMajor, CblasTrans, 3, 3, 1.0, m.element_, 3, element_, 1, 0.0, vr.element_, 1);
	return vr;
}

dmat3
dvec3::operator ^ (const dvec3 &v) const  // matrix = col vec x row vec
{
	dmat3 mr;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			mr.element_[i + 3 * j] = element_[i] * v.element_[j];
		}
	}

	return mr;
}

dvec3
dvec3::operator & (const dvec3 &v) const // cross product
{
	dvec3 vr;
	vr.element_[0] = element_[1] * v.element_[2] - element_[2] * v.element_[1];
	vr.element_[1] = element_[2] * v.element_[0] - element_[0] * v.element_[2];
	vr.element_[2] = element_[0] * v.element_[1] - element_[1] * v.element_[0];
	return vr;
}

dvec3
rmcp::operator * (double s, const dvec3 &v)
{
	dvec3 vr(v);
	cblas_dscal(3, s, vr.element_, 1);
	return vr;
}

dvec3
rmcp::cross(const dvec3 &x, const dvec3 &y)
{
	dvec3 vr;

	vr.element_[0] = x.element_[1] * y.element_[2] - x.element_[2] * y.element_[1];
	vr.element_[1] = x.element_[2] * y.element_[0] - x.element_[0] * y.element_[2];
	vr.element_[2] = x.element_[0] * y.element_[1] - x.element_[1] * y.element_[0];

	return vr;
}

dmat3
rmcp::cross(const dvec3 &v, const dmat3 &m)
{
	dmat3 mr;

	mr.element_[0] = v.element_[1] * m.element_[2] - v.element_[2] * m.element_[1];
	mr.element_[1] = v.element_[2] * m.element_[0] - v.element_[0] * m.element_[2];
	mr.element_[2] = v.element_[0] * m.element_[1] - v.element_[1] * m.element_[0];

	mr.element_[3] = v.element_[1] * m.element_[5] - v.element_[2] * m.element_[4];
	mr.element_[4] = v.element_[2] * m.element_[3] - v.element_[0] * m.element_[5];
	mr.element_[5] = v.element_[0] * m.element_[4] - v.element_[1] * m.element_[3];

	mr.element_[6] = v.element_[1] * m.element_[8] - v.element_[2] * m.element_[7];
	mr.element_[7] = v.element_[2] * m.element_[6] - v.element_[0] * m.element_[8];
	mr.element_[8] = v.element_[0] * m.element_[7] - v.element_[1] * m.element_[6];

	return mr;

}

dmat3
rmcp::cross(const dvec3 &v)
{
	dmat3 mr;
	
    mr.element_[0] = 0.0;
   	mr.element_[1] = v.element_[2];
	mr.element_[2] = -v.element_[1]; 
	mr.element_[3] = -v.element_[2];
	mr.element_[4] = 0.0; 
	mr.element_[5] = v.element_[0];
	mr.element_[6] = v.element_[1];
	mr.element_[7] = -v.element_[0];
	mr.element_[8] = 0.0; 

	return mr;	
}

dvec3
rmcp::inv_cross(const dmat3 &m)
{
	dvec3 vr;
	vr.element_[0] = m.element_[5];
	vr.element_[1] = m.element_[6];
	vr.element_[2] = m.element_[1];
	return vr;
}

dvec3 
rmcp::quaternion_rotate(const dvec3 &A, const dvec3 &P, double theta)	// axis A에 대한 point P의 theta만큼의 회전
{
	double Q[4];
	double c,s;
	dmat3 mr;
	dvec3 vr;

	c = cos(theta * 0.5);
	s = sin(theta * 0.5);

	Q[0] = c;
	Q[1] = s*A.element_[0];
	Q[2] = s*A.element_[1];
	Q[3] = s*A.element_[2];

	mr.element_[0] = 1.0 - 2.0*Q[2]*Q[2] - 2.0*Q[3]*Q[3];
	mr.element_[1] = 2.0*Q[1]*Q[2] + 2.0*Q[0]*Q[3];
	mr.element_[2] = 2.0*Q[1]*Q[3] - 2.0*Q[0]*Q[2];

	mr.element_[3] = 2.0*Q[1]*Q[2] - 2.0*Q[0]*Q[3];
	mr.element_[4] = 1.0 - 2.0*Q[1]*Q[1] - 2.0*Q[3]*Q[3];
	mr.element_[5] = 2.0*Q[2]*Q[3] + 2.0*Q[0]*Q[1];

	mr.element_[6] = 2.0*Q[1]*Q[3] + 2.0*Q[0]*Q[2];
	mr.element_[7] = 2.0*Q[2]*Q[3] - 2.0*Q[0]*Q[1];
	mr.element_[8] = 1.0 - 2.0*Q[1]*Q[1] - 2.0*Q[2]*Q[2];

	vr = mr * P;

	return vr;
}	
/* dmat3 */

dmat3::dmat3() : rmcp::dgqmatrix(3, 0.0)
{
}

dmat3::dmat3(const double element) : rmcp::dgqmatrix(3, element)
{
}

dmat3::dmat3(const double x, const double y, const double z) : rmcp::dgqmatrix(3, 0.0)
{
	/* x, y, z -> diagonal element */
	element_[0] = x;
	element_[4] = y;
	element_[8] = z;
}

dmat3::dmat3(const double element[]) : rmcp::dgqmatrix(3, element)
{
}

dmat3::dmat3(const dmat3 &m) : rmcp::dgqmatrix(3, m.element_)
{
}

const dmat3 &
dmat3::operator + (void) const
{
   	return *this;
}

dmat3 
dmat3::operator - (void) const
{
	dmat3 mr;

	for (i = 0; i < 9; i++) {
		mr.element_[i] = - element_[i];
	}

	return mr; 
}

dmat3 
dmat3::operator ~ (void) const
{ 
	dmat3 mr;

	for (i = 0; i < column_; i++) {
		for (j = 0; j < row_; j++) {
			mr.element_[i + j * column_] = element_[j + i * row_];
		}
	}

	return mr;
}

dmat3 &
dmat3::operator = (const dmat3 &m)
{
	cblas_dcopy(9, m.element_, 1, element_, 1);

	return *this;
}

dmat3 &
dmat3::operator = (const double element[])
{
	cblas_dcopy(9, element, 1, element_, 1);
	return *this;
}

dmat3 &
dmat3::operator += (const dmat3 &m)
{
	for (i = 0; i < 9; i++) {
		element_[i] += m.element_[i];
	}

	return *this;
}

dmat3 &
dmat3::operator -= (const dmat3 &m)
{
	assert (m.row_ == row_ &&
		   	m.column_ == column_ &&
		   "const dmat3 &dmat3::operator -= (const dmat3 &)");

	for (i = 0; i < 9; i++) {
		element_[i] -= m.element_[i];
	}

	return *this;
}

dmat3 &
dmat3::operator *= (double s)
{
	cblas_dscal(9, s, element_, 1);
	return *this;
}

dmat3 &
dmat3::operator *= (const dmat3 &m)
{
	rmcp::dmat3 mt(*this);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, mt.element_, 3, m.element_, 3, 0.0, element_, 3);
	return *this;
}

dmat3 &
dmat3::operator /= (double s)
{
	assert (s != 0.0 &&
			"const dmat3 &dmat3::operator /= (double)");

	cblas_dscal(9, 1.0 / s, element_, 1);
	return *this;
}

dmat3 
dmat3::operator + (const dmat3 &m) const
{
	dmat3 mr(*this);

	for (i = 0; i < 9; i++) {
		mr.element_[i] += m.element_[i];
	}

	return mr;
}

dmat3 
dmat3::operator - (const dmat3 &m) const
{
	dmat3 mr(*this);

	for (i = 0; i < 9; i++) {
		mr.element_[i] -= m.element_[i];
	}

	return mr;
}

dmat3 
dmat3::operator * (const double s) const
{
	dmat3 mr(*this);
	cblas_dscal(9, s, mr.element_, 1);
	return mr;
}

dmat3 
dmat3::operator / (const double s) const
{
	assert (s != 0.0 &&
			"dmat3 dmat3::operator / (double s) const");

	dmat3 mr(*this);
	cblas_dscal(9, 1.0 / s, element_, 1);
	return mr;
}

dmat3 
dmat3::operator * (const dmat3 &m) const
{
	dmat3 mr;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, element_, 3, m.element_, 3, 0.0, mr.element_, 3);
	return mr;
}

dvec3
dmat3::operator * (const dvec3 &v) const
{
	dvec3 vr;
	cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, element_, 3, v.element_, 1, 0.0, vr.element_, 1);
	return vr;
}

double
dmat3::det(void)
{
	return element_[0] * element_[4] * element_[8] + 
		element_[3] * element_[7] * element_[2] +
		element_[6] * element_[1] * element_[5] -
		element_[6] * element_[4] * element_[2] -
		element_[1] * element_[3] * element_[8] -
		element_[0] * element_[5] * element_[7];
}

void
dmat3::transpose(void)
{
	static double tmp[3];

	tmp[0] = element_[1];
	tmp[1] = element_[2];
	tmp[2] = element_[5];

	element_[1] = element_[3];
	element_[2] = element_[6];
	element_[5] = element_[7];

	element_[3] = tmp[0];
	element_[6] = tmp[1];
	element_[7] = tmp[2];
}

void
dmat3::inv(void) 
{
	static int *ipiv;
	static int info;
	static double *work;

	ipiv = new int[row_];
	work = new double[row_];

	dgetrf(&row_, &column_, element_, &row_, ipiv, &info);
	dgetri(&row_, element_, &row_, ipiv, work, &row_, &info);

 	delete[] ipiv;
	delete[] work;

	return;
}

dmat3
rmcp::operator * (double s, const dmat3 &m)
{
	dmat3 mr(m);
	cblas_dscal(9, s, mr.element_, 1);
	return mr;
}

dmat3
rmcp::inv(const dmat3 &m)
{
	dmat3 mr(m);

	static int ipiv[3];
	static int info;
	static double work[3];

	dgetrf(&mr.row_, &mr.column_, mr.element_, &mr.row_, ipiv, &info);
	dgetri(&mr.row_, mr.element_, &mr.row_, ipiv, work, &mr.row_, &info);

	return mr;
}

dmat3
rmcp::rotate(const dvec3 &ang, rmcp::ROT_TYPE type)
{
	return rmcp::rotate(ang.element_[0], ang.element_[1], ang.element_[2], type);
}

dmat3
rmcp::rotate(const double ang1, const double ang2, const double ang3, rmcp::ROT_TYPE type)
{
	static double sa, ca, sb, cb, sg, cg;

	dmat3 mr;

	sa = sin(ang1);
	ca = cos(ang1);
	sb = sin(ang2);
	cb = cos(ang2);
	sg = sin(ang3);
	cg = cos(ang3);

	switch (type) {
		case rmcp::ROT_RPY :
			mr.element_[0] = cg * cb;
			mr.element_[1] = sg * cb;
			mr.element_[2] = - sb;
			mr.element_[3] = cg * sb * sa - sg * ca;
			mr.element_[4] = sg * sb * sa + cg * ca;
			mr.element_[5] = cb * sa;
			mr.element_[6] = cg * sb * ca + sg * sa;
			mr.element_[7] = sg * sb * ca - cg * sa;
			mr.element_[8] = cb * ca;
			break;
		case rmcp::ROT_ZXZ:
			mr.element_[0] = ca * cg - sa * cb * sg;
			mr.element_[1] = sa * cg + ca * cb * sg;
			mr.element_[2] = sb * sg;
			mr.element_[3] = - ca * sg - sa * cb * cg;
			mr.element_[4] = - sa * sg + ca * cb * cg;
			mr.element_[5] = sb * cg;
			mr.element_[6] = sa * sb;
			mr.element_[7] = - ca * cb;
			mr.element_[8] = cb;
			break;
		case rmcp::ROT_ZYZ:
			mr.element_[0] = ca * cb * cg - sa * sg;
			mr.element_[1] = sa * cb * cg + ca * sg;
			mr.element_[2] = - sb * cg;
			mr.element_[3] = - ca * cb * sg - sa * cg;
			mr.element_[4] = - sa * cb * sg + ca * cg;
			mr.element_[5] = sb * sg;
			mr.element_[6] = ca * sb;
			mr.element_[7] = sa * sb;
			mr.element_[8] = cb;
			break;
	}

	return mr;
}

dmat3
rmcp::rotate(const double ang, dvec3 &axis)	// angle-axis
{
	static double st, ct, vt, t0, t1, t2;

	axis.normalize();

	st = sin(ang);
	ct = cos(ang);
	vt = 1.0 - ct;
	t0 = axis.element_[2] * st;
	t1 = axis.element_[1] * st;
	t2 = axis.element_[0] * st;

	rmcp::dmat3 mr;
	mr.element_[0] = axis.element_[0] * axis.element_[0] * vt + ct;
	mr.element_[1] = axis.element_[0] * axis.element_[1] * vt + t0;
	mr.element_[2] = axis.element_[0] * axis.element_[2] * vt - t1;
	mr.element_[3] = axis.element_[0] * axis.element_[1] * vt - t0;
	mr.element_[4] = axis.element_[1] * axis.element_[1] * vt + ct;
	mr.element_[5] = axis.element_[1] * axis.element_[2] * vt + t2;
	mr.element_[6] = axis.element_[0] * axis.element_[2] * vt + t1;
	mr.element_[7] = axis.element_[1] * axis.element_[2] * vt - t2;
	mr.element_[8] = axis.element_[2] * axis.element_[2] * vt + ct;
	return mr;
}

dvec3
rmcp::inv_rotate(const dmat3 &R, ROT_TYPE type)
{
	dvec3 vr;
	rmcp::inv_rotate(vr.element_[0], vr.element_[1], vr.element_[2], R, type);
	return vr;
}

void
rmcp::inv_rotate(double &ang1, double &ang2, double &ang3, const dmat3 &R, ROT_TYPE type)
{
	switch (type) {
		case rmcp::ROT_RPY :
			ang1 = atan2(R.element_[1], R.element_[0]);
			ang2 = asin(-R.element_[2]);
			ang3 = atan2(R.element_[5], R.element_[8]);
			break;
		case rmcp::ROT_ZXZ:
			ang1 = atan2(R.element_[6], -R.element_[7]);
			ang2 = acos(R.element_[8]);
			ang3 = atan2(R.element_[2], R.element_[5]);
			break;
		case rmcp::ROT_ZYZ:
			ang1 = atan2(R.element_[7], R.element_[6]);
			ang2 = acos(R.element_[8]);
			ang3 = atan2(R.element_[5], -R.element_[2]);
			break;
	}

}

void
rmcp::inv_rotate(double &ang, dvec3 &axis, const dmat3 &R)	// angle-axis
{
	ang = acos(0.5 * (R.element_[0] + R.element_[4] + R.element_[8] - 1.0));

	if (fabs(ang) < rmcp::TINY) {
		ang = 0.0;
		axis.element_[0] = 0.0;
		axis.element_[1] = 0.0;
		axis.element_[2] = 0.0;
	}

	double cof = ang / (2.0 * sin(ang));

	axis.element_[0] = cof * (R.element_[5] - R.element_[7]);
	axis.element_[1] = cof * (R.element_[6] - R.element_[2]);
	axis.element_[2] = cof * (R.element_[1] - R.element_[3]);
}

dhmat::dhmat() : rmcp::dgqmatrix(4)
{
	element_[0] = 1.0;
	element_[1] = 0.0;
	element_[2] = 0.0;
	element_[3] = 0.0;
	element_[4] = 0.0;
	element_[5] = 1.0;
	element_[6] = 0.0;
	element_[7] = 0.0;
	element_[8] = 0.0;
	element_[9] = 0.0;
	element_[10] = 1.0;
	element_[11] = 0.0;
	element_[12] = 0.0;
	element_[13] = 0.0;
	element_[14] = 0.0;
	element_[15] = 1.0;
}

dhmat::dhmat(const double element[]) : rmcp::dgqmatrix(4)
{
	element_[0] = element[0];
	element_[1] = element[1];
	element_[2] = element[2];
	element_[3] = 0.0;
	element_[4] = element[3];
	element_[5] = element[4];
	element_[6] = element[5];
	element_[7] = 0.0;
	element_[8] = element[6];
	element_[9] = element[7];
	element_[10] = element[8];
	element_[11] = 0.0;
	element_[12] = element[9];
	element_[13] = element[10];
	element_[14] = element[11];
	element_[15] = 1.0;
}

dhmat::dhmat(const rmcp::dmat3 &R, const rmcp::dvec3 &v)
	: rmcp::dgqmatrix(4)
{
	element_[0] = R.element_[0];
	element_[1] = R.element_[1];
	element_[2] = R.element_[2];
	element_[3] = 0.0;
	element_[4] = R.element_[3];
	element_[5] = R.element_[4];
	element_[6] = R.element_[5];
	element_[7] = 0.0;
	element_[8] = R.element_[6];
	element_[9] = R.element_[7];
	element_[10] = R.element_[8];
	element_[11] = 0.0;
	element_[12] = v.element_[0];
	element_[13] = v.element_[1];
	element_[14] = v.element_[2];
	element_[15] = 1.0;
}

dhmat::dhmat(const dhmat &H) : dgqmatrix(4, H.element_)
{
}

dhmat &
dhmat::operator = (const dhmat &m)
{
	cblas_dcopy(16, m.element_, 1, element_, 1);

	return *this;
}

dhmat &
dhmat::operator = (const double element[])
{
	element_[0] = element[0];
	element_[1] = element[1];
	element_[2] = element[2];
	element_[3] = 0.0;
	element_[4] = element[3];
	element_[5] = element[4];
	element_[6] = element[5];
	element_[7] = 0.0;
	element_[8] = element[6];
	element_[9] = element[7];
	element_[10] = element[8];
	element_[11] = 0.0;
	element_[12] = element[9];
	element_[13] = element[10];
	element_[14] = element[11];
	element_[15] = 1.0;

	return *this;
}

dhmat &
dhmat::operator *= (const dhmat &m)
{
	static double x0, x1, x2;

	element_[12] += 
		element_[0] * m.element_[12] + 
		element_[4] * m.element_[13] + 
		element_[8] * m.element_[14];
	element_[13] +=
	   	element_[1] * m.element_[12] + 
		element_[5] * m.element_[13] + 
		element_[9] * m.element_[14];
	element_[14] += 
		element_[2] * m.element_[12] + 
		element_[6] * m.element_[13] + 
		element_[10] * m.element_[14];

	x0 = 
		element_[0] * m.element_[0] + 
		element_[4] * m.element_[1] + 
		element_[8] * m.element_[2];
	x1 =
	   	element_[0] * m.element_[4] + 
		element_[4] * m.element_[5] + 
		element_[8] * m.element_[6];
	x2 =
	   	element_[0] * m.element_[8] + 
		element_[4] * m.element_[9] + 
		element_[8] * m.element_[10];

	element_[0] = x0;
   	element_[4] = x1;
   	element_[8] = x2;

	x0 = 
		element_[1] * m.element_[0] + 
		element_[5] * m.element_[1] + 
		element_[9] * m.element_[2];
	x1 = 
		element_[1] * m.element_[4] +
	   	element_[5] * m.element_[5] + 
		element_[9] * m.element_[6];
	x2 = 
		element_[1] * m.element_[8] + 
		element_[5] * m.element_[9] + 
		element_[9] * m.element_[10];

	element_[1] = x0;
   	element_[5] = x1;
   	element_[9] = x2;

	x0 = 
		element_[2] * m.element_[0] +
	   	element_[6] * m.element_[1] + 
		element_[10] * m.element_[2];
	x1 = 
		element_[2] * m.element_[4] + 
		element_[6] * m.element_[5] +
	   	element_[10] * m.element_[6];
	x2 = 
		element_[2] * m.element_[8] +
	   	element_[6] * m.element_[9] + 
		element_[10] * m.element_[10];

	element_[2] = x0;
   	element_[6] = x1;
   	element_[10] = x2;

	return  *this;
}

dhmat &
dhmat::operator /= (const dhmat &m)
{
	static double tmp[9];

	tmp[0] = 
		element_[0] * m.element_[0] + 
		element_[4] * m.element_[4] + 
		element_[8] * m.element_[8];
	tmp[1] =
	   	element_[1] * m.element_[0] + 
		element_[5] * m.element_[4] + 
		element_[9] * m.element_[8];
	tmp[2] = 
		element_[2] * m.element_[0] + 
		element_[6] * m.element_[4] + 
		element_[10] * m.element_[8];
	tmp[3] =
		element_[0] * m.element_[1] + 
		element_[4] * m.element_[5] + 
		element_[8] * m.element_[9];
	tmp[4] = 
		element_[1] * m.element_[1] + 
		element_[5] * m.element_[5] + 
		element_[9] * m.element_[9];
	tmp[5] = 
		element_[2] * m.element_[1] + 
		element_[6] * m.element_[5] + 
		element_[10] * m.element_[9];
	tmp[6] = 
		element_[0] * m.element_[2] + 
		element_[4] * m.element_[6] + 
		element_[8] * m.element_[10];
	tmp[7] = 
		element_[1] * m.element_[2] + 
		element_[5] * m.element_[6] + 
		element_[9] * m.element_[10];
	tmp[8] =
		element_[2] * m.element_[2] +
		element_[6] * m.element_[6] +
		element_[10] * m.element_[10];

	element_[0] = tmp[0];
   	element_[1] = tmp[1];
   	element_[2] = tmp[2];
	element_[4] = tmp[3];
   	element_[5] = tmp[4];
   	element_[6] = tmp[5];
	element_[8] = tmp[6];
   	element_[9] = tmp[7];
   	element_[10] = tmp[8]; 
	element_[12] -=
		tmp[0] * m.element_[12] + 
		tmp[3] * m.element_[13] + 
		tmp[6] * m.element_[14];
	element_[13] -=
		tmp[1] * m.element_[12] + 
		tmp[4] * m.element_[13] + 
		tmp[7] * m.element_[14];
	element_[14] -=
		tmp[2] * m.element_[12] + 
		tmp[5] * m.element_[13] + 
		tmp[8] * m.element_[14];

	return *this;
}

dhmat 
dhmat::operator * (const dhmat &m) const
{

	dhmat mr;
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 4, 4, 4, 1.0, element_, 4, m.element_, 4, 0.0, mr.element_, 4);
	return mr;
}

dhmat
dhmat::operator / (const dhmat& m) const
{
	static double tmp[9];

	tmp[0] = 
		element_[0] * m.element_[0] + 
		element_[4] * m.element_[4] + 
		element_[8] * m.element_[8];
	tmp[1] = 
		element_[1] * m.element_[0] + 
		element_[5] * m.element_[4] + 
		element_[9] * m.element_[8];
	tmp[2] = 
		element_[2] * m.element_[0] + 
		element_[6] * m.element_[4] + 
		element_[10] * m.element_[8];
	tmp[3] = 
		element_[0] * m.element_[1] + 
		element_[4] * m.element_[5] + 
		element_[8] * m.element_[9];
	tmp[4] = 
		element_[1] * m.element_[1] + 
		element_[5] * m.element_[5] + 
		element_[9] * m.element_[9];
	tmp[5] = 
		element_[2] * m.element_[1] + 
		element_[6] * m.element_[5] + 
		element_[10] * m.element_[9];
	tmp[6] = 
		element_[0] * m.element_[2] + 
		element_[4] * m.element_[6] + 
		element_[8] * m.element_[10];
	tmp[7] = 
		element_[1] * m.element_[2] + 
		element_[5] * m.element_[6] + 
		element_[9] * m.element_[10];
	tmp[8] = 
		element_[2] * m.element_[2] + 
		element_[6] * m.element_[6] + 
		element_[10] * m.element_[10];

	dhmat mr;

	mr.element_[0] = tmp[0];
	mr.element_[1] = tmp[1];
	mr.element_[2] = tmp[2];
	mr.element_[4] = tmp[3];
	mr.element_[5] = tmp[4];
	mr.element_[6] = tmp[5];
	mr.element_[8] = tmp[6];
	mr.element_[9] = tmp[7];
	mr.element_[10] = tmp[8];
	mr.element_[12] = element_[12] - tmp[0] * m.element_[12] - tmp[3] * m.element_[13] - tmp[6] * m.element_[14];
	mr.element_[13] = element_[13] - tmp[1] * m.element_[12] - tmp[4] * m.element_[13] - tmp[7] * m.element_[14];
	mr.element_[14] = element_[14] - tmp[2] * m.element_[12] - tmp[5] * m.element_[13] - tmp[8] * m.element_[14];

	return mr;
}

dvec3
dhmat::operator * (const dvec3 &v) const
{
	dvec3 vr;

	vr.element_[0] = 
		element_[0] * v.element_[0] + 
		element_[4] * v.element_[1] + 
		element_[8] * v.element_[3] + element_[12];
	vr.element_[1] = 
		element_[1] * v.element_[0] + 
		element_[5] * v.element_[1] + 
		element_[9] * v.element_[3] + element_[13];
	vr.element_[2] = 
		element_[2] * v.element_[0] + 
		element_[6] * v.element_[1] + 
		element_[10] * v.element_[3] + element_[14];

	return vr;
}

void
dhmat::set(const double element[])
{
	element_[0] = element[0];
	element_[1] = element[1];
	element_[2] = element[2];
	element_[3] = 0.0;
	element_[4] = element[3];
	element_[5] = element[4];
	element_[6] = element[5];
	element_[7] = 0.0;
	element_[8] = element[6];
	element_[9] = element[7];
	element_[10] = element[8];
	element_[11] = 0.0;
	element_[12] = element[9];
	element_[13] = element[10];
	element_[14] = element[11];
	element_[15] = 1.0;

	return;
}

void
dhmat::set(const dmat3 &R, const dvec3 &v)
{
	element_[0] = R.element_[0];
	element_[1] = R.element_[1];
	element_[2] = R.element_[2];
	element_[3] = 0.0;
	element_[4] = R.element_[3];
	element_[5] = R.element_[4];
	element_[6] = R.element_[5];
	element_[7] = 0.0;
	element_[8] = R.element_[6];
	element_[9] = R.element_[7];
	element_[10] = R.element_[8];
	element_[11] = 0.0;
	element_[12] = v.element_[0];
	element_[13] = v.element_[1];
	element_[14] = v.element_[2];
	element_[15] = 1.0;
}


void
dhmat::set(const dhmat &H)
{
	cblas_dcopy(16, H.element_, 1, element_, 1);
}

void
dhmat::set(const dmat3 &R)
{
	element_[0] = R.element_[0];
	element_[1] = R.element_[1];
	element_[2] = R.element_[2];
	element_[4] = R.element_[3];
	element_[5] = R.element_[4];
	element_[6] = R.element_[5];
	element_[8] = R.element_[6];
	element_[9] = R.element_[7];
	element_[10] = R.element_[8];
}

void
dhmat::set(const dvec3 &v)
{
	element_[12] = v.element_[0];
	element_[13] = v.element_[1];
	element_[14] = v.element_[2];
}

void
dhmat::inv(void) 
{
	static double *element;
	element = new double[16];

	element[0] = element_[0];
	element[1] = element_[4];
	element[2] = element_[8];
	element[3] = 0.0;
	element[4] = element_[1];
   	element[5] = element_[5];
   	element[6] = element_[9];
	element[7] = 0.0;
   	element[8] = element_[2];
   	element[9] = element_[6];
   	element[10] = element_[10];
	element[11] = 0.0;
	element[12] = - element_[0] * element_[12] - element_[1] * element_[13] - element_[2] * element_[14];
	element[13] = - element_[4] * element_[12] - element_[5] * element_[13] - element_[6] * element_[14];
	element[14] = - element_[8] * element_[12] - element_[9] * element_[13] - element_[10] * element_[14];
	element[15] = 1.0;
	delete[] this->element_;
	this->element_ = element;

	return;
}

dhmat
rmcp::inv(const dhmat &m)
{
	dhmat mr;

	mr.element_[1] = m.element_[4];
	mr.element_[2] = m.element_[8];
	mr.element_[4] = m.element_[1];
   	mr.element_[5] = m.element_[5];
   	mr.element_[6] = m.element_[9];
   	mr.element_[8] = m.element_[2];
   	mr.element_[9] = m.element_[6];
   	mr.element_[10] = m.element_[10];
	mr.element_[12] = - m.element_[0] * m.element_[12] - m.element_[1] * m.element_[13] - m.element_[2] * m.element_[14];
	mr.element_[13] = - m.element_[4] * m.element_[12] - m.element_[5] * m.element_[13] - m.element_[6] * m.element_[14];
	mr.element_[14] = - m.element_[8] * m.element_[12] - m.element_[9] * m.element_[13] - m.element_[10] * m.element_[14];
	
	return mr;
}


