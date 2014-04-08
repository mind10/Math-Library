
#include "rmcp/dso3.h"
#include "rmcp/dse3.h"

#include "rmcp/dvector.h"

#include <cmath>
#include <cassert>

using namespace rmcp;

static int i, j;


dso3::dso3(void) : dssmatrix(3, 0.0)
{

}

dso3::dso3(double d) : dssmatrix(3)
{
   	element_[0] = 0.0;
   	element_[1] = d;
	element_[2] = -d; 
	element_[3] = -d;
	element_[4] = 0.0; 
	element_[5] = d;
	element_[6] = d;
	element_[7] = -d;
	element_[8] = 0.0; 
}

dso3::dso3(const double v[]) : dssmatrix(3)
{
    element_[0] = 0.0;
   	element_[1] = v[2];
	element_[2] = -v[1]; 
	element_[3] = -v[2];
	element_[4] = 0.0; 
	element_[5] = v[0];
	element_[6] = v[1];
	element_[7] = -v[0];
	element_[8] = 0.0; 
}

dso3::dso3(double v0, double v1, double v2) : dssmatrix(3)
{ 
    element_[0] = 0.0;
   	element_[1] = v2;
	element_[2] = -v1; 
	element_[3] = -v2;
	element_[4] = 0.0; 
	element_[5] = v0;
	element_[6] = v1;
	element_[7] = -v0;
	element_[8] = 0.0; 
}

dso3::dso3(const rmcp::dvector &v) : dssmatrix(3)
{
	assert (v.size_ == 3 &&
			"dso3::dso3(const rmcp::dvector &)");

    element_[0] = 0.0;
   	element_[1] = v.element_[2];
	element_[2] = -v.element_[1]; 
	element_[3] = -v.element_[2];
	element_[4] = 0.0; 
	element_[5] = v.element_[0];
	element_[6] = v.element_[1];
	element_[7] = -v.element_[0];
	element_[8] = 0.0; 
}

dso3::dso3(const dso3 &w) : dssmatrix(3)
{
	cblas_dcopy(9, w.element_, 1, element_, 1);
}

const dso3&
dso3::operator + (void) const
{
   	return *this;
}

dso3 
dso3::operator - (void) const	
{
	dso3 wr;

   	wr.element_[1] = -element_[1];
	wr.element_[2] = -element_[2];
	wr.element_[3] = -element_[3];
	wr.element_[5] = -element_[5];
	wr.element_[6] = -element_[6];
	wr.element_[7] = -element_[7];

   	return wr;
} 

dso3 &
dso3::operator = (const dso3 &w)
{
   	element_[1] = w.element_[1];
	element_[2] = w.element_[2]; 
	element_[3] = w.element_[3];
	element_[5] = w.element_[5];
	element_[6] = w.element_[6];
	element_[7] = w.element_[7];

   	return *this;
}

dso3 &
dso3::operator = (const rmcp::dvector &v)
{
	assert (v.size_ == 3 &&
			"dse3 &dse3::operator = (const rmcp::dvector &)");

   	element_[1] = v.element_[2];
	element_[2] = -v.element_[1]; 
	element_[3] = -v.element_[2];
	element_[5] = v.element_[0];
	element_[6] = v.element_[1];
	element_[7] = -v.element_[0];

   	return *this;
}

dso3 &
dso3::operator = (const double v[])
{
   	element_[1] = v[2];
	element_[2] = -v[1]; 
	element_[3] = -v[2];
	element_[5] = v[0];
	element_[6] = v[1];
	element_[7] = -v[0];

   	return *this;
}

dso3 &
dso3::operator += (const dso3 &w)
{
   	element_[1] += w.element_[1];
	element_[2] += w.element_[2]; 
	element_[3] += w.element_[3];
	element_[5] += w.element_[5];
	element_[6] += w.element_[6];
	element_[7] += w.element_[7];

	return *this;
}

dso3&
dso3::operator -= (const dso3 &w)
{
   	element_[1] -= w.element_[1];
	element_[2] -= w.element_[2]; 
	element_[3] -= w.element_[3];
	element_[5] -= w.element_[5];
	element_[6] -= w.element_[6];
	element_[7] -= w.element_[7];

   	return *this;
}

dso3 & 
dso3::operator *= (double s)
{
  	element_[1] *= s;
	element_[2] *= s;
	element_[3] *= s;
	element_[5] *= s;
	element_[6] *= s;
	element_[7] *= s;

   	return *this;
}

dso3 &
dso3::operator /= (double s)
{
	assert ( s != 0.0 &&
		"dso3 &dso3::operator /= (double)");

	element_[1] /= s;
	element_[2] /= s;
	element_[3] /= s;
	element_[5] /= s;
	element_[6] /= s;
	element_[7] /= s;

	return *this;
}

dso3
dso3::operator + (const dso3 &w) const
{
	dso3 wr(*this);

  	wr.element_[1] += w.element_[1];
	wr.element_[2] += w.element_[2]; 
	wr.element_[3] += w.element_[3];
	wr.element_[5] += w.element_[5];
	wr.element_[6] += w.element_[6];
	wr.element_[7] += w.element_[7];

	return wr;
}

dso3
dso3::operator - (const dso3 &w) const
{
	dso3 wr(*this);

  	wr.element_[1] -= w.element_[1];
	wr.element_[2] -= w.element_[2]; 
	wr.element_[3] -= w.element_[3];
	wr.element_[5] -= w.element_[5];
	wr.element_[6] -= w.element_[6];
	wr.element_[7] -= w.element_[7];

	return wr;
}

dso3
dso3::operator * (double s) const
{
	dso3 wr(*this);

	wr.element_[1] *= s;
	wr.element_[2] *= s;
	wr.element_[3] *= s;
	wr.element_[5] *= s;
	wr.element_[6] *= s;
	wr.element_[7] *= s;

	return wr;
}

dso3
dso3::operator / (double s) const
{
	assert (s != 0.0 &&
		"dso3 dso3::operator / (double) const");

	dso3 wr(*this);

	wr.element_[1] /= s;
	wr.element_[2] /= s;
	wr.element_[3] /= s;
	wr.element_[5] /= s;
	wr.element_[6] /= s;
	wr.element_[7] /= s;

	return wr;
}

bool
dso3::operator == (const dso3 &w) const
{
   	return
	   	(element_[1] == w.element_[1]) &&
	   	(element_[2] == w.element_[2]) &&
	   	(element_[5] == w.element_[5]);
}

void
dso3::set(double v0, double v1, double v2)
{
   	element_[1] = v2;
	element_[2] = -v1; 
	element_[3] = -v2;
	element_[5] = v0;
	element_[6] = v1;
	element_[7] = -v0;
}

void
dso3::get(rmcp::dvector &w)		/* get [w], w : 3-dimensional vector */
{
	w.element_[0] = element_[5];
	w.element_[1] = element_[6];
	w.element_[2] = element_[1];
}

void
dso3::set(const double v[])
{
   	element_[1] = v[2];
	element_[2] = -v[1]; 
	element_[3] = -v[2];
	element_[5] = v[0];
	element_[6] = v[1];
	element_[7] = -v[0];
}

/* the norm of dso3 means the norm of 3-dimensional vector */
double
dso3::norm(void) const
{
	return std::sqrt(element_[5] * element_[5] +
		element_[6] * element_[6] +
		element_[3] * element_[3]);
}

void
dso3::log(const rmcp::dSO3& R)
{
	static double trR;
	static double w0, w1, w2;

	trR = R.trace();

	if (trR > 3 - rmcp::LOG_EPS) {
		// theta = 0
		// w = anything
		// so, [w]*theta = 0
		element_[1] = 0.0;
		element_[2] = 0.0; 
		element_[3] = 0.0;
		element_[5] = 0.0;
		element_[6] = 0.0;
		element_[7] = 0.0;
		return;
	}
	else if (trR < -1 + rmcp::LOG_EPS) {
		// theta = pi;
		w0 = std::sqrt((R.element_[0] + 1.0) * 0.5);
		w1= std::sqrt((R.element_[4] + 1.0) * 0.5);
		w2 = std::sqrt((R.element_[8] + 1.0) * 0.5);
		if (std::fabs(w0) > rmcp::LOG_EPS) {
			w1 = R.element_[3] / (2.0 * w0);
			w2 = R.element_[6] / (2.0 * w0);
		}
		else if (std::fabs(w1) > rmcp::LOG_EPS) {
			w0 = R.element_[3] / (2.0 * w1);
			w2 = R.element_[7] / (2.0 * w1);
		}
		else {
			w0 = R.element_[6] / (2.0 * w2);
			w1 = R.element_[7] / (2.0 * w2);
		}
		
		w0 *= rmcp::PI;
		w1 *= rmcp::PI;
		w2 *= rmcp::PI;

		element_[1] = w2;
		element_[2] = -w1; 
		element_[3] = -w2;
		element_[5] = w0;
		element_[6] = w1;
		element_[7] = -w0;

		return;
	}
	else {
		// 0 < theta < pi
		static double theta;
        theta = std::acos(0.5 * (trR - 1.0));
		theta /= (sin(theta) * 2.0);

		w0 = theta * (R.element_[5] - R.element_[7]);
		w1 = theta * (R.element_[6] - R.element_[2]);
		w2 = theta * (R.element_[1] - R.element_[3]);

		element_[1] = w2;
		element_[2] = -w1; 
		element_[3] = -w2;
		element_[5] = w0;
		element_[6] = w1;
		element_[7] = -w0;

		return;
	}
}

dSO3::dSO3(void) : donmatrix(3)
{
   	element_[0] = 1.0;
	element_[1] = 0.0;
	element_[2] = 0.0;
	element_[3] = 0.0;
	element_[4] = 1.0;
	element_[5] = 0.0;
	element_[6] = 0.0;
	element_[7] = 0.0;
	element_[8] = 1.0;
}

dSO3::dSO3(const dSO3 &R) : donmatrix(3, R.element_)
{
}

dSO3::dSO3(const double element[]) : donmatrix(3, element)
{
}

dSO3::dSO3(double R0, double R1, double R2, double R3, double R4,
		   double R5, double R6, double R7, double R8) : donmatrix(3)
{
   	element_[0] = R0;
   	element_[1] = R1;
   	element_[2] = R2;
   	element_[3] = R3;
   	element_[4] = R4;
   	element_[5] = R5;
   	element_[6] = R6; 
	element_[7] = R7;
   	element_[8] = R8;
}

dSO3 &
dSO3::operator = (const dSO3 &R)
{
	cblas_dcopy(9, R.element_, 1, element_, 1);
   	return *this;
}

dSO3
dSO3::operator * (const dSO3 &R) const
{
	return dSO3(
		element_[0] * R.element_[0] + element_[3] * R.element_[1] + element_[6] * R.element_[2],
	   	element_[1] * R.element_[0] + element_[4] * R.element_[1] + element_[7] * R.element_[2],
	   	element_[2] * R.element_[0] + element_[5] * R.element_[1] + element_[8] * R.element_[2],
		element_[0] * R.element_[3] + element_[3] * R.element_[4] + element_[6] * R.element_[5],
	   	element_[1] * R.element_[3] + element_[4] * R.element_[4] + element_[7] * R.element_[5],
	   	element_[2] * R.element_[3] + element_[5] * R.element_[4] + element_[8] * R.element_[5],
		element_[0] * R.element_[6] + element_[3] * R.element_[7] + element_[6] * R.element_[8],
	   	element_[1] * R.element_[6] + element_[4] * R.element_[7] + element_[7] * R.element_[8],
	   	element_[2] * R.element_[6] + element_[5] * R.element_[7] +	element_[8] * R.element_[8] );
}

dSO3 
dSO3::operator % (const dSO3 &R) const
{
	return dSO3(
		element_[0] * R.element_[0] + element_[1] * R.element_[1] + element_[2] * R.element_[2],
		element_[3] * R.element_[0] + element_[4] * R.element_[1] + element_[5] * R.element_[2],
		element_[6] * R.element_[0] + element_[7] * R.element_[1] + element_[8] * R.element_[2],
		element_[0] * R.element_[3] + element_[1] * R.element_[4] + element_[2] * R.element_[5],
		element_[3] * R.element_[3] + element_[4] * R.element_[4] + element_[5] * R.element_[5],
		element_[6] * R.element_[3] + element_[7] * R.element_[4] + element_[8] * R.element_[5],
		element_[0] * R.element_[6] + element_[1] * R.element_[7] + element_[2] * R.element_[8],
		element_[3] * R.element_[6] + element_[4] * R.element_[7] + element_[5] * R.element_[8],
		element_[6] * R.element_[6] + element_[7] * R.element_[7] +	element_[8] * R.element_[8] );
}

dSO3 &
dSO3::operator *= (const dSO3 &R)
{
	*this = *this * R;
	return *this;
}

rmcp::dvector
dSO3::operator * (const rmcp::dvector &p) const
{
	assert (p.size_ == 3 &&
			"rmcp::dvector dSO3::operator * (const rmcp::dvector &p) const");
	
	rmcp::dvector vr(3);

	vr.element_[0] = 
		element_[0] * p.element_[0] + 
		element_[3] * p.element_[1] + 
		element_[6] * p.element_[2];
	vr.element_[1] = 
		element_[1] * p.element_[0] + 
		element_[4] * p.element_[1] + 
		element_[7] * p.element_[2];
	vr.element_[2] =
		element_[2] * p.element_[0] +
	   	element_[5] * p.element_[1] +
	   	element_[8] * p.element_[2];

	return vr;
}

dSO3
dSO3::operator ~ (void) const	
{
   	return dSO3(element_[0], element_[3], element_[6], 
			element_[1], element_[4], element_[7],
		   	element_[2], element_[5], element_[8]);
}

void
dSO3::inv(void)
{
	this->transpose();
}

/* if the norm of w is less than rmcp::EXP_EPS, return identity matrix */
void
dSO3::exp(const dso3 &w)
{
	static double w0, w1, w2;
	static double theta;
	static double st, ct, vt, t0, t1, t2;

	theta = w.norm();

	if (std::fabs(theta) < rmcp::EXP_EPS) {
		element_[0] = 1.0;
		element_[1] = 0.0;
		element_[2] = 0.0;
		element_[3] = 0.0;
		element_[4] = 1.0;
		element_[5] = 0.0;
		element_[6] = 0.0;
		element_[7] = 0.0;
		element_[8] = 1.0;
		return;
	}
	else {
		w0 = w.element_[5] / theta;
		w1 = w.element_[6] / theta;
		w2 = -w.element_[3] / theta;

		st = sin(theta);
		ct = cos(theta);
		vt = 1.0 - ct;
		t0 = w2 * st;
		t1 = w1 * st;
		t2 = w0 * st;

		element_[0] = w0 * w0 * vt + ct;
		element_[1] = w0 * w1 * vt + t0;
		element_[2] = w0 * w2 * vt - t1;
		element_[3] = w0 * w1 * vt - t0;
		element_[4] = w1 * w1 * vt + ct;
		element_[5] = w1 * w2 * vt + t2;
		element_[6] = w0 * w2 * vt + t1;
		element_[7] = w1 * w2 * vt - t2;
		element_[8] = w2 * w2 * vt + ct;

		return;
	}
}

/* given normalized w */
void
dSO3::exp(const dso3 &w, double theta)
{
	static double w0, w1, w2;
	static double st, ct, vt, t0, t1, t2;

	w0 = w.element_[5];
	w1 = w.element_[6];
	w2 = -w.element_[3];

	st = sin(theta);
	ct = cos(theta);
	vt = 1.0 - ct;
	t0 = w2 * st;
	t1 = w1 * st;
	t2 = w0 * st;

	element_[0] = w0 * w0 * vt + ct;
	element_[1] = w0 * w1 * vt + t0;
	element_[2] = w0 * w2 * vt - t1;
	element_[3] = w0 * w1 * vt - t0;
	element_[4] = w1 * w1 * vt + ct;
	element_[5] = w1 * w2 * vt + t2;
	element_[6] = w0 * w2 * vt + t1;
	element_[7] = w1 * w2 * vt - t2;
	element_[8] = w2 * w2 * vt + ct;

	return;
}

void
dSO3::rot_x(double t)
{
	static double c, s;
	c = cos(t);
	s = sin(t);
	element_[0] = 1.0;
	element_[1] = 0.0;
	element_[2] = 0.0;
	element_[3] = 0.0;
	element_[4] = c;
	element_[5] = s;
	element_[6] = 0.0;
	element_[7] = -s;
	element_[8] = c;
}

void 
dSO3::rot_y(double t)
{ 
	static double c, s;
	c = cos(t);
	s = sin(t);
	element_[0] = c;
	element_[1] = 0.0;
	element_[2] = -s;
	element_[3] = 0.0;
	element_[4] = 1.0;
	element_[5] = 0.0;
	element_[6] = s;
	element_[7] = 0.0;
	element_[8] = c;
}

void
dSO3::rot_z(double t)
{
	static double c, s;
	c = cos(t);
	s = sin(t);
	element_[0] = c; 
	element_[1] = s;
	element_[2] = 0.0;
	element_[3] = -s;
	element_[4] = c;
	element_[5] = 0.0;
	element_[6] = 0.0;
	element_[7] = 0.0;
	element_[8] = 1.0;
}

void
dSO3::euler_zyx(double z, double y, double x)
{
	static double ca, sa, cb, sb, cg, sg;

	ca = cos(z);
	sa = sin(z);
	cb = cos(y);
	sb = sin(y);
	cg = cos(x);
	sg = sin(x);

	element_[0] = ca * cb;
	element_[1] = sa * cb;
	element_[2] = -sb;
	element_[3] = ca * sb * sg - sa * cg;
	element_[4] = sa * sb * sg + ca * cg;
	element_[5] = cb * sg;
	element_[6] = ca * sb * cg + sa * sg;
	element_[7] = sa * sb * cg - ca * sg;
	element_[8] = cb * cg;
}

void
dSO3::rollpitchyaw(double roll, double pitch, double yaw)
{
	static double sa, ca, sb, cb, sg, cg;

	sa = sin(roll);
	ca = cos(roll);
	sb = sin(pitch);
	cb = cos(pitch);
	sg = sin(yaw);
	cg = cos(yaw);

	element_[0] = cg * cb;
	element_[1] = sg * cb;
	element_[2] = - sb;
	element_[3] = cg * sb * sa - sg * ca;
	element_[4] = sg * sb * sa + cg * ca;
	element_[5] = cb * sa;
	element_[6] = cg * sb * ca + sg * sa;
	element_[7] = sg * sb * ca - cg * sa;
	element_[8] = cb * ca;

}

/* return  -PI/2 < beta <= PI/2, INV_EULER_EPS*/
int
dSO3::inv_euler_zyx(double &z, double &y, double &x)
{
	double cb;

	cb = std::sqrt(element_[0] * element_[0] + element_[1] * element_[1]);

	if (fabs(cb) < rmcp::INV_EULER_EPS) {
		z = 0.0;
		y = rmcp::PI_2;
		x = std::atan2(element_[3], element_[4]);  //atan2(r12, r22)
		return -1;
	}
	else {
		z = std::atan2(element_[1], element_[0]);
		y = std::atan2(-element_[2], cb);
		x = std::atan2(element_[5], element_[8]);
		return 0;
	}
}

void 
dSO3::inv_rollpitchyaw(double &roll, double &pitch, double &yaw)
{
	roll = std::atan2(element_[5], element_[8]); // roll
	pitch = std::asin(-element_[2]);	// pitch
	yaw = std::atan2(element_[1], element_[0]); // yaw
}



void
dSO3::euler_zyz(double z1, double y, double z2)
{
	double ca, sa, cb, sb, cg, sg;

	ca = cos(z1);
	sa = sin(z1);
	cb = cos(y);
	sb = sin(y);
	cg = cos(z2);
	sg = sin(z2);

	element_[0] = ca * cb * cg - sa * sg;
	element_[1] = sa * cb * cg + ca * sg;
	element_[2] = -sb * cg;
	element_[3] = -ca * cb * sg - sa * cg;
	element_[4] = ca * cg - sa * cb * sg;
	element_[5] = sb * sg;
	element_[6] = ca * sb;
	element_[7] = sa * sb;
	element_[8] = cb;

}

/* return 0 <= beta < PI, INV_EULER_EPS */
int
dSO3::inv_euler_zyz(double &z1, double &y, double &z2)
{
	double sb;

	sb = std::sqrt(element_[2] * element_[2] + element_[5] * element_[5]);

	if (fabs(sb) < rmcp::INV_EULER_EPS) {
		z1 = 0.0;
		y = 0.0;
		z2 = std::atan2(element_[1], element_[0]);  //atan2(r21, r11)
		return -1;
	}
	else {
		z1 = std::atan2(element_[7], element_[6]);
		y = std::atan2(sb, element_[8]);
		z2 = std::atan2(element_[5], -element_[2]);
		return 0;
	}
}

dSO3::EULER 
dSO3::inv_euler(double &a, double &b, double &c)
{
	/* inv_euler_zyx */
	double cb, sb;
	dSO3::EULER euler = dSO3::ZYZ;

	cb = std::sqrt(element_[0] * element_[0] + element_[1] * element_[1]);

	if (fabs(cb) > rmcp::INV_EULER_EPS) {
		double y = std::atan2(-element_[2], cb);
		if((y < rmcp::PI_4) && (y > -rmcp::PI_4)) {
			euler = dSO3::ZYX;
		}
	}

	switch (euler) {
		case ZYX :
			a = std::atan2(element_[1], element_[0]); // z
			b = std::atan2(-element_[2], cb); // y
			c = std::atan2(element_[5], element_[8]); // x
			break;
		case ZYZ : 
			sb = std::sqrt(element_[2] * element_[2] + element_[5] * element_[5]);
			a = std::atan2(element_[7], element_[6]); // z1
			b = std::atan2(sb, element_[8]); // y
			c = std::atan2(element_[5], -element_[2]); // z2
			break;
	}

	return euler;
}

///////////// friend 
//

dso3
rmcp::operator * (double s, const dso3 &w)
{
	dso3 wr(w);

  	wr.element_[1] *= s;
	wr.element_[2] *= s;
	wr.element_[3] *= s;
	wr.element_[5] *= s;
	wr.element_[6] *= s;
	wr.element_[7] *= s;

   	return wr;
}

dSO3
rmcp::inv(const dSO3 &R)
{
	rmcp::dSO3 mr(R);
	mr.transpose();
	return mr;
}

/* if the norm of w is less than rmcp::EXP_EPS, return identity matrix */
rmcp::dSO3
rmcp::exp(const dso3 &w)
{
	rmcp::dSO3 mr;
	mr.exp(w);
	return mr;
}

/* given normalized w */
rmcp::dSO3
rmcp::exp(const dso3 &w, double theta)
{
	rmcp::dSO3 mr;
	mr.exp(w, theta);
	return mr;
}

dso3
rmcp::log(const dSO3 &R)
{
	rmcp::dso3 mr;
	mr.log(R);
	return mr;
}

dSO3
rmcp::rot_x(double t)
{
	rmcp::dSO3 mr;
	mr.rot_x(t);
	return mr;
}

dSO3
rmcp::rot_y(double t)
{ 
	rmcp::dSO3 mr;
	mr.rot_y(t);
	return mr;
}

dSO3
rmcp::rot_z(double t) 
{ 
	rmcp::dSO3 mr;
	mr.rot_z(t);
	return mr;
}

dSO3
rmcp::euler_zyx(double z, double y, double x)
{
	rmcp::dSO3 mr;
	mr.euler_zyx(z, y, x);
	return mr;
}

/* return  -PI/2 < beta <= PI/2, INV_EULER_EPS*/
int
rmcp::inv_euler_zyx(dSO3 &R, double &z, double &y, double &x)
{
	return R.inv_euler_zyx(z, y, x);
}

dSO3
rmcp::euler_zyz(double z1, double y, double z2)
{
	rmcp::dSO3 mr;
	mr.euler_zyz(z1, y, z2);
	return mr;
}

/* return 0 <= beta < PI, INV_EULER_EPS */
int
rmcp::inv_euler_zyz(dSO3 &R, double &z1, double &y, double &z2)
{
	return R.inv_euler_zyz(z1, y, z2);
}


