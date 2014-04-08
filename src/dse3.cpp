
#include "rmcp/dse3.h"
#include "rmcp/dso3.h"

#include <cmath>
#include <cassert>

using namespace rmcp;

static int i, j;

dse3::dse3(void) : rmcp::dvector(6, 0.0)
{

}

dse3::dse3(double w0, double w1, double w2, double v0, double v1, double v2)
	: rmcp::dvector(6)
{
	element_[0] = w0;
	element_[1] = w1;
	element_[2] = w2;
   	element_[3] = v0;
	element_[4] = v1;
	element_[5] = v2;
}

dse3::dse3(const dse3 &S) : rmcp::dvector(6)
{
	for (i = 0; i < 6; i++) {
		element_[i] = S.element_[i];
	}
}

dse3::dse3(const rmcp::dvector &w, const rmcp::dvector &v) : rmcp::dvector(6)
{
	assert (w.size_ == 3 &&
			v.size_ == 3 &&
			"dse3::dse3(const rmcp::dvector &, const rmcp::dvector &)");

   	element_[0] = w.element_[0];
   	element_[1] = w.element_[1];
   	element_[2] = w.element_[2];
   	element_[3] = v.element_[0];
   	element_[4] = v.element_[1];
   	element_[5] = v.element_[2];
}

dse3::dse3(const double element[]) : rmcp::dvector(6, element)
{

}

const dse3 &
dse3::operator + (void) const
{
   	return *this;
}

dse3
dse3::operator - (void) const	
{
  	return dse3(-element_[0],
		   	-element_[1],
		   	-element_[2],
		   	-element_[3],
		   	-element_[4],
		   	-element_[5]);
}

dse3 &
dse3::operator = (const dse3 &S)	
{
   	element_[0] = S.element_[0];
   	element_[1] = S.element_[1];
   	element_[2] = S.element_[2];
   	element_[3] = S.element_[3];
   	element_[4] = S.element_[4];
   	element_[5] = S.element_[5];

   	return *this;
}

dse3 &
dse3::operator += (const dse3 &S)
{
   	element_[0] += S.element_[0];
   	element_[1] += S.element_[1];
   	element_[2] += S.element_[2];
   	element_[3] += S.element_[3];
   	element_[4] += S.element_[4];
   	element_[5] += S.element_[5];

   	return *this;
}

dse3 &
dse3::operator -= (const dse3 &S)	
{
   	element_[0] -= S.element_[0];
   	element_[1] -= S.element_[1];
   	element_[2] -= S.element_[2];
   	element_[3] -= S.element_[3];
   	element_[4] -= S.element_[4];
   	element_[5] -= S.element_[5];

	return *this;
}

dse3 &
dse3::operator *= (double s)
{
   	element_[0] *= s;
   	element_[1] *= s;
   	element_[2] *= s;
   	element_[3] *= s;
   	element_[4] *= s;
   	element_[5] *= s;
   
	return *this;
}

dse3 &
dse3::operator /= (double s)
{
	assert (s != 0.0 &&
			"dse3 &dse3::operator /= (double)");

   	s = 1.0 / s;
   	element_[0] *= s;
   	element_[1] *= s;
   	element_[2] *= s;
   	element_[3] *= s;
   	element_[4] *= s;
   	element_[5] *= s;
   
	return *this;
}

dse3
dse3::operator + (const dse3 &S) const
{
   	return dse3(element_[0] + S.element_[0],
		   	element_[1] + S.element_[1],
		   	element_[2] + S.element_[2],
		   	element_[3] + S.element_[3],
		   	element_[4] + S.element_[4],
		   	element_[5] + S.element_[5]);
}

dse3
dse3::operator - (const dse3 &S) const
{
   	return dse3(element_[0] - S.element_[0],
		   	element_[1] - S.element_[1],
		   	element_[2] - S.element_[2],
		   	element_[3] - S.element_[3],
		   	element_[4] - S.element_[4],
		   	element_[5] - S.element_[5]);
}

dse3
dse3::operator * (double s) const
{
   	return dse3(s * element_[0],
		   	s * element_[1],
		   	s * element_[2],
		   	s * element_[3],
		   	s * element_[4],
		   	s * element_[5]);
}

dse3
dse3::operator / (double s) const
{
	assert (s != 0.0 &&
		"dse3 dse3::operator / (double) const");

   	s = 1.0 / s;
   	return dse3(s * element_[0],
		   	s * element_[1],
		   	s * element_[2],
		   	s * element_[3],
		   	s * element_[4],
		   	s * element_[5]);
}

double
dse3::operator * (const dse3& s)
{
	return 
		element_[0] * s.element_[0] +
		element_[1] * s.element_[1] +
		element_[2] * s.element_[2] +
        element_[3] * s.element_[3] +
		element_[4] * s.element_[4] +
		element_[5] * s.element_[5];
}

dgqmatrix 
dse3::operator * (dSE3 &T) const
{
	double element[16];

	element[0] = -element_[2] * T.element_[1] + element_[1] * T.element_[2];
	element[1] =  element_[2] * T.element_[0] - element_[0] * T.element_[2];
	element[2] = -element_[1] * T.element_[0] + element_[0] * T.element_[1];
	element[3] = 0.0;
	element[4] = -element_[2] * T.element_[5] + element_[1] * T.element_[6];
	element[5] =  element_[2] * T.element_[4] - element_[0] * T.element_[6];
	element[6] = -element_[1] * T.element_[4] + element_[0] * T.element_[5];
	element[7] = 0.0;
	element[8] = -element_[2] * T.element_[9] + element_[1] * T.element_[10];
	element[9] =  element_[2] * T.element_[8] - element_[0] * T.element_[10];
	element[10] = -element_[1] * T.element_[8] + element_[0] * T.element_[9];
	element[11] = 0.0;
	element[12] = -element_[2] * T.element_[13] + element_[1] * T.element_[14] + element_[3];
	element[13] =  element_[2] * T.element_[12] - element_[0] * T.element_[14] + element_[4];
	element[14] = -element_[1] * T.element_[12] + element_[0] * T.element_[13] + element_[5];
	element[15] = 0.0;

	return dgqmatrix(4, element);
}

void
dse3::set(double w0, double w1, double w2, double v0, double v1, double v2)
{
   	element_[0] = w0;
   	element_[1] = w1;
   	element_[2] = w2;
   	element_[3] = v0;
   	element_[4] = v1;
   	element_[5] = v2;
}

void
dse3::set_angular(dvector &w)
{
	element_[0] = w.element_[0];
	element_[1] = w.element_[1];
	element_[2] = w.element_[2];
}

void
dse3::set_linear(dvector &v)
{
	element_[3] = v.element_[0];
	element_[4] = v.element_[1];
	element_[5] = v.element_[2];
}

dvector
dse3::get_angular(void)
{
	return dvector(3, element_);
}

dvector
dse3::get_linear(void)
{
	return dvector(3, element_+3);
}

void 
dse3::log(const dSE3 &T)
{
	// TODO
	double theta = 0.5 * (T.element_[0] + T.element_[5] + T.element_[10] - 1.0);
	theta = acos(theta);
	double cof = theta / (2.0 * sin(theta));
	double x = cof * (T.element_[6] - T.element_[9]);
	double y = cof * (T.element_[8] - T.element_[2]);
	double z = cof * (T.element_[1] - T.element_[4]);
	theta = sqrt(x * x + y * y + z * z);
	cof = (2.0 * sin(theta) - theta * (1.0 + cos(theta))) / (2.0 * theta * theta * sin(theta));
	/*
	return dse3( x, y, z,
		(1.0 - cof * (y * y + z * z)) * T.element_[12] + (0.5 * z + cof * x * y) * T.element_[13] + (cof * x * z - 0.5 * y) * T.element_[14],
		(1.0 - cof * (x * x + z * z)) * T.element_[13] + (0.5 * x + cof * z * y) * T.element_[14] + (cof * x * y - 0.5 * z) * T.element_[12],
		(1.0 - cof * (y * y + x * x)) * T.element_[14] + (0.5 * y + cof * x * z) * T.element_[12] + (cof * y * z - 0.5 * x) * T.element_[13]);
}
*/
}

void
dse3::Ad(const rmcp::dSE3 &T, const dse3 &S)
{
	element_[0] = 
		T.element_[0] * S.element_[0] +
	   	T.element_[4] * S.element_[1] +
	   	T.element_[8] * S.element_[2];
	element_[1] = 
		T.element_[1] * S.element_[0] +
	   	T.element_[5] * S.element_[1] +
	   	T.element_[9] * S.element_[2];
	element_[2] = 
		T.element_[2] * S.element_[0] +
	   	T.element_[6] * S.element_[1] +
	   	T.element_[10] * S.element_[2];
	element_[3] = 
		T.element_[13] * element_[2] -
	   	T.element_[14] * element_[1] +
	   	T.element_[0] * S.element_[3] +
	   	T.element_[4] * S.element_[4] +
	   	T.element_[8] * S.element_[5];
	element_[4] = 
		T.element_[14] * element_[0] -
	   	T.element_[12] * element_[2] +
	   	T.element_[1] * S.element_[3] +
	   	T.element_[5] * S.element_[4] +
	   	T.element_[9] * S.element_[5];
	element_[5] = 
		T.element_[12] * element_[1] -
	   	T.element_[13] * element_[0] +
	   	T.element_[2] * S.element_[3] +
	   	T.element_[6] * S.element_[4] +
	   	T.element_[10] * S.element_[5];
}

void
dse3::inv_Ad(const dSE3 &T, const dse3 &S)
{
	double tmp[3];

	tmp[0] =
	   	S.element_[3] +
	   	S.element_[1] * T.element_[14] -
	   	S.element_[2] * T.element_[13]; 
	tmp[1] =
	   	S.element_[4] +
	   	S.element_[2] * T.element_[12] -
	   	S.element_[0] * T.element_[14]; 
	tmp[2] =
	   	S.element_[5] +
	   	S.element_[0] * T.element_[13] -
	   	S.element_[1] * T.element_[12];

	element_[0] = 
		T.element_[0] * S.element_[0] +
	   	T.element_[1] * S.element_[1] + 
		T.element_[2] * S.element_[2];
	element_[1] =
   		T.element_[4] * S.element_[0] +
	   	T.element_[5] * S.element_[1] +
	   	T.element_[6] * S.element_[2];
	element_[2] = 
		T.element_[8] * S.element_[0] + 
		T.element_[9] * S.element_[1] + 
		T.element_[10] * S.element_[2];
	element_[3] =
   		T.element_[0] * tmp[0] + 
		T.element_[1] * tmp[1] +
	   	T.element_[2] * tmp[2];
	element_[4] =
   		T.element_[4] * tmp[0] +
	   	T.element_[5] * tmp[1] +
	   	T.element_[6] * tmp[2];
	element_[5] =
   		T.element_[8] * tmp[0] +
	   	T.element_[9] * tmp[1] +
	   	T.element_[10] * tmp[2];
}

void
dse3::ad(const dse3 &s1, const dse3 &s2)
{
	element_[0] =
		s1.element_[1] * s2.element_[2] -
	   	s1.element_[2] * s2.element_[1];
	element_[1] =
		s1.element_[2] * s2.element_[0] -
	   	s1.element_[0] * s2.element_[2];
	element_[2] =
		s1.element_[0] * s2.element_[1] -
	   	s1.element_[1] * s2.element_[0];
	element_[3] =
		s1.element_[1] * s2.element_[5] -
	   	s1.element_[2] * s2.element_[4] -
	   	s2.element_[1] * s1.element_[5] + 
		s2.element_[2] * s1.element_[4];
	element_[4] =
		s1.element_[2] * s2.element_[3] -
	   	s1.element_[0] * s2.element_[5] -
	   	s2.element_[2] * s1.element_[3] +
	   	s2.element_[0] * s1.element_[5];
	element_[5] =
		s1.element_[0] * s2.element_[4] -
	   	s1.element_[1] * s2.element_[3] -
	   	s2.element_[0] * s1.element_[4] +
	   	s2.element_[1] * s1.element_[3];
}

void
dse3::dAd(const rmcp::dSE3 &T, const dse3 &t)
{
	double tmp[3];

	tmp[0] = 
		t.element_[0] - 
		T.element_[13] * t.element_[5] + 
		T.element_[14] * t.element_[4];
	tmp[1] = 
		t.element_[1] - 
		T.element_[14] * t.element_[3] +
	   	T.element_[12] * t.element_[5]; 
	tmp[2] = 
		t.element_[2] - 
		T.element_[12] * t.element_[4] + 
		T.element_[13] * t.element_[3];

	element_[0] = 
		T.element_[0] * tmp[0] + 
		T.element_[1] * tmp[1] + 
		T.element_[2] * tmp[2];
	element_[1] = 
		T.element_[4] * tmp[0] + 
		T.element_[5] * tmp[1] +
	   	T.element_[6] * tmp[2];
	element_[2] = 
		T.element_[8] * tmp[0] +
	   	T.element_[9] * tmp[1] + 
		T.element_[10] * tmp[2];
	element_[3] = 
		T.element_[0] * t.element_[3] +
	   	T.element_[1] * t.element_[4] + 
		T.element_[2] * t.element_[5];
	element_[4] = 
		T.element_[4] * t.element_[3] +
	   	T.element_[5] * t.element_[4] +
	   	T.element_[6] * t.element_[5];
	element_[5] = 
		T.element_[8] * t.element_[3] +
	   	T.element_[9] * t.element_[4] +
	   	T.element_[10] * t.element_[5];
}

void
dse3::inv_dAd(const dSE3& T, const dse3& t)
{
	double tmp[3];
	tmp[0] = T.element_[0] * t.element_[3] + T.element_[4] * t.element_[4] + T.element_[8] * t.element_[5];
	tmp[1] = T.element_[1] * t.element_[3] + T.element_[5] * t.element_[4] + T.element_[9] * t.element_[5];
	tmp[2] = T.element_[2] * t.element_[3] + T.element_[6] * t.element_[4] + T.element_[10] * t.element_[5];

	element_[0] = T.element_[13] * tmp[2] - T.element_[14] * tmp[1] + T.element_[0] * t.element_[0] + T.element_[4] * t.element_[1] + T.element_[8] * t.element_[2];
	element_[1] = T.element_[14] * tmp[0] - T.element_[12] * tmp[2] + T.element_[1] * t.element_[0] + T.element_[5] * t.element_[1] + T.element_[9] * t.element_[2];
	element_[2] = T.element_[12] * tmp[1] - T.element_[13] * tmp[0] + T.element_[2] * t.element_[0] + T.element_[6] * t.element_[1] + T.element_[10] * t.element_[2];
	element_[3] = tmp[0];
	element_[4] = tmp[1];
	element_[5] = tmp[2];
}

void
dse3::dad(const dse3 &s, const dse3 &t)
{
	element_[0] =
		t.element_[1] * s.element_[2] - t.element_[2] * s.element_[1] +
	   	t.element_[4] * s.element_[5] - t.element_[5] * s.element_[4];
	element_[1] =
		t.element_[2] * s.element_[0] - t.element_[0] * s.element_[2] +
	   	t.element_[5] * s.element_[3] - t.element_[3] * s.element_[5];
	element_[2] =	
		t.element_[0] * s.element_[1] - t.element_[1] * s.element_[0] +
	   	t.element_[3] * s.element_[4] - t.element_[4] * s.element_[3];
	element_[3] =
		t.element_[4] * s.element_[2] - t.element_[5] * s.element_[1];
	element_[4] =
		t.element_[5] * s.element_[0] - t.element_[3] * s.element_[2];
	element_[5] =
		t.element_[3] * s.element_[1] - t.element_[4] * s.element_[0];
}

void
dse3::rotate(const dSE3 &T, const dse3 &S)
{
	element_[0] = T.element_[0] * S.element_[0] + T.element_[4] * S.element_[1] + T.element_[8] * S.element_[2];
	element_[1] = T.element_[1] * S.element_[0] + T.element_[5] * S.element_[1] + T.element_[9] * S.element_[2];
	element_[2] = T.element_[2] * S.element_[0] + T.element_[6] * S.element_[1] + T.element_[10] * S.element_[2];
	element_[3] = T.element_[0] * S.element_[3] + T.element_[4] * S.element_[4] + T.element_[8] * S.element_[5];
	element_[4] = T.element_[1] * S.element_[3] + T.element_[5] * S.element_[4] + T.element_[9] * S.element_[5];
	element_[5] = T.element_[2] * S.element_[3] + T.element_[6] * S.element_[4] + T.element_[10] * S.element_[5];

}
void
dse3::inv_rotate(const dSE3 &T, const dse3 &S)
{
	element_[0] = T.element_[0] * S.element_[0] + T.element_[1] * S.element_[1] + T.element_[2] * S.element_[2];
	element_[1] = T.element_[4] * S.element_[0] + T.element_[5] * S.element_[1] + T.element_[6] * S.element_[2];
	element_[2] = T.element_[8] * S.element_[0] + T.element_[9] * S.element_[1] + T.element_[10] * S.element_[2];
	element_[3] = T.element_[0] * S.element_[3] + T.element_[1] * S.element_[4] + T.element_[2] * S.element_[5];
	element_[4] = T.element_[4] * S.element_[3] + T.element_[5] * S.element_[4] + T.element_[6] * S.element_[5];
	element_[5] = T.element_[8] * S.element_[3] + T.element_[9] * S.element_[4] + T.element_[10] * S.element_[5];
}

rmcp::dmatrix 
dse3::export_dmatrix(void)
{
	rmcp::dmatrix M(4, 4, 0.0);

	M.element_[1] = element_[2];
	M.element_[2] = -element_[1];
	M.element_[4] = -element_[2];
	M.element_[6] = element_[0];
	M.element_[8] = element_[1];
	M.element_[9] = -element_[0];
	M.element_[12] = element_[3];
	M.element_[13] = element_[4];
	M.element_[14] = element_[5];

	return M;
}

/*
void 
dse3::multiply_T1sT2(const rmcp::dSE3 &T1, const dse3 &s, const dSE3 &T2)
{
	double tmp[16];

	tmp[ 0] =  T1.element_[ 4] * s.element_[2] - T1.element_[ 8] * s.element_[1];
	tmp[ 1] =  T1.element_[ 5] * s.element_[2] - T1.element_[ 9] * s.element_[1];
	tmp[ 2] =  T1.element_[ 6] * s.element_[2] - T1.element_[10] * s.element_[1];
	tmp[ 3] = 0.0 
	tmp[ 4] = -T1.element_[ 0] * s.element_[2] + T1.element_[ 8] * s.element_[0];
	tmp[ 5] = -T1.element_[ 1] * s.element_[2] + T1.element_[ 9] * s.element_[0];
	tmp[ 6] = -T1.element_[ 2] * s.element_[2] + T1.element_[10] * s.element_[0];
	tmp[ 7] = 0.0;
	tmp[ 8] =  T1.element_[ 0] * s.element_[1] - T1.element_[ 4] * s.element_[0];
	tmp[ 9] =  T1.element_[ 1] * s.element_[1] - T1.element_[ 5] * s.element_[0];
	tmp[10] =  T1.element_[ 2] * s.element_[1] - T1.element_[ 6] * s.element_[0];
	tmp[11] = 0.0;
	tmp[12] =  T1.element_[ 0] * s.element_[3] + T1.element_[ 4] * s.element_[4] + T1.element_[ 8] * s.element_[5];
	tmp[13] =  T1.element_[ 1] * s.element_[3] + T1.element_[ 5] * s.element_[4] + T1.element_[ 9] * s.element_[5];
	tmp[14] =  T1.element_[ 2] * s.element_[3] + T1.element_[ 6] * s.element_[4] + T1.element_[10] * s.element_[5];
	tmp[15] = 0.0;


	element_[ 0] = tmp[ 0] * T2.element_[ 0] + tmp[ 4] * T2.element_[ 1] + tmp[ 8] * T2.element_[ 2];
	element_[ 0] = tmp[ 0] * T2.element_[ 0] + tmp[ 4] * T2.element_[ 1] + tmp[ 8] * T2.element_[ 2];
	element_[ 0] = tmp[ 0] * T2.element_[ 0] + tmp[ 4] * T2.element_[ 1] + tmp[ 8] * T2.element_[ 2];



}
*/

dSE3::dSE3(void) : rmcp::dgqmatrix(4)
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

dSE3::dSE3(const rmcp::dSE3 &T) : rmcp::dgqmatrix(4)
{
   	element_[0] = T.element_[0];
   	element_[1] = T.element_[1];
   	element_[2] = T.element_[2];
   	element_[3] = 0.0;
	element_[4] = T.element_[4];
   	element_[5] = T.element_[5];
   	element_[6] = T.element_[6];
   	element_[7] = 0.0;
	element_[8] = T.element_[8];
   	element_[9] = T.element_[9];
   	element_[10] = T.element_[10];
   	element_[11] = 0.0;
    element_[12] = T.element_[12];
   	element_[13] = T.element_[13];
   	element_[14] = T.element_[14];
	element_[15] = 1.0;
}

dSE3::dSE3(double R0, double R1, double R2,
	   	double R3, double R4, double R5,
	   	double R6, double R7, double R8,
	   	double p0, double p1, double p2) : rmcp::dgqmatrix(4)
{
   	element_[0] = R0;
   	element_[1] = R1;
   	element_[2] = R2;
  	element_[3] = 0.0;
  	element_[4] = R3;
   	element_[5] = R4;
   	element_[6] = R5;
   	element_[7] = 0.0;
  	element_[8] = R6;
   	element_[9] = R7;
   	element_[10] = R8;
 	element_[11] = 0.0;
   	element_[12] = p0;
   	element_[13] = p1;
   	element_[14] = p2;
  	element_[15] = 1.0; 
}

dSE3::dSE3(const rmcp::dSO3 &R, const rmcp::dvector &p) : rmcp::dgqmatrix(4)
{
	assert (p.size_ == 3 &&
			"dSE3::dSE3(const rmcp::dSO3 &, const rmcp::dvector &)");

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
   	element_[12] = p.element_[0];
   	element_[13] = p.element_[1];
   	element_[14] = p.element_[2];
  	element_[15] = 1.0;
}

dSE3::dSE3(const rmcp::dvector &v0, 
		const rmcp::dvector &v1, 
		const rmcp::dvector &v2, 
		const rmcp::dvector &p) : rmcp::dgqmatrix(4)
{
	assert (v0.size_ == 3 &&
			v1.size_ == 3 &&
			v2.size_ == 3 &&
			p.size_ == 3 &&
			"dSE3::dSE3(const rmcp::dvector &, const rmcp::dvector &, const rmcp::dvector &, const rmcp::dvector &)");

	element_[0] = v0.element_[0];
	element_[1] = v0.element_[1];
	element_[2] = v0.element_[2];
	element_[3] = 0.0;
	element_[4] = v1.element_[0];
	element_[5] = v1.element_[1];
	element_[6] = v1.element_[2];
	element_[7] = 0.0;
	element_[8] = v2.element_[0];
	element_[9] = v2.element_[1];
	element_[10] = v2.element_[2];
	element_[11] = 0.0;
	element_[12] = p.element_[0];
	element_[13] = p.element_[1];
	element_[14] = p.element_[2];
	element_[15] = 1.0;	
}

dSE3::dSE3(const rmcp::dSO3 &R) : rmcp::dgqmatrix(4)
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
	element_[12] = 0.0;
	element_[13] = 0.0;
	element_[14] = 0.0;
   	element_[15] = 1.0;
}

dSE3::dSE3(const rmcp::dvector &p) : rmcp::dgqmatrix(4)
{
	assert (p.size_ == 3 &&
			"dSE3::dSE3(const rmcp::dvector &) : rmcp::dmatrix(4, 4)");

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
   	element_[12] = p.element_[0];
   	element_[13] = p.element_[1];
   	element_[14] = p.element_[2];
	element_[15] = 1.0;
}

dSE3::dSE3(const double s) : rmcp::dgqmatrix(4)
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
   	element_[12] = s;
	element_[13] = s;
	element_[14] = s;
	element_[15] = 1.0;
}

dSE3::dSE3(const double element[]) : rmcp::dgqmatrix(4)
{
   	element_[0] = element[0];
   	element_[1] = element[1];
   	element_[2] = element[2];
   	element_[3] = 0.0;
	element_[4] = element[4];
   	element_[5] = element[5];
   	element_[6] = element[6];
	element_[7] = 0.0;
	element_[8] = element[8];
   	element_[9] = element[9];
   	element_[10] = element[10];
	element_[11] = 0.0;
   	element_[12] = element[12];
   	element_[13] = element[13];
   	element_[14] = element[14];
   	element_[15] = 1.0;
}

dSE3&
dSE3::operator = (const dSE3& T)
{
   	element_[0] = T.element_[0];
   	element_[1] = T.element_[1];
   	element_[2] = T.element_[2];
   	element_[4] = T.element_[4];
   	element_[5] = T.element_[5];
   	element_[6] = T.element_[6];
   	element_[8] = T.element_[8];
   	element_[9] = T.element_[9];
   	element_[10] = T.element_[10];
   	element_[12] = T.element_[12];
   	element_[13] = T.element_[13];
   	element_[14] = T.element_[14]; 
	return *this;
}

bool
dSE3::operator == (const dSE3& T) const
{
	return (element_[0] == T.element_[0]) &&
	   	(element_[1] == T.element_[1]) &&
	   	(element_[2] == T.element_[2]) &&
	   	(element_[4] == T.element_[4]) &&
	   	(element_[5] == T.element_[5]) &&
	   	(element_[6] == T.element_[6]) &&
	   	(element_[8] == T.element_[8]) &&
	   	(element_[9] == T.element_[9]) &&
	   	(element_[10] == T.element_[10]) &&
		(element_[12] == T.element_[12]) &&
	   	(element_[13] == T.element_[13]) &&
	   	(element_[14] == T.element_[14]);
}

dSE3&
dSE3::operator *= (const dSE3& T)
{
	double x0, x1, x2;

	element_[12] += 
		element_[0] * T.element_[12] + 
		element_[4] * T.element_[13] + 
		element_[8] * T.element_[14];
	element_[13] +=
	   	element_[1] * T.element_[12] + 
		element_[5] * T.element_[13] + 
		element_[9] * T.element_[14];
	element_[14] += 
		element_[2] * T.element_[12] + 
		element_[6] * T.element_[13] + 
		element_[10] * T.element_[14];

	x0 = 
		element_[0] * T.element_[0] + 
		element_[4] * T.element_[1] + 
		element_[8] * T.element_[2];
	x1 =
	   	element_[0] * T.element_[4] + 
		element_[4] * T.element_[5] + 
		element_[8] * T.element_[6];
	x2 =
	   	element_[0] * T.element_[8] + 
		element_[4] * T.element_[9] + 
		element_[8] * T.element_[10];

	element_[0] = x0;
   	element_[4] = x1;
   	element_[8] = x2;

	x0 = 
		element_[1] * T.element_[0] + 
		element_[5] * T.element_[1] + 
		element_[9] * T.element_[2];
	x1 = 
		element_[1] * T.element_[4] +
	   	element_[5] * T.element_[5] + 
		element_[9] * T.element_[6];
	x2 = 
		element_[1] * T.element_[8] + 
		element_[5] * T.element_[9] + 
		element_[9] * T.element_[10];

	element_[1] = x0;
   	element_[5] = x1;
   	element_[9] = x2;

	x0 = 
		element_[2] * T.element_[0] +
	   	element_[6] * T.element_[1] + 
		element_[10] * T.element_[2];
	x1 = 
		element_[2] * T.element_[4] + 
		element_[6] * T.element_[5] +
	   	element_[10] * T.element_[6];
	x2 = 
		element_[2] * T.element_[8] +
	   	element_[6] * T.element_[9] + 
		element_[10] * T.element_[10];

	element_[2] = x0;
   	element_[6] = x1;
   	element_[10] = x2;

	return  *this;
}

dSE3 &
dSE3::operator /= (const dSE3 &T)
{
	double tmp[9];

	tmp[0] = 
		element_[0] * T.element_[0] + 
		element_[4] * T.element_[4] + 
		element_[8] * T.element_[8];
	tmp[1] =
	   	element_[1] * T.element_[0] + 
		element_[5] * T.element_[4] + 
		element_[9] * T.element_[8];
	tmp[2] = 
		element_[2] * T.element_[0] + 
		element_[6] * T.element_[4] + 
		element_[10] * T.element_[8];
	tmp[3] = element_[0] * T.element_[1] + 
		element_[4] * T.element_[5] + element_[8] * T.element_[9];
	tmp[4] = element_[1] * T.element_[1] + 
		element_[5] * T.element_[5] + element_[9] * T.element_[9];
	tmp[5] = element_[2] * T.element_[1] + 
		element_[6] * T.element_[5] + element_[10] * T.element_[9];
	tmp[6] = element_[0] * T.element_[2] + 
		element_[4] * T.element_[6] + element_[8] * T.element_[10];
	tmp[7] = element_[1] * T.element_[2] + 
		element_[5] * T.element_[6] + element_[9] * T.element_[10];
	tmp[8] = element_[2] * T.element_[2] + 
		element_[6] * T.element_[6] + element_[10] * T.element_[10];

	element_[0] = tmp[0];
   	element_[1] = tmp[1];
   	element_[2] = tmp[2];
	element_[4] = tmp[3];
   	element_[5] = tmp[4];
   	element_[6] = tmp[5];
	element_[8] = tmp[6];
   	element_[9] = tmp[7];
   	element_[10] = tmp[8]; 
	element_[12] -= tmp[0] * T.element_[12] + 
		tmp[3] * T.element_[13] + tmp[6] * T.element_[14];
	element_[13] -= tmp[1] * T.element_[12] + 
		tmp[4] * T.element_[13] + tmp[7] * T.element_[14];
	element_[14] -= tmp[2] * T.element_[12] + 
		tmp[5] * T.element_[13] + tmp[8] * T.element_[14];

	return *this;
}

dSE3&
dSE3::operator %= (const dSE3& T)
{
	double tmp[12];

	tmp[0] = element_[0];
	tmp[1] = element_[1];
	tmp[2] = element_[2];
	tmp[3] = element_[4];
	tmp[4] = element_[5];
	tmp[5] = element_[6];
	tmp[6] = element_[8];
	tmp[7] = element_[9];
	tmp[8] = element_[10];
	tmp[9] = T.element_[12] - element_[12];
	tmp[10] = T.element_[13] - element_[13]; 
	tmp[11] = T.element_[14] - element_[14];

	element_[0] = 
		tmp[0] * T.element_[0] + 
		tmp[1] * T.element_[1] + 
		tmp[2] * T.element_[2];
	element_[1] = 
		tmp[3] * T.element_[0] + 
		tmp[4] * T.element_[1] + tmp[5] * T.element_[2];
	element_[2] =
	   	tmp[6] * T.element_[0] +
	   	tmp[7] * T.element_[1] + tmp[8] * T.element_[2];
	element_[4] = 
		tmp[0] * T.element_[4] + 
		tmp[1] * T.element_[5] + tmp[2] * T.element_[6];
	element_[5] =
	   	tmp[3] * T.element_[4] + 
		tmp[4] * T.element_[5] + tmp[5] * T.element_[6];
	element_[6] =
	   	tmp[6] * T.element_[4] + 
		tmp[7] * T.element_[5] + tmp[8] * T.element_[6];
	element_[8] = 
		tmp[0] * T.element_[8] + 
		tmp[1] * T.element_[9] + tmp[2] * T.element_[10];
	element_[9] = 
		tmp[3] * T.element_[8] + 
		tmp[4] * T.element_[9] + tmp[5] * T.element_[10];
	element_[10] = 
		tmp[6] * T.element_[8] + 
		tmp[7] * T.element_[9] + tmp[8] * T.element_[10];
	element_[12] = tmp[0] * tmp[9] + tmp[1] * tmp[10] + tmp[2] * tmp[11];
	element_[13] = tmp[3] * tmp[9] + tmp[4] * tmp[10] + tmp[5] * tmp[11];
	element_[14] = tmp[6] * tmp[9] + tmp[7] * tmp[10] + tmp[8] * tmp[11];

	return *this;
}

dSE3
dSE3::operator * (const dSE3 &T) const
{
	return dSE3(
		element_[0] * T.element_[0] + element_[4] * T.element_[1] + element_[8] * T.element_[2],
		element_[1] * T.element_[0] + element_[5] * T.element_[1] + element_[9] * T.element_[2],
		element_[2] * T.element_[0] + element_[6] * T.element_[1] + element_[10] * T.element_[2],
		element_[0] * T.element_[4] + element_[4] * T.element_[5] + element_[8] * T.element_[6],
		element_[1] * T.element_[4] + element_[5] * T.element_[5] + element_[9] * T.element_[6],
		element_[2] * T.element_[4] + element_[6] * T.element_[5] + element_[10] * T.element_[6],
		element_[0] * T.element_[8] + element_[4] * T.element_[9] + element_[8] * T.element_[10],
		element_[1] * T.element_[8] + element_[5] * T.element_[9] + element_[9] * T.element_[10],
		element_[2] * T.element_[8] + element_[6] * T.element_[9] + element_[10] * T.element_[10],
		element_[12] + element_[0] * T.element_[12] + element_[4] * T.element_[13] + element_[8] * T.element_[14],
		element_[13] + element_[1] * T.element_[12] + element_[5] * T.element_[13] + element_[9] * T.element_[14],
		element_[14] + element_[2] * T.element_[12] + element_[6] * T.element_[13] + element_[10] * T.element_[14]
	   	);
}

rmcp::dvector
dSE3::operator * (const rmcp::dvector &p) const
{
	double tmp[3];

	tmp[0] = element_[12] + element_[0] * p.element_[0] + element_[4] * p.element_[1] + element_[8] * p.element_[2];
	tmp[1] = element_[13] + element_[1] * p.element_[0] + element_[5] * p.element_[1] + element_[9] * p.element_[2];
	tmp[2] = element_[14] + element_[2] * p.element_[0] + element_[6] * p.element_[1] + element_[10] * p.element_[2];

	return rmcp::dvector(3, tmp);
}

dSE3
dSE3::operator / (const dSE3& T) const
{
	double tmp[9] = { element_[0] * T.element_[0] + element_[4] * T.element_[4] + element_[8] * T.element_[8],
		element_[1] * T.element_[0] + element_[5] * T.element_[4] + element_[9] * T.element_[8],
		element_[2] * T.element_[0] + element_[6] * T.element_[4] + element_[10] * T.element_[8],
		element_[0] * T.element_[1] + element_[4] * T.element_[5] + element_[8] * T.element_[9],
		element_[1] * T.element_[1] + element_[5] * T.element_[5] + element_[9] * T.element_[9],
		element_[2] * T.element_[1] + element_[6] * T.element_[5] + element_[10] * T.element_[9],
		element_[0] * T.element_[2] + element_[4] * T.element_[6] + element_[8] * T.element_[10],
		element_[1] * T.element_[2] + element_[5] * T.element_[6] + element_[9] * T.element_[10],
		element_[2] * T.element_[2] + element_[6] * T.element_[6] + element_[10] * T.element_[10] };

	return dSE3(	tmp[0], tmp[1], tmp[2], tmp[3], tmp[4],
		   	tmp[5], tmp[6], tmp[7], tmp[8], 
		element_[12] - tmp[0] * T.element_[12] - tmp[3] * T.element_[13] - tmp[6] * T.element_[14],
		element_[13] - tmp[1] * T.element_[12] - tmp[4] * T.element_[13] - tmp[7] * T.element_[14],
		element_[14] - tmp[2] * T.element_[12] - tmp[5] * T.element_[13] - tmp[8] * T.element_[14] );
}

dSE3
dSE3::operator % (const dSE3& T) const
{
	return dSE3(	element_[0] * T.element_[0] + element_[1] * T.element_[1] + element_[2] * T.element_[2],
		element_[4] * T.element_[0] + element_[5] * T.element_[1] + element_[6] * T.element_[2],
		element_[8] * T.element_[0] + element_[9] * T.element_[1] + element_[10] * T.element_[2],
		element_[0] * T.element_[4] + element_[1] * T.element_[5] + element_[2] * T.element_[6],
		element_[4] * T.element_[4] + element_[5] * T.element_[5] + element_[6] * T.element_[6],
		element_[8] * T.element_[4] + element_[9] * T.element_[5] + element_[10] * T.element_[6],
		element_[0] * T.element_[8] + element_[1] * T.element_[9] + element_[2] * T.element_[10],
		element_[4] * T.element_[8] + element_[5] * T.element_[9] + element_[6] * T.element_[10],
		element_[8] * T.element_[8] + element_[9] * T.element_[9] + element_[10] * T.element_[10],
		element_[0] * (T.element_[12] - element_[12]) + element_[1] * (T.element_[13] - element_[13]) +
	   		element_[2] * (T.element_[14] - element_[14]),
		element_[4] * (T.element_[12] - element_[12]) + element_[5] * (T.element_[13] - element_[13]) +
	   		element_[6] * (T.element_[14] - element_[14]),
		element_[8] * (T.element_[12] - element_[12]) + element_[9] * (T.element_[13] - element_[13]) +
	   		element_[10] * (T.element_[14] - element_[14]) );
}

rmcp::dvector
dSE3::operator % (const rmcp::dvector& p) const
{
	assert (p.size_ == 3 &&
		"rmcp::dvector dSE3::operator % (const rmcp::dvector &) const");

	double tmp[] = { element_[0] * (p.element_[0] - element_[12]) + element_[1] * (p.element_[1] - element_[13]) +
		   			element_[2] * (p.element_[2] - element_[14]),
				 element_[4] * (p.element_[0] - element_[12]) + element_[5] * (p.element_[1] - element_[13]) +
					element_[6] * (p.element_[2] - element_[14]),
				 element_[8] * (p.element_[0] - element_[12]) + element_[9] * (p.element_[1] - element_[13]) +
				 element_[10] * (p.element_[2] - element_[14]) };

	return rmcp::dvector(3, tmp);
}

dgqmatrix 
dSE3::operator * (dse3 &A) const
{
	double element[16];

	element[0] = element_[4] * A.element_[2] - element_[8] * A.element_[1];
	element[1] = element_[5] * A.element_[2] - element_[9] * A.element_[1];
	element[2] = element_[6] * A.element_[2] - element_[10] * A.element_[1];
	element[3] = 0.0;
	element[4] = -element_[0] * A.element_[2] + element_[8] * A.element_[0];
	element[5] = -element_[1] * A.element_[2] + element_[9] * A.element_[0];
	element[6] = -element_[2] * A.element_[2] + element_[10] * A.element_[0];
	element[7] = 0.0;
	element[8] = element_[0] * A.element_[1] - element_[4] * A.element_[0];
	element[9] = element_[1] * A.element_[1] - element_[5] * A.element_[0];
	element[10] = element_[2] * A.element_[1] - element_[6] * A.element_[0];
	element[11] = 0.0;
	element[12] = element_[0] * A.element_[3] + element_[4] * A.element_[4] + element_[8] * A.element_[5];
	element[13] = element_[1] * A.element_[3] + element_[5] * A.element_[4] + element_[9] * A.element_[5];
	element[14] = element_[2] * A.element_[3] + element_[6] * A.element_[4] + element_[10] * A.element_[5];
	element[15] = 0.0;

	return dgqmatrix(4, element);
}

void
dSE3::set_identity(void)	
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

dSE3 &
dSE3::set_rotation(const dSO3 &R)
{
	element_[0] = R.element_[0];
   	element_[4] = R.element_[3];
   	element_[8] = R.element_[6];
	element_[1] = R.element_[1];
   	element_[5] = R.element_[4];
   	element_[9] = R.element_[7];
	element_[2] = R.element_[2];
   	element_[6] = R.element_[5];
   	element_[10] = R.element_[8];
	return  *this;
}

dSE3 &
dSE3::set_position(const rmcp::dvector &p)
{
	assert (p.size_ == 3 &&
			"dSE3 &dSE3::set_position(const rmcp::dvector &)");

	element_[12] = p.element_[0];
   	element_[13] = p.element_[1];
   	element_[14] = p.element_[2];

	return  *this;
}

dSE3 &
dSE3::translate(const rmcp::dvector &p)
{
	assert (p.size_ == 3 &&
			"dSE3 &dSE3::translate(const rmcp::dvector &)");

	element_[12] += p.element_[0];
   	element_[13] += p.element_[1];
   	element_[14] += p.element_[2];

	return  *this;
}

dSE3 &
dSE3::rotate(const dSO3 &R)
{
	double r1, r2, r3;

	r1 = element_[0] * R.element_[0] + 
		element_[4] * R.element_[1] + element_[8] * R.element_[2];
	r2 = element_[0] * R.element_[3] + 
		element_[4] * R.element_[4] + element_[8] * R.element_[5];
	r3 = element_[0] * R.element_[6] + 
		element_[4] * R.element_[7] + element_[8] * R.element_[8];

	element_[0] = r1;
   	element_[4] = r2;
   	element_[8] = r3;

	r1 = element_[1] * R.element_[0] + 
		element_[5] * R.element_[1] + element_[9] * R.element_[2];
	r2 = element_[1] * R.element_[3] + 
		element_[5] * R.element_[4] + element_[9] * R.element_[5];
	r3 = element_[1] * R.element_[6] + 
		element_[5] * R.element_[7] + element_[9] * R.element_[8];

	element_[1] = r1;
   	element_[5] = r2;
   	element_[9] = r3;

	r1 = element_[2] * R.element_[0] +
	   	element_[6] * R.element_[1] + element_[10] * R.element_[2];
	r2 = element_[2] * R.element_[3] +
	   	element_[6] * R.element_[4] + element_[10] * R.element_[5];
	r3 = element_[2] * R.element_[6] +
	   	element_[6] * R.element_[7] + element_[10] * R.element_[8];

	element_[2] = r1;
   	element_[6] = r2;
   	element_[10] = r3;

	return  *this;
}

rmcp::dvector
dSE3::position(void) const
{
	rmcp::dvector re(3);
	re.element_[0] = element_[12];
	re.element_[1] = element_[13];
	re.element_[2] = element_[14];

   	return re;
}

dSO3
dSE3::rotation(void) const
{
   	return dSO3(element_[0], element_[1], element_[2], element_[4], element_[5], element_[6], element_[8], element_[9], element_[10]); 
}

/* if the norm of w of S is less than rmcp::EXP_EPS, ... */
void
dSE3::exp(const dse3& S)
{
	double theta;
	theta = cblas_dnrm2(3, S.element_, 1);	// norm of w of S

	if (std::fabs(theta) < rmcp::EXP_EPS) {
		element_[0] = 1.0;
		element_[1] = 0.0;
		element_[2] = 0.0;
		element_[4] = 0.0;
		element_[5] = 1.0;
		element_[6] = 0.0;
		element_[8] = 0.0;
		element_[9] = 0.0;
		element_[10] = 1.0;
		element_[12] = S.element_[3];
		element_[13] = S.element_[4];
		element_[14] = S.element_[5];

		return;
	}
	else {
		static double s[6];
		static double st, ct, vt, ut, t0, t1, t2;

		for (i = 0; i < 6; i++) {
			s[i] = S.element_[i] / theta;
		}

		st = sin(theta);
		ct = cos(theta);
		vt = 1.0 - ct;
		ut = (theta - st) * (s[0] * s[3] + s[1] * s[4] + s[2] * s[5]);
		t0 = s[2] * st;
		t1 = s[1] * st;
		t2 = s[0] * st;

		element_[0] = s[0] * s[0] * vt + ct;
		element_[1] = s[0] * s[1] * vt + t0;
		element_[2] = s[0] * s[2] * vt - t1;
		element_[4] = s[0] * s[1] * vt - t0;
		element_[5] = s[1] * s[1] * vt + ct;
		element_[6] = s[1] * s[2] * vt + t2;
		element_[8] = s[0] * s[2] * vt + t1;
		element_[9] = s[1] * s[2] * vt - t2;
		element_[10] = s[2] * s[2] * vt + ct;
		element_[12] = st * s[3] + ut * s[0] + vt * (s[1] * s[5] - s[2] * s[4]);
		element_[13] = st * s[4] + ut * s[1] + vt * (s[2] * s[3] - s[0] * s[5]);
		element_[14] = st * s[5] + ut * s[2] + vt * (s[0] * s[4] - s[1] * s[3]);

		return;
	}
}

/* given w of S is normalized */
void
dSE3::exp(const dse3& S, double theta)
{
	static double st, ct, vt, ut, t0, t1, t2;

	st = sin(theta);
	ct = cos(theta);
	vt = 1.0 - ct;
	ut = (theta - st) * (S.element_[0] * S.element_[3] + S.element_[1] * S.element_[4] + S.element_[2] * S.element_[5]);
	t0 = S.element_[2] * st;
	t1 = S.element_[1] * st;
	t2 = S.element_[0] * st;

	element_[0] = S.element_[0] * S.element_[0] * vt + ct;
	element_[1] = S.element_[0] * S.element_[1] * vt + t0;
	element_[2] = S.element_[0] * S.element_[2] * vt - t1;
	element_[4] = S.element_[0] * S.element_[1] * vt - t0;
	element_[5] = S.element_[1] * S.element_[1] * vt + ct;
	element_[6] = S.element_[1] * S.element_[2] * vt + t2;
	element_[8] = S.element_[0] * S.element_[2] * vt + t1;
	element_[9] = S.element_[1] * S.element_[2] * vt - t2;
	element_[10] = S.element_[2] * S.element_[2] * vt + ct;
	element_[12] = st * S.element_[3] + ut * S.element_[0] + vt * (S.element_[1] * S.element_[5] - S.element_[2] * S.element_[4]);
	element_[13] = st * S.element_[4] + ut * S.element_[1] + vt * (S.element_[2] * S.element_[3] - S.element_[0] * S.element_[5]);
	element_[14] = st * S.element_[5] + ut * S.element_[2] + vt * (S.element_[0] * S.element_[4] - S.element_[1] * S.element_[3]);

	return;
}

void
dSE3::exp2(const dse3& S, double theta)
{
	static double st, ct, vt, ut, t0, t1, t2;

	st = sin(theta);
	ct = cos(theta);
	vt = 1.0 - ct;
	ut = (theta - st) * (S.element_[0] * S.element_[3] + S.element_[1] * S.element_[4] + S.element_[2] * S.element_[5]);
	t0 = S.element_[2] * st;
	t1 = S.element_[1] * st;
	t2 = S.element_[0] * st;

	element_[0] = 1 - (S.element_[1] * S.element_[1] + S.element_[2] * S.element_[2]) * vt;
	element_[1] = S.element_[0] * S.element_[1] * vt + t0;
	element_[2] = S.element_[0] * S.element_[2] * vt - t1;
	element_[4] = S.element_[0] * S.element_[1] * vt - t0;
	element_[5] = 1 - (S.element_[0] * S.element_[0] + S.element_[2] * S.element_[2]) * vt;
	element_[6] = S.element_[1] * S.element_[2] * vt + t2;
	element_[8] = S.element_[0] * S.element_[2] * vt + t1;
	element_[9] = S.element_[1] * S.element_[2] * vt - t2;
	element_[10] = S.element_[2] * S.element_[2] * vt + ct;
	element_[12] = st * S.element_[3] + ut * S.element_[0] + vt * (S.element_[1] * S.element_[5] - S.element_[2] * S.element_[4]);
	element_[13] = st * S.element_[4] + ut * S.element_[1] + vt * (S.element_[2] * S.element_[3] - S.element_[0] * S.element_[5]);
	element_[14] = st * S.element_[5] + ut * S.element_[2] + vt * (S.element_[0] * S.element_[4] - S.element_[1] * S.element_[3]);

	return;
}

void
dSE3::rot_x(double t)
{
	static double c, s;
	c = cos(t);
	s = sin(t);
	element_[0] = 1.0;
	element_[1] = 0.0;
	element_[2] = 0.0;
	element_[4] = 0.0;
	element_[5] = c;
	element_[6] = s;
	element_[8] = 0.0;
	element_[9] = -s;
	element_[10] = c;
	element_[12] = element_[13] = element_[14] = 0.0;
}

void 
dSE3::rot_y(double t)
{ 
	static double c, s;
	c = cos(t);
	s = sin(t);
	element_[0] = c;
	element_[1] = 0.0;
	element_[2] = -s;
	element_[4] = 0.0;
	element_[5] = 1.0;
	element_[6] = 0.0;
	element_[8] = s;
	element_[9] = 0.0;
	element_[10] = c;
	element_[12] = element_[13] = element_[14] = 0.0;
}

void
dSE3::rot_z(double t)
{
	static double c, s;
   	c = cos(t);
	s = sin(t);
   	element_[0] = c; 
	element_[1] = s;
   	element_[2] = 0.0;
   	element_[4] = -s;
   	element_[5] = c;
	element_[6] = 0.0;
	element_[8] = 0.0;
	element_[9] = 0.0;
	element_[10] = 1.0;
	element_[12] = element_[13] = element_[14] = 0.0;
}
	
void
dSE3::inv(void)
{
	static double element[16];

	cblas_dcopy(size_, element_, 1, element, 1);

	element_[1] = element[4];
	element_[2] = element[8];
	element_[4] = element[1];
	element_[6] = element[9];
	element_[8] = element[2];
	element_[9] = element[6];
	element_[12] = - element[0] * element[12] - element[1] * element[13] - element[2] * element[14];
	element_[13] = - element[4] * element[12] - element[5] * element[13] - element[6] * element[14];
	element_[14] = - element[8] * element[12] - element[9] * element[13] - element[10] * element[14];
}

rmcp::dmatrix 
dSE3::Ad(void)
{
	rmcp::dmatrix m(6, 6, 0.0);

	m.element_[0]  = element_[0];
	m.element_[1]  = element_[1];
	m.element_[2]  = element_[2];
	m.element_[6]  = element_[4];
	m.element_[7]  = element_[5];
	m.element_[8]  = element_[6];
	m.element_[12] = element_[8];
	m.element_[13] = element_[9];
	m.element_[14] = element_[10];

	m.element_[21] = element_[0];
	m.element_[22] = element_[1];
	m.element_[23] = element_[2];
	m.element_[27] = element_[4];
	m.element_[28] = element_[5];
	m.element_[29] = element_[6];
	m.element_[33] = element_[8];
	m.element_[34] = element_[9];
	m.element_[35] = element_[10];

	m.element_[3]  = -element_[14] * element_[1] + element_[13] * element_[2];
	m.element_[4]  =  element_[14] * element_[0] - element_[12] * element_[2];
	m.element_[5]  = -element_[13] * element_[0] + element_[12] * element_[1];
	m.element_[9]  = -element_[14] * element_[5] + element_[13] * element_[6];
	m.element_[10] =  element_[14] * element_[4] - element_[12] * element_[6];
	m.element_[11] = -element_[13] * element_[4] + element_[12] * element_[5];
	m.element_[15] = -element_[14] * element_[9] + element_[13] * element_[10];
	m.element_[16] =  element_[14] * element_[8] - element_[12] * element_[10];
	m.element_[17] = -element_[13] * element_[8] + element_[12] * element_[9];

	return m;
}


rmcp::dmatrix 
dSE3::inv_Ad(void)
{
	rmcp::dmatrix m(6, 6, 0.0);

	m.element_[0]  = element_[0];
	m.element_[1]  = element_[4];
	m.element_[2]  = element_[8];
	m.element_[6]  = element_[1];
	m.element_[7]  = element_[5];
	m.element_[8]  = element_[9];
	m.element_[12] = element_[2];
	m.element_[13] = element_[6];
	m.element_[14] = element_[10];

	m.element_[21] = element_[0];
	m.element_[22] = element_[4];
	m.element_[23] = element_[8];
	m.element_[27] = element_[1];
	m.element_[28] = element_[5];
	m.element_[29] = element_[9];
	m.element_[33] = element_[2];
	m.element_[34] = element_[6];
	m.element_[35] = element_[10];

	m.element_[3]  = -element_[14] * element_[1] + element_[13] * element_[2];
	m.element_[4]  = -element_[14] * element_[5] + element_[13] * element_[6];
	m.element_[5]  = -element_[14] * element_[9] + element_[13] * element_[10];
	m.element_[9]  =  element_[14] * element_[0] - element_[12] * element_[2];
	m.element_[10] =  element_[14] * element_[4] - element_[12] * element_[6];
	m.element_[11] =  element_[14] * element_[8] - element_[12] * element_[10];
	m.element_[15] = -element_[13] * element_[0] + element_[12] * element_[1];
	m.element_[16] = -element_[13] * element_[4] + element_[12] * element_[5];
	m.element_[17] = -element_[13] * element_[8] + element_[12] * element_[9];

	return m;
}



// friend

dse3
rmcp::operator * (double d, const dse3& s)
{
   	return dse3( d * s.element_[0],
		   	d * s.element_[1],
		   	d * s.element_[2],
		   	d * s.element_[3],
		   	d * s.element_[4],
		   	d * s.element_[5] );
}

/*
double
rmcp::square_sum(const dse3 &s)
{
   	return s.element_[0] * s.element_[0] + s.element_[1] * s.element_[1] + s.element_[2] * s.element_[2] +
	   	s.element_[3] * s.element_[3] + s.element_[4] * s.element_[4] + s.element_[5] * s.element_[5];
}
*/

rmcp::dSE3
rmcp::exp(const dse3& S)
{
	rmcp::dSE3 mr;
	mr.exp(S);
	return mr;

}

/* given normalized w of S */
rmcp::dSE3
rmcp::exp(const dse3& S, double theta)
{
	rmcp::dSE3 mr;
	mr.exp(S, theta);
	return mr;
}

dse3
rmcp::rotate(const dSE3 &T, const dse3 &S)
{
	rmcp::dse3 vr;
	vr.rotate(T, S);
	return vr;
}
dse3
rmcp::inv_rotate(const dSE3 &T, const dse3 &S)
{
	rmcp::dse3 vr;
	vr.inv_rotate(T, S);
	return vr;
}

/*
rmcp::dvector
rmcp::rotate(const dSE3 &T, const rmcp::dvector &v)
{
	rmcp::dvector vr(3);
	vr.element_[0] = 
		T.element_[0] * v.element_[0] +
		T.element_[4] * v.element_[1] +
		T.element_[8] * v.element_[2];
	vr.element_[1] = 
		T.element_[1] * v.element_[0] +
		T.element_[5] * v.element_[1] +
		T.element_[9] * v.element_[2];
	vr.element_[2] = 
		T.element_[2] * v.element_[0] +
		T.element_[6] * v.element_[1] +
		T.element_[10] * v.element_[2];
	return vr;
}

rmcp::dvector
rmcp::inv_rotate(const dSE3 &T, const rmcp::dvector &v)
{
	rmcp::dvector vr(3);
	vr.element_[0] = T.element_[0] * v.element_[0] + T.element_[1] * v.element_[1] + T.element_[2] * v.element_[2];
	vr.element_[1] = T.element_[4] * v.element_[0] + T.element_[5] * v.element_[1] + T.element_[6] * v.element_[2];
	vr.element_[2] = T.element_[8] * v.element_[0] + T.element_[9] * v.element_[1] + T.element_[10] * v.element_[2];
	return vr;
}
*/

dse3
rmcp::log(const dSE3 &T)
{
	rmcp::dse3 vr;
	vr.log(T);
	return vr;
}

dse3
rmcp::Ad(const dSE3& T, const dse3& s)
{
	rmcp::dse3 vr;
	vr.Ad(T, s);
	return vr;
}

dse3
rmcp::inv_Ad(const dSE3& T, const dse3& s)
{
	rmcp::dse3 vr;
	vr.inv_Ad(T, s);
	return vr;
}

/*
rmcp::dvector
rmcp::inv_Ad(const dSE3 &T, const rmcp::dvector &v)
{
	rmcp::dvector vr(3);

	vr.element_[0] = 
		T.element_[0] * v.element_[0] +
		T.element_[1] * v.element_[1] +
		T.element_[2] * v.element_[2];
	vr.element_[1] = 
		T.element_[4] * v.element_[0] +
		T.element_[5] * v.element_[1] +
		T.element_[6] * v.element_[2];
	vr.element_[2] = 
		T.element_[8] * v.element_[0] +
		T.element_[9] * v.element_[1] +
		T.element_[10] * v.element_[2];

	return vr;
}

dse3
rmcp::Ad(const rmcp::dvector& p, const dse3& s)
{
	return dse3(
				s.element_[0],
		   		s.element_[1],
			   	s.element_[2], 
				p.element_[1] * s.element_[2] - p.element_[2] * s.element_[1] + s.element_[3], 
				p.element_[2] * s.element_[0] - p.element_[0] * s.element_[2] + s.element_[4], 
				p.element_[0] * s.element_[1] - p.element_[1] * s.element_[0] + s.element_[5]	);
}
*/
/*
rmcp::dvector
rmcp::linear_Ad(const rmcp::dvector& p, const dse3& s)
{
	rmcp::dvector vr(3);
	vr.element_[0] = p.element_[1] * s.element_[2] - p.element_[2] * s.element_[1] + s.element_[3];
	vr.element_[1] = p.element_[2] * s.element_[0] - p.element_[0] * s.element_[2] + s.element_[4];
	vr.element_[2] = p.element_[0] * s.element_[1] - p.element_[1] * s.element_[0] + s.element_[5];

	return vr;
}
*/
dse3
rmcp::ad(const dse3& s1, const dse3& s2)
{
	rmcp::dse3 vr;
	vr.ad(s1, s2);
	return vr;
}
/*
rmcp::dvector
rmcp::ad(const rmcp::dvector& s1, const dse3& s2)
{
	rmcp::dvector vr(3);
	vr.element_[0] = s2.element_[2] * s1.element_[1] - s2.element_[1] * s1.element_[2];
	vr.element_[1] = s2.element_[0] * s1.element_[2] - s2.element_[2] * s1.element_[0];
	vr.element_[2] = s2.element_[1] * s1.element_[0] - s2.element_[0] * s1.element_[1];
	return vr;
}
*/
dse3
rmcp::dAd(const dSE3& T, const dse3& t)
{
	rmcp::dse3 vr;
	vr.dAd(T, t);
	return vr;
}

/*
dse3
rmcp::dAd(const rmcp::dvector& p, const dse3& t)
{
	return dse3(t.element_[0] - p.element_[1] * t.element_[5] + p.element_[2] * t.element_[4],
		   	t.element_[1] - p.element_[2] * t.element_[3] + p.element_[0] * t.element_[5],
		   	t.element_[2] - p.element_[0] * t.element_[4] + p.element_[1] * t.element_[3],
		   	t.element_[3],
		   	t.element_[4],
		   	t.element_[5]);
}
*/

dse3
rmcp::inv_dAd(const dSE3& T, const dse3& t)
{
	rmcp::dse3 vr;
	vr.inv_dAd(T, t);
	return vr;
}

dse3
rmcp::dad(const dse3& s, const dse3& t)
{
	rmcp::dse3 vr;
	vr.dad(s, t);
	return vr;
}

/*
dse3
rmcp::logR(const dSE3& T)
{
	double theta = 0.5 * (T.element_[0] + T.element_[5] + T.element_[10] - 1.0);
	if (fabs(theta) > 1.0 - rmcp::TINY) {
	   	return dse3(0.0, 0.0, 0.0, T.element_[12], T.element_[13], T.element_[14]);
	}

	theta = acos(theta);
	double cof = theta / (2.0 * sin(theta));
	return dse3( cof * (T.element_[6] - T.element_[9]),
		   	cof * (T.element_[8] - T.element_[2]),
		   	cof * (T.element_[1] - T.element_[4]),
		   	T.element_[12],
		   	T.element_[13],
		   	T.element_[14] );
}
*/
dSE3
rmcp::inv(const dSE3& T)
{
	rmcp::dSE3 mr(T);
	mr.inv();
	return mr;
}
/*
dse3
rmcp::linearize(const dSE3& T)
{
   	return dse3( 0.5 * (T.element_[6] - T.element_[9]),
		   	0.5 * (T.element_[8] - T.element_[2]),
		   	0.5 * (T.element_[1] - T.element_[4]),
		   	T.element_[12],
		   	T.element_[13],
		   	T.element_[14]);
}
*/	

