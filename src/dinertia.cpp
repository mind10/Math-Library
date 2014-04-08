
#include "rmcp/dinertia.h"
#include "rmcp/dvector.h"
#include "rmcp/dmatrix.h"
#include "rmcp/dse3.h"

#include <cassert>

using namespace rmcp;

dinertia::dinertia(void)
{
   	I_[0] = 1.0;
   	I_[1] = 1.0;
	I_[2] = 1.0;
	I_[3] = 0.0;
	I_[4] = 0.0;
	I_[5] = 0.0;
	
	r_[0] = 0.0;
	r_[1] = 0.0;
	r_[2] = 0.0;
	
	m_ = 1.0; 
}

dinertia::dinertia(double m)
{
   	I_[0] = m;
	I_[1] = m;
	I_[2] = m;
   	I_[3] = 0.0;
	I_[4] = 0.0;
	I_[5] = 0.0;
	
	r_[0] = 0.0;
	r_[1] = 0.0;
	r_[2] = 0.0;
	
	m_ = m;
}

dinertia::dinertia(double mass, double Ixx, double Iyy, double Izz)
{
	m_ = mass;
	
	I_[0] = Ixx;
   	I_[1] = Iyy;
   	I_[2] = Izz;
	I_[3] = 0.0;
	I_[4] = 0.0;
	I_[5] = 0.0;
	
	r_[0] = 0.0;
	r_[1] = 0.0;
	r_[2] = 0.0;
}

dinertia::dinertia(double mass, double Ixx, double Iyy, double Izz, double Ixy,
	   	double Ixz, double Iyz, double r0, double r1, double r2)
{
	m_ = mass;
	
	r_[0] = mass * r0;
	r_[1] = mass * r1;
	r_[2] = mass * r2;
	
	I_[0] = Ixx + mass * (r1 * r1 + r2 * r2);
	I_[1] = Iyy + mass * (r0 * r0 + r2 * r2);
	I_[2] = Izz + mass * (r0 * r0 + r1 * r1);
	I_[3] = Ixy - mass * r0 * r1;
	I_[4] = Ixz - mass * r2 * r0;
	I_[5] = Iyz - mass * r1 * r2;
}

dinertia::dinertia(const dinertia &I)
{
	I_[0] = I.I_[0];
	I_[1] = I.I_[1];
	I_[2] = I.I_[2];
	I_[3] = I.I_[3];
	I_[4] = I.I_[4];
	I_[5] = I.I_[5];
	
	r_[0] = I.r_[0];
	r_[1] = I.r_[1];
	r_[2] = I.r_[2];
	
	m_ = I.m_;	
}

dse3
dinertia::operator * (const dse3 &s) const
{
	return dse3(
		I_[0] * s.element_[0] + I_[3] * s.element_[1] + I_[4] * s.element_[2] +
		   	r_[1] * s.element_[5] -  r_[2] * s.element_[4],
		I_[3] * s.element_[0] + I_[1] * s.element_[1] + I_[5] * s.element_[2] +
			r_[2] * s.element_[3] - r_[0] * s.element_[5],
		 I_[4] * s.element_[0] + I_[5] * s.element_[1] + I_[2] * s.element_[2] +
			r_[0] * s.element_[4] - r_[1] * s.element_[3],
		 s.element_[1] * r_[2] - s.element_[2] * r_[1] + m_ * s.element_[3],
		 s.element_[2] * r_[0] - s.element_[0] * r_[2] + m_ * s.element_[4],
		 s.element_[0] * r_[1] - s.element_[1] * r_[0] + m_ * s.element_[5] );
}

dse3
dinertia::operator * (const rmcp::dvector &s) const
{
	assert (s.size_ == 3 &&
		"dse3 dinertia::operator * (const rmcp::dvector &) const");

	return dse3( r_[1] * s.element_[2] - r_[2] * s.element_[1],
		   		r_[2] * s.element_[0] - r_[0] * s.element_[2],
			   	r_[0] * s.element_[1] - r_[1] * s.element_[0],
			   	m_ * s.element_[0],
			   	m_ * s.element_[1],
			   	m_ * s.element_[2] );
}

dinertia &
dinertia::operator = (const dinertia &I)
{
	I_[0] = I.I_[0];
	I_[1] = I.I_[1];
	I_[2] = I.I_[2];
	I_[3] = I.I_[3];
	I_[4] = I.I_[4];
	I_[5] = I.I_[5];
	r_[0] = I.r_[0];
	r_[1] = I.r_[1];
	r_[2] = I.r_[2];
	m_ = I.m_;

	return *this;
}

dinertia
dinertia::operator + (const dinertia &I) const
{
	dinertia re(I);
	re.I_[0] += I_[0];
	re.I_[1] += I_[1];
	re.I_[2] += I_[2];
	re.I_[3] += I_[3];
	re.I_[4] += I_[4];
	re.I_[5] += I_[5];
	re.r_[0] += r_[0];
	re.r_[1] += r_[1];
	re.r_[2] += r_[2];
	re.m_ += m_;
	return re;
}

dinertia &
dinertia::operator *= (double x)
{
	I_[0] *= x;
	I_[1] *= x;
	I_[2] *= x;
	I_[3] *= x;
	I_[4] *= x;
	I_[5] *= x;
	r_[0] *= x;
	r_[1] *= x;
	r_[2] *= x;
	m_ *= x;

	return *this;
}

dainertia
dinertia::operator + (const dainertia &J) const
{
	return dainertia(J.J_[0] + I_[0], J.J_[6] + I_[3], J.J_[7] + I_[1],
		   	J.J_[12] + I_[4], J.J_[13] + I_[5], J.J_[14] + I_[2],
			J.J_[18], J.J_[19] + r_[2], J.J_[20]-r_[1], J.J_[24]-r_[2],
		   	J.J_[25], J.J_[26] + r_[0], J.J_[30] + r_[1], J.J_[31]-r_[0],
		   	J.J_[32], J.J_[21] + m_, J.J_[27], J.J_[28] + m_, J.J_[33],
		   	J.J_[34], J.J_[35] + m_);
}

dainertia
dinertia::operator - (const dainertia &J) const
{
	return dainertia( -J.J_[0] + I_[0], -J.J_[6] + I_[3], -J.J_[7] + I_[1],
		   	-J.J_[12] + I_[4], -J.J_[13] + I_[5], -J.J_[14] + I_[2],
		   	-J.J_[18], -J.J_[19] + r_[2], -J.J_[20] - r_[1],
		   	-J.J_[24] - r_[2], -J.J_[25], -J.J_[26] + r_[0],
		   	-J.J_[30] + r_[1], -J.J_[31] - r_[0], -J.J_[32],
			-J.J_[21] + m_, -J.J_[27], -J.J_[28] + m_,
		   	-J.J_[33], -J.J_[34], -J.J_[35] + m_ ); 
}

bool
dinertia::operator == (const dinertia &I) const
{
   	return  (I_[0] == I.I_[0]) &&
	   	(I_[1] == I.I_[1]) &&
	   	(I_[2] == I.I_[2]) &&
	   	(I_[3] == I.I_[3]) &&
	   	(I_[4] == I.I_[4]) &&
	   	(I_[5] == I.I_[5]) &&
	   	(r_[0] == I.r_[0]) &&
	   	(r_[1] == I.r_[1]) &&
	   	(r_[2] == I.r_[2]) &&
	   	(m_ == I.m_);
}

void
dinertia::set_mass(double mass)
{
	assert (mass != 0.0 &&
			"void dinertia::set_mass(double)");

	// firstly, remove the effect of old mass(m_)
 	r_[0] /= m_;
	r_[1] /= m_;
	r_[2] /= m_;

	I_[0] -= m_ * (r_[1] * r_[1] + r_[2] * r_[2]); // Ixx
   	I_[1] -= m_ * (r_[0] * r_[0] + r_[2] * r_[2]); // Iyy
   	I_[2] -= m_ * (r_[0] * r_[0] + r_[1] * r_[1]); // Izz
   	I_[3] += m_ * r_[0] * r_[1]; // Ixy
   	I_[4] += m_ * r_[2] * r_[0]; // Ixz
   	I_[5] += m_ * r_[1] * r_[2]; // Iyz

	// secondly, update mass and center of mass and inertia
	I_[0] += mass * (r_[1] * r_[1] + r_[2] * r_[2]);
	I_[1] += mass * (r_[0] * r_[0] + r_[2] * r_[2]);
	I_[2] += mass * (r_[0] * r_[0] + r_[1] * r_[1]);
	I_[3] -= mass * r_[0] * r_[1];
	I_[4] -= mass * r_[2] * r_[0];
	I_[5] -= mass * r_[1] * r_[2];
	
	r_[0] *= mass;
	r_[1] *= mass;
	r_[2] *= mass;

	m_ = mass;
}

double
dinertia::mass(void) const
{
   	return m_;
}

void
dinertia::set_inertia(double Ixx, double Iyy, double Izz, double Ixy,
	   	double Ixz, double Iyz)
{
	static double r0, r1, r2;
	
 	r0 = r_[0] / m_;
	r1 = r_[1] / m_;
	r2 = r_[2] / m_;

	I_[0] = Ixx + m_ * (r1 * r1 + r2 * r2);
	I_[1] = Iyy + m_ * (r0 * r0 + r2 * r2);
	I_[2] = Izz + m_ * (r0 * r0 + r1 * r1);
	I_[3] = Ixy - m_ * r0 * r1;
	I_[4] = Ixz - m_ * r2 * r0;
	I_[5] = Iyz - m_ * r1 * r2;
}

void
dinertia::get_inertia(double I[])
{
	static double r0, r1, r2;
	
 	r0 = r_[0] / m_;
	r1 = r_[1] / m_;
	r2 = r_[2] / m_;

   	I[0] = I_[0] - m_ * (r1 * r1 + r2 * r2); // Ixx
   	I[1] = I_[1] - m_ * (r0 * r0 + r2 * r2); // Iyy
   	I[2] = I_[2] - m_ * (r0 * r0 + r1 * r1); // Izz
   	I[3] = I_[3] + m_ * r0 * r1; // Ixy
   	I[4] = I_[4] + m_ * r2 * r0; // Ixz
   	I[5] = I_[5] + m_ * r1 * r2; // Iyz
} 

void
dinertia::set_com(const rmcp::dvector &v)
{
	assert (v.size_ == 3 &&
		"void dinertia::set_com(const rmcp::dvector &)");

	// firstly, remove the effect of old mass(m_)
 	r_[0] /= m_;
	r_[1] /= m_;
	r_[2] /= m_;

	I_[0] -= m_ * (r_[1] * r_[1] + r_[2] * r_[2]); // Ixx
   	I_[1] -= m_ * (r_[0] * r_[0] + r_[2] * r_[2]); // Iyy
   	I_[2] -= m_ * (r_[0] * r_[0] + r_[1] * r_[1]); // Izz
   	I_[3] += m_ * r_[0] * r_[1]; // Ixy
   	I_[4] += m_ * r_[2] * r_[0]; // Ixz
   	I_[5] += m_ * r_[1] * r_[2]; // Iyz

	// secondly, update mass and center of mass and inertia
	I_[0] += m_ * (v.element_[1] * v.element_[1] + v.element_[2] * v.element_[2]);
	I_[1] += m_ * (v.element_[0] * v.element_[0] + v.element_[2] * v.element_[2]);
	I_[2] += m_ * (v.element_[0] * v.element_[0] + v.element_[1] * v.element_[1]);
	I_[3] -= m_ * v.element_[0] * v.element_[1];
	I_[4] -= m_ * v.element_[2] * v.element_[0];
	I_[5] -= m_ * v.element_[1] * v.element_[2];
	
	r_[0] = m_ * v.element_[0];
	r_[1] = m_ * v.element_[1];
	r_[2] = m_ * v.element_[2];
}

rmcp::dvector
dinertia::com(void) const
{
	static double tmp[3];
	
	tmp[0] = r_[0] / m_;
	tmp[1] = r_[1] / m_;
	tmp[2] = r_[2] / m_;
	
   	return rmcp::dvector(3, tmp);
}

dinertia
dinertia::transform(const dSE3 &T) const
{
	dinertia J(m_);

	J.I_[0] = I_[0] + m_ * T.element_[14] * T.element_[14] + m_ * T.element_[13] * T.element_[13] -
	   	2.0 * r_[2] * T.element_[14] - 2.0 * r_[1] * T.element_[13];
	J.I_[3] = I_[3] + T.element_[13] * r_[0] + T.element_[12] * r_[1] - 
		m_ * T.element_[13] * T.element_[12];
	J.I_[4] = I_[4] + T.element_[14] * r_[0] + T.element_[12] * r_[2] - 
		m_ * T.element_[14] * T.element_[12];
	J.I_[1] = I_[1] + m_ * T.element_[14] * T.element_[14] + m_ * T.element_[12] * T.element_[12] - 
		2.0 * r_[2] * T.element_[14] - 2.0 * r_[0] * T.element_[12];
	J.I_[5] = I_[5] + T.element_[14] * r_[1] + T.element_[13] * r_[2] - 
		m_ * T.element_[14] * T.element_[13];
	J.I_[2] = I_[2] + m_ * T.element_[13] * T.element_[13] + m_ * T.element_[12] * T.element_[12] - 
		2.0 * r_[1] * T.element_[13] - 2.0 * r_[0] * T.element_[12];

	// _tmp = Transpose of Rotation Part of T * J.I_
	double tmp[9] = { T.element_[0] * J.I_[0] + T.element_[1] * J.I_[3] + T.element_[2] * J.I_[4],
		T.element_[4] * J.I_[0] + T.element_[5] * J.I_[3] + T.element_[6] * J.I_[4],
		T.element_[8] * J.I_[0] + T.element_[9] * J.I_[3] + T.element_[10] * J.I_[4],
		T.element_[0] * J.I_[3] + T.element_[1] * J.I_[1] + T.element_[2] * J.I_[5],
		T.element_[4] * J.I_[3] + T.element_[5] * J.I_[1] + T.element_[6] * J.I_[5],
		T.element_[8] * J.I_[3] + T.element_[9] * J.I_[1] + T.element_[10] * J.I_[5],
		T.element_[0] * J.I_[4] + T.element_[1] * J.I_[5] + T.element_[2] * J.I_[2],
		T.element_[4] * J.I_[4] + T.element_[5] * J.I_[5] + T.element_[6] * J.I_[2],
		T.element_[8] * J.I_[4] + T.element_[9] * J.I_[5] + T.element_[10] * J.I_[2] };
	// J.I_ = tmp * Rotation Part of T
	J.I_[0] = tmp[0] * T.element_[0] + tmp[3] * T.element_[1] + tmp[6] * T.element_[2];
	J.I_[3] = tmp[1] * T.element_[0] + tmp[4] * T.element_[1] + tmp[7] * T.element_[2];
	J.I_[4] = tmp[2] * T.element_[0] + tmp[5] * T.element_[1] + tmp[8] * T.element_[2];
	J.I_[1] = tmp[1] * T.element_[4] + tmp[4] * T.element_[5] + tmp[7] * T.element_[6];
	J.I_[5] = tmp[2] * T.element_[4] + tmp[5] * T.element_[5] + tmp[8] * T.element_[6];
	J.I_[2] = tmp[2] * T.element_[8] + tmp[5] * T.element_[9] + tmp[8] * T.element_[10];

	J.r_[0] = T.element_[0] * (r_[0] - m_ * T.element_[12]) +
	   	T.element_[1] * (r_[1] - m_ * T.element_[13]) + T.element_[2] * (r_[2] - m_ * T.element_[14]);
	J.r_[1] = T.element_[4] * (r_[0] - m_ * T.element_[12]) +
	   	T.element_[5] * (r_[1] - m_ * T.element_[13]) + T.element_[6] * (r_[2] - m_ * T.element_[14]);
	J.r_[2] = T.element_[8] * (r_[0] - m_ * T.element_[12]) +
	   	T.element_[9] * (r_[1] - m_ * T.element_[13]) + T.element_[10] * (r_[2] - m_ * T.element_[14]);

	return J;
}

void
dinertia::to_array(double I[]) const
{
	I[0] = I_[0];
	I[6] = I_[3];
	I[12] = I_[4];
	I[18] = 0.0;
	I[24] = -r_[2];
	I[30] = r_[1];
	I[1] = I_[3];
	I[7] = I_[1];
	I[13] = I_[5];
	I[19] = r_[2];
	I[25] = 0.0;
	I[31] = -r_[0];
	I[2] = I_[4];
	I[8] = I_[5];
	I[14] = I_[2];
	I[20] = -r_[1];
	I[26] = r_[0];
	I[32] = 0.0;
	I[3] = 0.0;
	I[9] = r_[2];
	I[15] = -r_[1];
	I[21] = m_;
	I[27] = 0.0;
	I[33] = 0.0;
	I[4] = -r_[2];
	I[10] = 0.0;
	I[16] = r_[0];
	I[22] = 0.0;
	I[28] = m_;
	I[34] = 0.0;
	I[5] =  r_[1];
	I[11] = -r_[0];
	I[17] = 0.0;
	I[23] = 0.0;
	I[29] = 0.0;
	I[35] = m_;
}

dainertia::dainertia(void)
{
	J_[0] = 0.0;
	J_[6] = 0.0;
   	J_[7] = 0.0;
	J_[12] = 0.0;
	J_[13] = 0.0;
	J_[14] = 0.0;
	J_[18] = 0.0;
	J_[19] = 0.0;
	J_[20] = 0.0;
	J_[24] = 0.0;
	J_[25] = 0.0;
	J_[26] = 0.0;
	J_[30] = 0.0;
	J_[31] = 0.0;
	J_[32] = 0.0;
	J_[21] = 0.0;
	J_[27] = 0.0;
	J_[28] = 0.0;
	J_[33] = 0.0;
	J_[34] = 0.0;
	J_[35] = 0.0;	
}

dainertia::dainertia(double d)
{
	J_[0] = d;
	J_[6] = d;
	J_[7] = d;
	J_[12] = d;
   	J_[13] = d;
	J_[14] = d;
	J_[18] = d;
	J_[19] = d;
	J_[20] = d;
	J_[24] = d;
	J_[25] = d;
	J_[26] = d;
	J_[30] = d;
	J_[31] = d;
	J_[32] = d;
	J_[21] = d;
	J_[27] = d;
	J_[28] = d;
	J_[33] = d;
	J_[34] = d;
	J_[35] = d;
}

dainertia::dainertia(const dinertia& I)
{
	J_[0] = I.I_[0];
	J_[6] = I.I_[3];
	J_[7] = I.I_[1];
	J_[12] = I.I_[4];
	J_[13] = I.I_[5];
	J_[14] = I.I_[2];
	J_[18] = 0.0;
	J_[25] = 0.0;
	J_[32] = 0.0;
	J_[27] = 0.0;
	J_[33] = 0.0;
	J_[34] = 0.0;
	J_[26] = I.r_[0];
	J_[31] = -I.r_[0];
	J_[30] = I.r_[1];
	J_[20] = -I.r_[1];
	J_[19] = I.r_[2];
	J_[24] = -I.r_[2];
	J_[21] = I.m_;
	J_[28] = I.m_;
	J_[35] = I.m_;
}

dainertia::dainertia(double a0, double a3, double a4, double a6, double a7,
	   	double a8, double b0, double b1, double b2, double b3, double b4,
	   	double b5, double b6, double b7, double b8, double c0, double c3,
	   	double c4, double c6, double c7, double c8)
{
	J_[0] = a0;
	J_[6] = a3;
	J_[7] = a4;
	J_[12] = a6;
	J_[13] = a7;
	J_[14] = a8;
	J_[18] = b0;
	J_[19] = b1;
	J_[20] = b2;
	J_[24] = b3;
	J_[25] = b4;
	J_[26] = b5;
	J_[30] = b6;
	J_[31] = b7;
	J_[32] = b8;
	J_[21] = c0;
	J_[27] = c3;
	J_[28] = c4;
	J_[33] = c6;
	J_[34] = c7;
	J_[35] = c8;
}

dainertia::dainertia(const dainertia& J)
{
	J_[0] = J.J_[0];
	J_[6] = J.J_[6];
	J_[7] = J.J_[7];
	J_[12] = J.J_[12];
	J_[13] = J.J_[13];
	J_[14] = J.J_[14];
	J_[18] = J.J_[18];
	J_[19] = J.J_[19];
	J_[20] = J.J_[20];
	J_[24] = J.J_[24];
	J_[25] = J.J_[25];
	J_[26] = J.J_[26];
	J_[30] = J.J_[30];
	J_[31] = J.J_[31];
	J_[32] = J.J_[32];
	J_[21] = J.J_[21];
	J_[27] = J.J_[27];
	J_[28] = J.J_[28];
	J_[33] = J.J_[33];
	J_[34] = J.J_[34];
	J_[35] = J.J_[35];
}

dainertia 
dainertia::operator - (void) const
{
	return dainertia( -J_[0], -J_[6], -J_[7], -J_[12], -J_[13],
		   	-J_[14], -J_[18], -J_[19], -J_[20], -J_[24],
		   	-J_[25], -J_[26], -J_[30], -J_[31], -J_[32],
		   	-J_[21], -J_[27], -J_[28], -J_[33], -J_[34],
		   	-J_[35] );
}

dse3
dainertia::operator * (const dse3& a) const
{
	return dse3(J_[0] * a.element_[0] + J_[6] * a.element_[1] + J_[12] * a.element_[2] + J_[18] * a.element_[3] + J_[24] * a.element_[4] + J_[30] * a.element_[5],
		J_[6] * a.element_[0] + J_[7] * a.element_[1] + J_[13] * a.element_[2] + J_[19] * a.element_[3] + J_[25] * a.element_[4] + J_[31] * a.element_[5],
		J_[12] * a.element_[0] + J_[13] * a.element_[1] + J_[14] * a.element_[2] + J_[20] * a.element_[3] + J_[26] * a.element_[4] + J_[32] * a.element_[5],
		J_[18] * a.element_[0] + J_[19] * a.element_[1] + J_[20] * a.element_[2] + J_[21] * a.element_[3] + J_[27] * a.element_[4] + J_[33] * a.element_[5],
		J_[24] * a.element_[0] + J_[25] * a.element_[1] + J_[26] * a.element_[2] + J_[27] * a.element_[3] + J_[28] * a.element_[4] + J_[34] * a.element_[5],
		J_[30] * a.element_[0] + J_[31] * a.element_[1] + J_[32] * a.element_[2] + J_[33] * a.element_[3] + J_[34] * a.element_[4] + J_[35] * a.element_[5]);
}

dainertia
dainertia::operator + (const dainertia& J) const
{
	return dainertia(J_[0] + J.J_[0], J_[6] + J.J_[6], J_[7] + J.J_[7], J_[12] + J.J_[12], J_[13] + J.J_[13], J_[14] + J.J_[14], 
		J_[18] + J.J_[18], J_[19] + J.J_[19], J_[20] + J.J_[20], J_[24] + J.J_[24], J_[25] + J.J_[25], J_[26] + J.J_[26], J_[30] + J.J_[30], J_[31] + J.J_[31], J_[32] + J.J_[32],
		J_[21] + J.J_[21], J_[27] + J.J_[27], J_[28] + J.J_[28], J_[33] + J.J_[33], J_[34] + J.J_[34], J_[35] + J.J_[35]);
}

dainertia
dainertia::operator + (const dinertia& J) const
{
	return dainertia(J_[0] + J.I_[0], J_[6] + J.I_[3], J_[7] + J.I_[1], J_[12] + J.I_[4], J_[13] + J.I_[5], J_[14] + J.I_[2],
		J_[18], J_[19] + J.r_[2], J_[20]-J.r_[1], J_[24]-J.r_[2], J_[25], J_[26] + J.r_[0], J_[30] + J.r_[1], J_[31]-J.r_[0], J_[32],
		J_[21] + J.m_, J_[27], J_[28] + J.m_, J_[33], J_[34], J_[35] + J.m_);
}

dainertia
dainertia::operator - (const dainertia& J) const
{
	return dainertia(J_[0] - J.J_[0], J_[6] - J.J_[6], J_[7] - J.J_[7], J_[12] - J.J_[12], J_[13] - J.J_[13], J_[14] - J.J_[14],
		J_[18] - J.J_[18], J_[19] - J.J_[19], J_[20] - J.J_[20], J_[24] - J.J_[24], J_[25] - J.J_[25], J_[26] - J.J_[26], J_[30] - J.J_[30], J_[31] - J.J_[31], J_[32] - J.J_[32],
		J_[21] - J.J_[21], J_[27] - J.J_[27], J_[28] - J.J_[28], J_[33] - J.J_[33], J_[34] - J.J_[34], J_[35] - J.J_[35]);
}

dainertia
dainertia::operator - (const dinertia& J) const
{
	return dainertia(J_[0] - J.I_[0], J_[6] - J.I_[3], J_[7] - J.I_[1], J_[12] - J.I_[4], J_[13] - J.I_[5], J_[14] - J.I_[2],
		J_[18], J_[19] - J.r_[2], J_[20]+J.r_[1], J_[24]+J.r_[2], J_[25], J_[26] - J.r_[0], J_[30] - J.r_[1], J_[31]+J.r_[0], J_[32],
		J_[21] - J.m_, J_[27], J_[28] - J.m_, J_[33], J_[34], J_[35] - J.m_);
}

dainertia&
dainertia::operator += (const dainertia& J)
{
	J_[0] += J.J_[0];
	J_[6] += J.J_[6];
	J_[7] += J.J_[7];
	J_[12] += J.J_[12];
	J_[13] += J.J_[13];
	J_[14] += J.J_[14];
	J_[18] += J.J_[18];
	J_[19] += J.J_[19];
	J_[20] += J.J_[20];
	J_[24] += J.J_[24];	
	J_[25] += J.J_[25];
	J_[26] += J.J_[26];
	J_[30] += J.J_[30];
	J_[31] += J.J_[31];
	J_[32] += J.J_[32];
	J_[21] += J.J_[21];
	J_[27] += J.J_[27];
	J_[28] += J.J_[28];
	J_[33] += J.J_[33];
	J_[34] += J.J_[34];
	J_[35] += J.J_[35];

	return *this;
}

dainertia&
dainertia::operator -= (const dainertia& J)
{
	J_[0] -= J.J_[0];
	J_[6] -= J.J_[6];
	J_[7] -= J.J_[7];
	J_[12] -= J.J_[12];
	J_[13] -= J.J_[13];
	J_[14] -= J.J_[14];
	J_[18] -= J.J_[18];
	J_[19] -= J.J_[19];
	J_[20] -= J.J_[20];
	J_[24] -= J.J_[24];
	J_[25] -= J.J_[25];
	J_[26] -= J.J_[26];
	J_[30] -= J.J_[30];
	J_[31] -= J.J_[31];
	J_[32] -= J.J_[32];
	J_[21] -= J.J_[21];
	J_[27] -= J.J_[27];	
	J_[28] -= J.J_[28];
	J_[33] -= J.J_[33];
	J_[34] -= J.J_[34];
	J_[35] -= J.J_[35];

	return *this;
}

dainertia&
dainertia::operator += (const dinertia& J)
{
	J_[0] += J.I_[0];
	J_[6] += J.I_[3];
	J_[7] += J.I_[1];
	J_[12] += J.I_[4];
	J_[13] += J.I_[5];
	J_[14] += J.I_[2];
	J_[19] += J.r_[2];
	J_[20] -= J.r_[1];
	J_[24] -= J.r_[2];
	J_[26] += J.r_[0];
	J_[30] += J.r_[1];
	J_[31] -= J.r_[0];
	J_[21] += J.m_;
	J_[28] += J.m_;
	J_[35] += J.m_;

	return *this;
}

dainertia&
dainertia::operator -= (const dinertia& J)
{
	J_[0] -= J.I_[0];
	J_[6] -= J.I_[3];
	J_[7] -= J.I_[1];
	J_[12] -= J.I_[4];
	J_[13] -= J.I_[5];
	J_[14] -= J.I_[2];
	J_[19] -= J.r_[2];
	J_[20] += J.r_[1];
	J_[24] += J.r_[2];
	J_[26] -= J.r_[0];
	J_[30] -= J.r_[1];
	J_[31] += J.r_[0];
	J_[21] -= J.m_;
	J_[28] -= J.m_;
	J_[35] -= J.m_;

	return *this;
}

void
dainertia::to_array(double I[]) const
{
	I[0] = J_[0];
   	I[6] = J_[6];
   	I[12] = J_[12];
   	I[18] = J_[18];
   	I[24] = J_[24];
   	I[30] = J_[30];
	I[1] = J_[6];
   	I[7] = J_[7];
   	I[13] = J_[13];
   	I[19] = J_[19];
   	I[25] = J_[25];
   	I[31] = J_[31];
	I[2] = J_[12];
   	I[8] = J_[13];
   	I[14] = J_[14];
   	I[20] = J_[20];
   	I[26] = J_[26];
   	I[32] = J_[32];
	I[3] = J_[18];
   	I[9] = J_[19];
   	I[15] = J_[20];
   	I[21] = J_[21];
   	I[27] = J_[27];
   	I[33] = J_[33];
	I[4] = J_[24];
   	I[10] = J_[25];
   	I[16] = J_[26];
   	I[22] = J_[27];
   	I[28] = J_[28];
   	I[34] = J_[34];
	I[5] = J_[30];
   	I[11] = J_[31];
   	I[17] = J_[32];
   	I[23] = J_[33];
   	I[29] = J_[34];
   	I[35] = J_[35];
}

dainertia
dainertia::transform(const dSE3& T) const
{
	double R_[] = {	T.element_[0], T.element_[1], T.element_[2], T.element_[4], T.element_[5], T.element_[6], T.element_[8], T.element_[9], T.element_[10] };

	double _D[] = {	J_[18] + T.element_[14] * J_[27] - T.element_[13] * J_[33], J_[19] - T.element_[14] * J_[21] + T.element_[12] * J_[33], J_[20] + T.element_[13] * J_[21] - T.element_[12] * J_[27],
		J_[24] + T.element_[14] * J_[28] - T.element_[13] * J_[34], J_[25] - T.element_[14] * J_[27] + T.element_[12] * J_[34], J_[26] + T.element_[13] * J_[27] - T.element_[12] * J_[28],
		J_[30] + T.element_[14] * J_[34] - T.element_[13] * J_[35], J_[31] - T.element_[14] * J_[33] + T.element_[12] * J_[35], J_[32] + T.element_[13] * J_[33] - T.element_[12] * J_[34] };

	double _E[] = {	J_[0] + T.element_[14] * J_[24] - T.element_[13] * J_[30] + _D[3] * T.element_[14] - _D[6] * T.element_[13], 0.0, 0.0,
		J_[6] + T.element_[14] * J_[25] - T.element_[13] * J_[31] - _D[0] * T.element_[14] + _D[6] * T.element_[12], J_[7] - T.element_[14] * J_[19] + T.element_[12] * J_[31] - _D[1] * T.element_[14] + _D[7] * T.element_[12], 0.0,
		J_[12] + T.element_[14] * J_[26] - T.element_[13] * J_[32] + _D[0] * T.element_[13] - _D[3] * T.element_[12], J_[13] - T.element_[14] * J_[20] + T.element_[12] * J_[32] + _D[1] * T.element_[13] - _D[4] * T.element_[12], J_[14] + T.element_[13] * J_[20] - T.element_[12] * J_[26] + _D[2] * T.element_[13] - _D[5] * T.element_[12]	};

	return dainertia((R_[0] * _E[0] + R_[1] * _E[3] + R_[2] * _E[6]) * R_[0] + (R_[0] * _E[3] + R_[1] * _E[4] + R_[2] * _E[7]) * R_[1] + (R_[0] * _E[6] + R_[1] * _E[7] + R_[2] * _E[8]) * R_[2],
		(R_[0] * _E[0] + R_[1] * _E[3] + R_[2] * _E[6]) * R_[3] + (R_[0] * _E[3] + R_[1] * _E[4] + R_[2] * _E[7]) * R_[4] + (R_[0] * _E[6] + R_[1] * _E[7] + R_[2] * _E[8]) * R_[5],
		(R_[3] * _E[0] + R_[4] * _E[3] + R_[5] * _E[6]) * R_[3] + (R_[3] * _E[3] + R_[4] * _E[4] + R_[5] * _E[7]) * R_[4] + (R_[3] * _E[6] + R_[4] * _E[7] + R_[5] * _E[8]) * R_[5],
		(R_[0] * _E[0] + R_[1] * _E[3] + R_[2] * _E[6]) * R_[6] + (R_[0] * _E[3] + R_[1] * _E[4] + R_[2] * _E[7]) * R_[7] + (R_[0] * _E[6] + R_[1] * _E[7] + R_[2] * _E[8]) * R_[8],
		(R_[3] * _E[0] + R_[4] * _E[3] + R_[5] * _E[6]) * R_[6] + (R_[3] * _E[3] + R_[4] * _E[4] + R_[5] * _E[7]) * R_[7] + (R_[3] * _E[6] + R_[4] * _E[7] + R_[5] * _E[8]) * R_[8],
		(R_[6] * _E[0] + R_[7] * _E[3] + R_[8] * _E[6]) * R_[6] + (R_[6] * _E[3] + R_[7] * _E[4] + R_[8] * _E[7]) * R_[7] + (R_[6] * _E[6] + R_[7] * _E[7] + R_[8] * _E[8]) * R_[8],
		(R_[0] * _D[0] + R_[1] * _D[1] + R_[2] * _D[2]) * R_[0] + (R_[0] * _D[3] + R_[1] * _D[4] + R_[2] * _D[5]) * R_[1] + (R_[0] * _D[6] + R_[1] * _D[7] + R_[2] * _D[8]) * R_[2],
		(R_[3] * _D[0] + R_[4] * _D[1] + R_[5] * _D[2]) * R_[0] + (R_[3] * _D[3] + R_[4] * _D[4] + R_[5] * _D[5]) * R_[1] + (R_[3] * _D[6] + R_[4] * _D[7] + R_[5] * _D[8]) * R_[2],
		(R_[6] * _D[0] + R_[7] * _D[1] + R_[8] * _D[2]) * R_[0] + (R_[6] * _D[3] + R_[7] * _D[4] + R_[8] * _D[5]) * R_[1] + (R_[6] * _D[6] + R_[7] * _D[7] + R_[8] * _D[8]) * R_[2],
		(R_[0] * _D[0] + R_[1] * _D[1] + R_[2] * _D[2]) * R_[3] + (R_[0] * _D[3] + R_[1] * _D[4] + R_[2] * _D[5]) * R_[4] + (R_[0] * _D[6] + R_[1] * _D[7] + R_[2] * _D[8]) * R_[5],
		(R_[3] * _D[0] + R_[4] * _D[1] + R_[5] * _D[2]) * R_[3] + (R_[3] * _D[3] + R_[4] * _D[4] + R_[5] * _D[5]) * R_[4] + (R_[3] * _D[6] + R_[4] * _D[7] + R_[5] * _D[8]) * R_[5],
		(R_[6] * _D[0] + R_[7] * _D[1] + R_[8] * _D[2]) * R_[3] + (R_[6] * _D[3] + R_[7] * _D[4] + R_[8] * _D[5]) * R_[4] + (R_[6] * _D[6] + R_[7] * _D[7] + R_[8] * _D[8]) * R_[5],
		(R_[0] * _D[0] + R_[1] * _D[1] + R_[2] * _D[2]) * R_[6] + (R_[0] * _D[3] + R_[1] * _D[4] + R_[2] * _D[5]) * R_[7] + (R_[0] * _D[6] + R_[1] * _D[7] + R_[2] * _D[8]) * R_[8],
		(R_[3] * _D[0] + R_[4] * _D[1] + R_[5] * _D[2]) * R_[6] + (R_[3] * _D[3] + R_[4] * _D[4] + R_[5] * _D[5]) * R_[7] + (R_[3] * _D[6] + R_[4] * _D[7] + R_[5] * _D[8]) * R_[8],
		(R_[6] * _D[0] + R_[7] * _D[1] + R_[8] * _D[2]) * R_[6] + (R_[6] * _D[3] + R_[7] * _D[4] + R_[8] * _D[5]) * R_[7] + (R_[6] * _D[6] + R_[7] * _D[7] + R_[8] * _D[8]) * R_[8],
		(R_[0] * J_[21] + R_[1] * J_[27] + R_[2] * J_[33]) * R_[0] + (R_[0] * J_[27] + R_[1] * J_[28] + R_[2] * J_[34]) * R_[1] + (R_[0] * J_[33] + R_[1] * J_[34] + R_[2] * J_[35]) * R_[2],
		(R_[0] * J_[21] + R_[1] * J_[27] + R_[2] * J_[33]) * R_[3] + (R_[0] * J_[27] + R_[1] * J_[28] + R_[2] * J_[34]) * R_[4] + (R_[0] * J_[33] + R_[1] * J_[34] + R_[2] * J_[35]) * R_[5],
		(R_[3] * J_[21] + R_[4] * J_[27] + R_[5] * J_[33]) * R_[3] + (R_[3] * J_[27] + R_[4] * J_[28] + R_[5] * J_[34]) * R_[4] + (R_[3] * J_[33] + R_[4] * J_[34] + R_[5] * J_[35]) * R_[5],
		(R_[0] * J_[21] + R_[1] * J_[27] + R_[2] * J_[33]) * R_[6] + (R_[0] * J_[27] + R_[1] * J_[28] + R_[2] * J_[34]) * R_[7] + (R_[0] * J_[33] + R_[1] * J_[34] + R_[2] * J_[35]) * R_[8],
		(R_[3] * J_[21] + R_[4] * J_[27] + R_[5] * J_[33]) * R_[6] + (R_[3] * J_[27] + R_[4] * J_[28] + R_[5] * J_[34]) * R_[7] + (R_[3] * J_[33] + R_[4] * J_[34] + R_[5] * J_[35]) * R_[8],
		(R_[6] * J_[21] + R_[7] * J_[27] + R_[8] * J_[33]) * R_[6] + (R_[6] * J_[27] + R_[7] * J_[28] + R_[8] * J_[34]) * R_[7] + (R_[6] * J_[33] + R_[7] * J_[34] + R_[8] * J_[35]) * R_[8]);
}

void
dainertia::add_transform(const dainertia& J, const dSE3& T)
{
	double R_[] = {	T.element_[0], T.element_[1], T.element_[2], T.element_[4], T.element_[5], T.element_[6], T.element_[8], T.element_[9], T.element_[10] };

	double _D[] = {	J.J_[18] + T.element_[14] * J.J_[27] - T.element_[13] * J.J_[33], J.J_[19] - T.element_[14] * J.J_[21] + T.element_[12] * J.J_[33], J.J_[20] + T.element_[13] * J.J_[21] - T.element_[12] * J.J_[27],
		J.J_[24] + T.element_[14] * J.J_[28] - T.element_[13] * J.J_[34], J.J_[25] - T.element_[14] * J.J_[27] + T.element_[12] * J.J_[34], J.J_[26] + T.element_[13] * J.J_[27] - T.element_[12] * J.J_[28],
		J.J_[30] + T.element_[14] * J.J_[34] - T.element_[13] * J.J_[35], J.J_[31] - T.element_[14] * J.J_[33] + T.element_[12] * J.J_[35], J.J_[32] + T.element_[13] * J.J_[33] - T.element_[12] * J.J_[34] };

	double _E[] = {	J.J_[0] + T.element_[14] * J.J_[24] - T.element_[13] * J.J_[30] + _D[3] * T.element_[14] - _D[6] * T.element_[13], 0.0, 0.0,
		J.J_[6] + T.element_[14] * J.J_[25] - T.element_[13] * J.J_[31] - _D[0] * T.element_[14] + _D[6] * T.element_[12], J.J_[7] - T.element_[14] * J.J_[19] + T.element_[12] * J.J_[31] - _D[1] * T.element_[14] + _D[7] * T.element_[12], 0.0,
		J.J_[12] + T.element_[14] * J.J_[26] - T.element_[13] * J.J_[32] + _D[0] * T.element_[13] - _D[3] * T.element_[12], J.J_[13] - T.element_[14] * J.J_[20] + T.element_[12] * J.J_[32] + _D[1] * T.element_[13] - _D[4] * T.element_[12], J.J_[14] + T.element_[13] * J.J_[20] - T.element_[12] * J.J_[26] + _D[2] * T.element_[13] - _D[5] * T.element_[12]	};

		J_[0] += (R_[0] * _E[0] + R_[1] * _E[3] + R_[2] * _E[6]) * R_[0] + (R_[0] * _E[3] + R_[1] * _E[4] + R_[2] * _E[7]) * R_[1] + (R_[0] * _E[6] + R_[1] * _E[7] + R_[2] * _E[8]) * R_[2];
		J_[6] += (R_[0] * _E[0] + R_[1] * _E[3] + R_[2] * _E[6]) * R_[3] + (R_[0] * _E[3] + R_[1] * _E[4] + R_[2] * _E[7]) * R_[4] + (R_[0] * _E[6] + R_[1] * _E[7] + R_[2] * _E[8]) * R_[5];
		J_[7] += (R_[3] * _E[0] + R_[4] * _E[3] + R_[5] * _E[6]) * R_[3] + (R_[3] * _E[3] + R_[4] * _E[4] + R_[5] * _E[7]) * R_[4] + (R_[3] * _E[6] + R_[4] * _E[7] + R_[5] * _E[8]) * R_[5];
		J_[12] += (R_[0] * _E[0] + R_[1] * _E[3] + R_[2] * _E[6]) * R_[6] + (R_[0] * _E[3] + R_[1] * _E[4] + R_[2] * _E[7]) * R_[7] + (R_[0] * _E[6] + R_[1] * _E[7] + R_[2] * _E[8]) * R_[8];
		J_[13] += (R_[3] * _E[0] + R_[4] * _E[3] + R_[5] * _E[6]) * R_[6] + (R_[3] * _E[3] + R_[4] * _E[4] + R_[5] * _E[7]) * R_[7] + (R_[3] * _E[6] + R_[4] * _E[7] + R_[5] * _E[8]) * R_[8];
		J_[14] += (R_[6] * _E[0] + R_[7] * _E[3] + R_[8] * _E[6]) * R_[6] + (R_[6] * _E[3] + R_[7] * _E[4] + R_[8] * _E[7]) * R_[7] + (R_[6] * _E[6] + R_[7] * _E[7] + R_[8] * _E[8]) * R_[8];
		J_[18] += (R_[0] * _D[0] + R_[1] * _D[1] + R_[2] * _D[2]) * R_[0] + (R_[0] * _D[3] + R_[1] * _D[4] + R_[2] * _D[5]) * R_[1] + (R_[0] * _D[6] + R_[1] * _D[7] + R_[2] * _D[8]) * R_[2];
		J_[19] += (R_[3] * _D[0] + R_[4] * _D[1] + R_[5] * _D[2]) * R_[0] + (R_[3] * _D[3] + R_[4] * _D[4] + R_[5] * _D[5]) * R_[1] + (R_[3] * _D[6] + R_[4] * _D[7] + R_[5] * _D[8]) * R_[2];
		J_[20] += (R_[6] * _D[0] + R_[7] * _D[1] + R_[8] * _D[2]) * R_[0] + (R_[6] * _D[3] + R_[7] * _D[4] + R_[8] * _D[5]) * R_[1] + (R_[6] * _D[6] + R_[7] * _D[7] + R_[8] * _D[8]) * R_[2];
		J_[24] += (R_[0] * _D[0] + R_[1] * _D[1] + R_[2] * _D[2]) * R_[3] + (R_[0] * _D[3] + R_[1] * _D[4] + R_[2] * _D[5]) * R_[4] + (R_[0] * _D[6] + R_[1] * _D[7] + R_[2] * _D[8]) * R_[5];
		J_[25] += (R_[3] * _D[0] + R_[4] * _D[1] + R_[5] * _D[2]) * R_[3] + (R_[3] * _D[3] + R_[4] * _D[4] + R_[5] * _D[5]) * R_[4] + (R_[3] * _D[6] + R_[4] * _D[7] + R_[5] * _D[8]) * R_[5];
		J_[26] += (R_[6] * _D[0] + R_[7] * _D[1] + R_[8] * _D[2]) * R_[3] + (R_[6] * _D[3] + R_[7] * _D[4] + R_[8] * _D[5]) * R_[4] + (R_[6] * _D[6] + R_[7] * _D[7] + R_[8] * _D[8]) * R_[5];
		J_[30] += (R_[0] * _D[0] + R_[1] * _D[1] + R_[2] * _D[2]) * R_[6] + (R_[0] * _D[3] + R_[1] * _D[4] + R_[2] * _D[5]) * R_[7] + (R_[0] * _D[6] + R_[1] * _D[7] + R_[2] * _D[8]) * R_[8];
		J_[31] += (R_[3] * _D[0] + R_[4] * _D[1] + R_[5] * _D[2]) * R_[6] + (R_[3] * _D[3] + R_[4] * _D[4] + R_[5] * _D[5]) * R_[7] + (R_[3] * _D[6] + R_[4] * _D[7] + R_[5] * _D[8]) * R_[8];
		J_[32] += (R_[6] * _D[0] + R_[7] * _D[1] + R_[8] * _D[2]) * R_[6] + (R_[6] * _D[3] + R_[7] * _D[4] + R_[8] * _D[5]) * R_[7] + (R_[6] * _D[6] + R_[7] * _D[7] + R_[8] * _D[8]) * R_[8];
		J_[21] += (R_[0] * J.J_[21] + R_[1] * J.J_[27] + R_[2] * J.J_[33]) * R_[0] + (R_[0] * J.J_[27] + R_[1] * J.J_[28] + R_[2] * J.J_[34]) * R_[1] + (R_[0] * J.J_[33] + R_[1] * J.J_[34] + R_[2] * J.J_[35]) * R_[2];
		J_[27] += (R_[0] * J.J_[21] + R_[1] * J.J_[27] + R_[2] * J.J_[33]) * R_[3] + (R_[0] * J.J_[27] + R_[1] * J.J_[28] + R_[2] * J.J_[34]) * R_[4] + (R_[0] * J.J_[33] + R_[1] * J.J_[34] + R_[2] * J.J_[35]) * R_[5];
		J_[28] += (R_[3] * J.J_[21] + R_[4] * J.J_[27] + R_[5] * J.J_[33]) * R_[3] + (R_[3] * J.J_[27] + R_[4] * J.J_[28] + R_[5] * J.J_[34]) * R_[4] + (R_[3] * J.J_[33] + R_[4] * J.J_[34] + R_[5] * J.J_[35]) * R_[5];
		J_[33] += (R_[0] * J.J_[21] + R_[1] * J.J_[27] + R_[2] * J.J_[33]) * R_[6] + (R_[0] * J.J_[27] + R_[1] * J.J_[28] + R_[2] * J.J_[34]) * R_[7] + (R_[0] * J.J_[33] + R_[1] * J.J_[34] + R_[2] * J.J_[35]) * R_[8];
		J_[34] += (R_[3] * J.J_[21] + R_[4] * J.J_[27] + R_[5] * J.J_[33]) * R_[6] + (R_[3] * J.J_[27] + R_[4] * J.J_[28] + R_[5] * J.J_[34]) * R_[7] + (R_[3] * J.J_[33] + R_[4] * J.J_[34] + R_[5] * J.J_[35]) * R_[8];
		J_[35] += (R_[6] * J.J_[21] + R_[7] * J.J_[27] + R_[8] * J.J_[33]) * R_[6] + (R_[6] * J.J_[27] + R_[7] * J.J_[28] + R_[8] * J.J_[34]) * R_[7] + (R_[6] * J.J_[33] + R_[7] * J.J_[34] + R_[8] * J.J_[35]) * R_[8];
}

void 
dainertia::subtract_kronecker_product(const dse3& x, const dse3& y)
{
	J_[0] -= x.element_[0] * y.element_[0];
	J_[6] -= 0.5 * (x.element_[0] * y.element_[1] + x.element_[1] * y.element_[0]);
	J_[7] -= x.element_[1] * y.element_[1];
	J_[12] -= 0.5 * (x.element_[0] * y.element_[2] + x.element_[2] * y.element_[0]);
	J_[13] -= 0.5 * (x.element_[1] * y.element_[2] + x.element_[2] * y.element_[1]);
	J_[14] -= x.element_[2] * y.element_[2];
	J_[18] -= 0.5 * (x.element_[0] * y.element_[3] + x.element_[3] * y.element_[0]);
	J_[19] -= 0.5 * (x.element_[1] * y.element_[3] + x.element_[3] * y.element_[1]);
	J_[20] -= 0.5 * (x.element_[2] * y.element_[3] + x.element_[3] * y.element_[2]);
	J_[24] -= 0.5 * (x.element_[0] * y.element_[4] + x.element_[4] * y.element_[0]);
	J_[25] -= 0.5 * (x.element_[1] * y.element_[4] + x.element_[4] * y.element_[1]);
	J_[26] -= 0.5 * (x.element_[2] * y.element_[4] + x.element_[4] * y.element_[2]);
	J_[30] -= 0.5 * (x.element_[0] * y.element_[5] + x.element_[5] * y.element_[0]);
	J_[31] -= 0.5 * (x.element_[1] * y.element_[5] + x.element_[5] * y.element_[1]);
	J_[32] -= 0.5 * (x.element_[2] * y.element_[5] + x.element_[5] * y.element_[2]);
	J_[21] -= x.element_[3] * y.element_[3];
	J_[27] -= 0.5 * (x.element_[3] * y.element_[4] + x.element_[4] * y.element_[3]);
	J_[28] -= x.element_[4] * y.element_[4];
    J_[33] -= 0.5 * (x.element_[3] * y.element_[5] + x.element_[5] * y.element_[3]);
	J_[34] -= 0.5 * (x.element_[4] * y.element_[5] + x.element_[5] * y.element_[4]);
	J_[35] -= x.element_[5] * y.element_[5];
}
	
dse3
dainertia::operator % (const dse3& t) const
{
	// symmetric matrix ????
/* TODO	

	rmcp::dgematrix I(6, 6, J_);
	rmcp::dgematrix x(4, 4);
	rmcp::dgematrix S(6, 1, t.element_);
	bool re = solve_poAxeB(I, x, S);
	assert(re);

	return dse3(x[0], x[1], x[2], x[3], x[4], x[5]);

	*/

	return dse3();
}

bool
dainertia::operator == (const dainertia& J) const
{
	return (J_[0] == J.J_[0]) && (J_[1] == J.J_[1]) && (J_[2] == J.J_[2]) &&
	   	(J_[3] == J.J_[3]) && (J_[4] == J.J_[4]) && (J_[5] == J.J_[5]) &&
	   	(J_[6] == J.J_[6]) && (J_[7] == J.J_[7]) && (J_[8] == J.J_[8]) &&
	   	(J_[9] == J.J_[9]) && (J_[10] == J.J_[10]) && (J_[11] == J.J_[11]) &&
	   	(J_[12] == J.J_[12]) && (J_[13] == J.J_[13]) && (J_[14] == J.J_[14]) &&
	   	(J_[15] == J.J_[15]) && (J_[16] == J.J_[16]) && (J_[17] == J.J_[17]) &&
	   	(J_[18] == J.J_[18]) && (J_[19] == J.J_[19]) && (J_[20] == J.J_[20]) &&
	   	(J_[21] == J.J_[21]) && (J_[22] == J.J_[22]) && (J_[23] == J.J_[23]) &&
	   	(J_[24] == J.J_[24]) && (J_[25] == J.J_[25]) && (J_[26] == J.J_[26]) &&
	   	(J_[27] == J.J_[27]) && (J_[28] == J.J_[28]) && (J_[29] == J.J_[29]) &&
		(J_[30] == J.J_[30]) && (J_[31] == J.J_[31]) && (J_[32] == J.J_[32]) &&
	   	(J_[33] == J.J_[33]) && (J_[34] == J.J_[34]) && (J_[35] == J.J_[35]); 
}

rmcp::dgematrix
dainertia::to_dgematrix(void) const
{
	double d[] = {  J_[0], J_[6], J_[12], J_[18], J_[24], J_[30],
					J_[6], J_[7], J_[13], J_[19], J_[25], J_[31],
					J_[12], J_[13], J_[14], J_[20], J_[26], J_[32],
					J_[18], J_[19], J_[20], J_[21], J_[27], J_[33], 
					J_[24], J_[25], J_[26], J_[27], J_[28], J_[34],
					J_[30], J_[31], J_[32], J_[33], J_[34], J_[35] };
	return rmcp::dgematrix(6, 6, d);
}


dinertia
rmcp::operator * (double x, const dinertia& I)
{
	dinertia re;

	re.I_[0] = x * I.I_[0];
	re.I_[1] = x * I.I_[1];
	re.I_[2] = x * I.I_[2];
	re.I_[3] = x * I.I_[3];
	re.I_[4] = x * I.I_[4];
	re.I_[5] = x * I.I_[5];
	re.r_[0] = x * I.r_[0];
	re.r_[1] = x * I.r_[1];
	re.r_[2] = x * I.r_[2];
	re.m_ = x * I.m_;

	return re;
}


dainertia
rmcp::inv(const dinertia& I)
{
	/* TODO

	rmcp::dgematrix mI(6, 6, 0.0);
	
	mI(0, 0) = I.I_[0];
	mI(1, 1) = I.I_[1];
	mI(2, 2) = I.I_[2];
	mI(0, 1) = mI(1,0) = I.I_[3];
	mI(0, 2) = mI(2,0) = I.I_[4];
	mI(1, 2) = mI(2,1) = I.I_[5];

	mI(3, 3) = mI(4,4) = mI(5,5) = I.m_;
	mI(3, 4) = mI(3,5) = mI(4,3) = mI(4,5) = 0.0;
	mI(5,3) = mI(5,4) = mI(0,3) = mI(1,4) = 0.0;
   	mI(2,5) = mI(3,0) = mI(4,1) = mI(5,2) = 0.0;

	mI(2,4) = mI(4,2) = I.r_[0];
	mI(1,5) = mI(5,1) = - I.r_[0];
	mI(0,5) = mI(5,0) = I.r_[1];
	mI(2,3) = mI(3,2) = - I.r_[1];
	mI(1,3) = mI(3,1) = I.r_[2];
	mI(0,4) = mI(4,0) = - I.r_[2];	

	mI = rmcp::inv(mI);

	return dainertia(mI(0,0), mI(0,1), mI(1,1), mI(2,0), mI(2,1), mI(2,2), 
		mI(0,3), mI(1,3), mI(2,3), mI(0,4), mI(1,4), mI(2,4), mI(0,5),
	   	mI(1,5), mI(2,5), 
		mI(3,3), mI(3,4), mI(4,4), mI(5,3), mI(5,4), mI(5,5));
		*/

	return dainertia();
}

dainertia
rmcp::kronecker_product(const dse3& x, const dse3& y)
{
	return dainertia( x.element_[0] * y.element_[0],
		   	0.5 * (x.element_[0] * y.element_[1] + x.element_[1] * y.element_[0]),
		   	x.element_[1] * y.element_[1],
		   	0.5 * (x.element_[0] * y.element_[2] + x.element_[2] * y.element_[0]),
		   	0.5 * (x.element_[1] * y.element_[2] + x.element_[2] * y.element_[1]),
		   	x.element_[2] * y.element_[2],
			0.5 * (x.element_[0] * y.element_[3] + x.element_[3] * y.element_[0]),
		   	0.5 * (x.element_[1] * y.element_[3] + x.element_[3] * y.element_[1]),
		   	0.5 * (x.element_[2] * y.element_[3] + x.element_[3] * y.element_[2]),
		   	0.5 * (x.element_[0] * y.element_[4] + x.element_[4] * y.element_[0]),
		   	0.5 * (x.element_[1] * y.element_[4] + x.element_[4] * y.element_[1]),
		   	0.5 * (x.element_[2] * y.element_[4] + x.element_[4] * y.element_[2]),
		   	0.5 * (x.element_[0] * y.element_[5] + x.element_[5] * y.element_[0]),
		   	0.5 * (x.element_[1] * y.element_[5] + x.element_[5] * y.element_[1]),
		   	0.5 * (x.element_[2] * y.element_[5] + x.element_[5] * y.element_[2]),
			x.element_[3] * y.element_[3],
		   	0.5 * (x.element_[3] * y.element_[4] + x.element_[4] * y.element_[3]),
		   	x.element_[4] * y.element_[4],
		   	0.5 * (x.element_[3] * y.element_[5] + x.element_[5] * y.element_[3]),
		   	0.5 * (x.element_[4] * y.element_[5] + x.element_[5] * y.element_[4]),
		   	x.element_[5] * y.element_[5] );
}

dainertia
rmcp::kronecker_product(const dse3& x, const dse3& y, double d)
{
	return dainertia( d * x.element_[0] * y.element_[0],
		   	d * 0.5 * (x.element_[0] * y.element_[1] + x.element_[1] * y.element_[0]),
		   	d * x.element_[1] * y.element_[1],
		   	d * 0.5 * (x.element_[0] * y.element_[2] + x.element_[2] * y.element_[0]),
		   	d * 0.5 * (x.element_[1] * y.element_[2] + x.element_[2] * y.element_[1]),
		   	d * x.element_[2] * y.element_[2],
			d * 0.5 * (x.element_[0] * y.element_[3] + x.element_[3] * y.element_[0]),
		   	d * 0.5 * (x.element_[1] * y.element_[3] + x.element_[3] * y.element_[1]),
		   	d * 0.5 * (x.element_[2] * y.element_[3] + x.element_[3] * y.element_[2]),
		   	d * 0.5 * (x.element_[0] * y.element_[4] + x.element_[4] * y.element_[0]),
		   	d * 0.5 * (x.element_[1] * y.element_[4] + x.element_[4] * y.element_[1]),
		   	d * 0.5 * (x.element_[2] * y.element_[4] + x.element_[4] * y.element_[2]),
		   	d * 0.5 * (x.element_[0] * y.element_[5] + x.element_[5] * y.element_[0]),
		   	d * 0.5 * (x.element_[1] * y.element_[5] + x.element_[5] * y.element_[1]),
		   	d * 0.5 * (x.element_[2] * y.element_[5] + x.element_[5] * y.element_[2]),
			d * x.element_[3] * y.element_[3],
		   	d * 0.5 * (x.element_[3] * y.element_[4] + x.element_[4] * y.element_[3]),
		   	d * x.element_[4] * y.element_[4],
		   	d * 0.5 * (x.element_[3] * y.element_[5] + x.element_[5] * y.element_[3]),
		   	d * 0.5 * (x.element_[4] * y.element_[5] + x.element_[5] * y.element_[4]),
		   	d * x.element_[5] * y.element_[5] );
}
