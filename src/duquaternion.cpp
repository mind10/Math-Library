
#include <rmcp/duquaternion.h>
#include <rmcp/dvector.h>
#include <rmcp/dso3.h>
#include <cmath>

using namespace rmcp;

duquaternion::duquaternion(void) : rmcp::dquaternion()
{

}

duquaternion::duquaternion(double s, double x, double y, double z)
: rmcp::dquaternion(s, x, y, z)
{
	this->normalize();
}

void
duquaternion::euler2quaternion(rmcp::dvector euler_rpy)
{
	element_[0] = cos( (euler_rpy[0] + euler_rpy[2]) * 0.5 ) * cos(euler_rpy[1] * 0.5);
	element_[1] = cos( (euler_rpy[0] - euler_rpy[2]) * 0.5 ) * sin(euler_rpy[1] * 0.5);
	element_[2] = sin( (euler_rpy[0] - euler_rpy[2]) * 0.5 ) * sin(euler_rpy[1] * 0.5);
	element_[3] = sin( (euler_rpy[0] + euler_rpy[2]) * 0.5 ) * cos(euler_rpy[1] * 0.5); 

}

void
duquaternion::quaternion2euler(rmcp::dvector &euler_rpy)  
{
	euler_rpy.element_[0] = atan2( (element_[1]*element_[3] + element_[2]*element_[0]) , -(element_[2]*element_[3] - element_[1]*element_[0]) );
	euler_rpy.element_[1] = acos(-element_[1]*element_[1]-element_[2]*element_[2]+element_[3]*element_[3]+element_[0]*element_[0]);
	euler_rpy.element_[2] = atan2( (element_[1]*element_[3] - element_[2]*element_[0]) , (element_[2]*element_[3]+element_[1]*element_[0]) );
}

void
duquaternion::euler2quaternion(rmcp::dvector axis, double angle)
{
	rmcp::dvector norm_axis = axis.normalize();
	element_[0] = cos(angle * 0.5);
	element_[1] = norm_axis[0]*sin(angle * 0.5);
	element_[2] = norm_axis[1]*sin(angle * 0.5);
	element_[3] = norm_axis[2]*sin(angle * 0.5);
}

void
duquaternion::quaternion2euler(rmcp::dvector axis, double angle)
{
	double sinq;
	angle = 2.0 * acos(element_[0]);
	sinq = sin(angle * 0.5);

	axis.element_[0] = element_[1] * sinq;
	axis.element_[1] = element_[2] * sinq;
	axis.element_[2] = element_[3] * sinq;
}


void 
duquaternion::get_rotation(rmcp::dvector axis, double angle)
{
	double sinq;
	angle = 2.0 * acos(element_[0]);
	sinq = sin(angle * 0.5);

	axis.element_[0] = element_[1] * sinq;
	axis.element_[1] = element_[2] * sinq;
	axis.element_[2] = element_[3] * sinq;
}

void
duquaternion::convert(dSO3 &R)
{
	R.element_[0] = 1.0 - 2.0 * (element_[2] * element_[2] + element_[3] * element_[3]);
	R.element_[1] = 2.0 * element_[1] * element_[2] + 2.0 * element_[0] * element_[3];
	R.element_[2] = -2.0 * element_[0] * element_[2] + 2.0 * element_[1] * element_[3];
	R.element_[3] = 2.0 * element_[1] * element_[2] - 2.0 * element_[0] * element_[3];
	R.element_[4] = 1.0 - 2.0 * (element_[1] * element_[1] + element_[3] * element_[3]);
	R.element_[5] = 2.0 * element_[0] * element_[1] + 2.0 * element_[2] * element_[3];
	R.element_[6] = 2.0 * element_[0] * element_[2] + 2.0 * element_[1] * element_[3];
	R.element_[7] = -2.0 * element_[0] * element_[1] + 2.0 * element_[2] * element_[3];
	R.element_[8] = 1.0 - 2.0 * (element_[1] * element_[1] + element_[2] * element_[2]);
}

