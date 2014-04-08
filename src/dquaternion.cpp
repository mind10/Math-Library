
#include <rmcp/dquaternion.h>

using namespace rmcp;

dquaternion::dquaternion(void) : rmcp::dvector(4, 0.0)
{
	element_[3] = 1.0;
}

dquaternion::dquaternion(double s, double x, double y, double z)
: rmcp::dvector(4)
{
	element_[0] = s;
	element_[1] = x;
	element_[2] = y;
	element_[3] = z;
}

void 
dquaternion::set(const double s, const double x, const double y, const double z)
{
	element_[0] = s;
	element_[1] = x;
	element_[2] = y;
	element_[3] = z;
}

void
dquaternion::get(double &s, double &x, double &y, double &z)
{
	s = element_[0];
	x = element_[1];
	y = element_[2];
	z = element_[3];
}

