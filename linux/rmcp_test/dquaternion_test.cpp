


#include "rmcp/dquaternion.h"
#include "utility.h"

void dquaternion_constructor_test(void)
{
	std::cout << " === quaternion constructor === " << std::endl;

	rmcp::dquaternion q1(1, 2, 3, 4);

	util_cout("q1 = ", q1);

	util_pause();
}

void dquaternion_test(void)
{
	dquaternion_constructor_test();
}


