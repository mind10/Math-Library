
#include "rmcp/duquaternion.h"
#include "rmcp/dso3.h"
#include "utility.h"
#include <math.h>

void duquaternion_constructor_test(void)
{
	std::cout << " === duquaternion constructor === " << std::endl;

	rmcp::duquaternion q1(1, 2, 3, 4);

	util_cout("q1 = ", q1);

	util_pause();
}

void duquaternion_convert_test(void)
{
	double n[3] = { 1, 0, 0};
	double theta = 90 * rmcp::PI /180.0;

	rmcp::duquaternion q1(cos(theta/2.0), sin(theta/2.0)*n[0], sin(theta/2.0)*n[1], sin(theta/2.0)*n[2]);
	util_cout("q1 = ", q1);
	
	rmcp::dSO3 R1;
	q1.convert(R1);
	util_cout("R1 = ", R1);

	util_pause();

}

void duquaternion_test(void)
{
	duquaternion_constructor_test();
	duquaternion_convert_test();
}
