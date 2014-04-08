
#include <iostream>
#include <rmcp/dvector.h>
#include <rmcp/dmatrix.h>

int main(void)
{
	rmcp::dvector v1(4, 0.0);

	v1[0] = 1.0;
	v1[1] = 2.0;
	v1[2] = 3.0;
	v1[3] = 4.0;

	double element[16] = { 
		0.1, 0.2, 0.3, 0.4, 
		0.5, 0.6, 0.7, 0.8, 
		0.9, 1.0, 1.1, 1.2, 
		1.3, 1.4, 1.5, 1.6
	};
	
	rmcp::dmatrix M(4, 4, element);
	rmcp::dvector v2;

	v2 = M * v1;

	std::cout << "v1 = " << v1 << std::endl;
	std::cout << "M = " << M << std::endl;
	std::cout << "v2 = " << v2 << std::endl;

	return 0;
}

