#include "rmcp/dgqmatrix.h"

#include <iostream>

void dgqmatrix_test()
{
	rmcp::dgqmatrix a(5, 3.0);

	std::cout << a;

	a[3] = 1.0;

	a[5] = a[3];

	a(0,0) = 8.0;
	a(3,3) = a(0,0);

	std::cout << a;

	rmcp::dgqmatrix b(a);
	b *= 2.0;

	b.resize(3, 3, 1.5);
	std::cout << b.normalize() << b;

	a.swap_column(0, 3);
	std::cout << a;

	a.swap_row(0, 3);
	std::cout << a;


}

