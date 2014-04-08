
#include "rmcp/dcspline.h"
#include "utility.h"
#include <fstream>
#include <cmath>

#define NUM 301
#define OFFSET 0.45
#define H 0.05

void dcspline_test(void)
{
	rmcp::dvector x(NUM);


	for (int i = 0; i < NUM; i++) {
		x[i] = (double)i*H + OFFSET;
	}

	rmcp::dvector y(NUM);

	for (int i = 0; i < NUM; i++) {
		y[i] = std::sin(x[i]);
	}

	double yp1 = std::cos(x[0]);
	double ypn = std::cos(x[NUM-1]);

//	double yp1 = (y[1] - y[0]) / H;
//	double ypn = (y[NUM-1] - y[NUM-2]) / H;

	rmcp::dvector xx(NUM);
	for (int i = 0; i < NUM; i++) {
		xx[i] = (double)i*H + OFFSET;
	}

	rmcp::dcspline spl;

	spl.init(rmcp::dcspline::NOTAKNOT, x, y);
//	spl.init(rmcp::dcspline::CLAMPED, x, y, yp1, ypn);
//	spl.init(rmcp::dcspline::NATURAL, x, y);

	rmcp::dvector yy(NUM);
	rmcp::dvector dyy(NUM);
	rmcp::dvector ddyy(NUM);

	
	for (int i = 0; i < NUM; i++) {
		spl.interp(x, xx[i], yy[i], dyy[i], ddyy[i]);
	}
	
//	spl.interp(x, xx, yy, dyy, ddyy);

	std::ofstream fout;
	fout.open("dcspline_test.txt", std::ios::out);

	fout << "xx\tyy\tdyy\tddyy" << std::endl;

	for (int i = 0; i < NUM; i++) {
		fout << xx[i] << "\t" << yy[i] << "\t" << dyy[i] << "\t" << ddyy[i] << "\t" << y[i] << std::endl;
	}
	fout.close();

}
