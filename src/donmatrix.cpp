#include "rmcp/donmatrix.h"

using namespace rmcp;

static int i, j;

donmatrix::donmatrix(const int row) : dormatrix(row)
{
}

donmatrix::donmatrix(const int row, const double element)
	: dormatrix(row, element)
{
}

donmatrix::donmatrix(const int row, const double element[])
	: dormatrix(row, element)
{
}

donmatrix::donmatrix(const donmatrix &m) 
	: dormatrix(m.row_, m.element_)
{
}


void
donmatrix::inv(void)
{
	// inverse is equal to transpose
	
	static double *element;
	element = new double[size_];

	cblas_dcopy(size_, element_, 1, element, 1);

	for (i = 0; i < row_; i++) {
		for (j = 0; j < row_; j++) {
			element_[i + j * row_] = element[j + i * row_];
		}
	}

	delete element;

}


donmatrix
rmcp::inv(const donmatrix &m)
{
	donmatrix mr(m.row_);

	for (i = 0; i < m.row_; i++) {
		for (j = 0; j < m.row_; j++) {
			mr.element_[i + j * m.row_] = m.element_[j + i * m.row_];
		}
	}

	return mr;
}
