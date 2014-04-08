
#include "rmcp/dgqmatrix.h"
#include <cassert>

using namespace rmcp;

static int i, j;

dgqmatrix::dgqmatrix(const int row) : dmatrix(row, row)
{
}

dgqmatrix::dgqmatrix(const int row, const double element)
	: dmatrix(row, row, element)
{
}

dgqmatrix::dgqmatrix(const int row, const double element[])
	: dmatrix(row, row, element)
{
}

dgqmatrix::dgqmatrix(const dgqmatrix &m) 
	: dmatrix(m.row_, m.row_, m.element_)
{
	assert (m.row_ == m.column_ &&
			"dgqmatrix::dgqmatrix(const dgqmatrix &)");
}

bool
dgqmatrix::identity(void)
{
	static bool re;
	re = true;

	for (i = 0; i < row_; i++) {
		for (j = 0; j < i; j++) {
			re = re && (element_[i + j * row_] == 0.0);
		}
		re = re && (element_[i + i * row_] == 1.0);
		for (j = i + 1; j < column_; j++) {
			re = re && (element_[i + j * row_] == 0.0);
		}
	}

	return re;
}

void
dgqmatrix::set_identity(void)
{
	for (i = 0; i < row_; i++) {
		for (j = 0; j < i; j++) {
			element_[i + j * row_] = 0.0;
		}
		element_[i + i * row_] = 1.0;
		for (j = i + 1; j < column_; j++) {
			element_[i + j * row_] = 0.0;
		}
	}
}

void
dgqmatrix::inv(void) 
{
	static int *ipiv;
	static int info;
	static double *work;

	ipiv = new int[row_];
	work = new double[row_];
	
	dgetrf(&row_, &column_, element_, &row_, ipiv, &info);
	dgetri(&row_, element_, &row_, ipiv, work, &row_, &info);

 	delete[] ipiv;
	delete[] work;

	return;
}
