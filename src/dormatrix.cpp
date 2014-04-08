
#include "rmcp/dormatrix.h"
#include <cassert>

using namespace rmcp;

dormatrix::dormatrix(const int row) : dgqmatrix(row)
{
}

dormatrix::dormatrix(const int row, const double element)
	: dgqmatrix(row, element)
{
}

dormatrix::dormatrix(const int row, const double element[])
	: dgqmatrix(row, element)
{
}

dormatrix::dormatrix(const dormatrix &m) 
	: dgqmatrix(m.row_, m.element_)
{
}

