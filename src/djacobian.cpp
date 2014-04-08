
#include "rmcp/djacobian.h"

#include <cassert>

using namespace rmcp;

djacobian::djacobian(int dof) : rmcp::dmatrix(6, dof)
{
}

djacobian::djacobian(const djacobian &J)
	: rmcp::dmatrix(J.row_, J.column_, J.element_)
{
}


rmcp::dse3
djacobian::get(const int i)
{
	assert (i > -1 &&
			i < column_ &&
			"rmcp::dse3 djacobian::get(const int)");

	return rmcp::dse3(element_ + 6 * i);

}

void
djacobian::set(const int i, const rmcp::dse3 &S)
{
	assert (i > -1 &&
			i < column_ &&
			"void djacobian::set(const int, const rmcp::dse3)");

	// overwrite S(6 dimensional vector) to ith column of djacobian matrix
	cblas_dcopy(6, S.element_, 1, element_ + 6 * i, 1);
}
