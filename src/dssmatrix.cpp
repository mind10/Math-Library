
#include "rmcp/dssmatrix.h"
#include <cassert>

using namespace rmcp;

static int i, j;

dssmatrix::dssmatrix(const int row) : dgqmatrix(row)
{
}

dssmatrix::dssmatrix(const int row, const double element)
	: dgqmatrix(row)
{
	for (i = 0; i < row_; i++) {
		element_[i + i * row_] = element;
		for (j = 0; j < i; j++) {
			element_[i + j * row_] = element;
		}
		element_[i + i * row_] = element;
		for (j = i + 1; j < row_; j++) {
			element_[i + j * row_] = - element;
		}	
	}
}

dssmatrix::dssmatrix(const int row, const double element[])
	: dgqmatrix(row)
{
	for (i = 0; i < row_; i++ ) {
		for (j = 0; j <= i; j++) {
			element_[i + j * row_] = element[i + i * j];
		}
	}
			
}

dssmatrix::dssmatrix(const dssmatrix &m) 
	: dgqmatrix(m.row_, m.element_)
{
}

