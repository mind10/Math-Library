
#include <rmcp/dpolynomial.h>
#include <rmcp/dgqmatrix.h>
#include <rmcp/dvector.h>

rmcp::dvector 
rmcp::poly3coeff(double x0, double xf, double y0, double yf, double dy0, double dyf)
{
	double mat[16];

    mat[0] = 1.0;
	mat[1] = 1.0;
	mat[2] = 0.0;
	mat[3] = 0.0;

    mat[4] = x0;
	mat[5] = xf;
	mat[6] = 1.0;
	mat[7] = 1.0;

	mat[8] = mat[4] * x0;
	mat[9] = mat[5] * xf;
	mat[10] = 2.0 * mat[4];
	mat[11] = 2.0 * mat[5];

	mat[12] = mat[8] * x0;
    mat[13] = mat[9] * xf;
    mat[14] = 3.0 * mat[8];
    mat[15] = 3.0 * mat[9];

	rmcp::dgqmatrix M(4, mat);
	rmcp::dvector y(4);

	M.inv();

	y[0] = y0;
    y[1] = yf;
    y[2] = dy0;
    y[3] = dyf;

	return M * y;
}

double
rmcp::poly3val(rmcp::dvector coeff, double x)
{
	return 
		coeff[0] +
	   	coeff[1] * x +
	   	coeff[2] * x * x +
		coeff[3] * x * x * x; 
}

double 
rmcp::poly3der(rmcp::dvector coeff, double x)
{
	return
	   	coeff[1] +
		coeff[2] * 2.0 * x +
		coeff[3] * 3.0 * x * x;
}

rmcp::dvector
rmcp::poly5coeff(double x0, double xf, double y0, double yf, double dy0, double dyf, double ddy0, double ddyf)
{
	double mat[36];

	mat[0] = 1.0;
	mat[1] = 1.0;
	mat[2] = 0.0;
	mat[3] = 0.0;
	mat[4] = 0.0;
	mat[5] = 0.0;

	mat[6] = x0;
	mat[7] = xf;
	mat[8] = 1.0;
	mat[9] = 1.0;
	mat[10] = 0.0;
	mat[11] = 0.0;

	mat[12] = mat[6] * x0;
	mat[13] = mat[7] * xf;
	mat[14] = 2.0 * mat[6];
	mat[15] = 2.0 * mat[7];
	mat[16] = 2.0;
	mat[17] = 2.0;

	mat[18] = mat[12] * x0;
	mat[19] = mat[13] * xf;
	mat[20] = 3.0 * mat[12];
	mat[21] = 3.0 * mat[13];
	mat[22] = 6.0 * x0;
	mat[23] = 6.0 * xf;

	mat[24] = mat[18] * x0;
	mat[25] = mat[19] * xf;
	mat[26] = 4.0 * mat[18];
	mat[27] = 4.0 * mat[19];
	mat[28] = 12.0 * mat[12];
	mat[29] = 12.0 * mat[13];

	mat[30] = mat[24] * x0;
	mat[31] = mat[25] * xf;
	mat[32] = 5.0 * mat[24];
	mat[33] = 5.0 * mat[25];
	mat[34] = 20.0 * mat[18];
	mat[35] = 20.0 * mat[19];

	rmcp::dgqmatrix M(6, mat);
	rmcp::dvector y(6);

	M.inv();

	y[0] = y0;
	y[1] = yf;
	y[2] = dy0;
	y[3] = dyf;
	y[4] = ddy0;
	y[5] = ddyf;

	return M * y;
}

double
rmcp::poly5val(rmcp::dvector coeff, double x)
{
	return 
		coeff[0] +
		coeff[1] * x +
		coeff[2] * x * x +
		coeff[3] * x * x * x + 
		coeff[4] * x * x * x * x +
		coeff[5] * x * x * x * x * x;
}

double 
rmcp::poly5der(rmcp::dvector coeff, double x)
{
	return
		coeff[1] +
		coeff[2] * 2.0 * x +
		coeff[3] * 3.0 * x * x +
		coeff[4] * 4.0 * x * x * x +
		coeff[5] * 5.0 * x * x * x * x;
}

double 
rmcp::poly5secder(rmcp::dvector coeff, double x)
{
	return
		coeff[2] * 2.0 +
		coeff[3] * 6.0 * x +
		coeff[4] * 12.0 * x * x +
		coeff[5] * 20.0 * x * x * x;
}
