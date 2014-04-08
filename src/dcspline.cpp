
#include <rmcp/dcspline.h>
#include <rmcp/dgqmatrix.h>
#include <cassert>
#include <cmath>
#include <fstream>

using namespace rmcp;

dcspline::dcspline() : coeff_(1)
{

}

bool 
dcspline::init(rmcp::dcspline::TYPE type, dvector &x, dvector &y, double yp1, double ypn)
{
	assert ( (x.size() == y.size()) && "dcspline::init()");

	switch (type) {
		case rmcp::dcspline::CLAMPED :
			{
				int size = x.size();

				rmcp::dgqmatrix m((size-1)*4, 0.0);
				rmcp::dvector n((size-1)*4, 0.0);

				coeff_.resize((size-1)*4);
				m(0, 0) = 1.0;
				m(0, 1) = x[0];
				m(0, 2) = x[0] * x[0];
				m(0, 3) = x[0] * x[0] * x[0];

				m(1, 1) = 1.0;
				m(1, 2) = 2.0 * x[0];
				m(1, 3) = 3.0 * x[0] * x[0];

				m(2, (size-1)*4-4) = 1.0;
				m(2, (size-1)*4-3) = x[size-1];
				m(2, (size-1)*4-2) = x[size-1] * x[size-1];
				m(2, (size-1)*4-1) = x[size-1] * x[size-1] * x[size-1];

				m(3, (size-1)*4-3) = 1.0;
				m(3, (size-1)*4-2) = 2.0 * x[size-1];
				m(3, (size-1)*4-1) = 3.0 * x[size-1] * x[size-1];

				int row, column;
				for (int i = 1; i < size-1; i++) {
					row = (i-1)*4+4;
					column = (i-1)*4;
					m(row+0, column+0) = 1.0;
					m(row+0, column+1) = x[i];
					m(row+0, column+2) = x[i] * x[i];
					m(row+0, column+3) = x[i] * x[i] * x[i];

					m(row+1, column+1) = 1.0;
					m(row+1, column+2) = 2.0 * x[i];
					m(row+1, column+3) = 3.0 * x[i] * x[i];
					m(row+1, column+5) = -1.0;
					m(row+1, column+6) = -2.0 * x[i];
					m(row+1, column+7) = -3.0 * x[i] * x[i];

					m(row+2, column+2) = 2.0;
					m(row+2, column+3) = 6.0 * x[i];
					m(row+2, column+6) = -2.0;
					m(row+2, column+7) = -6.0 * x[i];

					m(row+3, column+4) = 1.0;
					m(row+3, column+5) = x[i];
					m(row+3, column+6) = x[i] * x[i];
					m(row+3, column+7) = x[i] * x[i] * x[i];
				}

				n[0] = y[0];
				n[1] = yp1;
				n[2] = y[size-1];
				n[3] = ypn;

				for (int i = 1; i < size-1; i++) {
					row = (i-1)*4+4;
					n[row] = y[i];
					n[row+3] = y[i];
				}
				/*
				std::ofstream fout;
				fout.open("dcspline_m.txt", std::ios::out);
				fout << m << std::endl;
				fout << n << std::endl;
				*/
				m.inv();
				coeff_ = m * n;
				/*
				fout << m << std::endl;
				fout << coeff_ << std::endl;
				fout.close();
				*/
			}
			break;
		case rmcp::dcspline::NATURAL :
			{
				int size = x.size();
				rmcp::dgqmatrix m((size-1)*4, 0.0);
				rmcp::dvector n((size-1)*4, 0.0);
				coeff_.resize((size-1)*4);

				m(0, 0) = 1.0;
				m(0, 1) = x[0];
				m(0, 2) = x[0] * x[0];
				m(0, 3) = x[0] * x[0] * x[0];

				m(1, 2) = 2.0;
				m(1, 3) = 6.0 * x[0];

				m(2, (size-1)*4-4) = 1.0;
				m(2, (size-1)*4-3) = x[size-1];
				m(2, (size-1)*4-2) = x[size-1] * x[size-1];
				m(2, (size-1)*4-1) = x[size-1] * x[size-1] * x[size-1];

				m(3, (size-1)*4-2) = 2.0;
				m(3, (size-1)*4-1) = 6.0 * x[size-1];

				int row, column;
				for (int i = 1; i < size-1; i++) {
					row = (i-1)*4+4;
					column = (i-1)*4;
					m(row+0, column+0) = 1.0;
					m(row+0, column+1) = x[i];
					m(row+0, column+2) = x[i] * x[i];
					m(row+0, column+3) = x[i] * x[i] * x[i];

					m(row+1, column+1) = 1.0;
					m(row+1, column+2) = 2.0 * x[i];
					m(row+1, column+3) = 3.0 * x[i] * x[i];
					m(row+1, column+5) = -1.0;
					m(row+1, column+6) = -2.0 * x[i];
					m(row+1, column+7) = -3.0 * x[i] * x[i];

					m(row+2, column+2) = 2.0;
					m(row+2, column+3) = 6.0 * x[i];
					m(row+2, column+6) = -2.0;
					m(row+2, column+7) = -6.0 * x[i];

					m(row+3, column+4) = 1.0;
					m(row+3, column+5) = x[i];
					m(row+3, column+6) = x[i] * x[i];
					m(row+3, column+7) = x[i] * x[i] * x[i];
				}

				n[0] = y[0];
				n[1] = 0.0;
				n[2] = y[size-1];
				n[3] = 0.0;

				for (int i = 1; i < size-1; i++) {
					row = (i-1)*4+4;
					n[row] = y[i];
					n[row+3] = y[i];
				}

				/*
				std::ofstream fout;
				fout.open("dcspline_m.txt", std::ios::out);
				fout << m << std::endl;
				fout << n << std::endl;
				*/

				m.inv();
				coeff_ = m * n;

				/*
				fout << m << std::endl;
				fout << coeff_ << std::endl;
				fout.close();
				*/
			}
			break;
		case rmcp::dcspline::NOTAKNOT :
			{
				int size = x.size();
				rmcp::dgqmatrix m((size-1)*4, 0.0);
				rmcp::dvector n((size-1)*4, 0.0);
				coeff_.resize((size-1)*4);

				m(0, 0) = 1.0;
				m(0, 1) = x[0];
				m(0, 2) = x[0] * x[0];
				m(0, 3) = x[0] * x[0] * x[0];

				m(1, 3) = 6.0;
				m(1, 7) = -6.0;

				m(2, (size-1)*4-4) = 1.0;
				m(2, (size-1)*4-3) = x[size-1];
				m(2, (size-1)*4-2) = x[size-1] * x[size-1];
				m(2, (size-1)*4-1) = x[size-1] * x[size-1] * x[size-1];

				m(3, (size-1)*4-5) = 6.0;
				m(3, (size-1)*4-1) = -6.0;

				int row, column;
				for (int i = 1; i < size-1; i++) {
					row = (i-1)*4+4;
					column = (i-1)*4;
					m(row+0, column+0) = 1.0;
					m(row+0, column+1) = x[i];
					m(row+0, column+2) = x[i] * x[i];
					m(row+0, column+3) = x[i] * x[i] * x[i];

					m(row+1, column+1) = 1.0;
					m(row+1, column+2) = 2.0 * x[i];
					m(row+1, column+3) = 3.0 * x[i] * x[i];
					m(row+1, column+5) = -1.0;
					m(row+1, column+6) = -2.0 * x[i];
					m(row+1, column+7) = -3.0 * x[i] * x[i];

					m(row+2, column+2) = 2.0;
					m(row+2, column+3) = 6.0 * x[i];
					m(row+2, column+6) = -2.0;
					m(row+2, column+7) = -6.0 * x[i];

					m(row+3, column+4) = 1.0;
					m(row+3, column+5) = x[i];
					m(row+3, column+6) = x[i] * x[i];
					m(row+3, column+7) = x[i] * x[i] * x[i];
				}

				n[0] = y[0];
				n[1] = 0.0;
				n[2] = y[size-1];
				n[3] = 0.0;

				for (int i = 1; i < size-1; i++) {
					row = (i-1)*4+4;
					n[row] = y[i];
					n[row+3] = y[i];
				}

				/*
				std::ofstream fout;
				fout.open("dcspline_m.txt", std::ios::out);
				fout << m << std::endl;
				fout << n << std::endl;
				*/

				m.inv();
				coeff_ = m * n;

				/*
				fout << m << std::endl;
				fout << coeff_ << std::endl;
				fout.close();
				*/
			}
			break;
		default :
			return false;
			break;
	}

	return true;
}

bool
dcspline::interp(dvector &x, double new_x, double &new_y, double &new_dy, double &new_ddy)
{
	int klo = 0;
	int khi = x.size() - 1;
	int k;

	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		if (x[k] > new_x) khi = k;
		else klo = k;
	}

	double A, B, C, D;
	D = coeff_[klo*4+0];
	C = coeff_[klo*4+1];
	B = coeff_[klo*4+2];
	A = coeff_[klo*4+3];

	new_y = D + C * new_x + B * new_x * new_x + A * new_x * new_x * new_x;
	new_dy = C + 2.0 * B * new_x + 3.0 * A * new_x * new_x; 
	new_ddy = 2.0 * B + 6.0 * A * new_x;

	return true;
}

bool 
dcspline::interp(dvector &x, dvector &new_x, dvector &new_y, dvector &new_dy, dvector &new_ddy)
{
	int k, klo, khi;
	double A, B, C, D;
	double xi;

	new_y.resize(new_x.size());
	new_dy.resize(new_x.size());
	new_ddy.resize(new_x.size());

	for (int i = 0; i < new_x.size(); i++) {
		klo = 0;
		khi = x.size() - 1;

		xi = new_x[i]; 
		while (khi - klo > 1) {
			k = (khi + klo) >> 1;
			if (x[k] > xi) khi = k;
			else klo = k;
		}

		D = coeff_[klo*4+0];
		C = coeff_[klo*4+1];
		B = coeff_[klo*4+2];
		A = coeff_[klo*4+3];
		
		new_y[i] = D + C * xi + B * pow(xi, 2) + A * pow(xi, 3);
		new_dy[i] = C + 2.0 * B * xi + 3.0 * A * std::pow(xi, 2); 
		new_ddy[i] = 2.0 * B + 6.0 * A * xi;
	}

	return true;
}
