#include "rmcp/dsqobjectfn.h"
#include "rmcp/dsqoptimize.h"
#include <cmath>

class fmincon : public rmcp::dsqobjectfn {
	public:
		double evaluate(double x[]);
		void gradient(double x[], double grad[]) { }
		void nonlinear_constraint(double x[], double con[], int err[]) { }
		void nonlinear_constraint(double x[], int i, double &con, int &err) { }
		void nonlinear_gradient(double x[], int i, double grad[]) { }
};

double
fmincon::evaluate(double x[])
{
	return 8.0*pow( (120.0-x[0])/120.0 ,2) +
	7.75* (1.64*x[0]+360.0 ) * pow( (x[0]+360.0-x[1])/(x[0]+360.0) ,4)  +
	3.0* (2.69*x[0]+2.98*x[1]+800) * pow((x[0]+x[1]+800-x[2])/(x[0]+x[1]+800),4)
	+ 5.2*(4.41*x[0]+8.88*x[1]+2.94*x[2])*
	pow( (x[0]+x[1]+x[2]-x[3])/(x[0]+x[1]+x[2]) , 5)
	+ 10.6*(7.23*x[0]+26.42*x[3]+8.64*x[2]+x[3])*
	pow( (x[0]+x[1]+x[2]+x[3]-x[4])/(x[0]+x[1]+x[2]+x[3]) ,5)  ;
}

class portfol1 : public rmcp::dsqobjectfn {
public:
	double evaluate(double x[]);
	void gradient(double x[], double grad[]) { }
	void nonlinear_constraint(double x[], double con[], int err[]) { }
	void nonlinear_constraint(double x[], int i, double &con, int &err) { }
	void nonlinear_gradient(double x[], int i, double grad[]) { }
};

double
portfol1::evaluate(double x[])
{
	static double aloc[13] = { 0.05,  0.1, 0.15, 0.08, 0.07, 0.1,
		0.13, 0.12, 0.056, 0.021, 0.063,0.036, 0.024};
	static int i; 
	static double sum;

	sum = 0.0;
	for (i = 0; i < 13; i++) {
		sum += pow(((x[i]-aloc[i])/aloc[i]), 2);
	}
	return sum;
}

class simple : public rmcp::dsqobjectfn {
public:
	double evaluate(double x[]);
	void gradient(double x[], double grad[]);
	void nonlinear_constraint(double x[], double con[], int err[]);
	void nonlinear_constraint(double x[], int i, double &con, int &err);
	void nonlinear_gradient(double x[], int i, double grad[]);

};

double
simple::evaluate(double x[])
{
	return x[0] * x[0] + x[1] * x[1];
}

void 
simple::gradient(double x[], double grad[])
{
	grad[0] = 2.0 * x[0];
	grad[1] = 2.0 * x[1];
}

void
simple::nonlinear_constraint(double x[], double con[], int err[])
{
//	err[0] = 0;
	con[0] = x[0] * x[1];
}

void
simple::nonlinear_constraint(double x[], int i, double &con, int &err)
{
	if (i == 0) {
		con = x[0] * x[1];
	}
}

void
simple::nonlinear_gradient(double x[], int i, double grad[])
{
	if (i == 0)	 {
		grad[0] = x[1];
		grad[1] = x[0];
	}
}

class hs114_new : public rmcp::dsqobjectfn {
	public:
		double evaluate(double x[]);
		void gradient(double x[], double grad[]);
		void nonlinear_constraint(double x[], double con[], int err[]);
		void nonlinear_constraint(double x[], int i, double &con, int &err);
		void nonlinear_gradient(double x[], int i, double grad[]);

};

double
hs114_new::evaluate(double x[])
{
	return 5.04e0 * x[0]
	+ .035e0 * x[1]
	+ 10.e0 * x[2]
	+ 3.36e0 * x[4]
	- .063e0 * x[3] * x[6];
}

void
hs114_new::gradient(double x[], double grad[])
{
	grad[0] = 5.04;
	grad[1] = 0.035;
	grad[2] = 10.0;
	grad[3] = -0.063 * x[6];
	grad[4] = 3.36;
	grad[5] = 0.0;
	grad[6] = -0.063 * x[3];
	grad[7] = 0.0;
	grad[8] = 0.0;
	grad[9] = 0.0;
}

void
hs114_new::nonlinear_constraint(double x[], double con[], int err[])
{
	const double a = 0.99;
	const double b = 0.9;

	con[0] = 1.12 * x[0] + 0.13167 * x[0] * x[7] - 0.00667 * x[0] * x[7] * x[7] - a * x[3];
	con[1] = -1.12 * x[0] - 0.13167 * x[0] * x[7] + 0.00667 * x[0] * x[7] * x[7] + 1.0 / a * x[3];
	con[2] = 1.098 * x[7] - 0.038 * x[7] * x[7] + 0.325 * x[5] - a * x[6];
	con[3] = -1.098 * x[7] + 0.038 * x[7] * x[7] - 0.325 * x[5] + 1.0 / a * x[6];
	con[4] = 98000.0 * x[2] / (x[3] * x[8] + 1000.0 * x[2]) - x[5];
	con[5] = (x[1] + x[4]) / x[0] - x[7];
}

void
hs114_new::nonlinear_constraint(double x[], int i, double &con, int &err)
{
	const double a = 0.99;
	const double b = 0.9;

	switch (i) {
		case 0 :
			con = 1.12 * x[0] + 0.13167 * x[0] * x[7] - 0.00667 * x[0] * x[7] * x[7] - a * x[3];
			break;
		case 1:
			con = -1.12 * x[0] - 0.13167 * x[0] * x[7] + 0.00667 * x[0] * x[7] * x[7] + 1.0 / a * x[3];
			break;
		case 2:
			con = 1.098 * x[7] - 0.038 * x[7] * x[7] + 0.325 * x[5] - a * x[6];
			break;
		case 3:
			con = -1.098 * x[7] + 0.038 * x[7] * x[7] - 0.325 * x[5] + 1.0 / a * x[6];
			break;
		case 4:
			con = 98000.0 * x[2] / (x[3] * x[8] + 1000.0 * x[2]) - x[5];
			break;
		case 5:
			con = (x[1] + x[4]) / x[0] - x[7];
			break;
	}
}

void
hs114_new::nonlinear_gradient(double x[], int i, double grad[])
{
	const double a = 0.99;
	const double b = 0.9;

	double temp;
	double temp2;

	switch (i) {
		case 0:
			grad[0] = 1.12 + 0.13167 * x[7] - 0.00667 * x[7] * x[7];
			grad[1] = 0.0;
			grad[2] = 0.0;
			grad[3] = - a;
			grad[4] = 0.0;
			grad[5] = 0.0;
			grad[6] = 0.0;
			grad[7] = 0.13167 * x[0] - 0.01334 * x[0] * x[7];
			grad[8] = 0.0;
			grad[9] = 0.0;
			break;
		case 1:
			grad[0] = -1.12 - 0.13167 * x[7] + 0.00667 * x[7] * x[7];
			grad[1] = 0.0;
			grad[2] = 0.0;
			grad[3] = 1.0 / a;
			grad[4] = 0.0;
			grad[5] = 0.0;
			grad[6] = 0.0;
			grad[7] = -0.13167 * x[0] + 0.01334 * x[0] * x[7];
			grad[8] = 0.0;
			grad[9] = 0.0;
			break;
		case 2:
			grad[0] = 0.0;
			grad[1] = 0.0;
			grad[2] = 0.0;
			grad[3] = 0.0;
			grad[4] = 0.0;
			grad[5] = 0.325;
			grad[6] = -a;
			grad[7] = 1.098 - 0.076 * x[7];
			grad[8] = 0.0;
			grad[9] = 0.0;
			break;
		case 3:
			grad[0] = 0.0;
			grad[1] = 0.0;
			grad[2] = 0.0;
			grad[3] = 0.0;
			grad[4] = 0.0;
			grad[5] = -0.325;
			grad[6] = 1.0 / a;
			grad[7] = -1.098 + 0.076 * x[7];
			grad[8] = 0.0;
			grad[9] = 0.0;
			break;
		case 4:
			temp = 98000.0 / (x[3] * x[8] + 1000.0 * x[2]);
			temp2 = temp / (x[3] * x[8] + 1.e3 * x[2]) * x[2];
			grad[0] = 0.0;
			grad[1] = 0.0;
			grad[2] = temp - 1000 * temp2;
			grad[3] = - x[8] * temp2;
			grad[4] = 0.0;
			grad[5] = -1.0;
			grad[6] = 0.0;
			grad[7] = 0.0;
			grad[8] = - x[3]  * temp2;
			grad[9] = 0.0;
			break;
		case 5:
			grad[0] = - (x[1] + x[4]) / (x[0] * x[0]);
			grad[1] = 1.0 / x[0];
			grad[2] = 0.0;
			grad[3] = 0.0;
			grad[4] = 1.0 / x[0];
			grad[5] = 0.0;
			grad[6] = 0.0;
			grad[7] = -1.0;
			grad[8] = 0.0;
			grad[9] = 0.0;
			break;
	}
			
}

void hs114_new_test()
{
	const double a = 0.99;
	const double b = 0.9;

	hs114_new fn;

	rmcp::dsqoptimize OPT(&fn);
	OPT.init(10, 5, 6, 4000, 20);

	OPT.set_upper(0, 2000.0);
	OPT.set_upper(1, 16000.0);
    OPT.set_upper(2, 120.0);
	OPT.set_upper(3, 5000.0);
	OPT.set_upper(4, 2000.0);
	OPT.set_upper(5, 93.0);
	OPT.set_upper(6, 95.0);
	OPT.set_upper(7, 12.0);
	OPT.set_upper(8, 4.0);
	OPT.set_upper(9, 162.0);

	OPT.set_lower(0, 0.00001);
	OPT.set_lower(1, 0.00001);
	OPT.set_lower(2, 0.00001);
	OPT.set_lower(3, 0.00001);
	OPT.set_lower(4, 0.00001);
	OPT.set_lower(5, 85.0);
	OPT.set_lower(6, 90.0);
	OPT.set_lower(7, 3.0);
	OPT.set_lower(8, 1.2);
	OPT.set_lower(9, 145.0);

	OPT.set_lincon_upper(4, 0.0);

	OPT.set_lincon_lower(0, -35.82);
	OPT.set_lincon_lower(1, 133.0);
	OPT.set_lincon_lower(2, 35.82);
	OPT.set_lincon_lower(3, -133.0);
	OPT.set_lincon_lower(4, 0.0);

	OPT.set_nonlincon_upper(4, 0.0);
	OPT.set_nonlincon_upper(5, 0.0);

	OPT.set_nonlincon_lower(0, 0.0);
	OPT.set_nonlincon_lower(1, 0.0);
	OPT.set_nonlincon_lower(2, -57.425);
	OPT.set_nonlincon_lower(3, 57.425);
	OPT.set_nonlincon_lower(4, 0.0);
	OPT.set_nonlincon_lower(5, 0.0);

	rmcp::dmatrix lincon(5, 10, 0.0);

	lincon(0, 8) = -b;
	lincon(0, 9) = -0.222;
	lincon(1, 6) = 3.0;
	lincon(1, 9) = -a;
	lincon(2, 8) = 1.0 / b;
	lincon(2, 9) = 0.222;
	lincon(3, 6) = -3.0;
	lincon(3, 9) = 1.0 / a;
	lincon(4, 0) = -1.0;
	lincon(4, 3) = 1.22;
	lincon(4, 4) = -1.0;

	OPT.set_analytic(true);
	OPT.set_silent(false);
	OPT.set_logfile("hs114_new");
	OPT.set_linear_constraint(lincon);
	OPT.set_epsdif(1.e-16);
	OPT.set_del0(0.2e0);
	OPT.set_tau0(1.e0);
	OPT.set_te0(true);

	double initial_x_array[] = { 1745.e0, 12.e3, 11.e1, 3048.e0, 1974.e0, 89.2e0, 92.8e0, 8.e0, 3.6e0, 149.e0 };
	rmcp::dvector initial_x(10, initial_x_array);

	OPT.solve(initial_x);

	rmcp::dvector optimal(10);

	optimal = OPT.optimal_solution();

	std::cout << optimal << std::endl;
	std::cout << OPT.optimal_value() << std::endl;

	std::cout << OPT.elapsed_time() << std::endl;
	std::cout << "optite = " << OPT.get_optite() << std::endl;

}

void simple_test()
{
	simple fn;

	rmcp::dsqoptimize OPT(&fn);
	OPT.init(2, 0, 1, 4000, 20);

	OPT.set_upper(0, 100.0);
	OPT.set_upper(1, 100.0);
	OPT.set_lower(0, -100.0);
	OPT.set_lower(1, -100.0);
	OPT.set_nonlincon_upper(0, 5.0);
	OPT.set_nonlincon_lower(0, 5.0);

	OPT.set_analytic(true);
	OPT.set_silent(false);
	OPT.set_logfile("simple2");
	OPT.set_epsdif(1.e-16);
	OPT.set_epsfcn(1.e-16);
	OPT.set_taubnd(1.0);
	OPT.set_del0(1.0e0);
	OPT.set_tau0(1.e1);
	OPT.set_te2(true);
	OPT.set_te3(true);

	rmcp::dvector initial_x(2, 10.0);

	OPT.solve(initial_x);

	rmcp::dvector optimal(2);

	optimal = OPT.optimal_solution();

	std::cout << optimal << std::endl;
	std::cout << OPT.optimal_value() << std::endl;

	std::cout << OPT.elapsed_time() << std::endl;
	std::cout << "optite = " << OPT.get_optite() << std::endl;

}

void fmincon_test(void)
{

	fmincon fn;

	rmcp::dsqoptimize OPT(&fn);
	OPT.init(5, 5, 0, 4000, 20);

	OPT.set_lower(0, 1.e-5);
	OPT.set_lower(1, 1.e-5);
	OPT.set_lower(2, 1.e-5);
	OPT.set_lower(3, 1.e-5);
	OPT.set_lower(4, 1.e-5);

	OPT.set_lincon_upper(0, 19.0);
	OPT.set_lincon_upper(1, 360.0);
	OPT.set_lincon_upper(2, 800.0);
	OPT.set_lincon_upper(3, -1.e-4);
	OPT.set_lincon_upper(4, -1.e-4);

	double lincon_array[] = { 1.0, -1.0, -1.0, -1.0, -1.0,
		1.0, 1.0, -1.0, -1.0, -1.0,
		1.0, 0.0, 1.0, -1.0, -1.0,
		1.0, 0.0, 0.0, 1.0, -1.0,
		1.0, 0.0, 0.0, 0.0, 1.0 };  
	rmcp::dmatrix lincon(5, 5, lincon_array);

	OPT.set_analytic(false);
	OPT.set_silent(false);
	OPT.set_logfile("fmincon");
	OPT.set_linear_constraint(lincon);
	OPT.set_epsdif(1.e-16);
	OPT.set_epsfcn(1.e-16);
	OPT.set_taubnd(1.0);
	OPT.set_del0(0.2e0);
	OPT.set_tau0(1.e-5);

	rmcp::dvector initial_x(5, 0.1);

	OPT.solve(initial_x);

	rmcp::dvector optimal(5);

	optimal = OPT.optimal_solution();

	std::cout << optimal << std::endl;

	std::cout << "optimal value: " << OPT.optimal_value() << std::endl;

	std::cout << "time : " << OPT.elapsed_time() << std::endl;
	std::cout << "optite = " << OPT.get_optite() << std::endl;

}

void portfol1_test(void)
{

	portfol1 fn;

	rmcp::dsqoptimize OPT(&fn);
	OPT.init(13, 5, 0, 4000, 20);

	for (int j = 0; j < 13; j++) {
		OPT.set_upper(j, 1.0);
	}
	for (int j = 0; j < 13; j++) {
		OPT.set_lower(j, 0.0);
	}
	OPT.set_lower(2, 0.2);
	OPT.set_lower(5, 0.1);

	OPT.set_lincon_upper(0, 1.0);
	OPT.set_lincon_upper(1, 0.35);
	OPT.set_lincon_upper(2, 0.3);
	OPT.set_lincon_upper(3, 0.35);

	OPT.set_lincon_lower(0, 1.0);
	OPT.set_lincon_lower(1, 0.35);
	OPT.set_lincon_lower(2, 0.3);
	OPT.set_lincon_lower(3, 0.35);
	OPT.set_lincon_lower(4, 0.3);

	rmcp::dmatrix lincon(5, 13, 0.0);
	int i;
	for (i = 0; i < 13; i++) lincon(0, i) = 1.0;
	lincon(1, 0) = 1.0;
	lincon(1, 1) = 1.0;
	for (i = 8; i < 13; i++) { lincon(1, i) = 1.0; lincon(4, i) = 1.0; }
	lincon(2, 2) = 1.0;
	lincon(2, 3) = 1.0;
	lincon(2, 4) = 1.0;
	lincon(3, 5) = 1.0;
	lincon(3, 6) = 1.0;
	lincon(3, 7) = 1.0;

	OPT.set_linear_constraint(lincon);

	OPT.set_analytic(false);
	OPT.set_silent(false);
	OPT.set_logfile("portfol1");
	OPT.set_difftype(2);
	OPT.set_epsdif(1.e-16);
	OPT.set_epsfcn(1.e-16);
	OPT.set_taubnd(1.0);
	OPT.set_del0(0.2e0);
	OPT.set_tau0(1.e0);
	OPT.set_te0(true);

	rmcp::dvector initial_x(13, 10.e0);

	OPT.solve(initial_x);

	rmcp::dvector optimal(13);

	optimal = OPT.optimal_solution();

	std::cout << optimal << std::endl;

	std::cout << "optimal value: " << OPT.optimal_value() << std::endl;

	std::cout << "time : " << OPT.elapsed_time() << std::endl;

}


void dsqoptimize_test(void)
{
	simple_test();
//	portfol1_test();
//	fmincon_test();
//	hs114_new_test();
}
