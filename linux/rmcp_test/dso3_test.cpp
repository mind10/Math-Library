
#include "rmcp/dso3.h"
#include "rmcp/dvector.h"
#include "utility.h"

#include <iostream>

void dso3_constructor_test(void)
{
	std::cout << " === dso3 constructor ===" << std::endl;

	rmcp::dso3 a;
	util_cout("a = ", a);
	
	rmcp::dso3 b(2.0);
	util_cout("b = ", b);

	double element[] = { 1.0, 2.0, 3.0 };

	rmcp::dso3 c(element);
	util_cout("c = ", c);

	rmcp::dso3 d(1.0, 2.0, 3.0);
	util_cout("d = ", d);

	rmcp::dvector v(3, element);
	rmcp::dso3 e(v);
	util_cout("e = ", e);

	rmcp::dso3 f(e);

	util_cout("f = ", f);

}
void dso3_operator_test(void)
{
	std::cout << " === dso3 operator ===" << std::endl;

	double element[] = { 1.1, 2.2, 3.3};

	rmcp::dso3 a(element);
	util_cout("a = ", a);

	util_cout("+a = ", +a);
	util_cout("-a = ", -a);

	rmcp::dso3 b;
	b = a;
	util_cout("b = ", b);

	rmcp::dvector v(3, element);

	b = v;
	util_cout("b = ", b);

	b = element;
	util_cout("b = ", b);

	a += b;
	util_cout("a = ", a);
	a -= b;
	util_cout("a = ", a);

	a *= 3.0;
	util_cout("a = ", a);

	a /= 3.0;
	util_cout("a = ", a);

	util_cout(" a + b = ", a+b);
	util_cout(" a - b = ", a-b);

	util_cout(" a * 2.0 = ", a * 2.0);
	util_cout(" a / 2.0 = ", a / 2.0);

	a = b;
	util_cout("a = ", a);
	util_cout("b = ", b);

	std::cout << (a == b) << std::endl;
}

void dso3_member_test(void)
{
	std::cout << " === dso3 member ===" << std::endl;

	rmcp::dso3 a;
	a.set(1.0, 2.0, 3.0);
	util_cout("a = ", a);
	double element[] = { 4.0, 5.0, 6.0};
	a.set(element);
	util_cout("a = ", a);

	std::cout << "norm of a = " << a.norm() << std::endl;
}

void dso3_friend_test(void)
{
	std::cout << " === dso3 friend ===" << std::endl;

	rmcp::dso3 a;
	a.set(1.0, 2.0, 3.0);
	util_cout("a = ", a);
	rmcp::dSO3 R = exp(a);
	util_cout("R = ", R);

	a.normalize();
	R = exp(a, 2.1);
	util_cout("R = ", R);

}	

void dSO3_constructor_test(void)
{
	std::cout << " === dSO3 constructor ===" << std::endl;

	rmcp::dSO3 R;
	util_cout("R = ", R);

	double element[] = { 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0 };
	rmcp::dSO3 M(element);
	util_cout("M = ", M);

	rmcp::dSO3 Y(M);
	util_cout("Y = ", Y);

	rmcp::dSO3 W(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	util_cout("W = ", W);

}
void dSO3_operator_test(void)
{
	std::cout << " === dSO3 operator ===" << std::endl;

	rmcp::dSO3 W(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	rmcp::dSO3 R = W;
	util_cout("R = ", R);
	util_cout("~R = ", ~R);

	R *= W;
	util_cout("R = ", R);

	rmcp::dSO3 M = R * W;
	util_cout("M = ", M);


	rmcp::dvector v(3);
	v[0] = 1.2;
	v[1] = 2.2;
	v[2] = 3.3;

	rmcp::dvector new_v = R * v;
	util_cout("new_v = ", new_v);


}
void dSO3_member_test(void)
{
	std::cout << " === dSO3 member ===" << std::endl;

}

void dSO3_inv_euler_test(void)
{
	rmcp::dSO3 R1, R2;
	double a, b, c, y;

	rmcp::dSO3::EULER euler;

	for (int i = 0; i <= 360; i++) {
		y = i * rmcp::PI / 180.0;
		R1.rot_y(y);
		euler = R1.inv_euler(a, b, c);
		if (euler == rmcp::dSO3::ZYX) {
			R2.euler_zyx(a, b, c);
		}
		else if (euler == rmcp::dSO3::ZYZ) {
			R2.euler_zyz(a, b, c);
		}

		std::cout << "i = " << i << " deg ," ;
		std::cout << "euler = " << euler << std::endl;
		std::cout << "R1 = " << R2;
		std::cout << "R1^(T) * R2 = " << R1 % R2;
		util_pause();
	}
}
void dSO3_friend_test(void)
{
	std::cout << " === dSO3 friend ===" << std::endl;
/*
	rmcp::dSO3 W(0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
	rmcp::dSO3 R = rmcp::inv(W);
	util_cout("R = ", R);

	rmcp::dso3 w = rmcp::log(R);
	util_cout("w = ", w);
	rmcp::dangles angle = { rmcp::PI/4.0, rmcp::PI/7.0, rmcp::PI/3.0 };

	R = rmcp::euler_zyx(angle);
	util_cout("R = ", R);

	rmcp::dangles angle_r = rmcp::inv_euler_zyx(R);

	std::cout << "alpha = " << angle_r.alpha << " beta = " << angle_r.beta << " gamma = " << angle_r.gamma << std::endl;

	R = rmcp::euler_zyz(angle);
	util_cout("R = ", R);

	angle_r = rmcp::inv_euler_zyz(R);

	std::cout << "alpha = " << angle_r.alpha << " beta = " << angle_r.beta << " gamma = " << angle_r.gamma << std::endl;

	util_cout("log R = ", rmcp::log(R));
	util_cout("exp(log R) = ", rmcp::exp(rmcp::log(R)));
	*/
}

void dso3_test(void)
{
	dso3_constructor_test();
	util_pause();
	dso3_operator_test();
	util_pause();
	dso3_member_test();
	util_pause();
	dso3_friend_test();
	util_pause();
}

void dSO3_test(void)
{
/*	dSO3_constructor_test();
	util_pause();
	dSO3_operator_test();
	util_pause();
	dSO3_member_test();
	util_pause();
	dSO3_friend_test();
	util_pause();*/

	dSO3_inv_euler_test();
	util_pause();
}
