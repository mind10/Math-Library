
#include "rmcp/dse3.h"
#include "utility.h"

#include <iostream>

void dse3_constructor_test(void)
{
	std::cout << " === dse3 constructor ===" << std::endl;

	rmcp::dse3 S;
	util_cout("S = ", S);

	rmcp::dse3 L(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
	util_cout("L = ", L);

	rmcp::dse3 Q(L);
	util_cout("Q = ", Q);

	double element[] = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
	rmcp::dvector w(3, element);
	rmcp::dvector v(3, element+3);

	rmcp::dse3 U(w, v);
	util_cout("U = ", U);

	rmcp::dse3 Y(element);
	util_cout("Y = ", Y);

}

void dse3_operator_test(void)
{
	std::cout << " === dse3 operator ===" << std::endl;

	double element[] = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6};
	rmcp::dse3 Y(element);
	util_cout("Y = ", Y);
	util_cout("+Y = ", +Y);
	util_cout("-Y = ", -Y);

	rmcp::dse3 S = Y;
	util_cout("S = ", S);

	S += Y;
	util_cout("S = ", S);

	S -= Y;
	util_cout("S = ", S);

	Y *= 2.0;
	util_cout("Y = ", Y);

	Y /= 2.0;
	util_cout("Y = ", Y);

	util_cout("S + Y = ", S + Y);
	util_cout("S - Y = ", S - Y);

	util_cout("Y * 2.0 = ", Y * 2.0);
	util_cout("2.0 * Y = ", 2.0 * Y);
	util_cout("Y / 2.0 = ", Y / 2.0);

	std::cout << S * Y << std::endl;

}
void dse3_member_test(void)
{
	std::cout << " === dse3 member ===" << std::endl;
	rmcp::dse3 A;
	A.set(1.1, 2.2, 3.3, 4.4, 5.5, 6.6);
	util_cout("A = ", A);

}
void dse3_friend_test(void)
{
	std::cout << " === dse3 friend ===" << std::endl;

}
void dSE3_constructor_test(void)
{
	std::cout << " === dSE3 constructor ===" << std::endl;

}

void dSE3_operator_test(void)
{
	std::cout << " === dSE3 operator ===" << std::endl;

	rmcp::dSE3 T(1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 12.12, 13.13, 14.14);
	rmcp::dvector a(0.1, 1.1, 2.2);

	util_cout("T = ", T);
	util_cout("a = ", a);

	util_cout("T % a = T^(-1) * a ", T % a);

	util_cout("T * a = T * a ", T * a);

}
void dSE3_member_test(void)
{
	std::cout << " === dSE3 member ===" << std::endl;

}
void dSE3_friend_test(void)
{
	std::cout << " === dSE3 friend ===" << std::endl;

}


void dse3_test(void)
{
	dse3_constructor_test();
	util_pause();
	dse3_operator_test();
	util_pause();
	dse3_member_test();
	util_pause();
	dse3_friend_test();
	util_pause();
}

void dSE3_test(void)
{
//	dSE3_constructor_test();
//	util_pause();
	dSE3_operator_test();
//	util_pause();
//	dSE3_member_test();
//	util_pause();
//	dSE3_friend_test();
//	util_pause();
	

//	rmcp::dSE3 tt(0.5, 1.1, 2.2, 5.4, 5.5, 6.6, 8.8, 9.9, 10.10, 12.12, 13.13, 14.14);
//	rmcp::dse3 aa(rmcp::dvector(0.1, 1.1, 2.2), rmcp::dvector(3.3, 4.4, 5.5));

//	std::cout << "aa" <<  tt * aa << std::endl;

}

