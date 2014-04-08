#include "rmcp/drobotics.h"
#include "utility.h"

#include <iostream>

/*
 * test function for rmcp::dvec3, rmcp::dmat3, rmcp::dhmat class
 */

void dvec3_constructor_test(void)
{
	std::cout << " ========== dvec3 constructor ==========" << std::endl;
	
	std::cout << " --- default constructor --- " << std::endl;
	rmcp::dvec3 a;
	util_cout("a = ", a);
	std::cout << " --- constructor with one number --- " << std::endl;
	rmcp::dvec3 b(1.0);
	std::cout << ">>> rmcp::dvec3 b(1.0);" << std::endl;
	std::cout << "b = " << b << std::endl;

	std::cout << " --- constructor with three numbers --- " << std::endl;
	rmcp::dvec3 c(1.0, 2.0, 3.0);
	std::cout << ">>> rmcp::dvec3 c(1.0, 2.0, 3.0);" << std::endl;
	std::cout << "c = " << c << std::endl;

	std::cout << " --- constructor with array --- " << std::endl;
	double element[] = { 1.1, 1.2, 1.3 };
	rmcp::dvec3 d(element);
	std::cout << ">>> double element[] = { 1.1, 1.2, 1.3 }; rmcp::dvec3 d(element);  " << std::endl;
	std::cout << "d = " << d << std::endl;

	std::cout << " --- copy constructor --- " << std::endl;
	rmcp::dvec3 e(d);
	std::cout << ">>> rmcp::dvec3 e(d);" << std::endl;
	std::cout << "e = " << e << std::endl;
}

void dvec3_operator_test(void)
{
	std::cout << " ========== dvec3 operator ==========" << std::endl;
	
	std::cout << " --- access the part of a vector --- " << std::endl;
	rmcp::dvec3 a(1.0);
	std::cout << "a = " << a << std::endl;
	a(1) = 3.4;
	std::cout << ">>> a(1) = 3.4;  " << std::endl;
	std::cout << "a = " << a << std::endl;
	double a1 = a(1);
	std::cout << ">>> double a1 = a(1);  " << std::endl;
	std::cout << "a1 = " <<a1 << std::endl << std::endl;
	a[2] = 2.1;
	std::cout << ">>> a[2] = 2.1;  " << std::endl;
	std::cout << "a = " << a << std::endl;
	double a2 = a[2];
	std::cout << ">>> double a2 = a[2];  " << std::endl;
	std::cout << "a2 = " << a2 << std::endl << std::endl;
	std::cout << " On vector, () operator is the same as [] operator " << std::endl;

	util_pause();

	std::cout << " --- '+', '-' unary operator --- " << std::endl;
	std::cout << "a = " << a << std::endl;
	std::cout << "+a = " << +a << std::endl;
	std::cout << "-a = " << -a << std::endl;

	util_pause();

	std::cout << " --- '=' operator --- " << std::endl;
	std::cout << "a = " << a << std::endl;
	rmcp::dvec3 k = a;
	std::cout << ">>> rmcp::dvec3 k = a;" << std::endl;
	std::cout << "k = " << k << std::endl;
	k = 2.345;
	std::cout << ">>> k = 2.345; " << std::endl;
	std::cout << "k = " << k << std::endl;
	double data[] = { 7.6, 5.4, 9.2 };
	k = data;
	std::cout << ">>> double data[] = { 7.6, 5.4, 9.2 }; k = data;" << std::endl;
	std::cout << "k = " << k << std::endl;

	util_pause();

	std::cout << " ---- '+', '-', '*', '/', '&' operator ---" << std::endl;
	rmcp::dvec3 b(a);
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << a << std::endl;
	std::cout << "a + b = " << a + b << std::endl;
	std::cout << "a - b = " << a - b << std::endl;
	std::cout << "a * 2.0 = " << a * 2.0 << std::endl;
	std::cout << "2.0 * a = " << 2.0 * b << std::endl;
	std::cout << "a / 2.0 = " << a / 2.0 << std::endl;
	std::cout << "a * b (i.e. inner product) = " << a * b << std::endl << std::endl;
	std::cout << "a & b (i.e. cross product) = " << (a & b) << std::endl;
	std::cout << "a ^ b (i.e. matrix = column vector x row vector) = " << (a ^ b) << std::endl;
	rmcp::dmat3 c(3.4);
	std::cout << "c = " << c << std::endl;
	std::cout << "a * c (i.e. vector = row vector x matrix) = " << a * c << std::endl;
	
	std::cout << " --- '+=', '-=', '*=', '/=' operator --- " << std::endl;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << a << std::endl;
	a += b;
	std::cout << ">>> a += b;" << std::endl;
	std::cout << "a = " << a << std::endl;
	a -= b;
	std::cout << ">>> a -= b;" << std::endl;
	std::cout << "a = " << a << std::endl;
	a *= 2.0;
	std::cout << ">>> a *= 2.0;" << std::endl;
	std::cout << "a = " << a << std::endl;
	a /= 2.0;
	std::cout << ">>> a /= 2.0;" << std::endl;
	std::cout << "a = " << a << std::endl;

}

void dvec3_friend_test(void)
{
	std::cout << " ========== dvec3 friend function ==========" << std::endl;

	std::cout << " --- cross(...) : cross product --- " << std::endl;
	rmcp::dvec3 x;
	x[0] = 1.0; x[1] = 0.0; x[2] = 0.0;
	rmcp::dvec3 y;
	y[0] = 0.0; y[1] = 1.0; y[2] = 0.0;
	std::cout << "x = " << x << std::endl;
	std::cout << "y = " << y << std::endl;
	rmcp::dvec3 z = cross(x, y);
	std::cout << ">>> rmcp::dvec3 z = cross(x, y);" << std::endl;
	std::cout << "z = " << z << std::endl;
	rmcp::dmat3 m;
	m.set_column(x, 0);
	m.set_column(y, 1);
	m.set_column(z, 2);
	std::cout << "m = " << m << std::endl;
	rmcp::dmat3 n = cross(z, m);
	std::cout << ">>> rmcp::dmat3 n = cross(z, m);" << std::endl;
	std::cout << "n = " << n << std::endl;

	util_pause();

	std::cout << " --- cross(...) and inv_cross: get skew-symmetric matrix --- " << std::endl;

	rmcp::dvec3 w(1.2, 2.3, 3.4);
	n = cross(w);
	std::cout << "w = " << w << std::endl;
	std::cout << ">>> n = cross(w);" << std::endl;
	std::cout << "n = " << n << std::endl;
	w = inv_cross(n);	
	std::cout << ">>> w = inv_cross(n);" << std::endl;
	std::cout << "w = " << w << std::endl;
}

void dmat3_constructor_test(void)
{
	std::cout << " ========== dmat3 constructor ==========" << std::endl;

	std::cout << " --- default constructor --- " << std::endl;
	rmcp::dmat3 a;
	std::cout << ">>> rmcp::dmat3 a;" << std::endl;
	std::cout << "a = " << a << std::endl;

	std::cout << " --- constructor with one number --- " << std::endl;
	rmcp::dmat3 b(1.0);
	std::cout << ">>> rmcp::dmat3 b(1.0);" << std::endl;
	std::cout << "b = " << b << std::endl;

	std::cout << " --- constructor with three numbers --- " << std::endl;
	rmcp::dmat3 c(1.0, 2.0, 3.0);
	std::cout << ">>> rmcp::dmat3 c(1.0, 2.0, 3.0);" << std::endl;
	std::cout << "c = " << c << std::endl;

	std::cout << " --- constructor with array --- " << std::endl;
	double element[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3};
	rmcp::dmat3 d(element);
	std::cout << ">>> double element[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3}; rmcp::dmat3 d(element);  " << std::endl;
	std::cout << "d = " << d << std::endl;

	std::cout << " --- copy constructor --- " << std::endl;
	rmcp::dmat3 e(d);
	std::cout << ">>> rmcp::dmat3 e(d);" << std::endl;
	std::cout << "e = " << e << std::endl;

}

void dmat3_operator_test(void)
{
	std::cout << " ========== dmat3 operator ==========" << std::endl;

	std::cout << " --- access the part of a matrix --- " << std::endl;
	double element[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3};
	rmcp::dmat3 a(element);
	util_cout("a = ", a);
	a(0, 0) = 5.0;
	std::cout << ">>> a(0, 0) = 5.0;  /* first row, first column */" << std::endl;
	util_cout("a = ", a);

	double a12 = a(1, 2);
	std::cout << ">>> double a12 = a(1, 2);  /* second row, third column */" << std::endl;
	std::cout << "a12 = " << a12 << std::endl << std::endl;
	a[2] = 7.0;
	std::cout << ">>> a[2] = 7.0; /* third element, in column by column order */" << std::endl;
	util_cout("a = ", a);

	double a6 = a[6];
	std::cout << ">>> double a6 = a[6];  /* seventh element, in column by column order */" << std::endl;
	std::cout << "a6 = " << a6 << std::endl << std::endl;
	std::cout << " On matrix, () operator is not the same as [] operator " << std::endl;

	util_pause();

	std::cout << " --- '+', '-', '~' unary operator --- " << std::endl;
	util_cout("a = ", a);
	util_cout("+a = ", +a);
	util_cout("-a = ", -a);
	util_cout("~a (i.e. transpose) = ", ~a);

	util_pause();

	std::cout << " --- '=' operator --- " << std::endl;
	util_cout("a = ", a);
	rmcp::dmat3 k = a;
	std::cout << ">>> rmcp::dmat3 k = a;" << std::endl;
	util_cout("k = ", k);
	k = 2.345;
	std::cout << ">>> k = 2.345; " << std::endl;
	std::cout << "k = " << k << std::endl;
 	double data[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3};
	k = data;
	std::cout << ">>> double data[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3}; k = data;" << std::endl;
	util_cout("k = ", k);

	util_pause();

	std::cout << " ---- '+', '-', '*', '/', operator ---" << std::endl;
	rmcp::dmat3 b(a);
	util_cout("a = ", a);
	util_cout("b = ", b);
	util_cout("a + b = ", a + b);
	util_cout("a - b = ", a - b);
	util_cout("a * 2.0 = ", a * 2.0);
	util_cout("2.0 * a = ", 2.0 * b);
	util_cout("a / 2.0 = ", a / 2.0);
	util_cout("a * b = ", a * b);

	rmcp::dvec3 c(1.0, 2.0, 3.0);
	util_cout("c = ", c);
	util_cout("a * c = ", a * c);

	std::cout << " --- '+=', '-=', '*=', '/=' operator --- " << std::endl;
	util_cout("a = ", a);
	util_cout("b = ", b);

	a += b;
	std::cout << ">>> a += b;" << std::endl;
	util_cout("a = ", a);
	a -= b;
	std::cout << ">>> a -= b;" << std::endl;
	util_cout("a = ", a);
	a *= 2.0;
	std::cout << ">>> a *= 2.0;" << std::endl;
	util_cout("a = ", a);
	a /= 2.0;
	std::cout << ">>> a /= 2.0;" << std::endl;
	util_cout("a = ", a);
	a *= b;
	std::cout << ">>> a *= b;" << std::endl;
	util_cout("a = ", a);

}

void dmat3_member_test(void)
{
	std::cout << " ========== dmat3 member functions ==========" << std::endl;

	std::cout << " --- double det() : get a determinant of matrix --- " << std::endl;
	double data[] = { 1.9, 2.1, 3.1, 1.2, 2.5, 3.2, 1.3, 2.3, 3.3};
	rmcp::dmat3 a(data);
	std::cout << "a = " << a << std::endl;
	std::cout << "a.det() = " << a.det() << std::endl;

	util_pause();

	std::cout << " --- void transpose() : get transposed matrix --- " << std::endl;
	a.transpose();
	std::cout << ">>> a.transpose();" << std::endl;
	std::cout << "a = " << a << std::endl;

	util_pause();

	std::cout << " --- void inv() : get inverse matrix --- " << std::endl;
	a.inv();
	std::cout << ">>> a.inv();" << std::endl;
	std::cout << "a = " << a << std::endl;

}

void dmat3_friend_test(void)
{
	std::cout << " ========== dmat3 friend functions ==========" << std::endl;

	std::cout << " --- dmat3 rotate(...) : rotation, e.g. Roll-Pitch-Yaw, Euler Angles --- " << std::endl;

	rmcp::dvec3 angle(rmcp::PI/4.0, rmcp::PI/4.0, 0.0);
	std::cout << "angle = " << angle << std::endl;

	rmcp::dmat3 R = rotate(angle, rmcp::ROT_RPY);
	std::cout << ">>> rmcp::dmat3 R = rotate(angle, rmcp::ROT_RPY);" << std::endl;
	std::cout << "R = " << R << std::endl;

	R = rotate(angle, rmcp::ROT_ZYZ);
	std::cout << ">>> R = rotate(angle, rmcp::ROT_ZYZ);" << std::endl;
	std::cout << "R = " << R << std::endl;

	R = rotate(angle, rmcp::ROT_ZXZ);
	std::cout << ">>> R = rotate(angle, rmcp::ROT_ZXZ);" << std::endl;
	std::cout << "R = " << R << std::endl;

	R = rotate(rmcp::PI/4.0, rmcp::PI/4.0, 0.0, rmcp::ROT_RPY);
	std::cout << ">>> R = rotate(rmcp::PI/4.0, rmcp::PI/4.0, 0.0, rmcp::ROT_RPY);" << std::endl;
	std::cout << "R = " << R << std::endl;
	R = rotate(rmcp::PI/4.0, rmcp::PI/4.0, 0.0, rmcp::ROT_ZYZ);
	std::cout << ">>> 	R = rotate(rmcp::PI/4.0, rmcp::PI/4.0, 0.0, rmcp::ROT_ZYZ);" << std::endl;
	std::cout << "R = " << R << std::endl;
	R = rotate(rmcp::PI/4.0, rmcp::PI/4.0, 0.0, rmcp::ROT_ZXZ);
	std::cout << ">>> 	R = rotate(rmcp::PI/4.0, rmcp::PI/4.0, 0.0, rmcp::ROT_ZXZ);" << std::endl;
	std::cout << "R = " << R << std::endl;

	rmcp::dvec3 axis(1.0, 1.0, 0.0);
	std::cout << "axis = " << axis << std::endl;
	R = rotate(rmcp::PI/4.0, axis);

	std::cout << ">>> R = rotate(rmcp::PI/4.0, axis);" << std::endl;
	std::cout << "R = " << R << std::endl;

	std::cout << " --- inv_rotate(...) : get angles from a rotation matrix --- " << std::endl;
	std::cout << "R = " << R << std::endl;
	angle = inv_rotate(R, rmcp::ROT_RPY);
	std::cout << ">>> angle = inv_rotate(R, rmcp::ROT_RPY); " << std::endl;
	std::cout << "angle = " << angle << std::endl;
	angle = inv_rotate(R, rmcp::ROT_ZYZ);
	std::cout << ">>> angle = inv_rotate(R, rmcp::ROT_ZYZ); " << std::endl;
	std::cout << "angle = " << angle << std::endl;
	angle = inv_rotate(R, rmcp::ROT_ZXZ);
	std::cout << ">>> angle = inv_rotate(R, rmcp::ROT_ZXZ); " << std::endl;
	std::cout << "angle = " << angle << std::endl;

	double ang1, ang2, ang3;
	inv_rotate(ang1, ang2, ang3, R, rmcp::ROT_RPY);
	std::cout << ">>> double ang1, ang2, ang3;" << std::endl;
	std::cout << ">>> inv_rotate(ang1, ang2, ang3, R, rmcp::ROT_RPY);" << std::endl;
	std::cout << "ang1 = " << ang1 << '\t' << "ang2 = " << ang2 << '\t' << "ang3 = " << ang3 << std::endl;

	inv_rotate(ang1, ang2, ang3, R, rmcp::ROT_ZYZ);
	std::cout << ">>> inv_rotate(ang1, ang2, ang3, R, rmcp::ROT_ZYZ);" << std::endl;
	std::cout << "ang1 = " << ang1 << '\t' << "ang2 = " << ang2 << '\t' << "ang3 = " << ang3 << std::endl;

	inv_rotate(ang1, ang2, ang3, R, rmcp::ROT_ZXZ);
	std::cout << ">>> inv_rotate(ang1, ang2, ang3, R, rmcp::ROT_ZXZ);" << std::endl;
	std::cout << "ang1 = " << ang1 << '\t' << "ang2 = " << ang2 << '\t' << "ang3 = " << ang3 << std::endl;

	inv_rotate(ang1, axis, R);
	std::cout << ">>> inv_rotate(ang1, axis, R);" << std::endl;
	std::cout << "ang1 = " << ang1 << std::endl;
	std::cout << "axis = " << axis << std::endl;

}

void dhmat_constructor_test(void)
{
	std::cout << " ========== dhmat constructor ==========" << std::endl;

	std::cout << " --- default constructor --- " << std::endl;
	rmcp::dhmat a;
	std::cout << ">>> rmcp::dhmat a;" << std::endl;
	std::cout << "a = " << a << std::endl;

	std::cout << " --- constructor with array --- " << std::endl;
	double element[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3, 1.4, 2.4, 3.4};
	rmcp::dhmat d(element);
	std::cout << ">>> double element[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3, 1.4, 2.4, 3.4};" << std::endl;
	std::cout << ">>> rmcp::dhmat d(element);  " << std::endl;
	std::cout << "d = " << d << std::endl;

	std::cout << " --- constructor with dmat3 and dvec3 --- " << std::endl;
	rmcp::dmat3 R(element);
	rmcp::dvec3 v(element+9);
	std::cout << "R = " << R << std::endl;
	std::cout << "v = " << v << std::endl;
	rmcp::dhmat H(R, v);
	std::cout << ">>> rmcp::dhmat H(R, v);" << std::endl;
	std::cout << "H = " << H << std::endl;

	std::cout << " --- copy constructor --- " << std::endl;
	rmcp::dhmat e(d);
	std::cout << ">>> rmcp::dhmat e(d);" << std::endl;
	std::cout << "e = " << e << std::endl;

}

void dhmat_operator_test(void)
{
	std::cout << " ========== dhmat operator ==========" << std::endl;

	std::cout << " --- '=' operator --- " << std::endl;
	double element[] = { 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0, 0, -1.0, 1.0, 2.0, 3.0};
	rmcp::dhmat a = element;

	std::cout << ">>> double element[] = { 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0, 0, -1.0, 1.0, 2.0, 3.0};" << std::endl;
	std::cout << ">>> rmcp::dhmat a = element;" << std::endl;
	std::cout << "a = " << a << std::endl;

	rmcp::dhmat b = a;
	std::cout << ">>> rmcp::dhmat b = a;" << std::endl;
	std::cout << "b = " << b << std::endl;

	util_pause();

	std::cout << " --- '*=', '/=' operator --- " << std::endl;

	a *= b;
	std::cout << ">>> a *= b;" << std::endl;
	std::cout << "a = " << a << std::endl;

	a /= b;
	std::cout << ">>> a /= b; (i.e. a = a x b^(-1))" << std::endl;
	std::cout << "a = " << a << std::endl;

	std::cout << " --- '*', '/' operator --- " << std::endl;
	rmcp::dhmat c;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;

	c = a * b;
	std::cout << ">>> c = a * b;" << std::endl;
	std::cout << "c = " << c << std::endl;

	c = a / b;
	std::cout << ">>> c = a / b (i.e. c = a x b^(-1));" << std::endl;
	std::cout << "c = " << c << std::endl;
    
}

void dhmat_member_test(void)
{
	std::cout << " ========== dhmat member function ==========" << std::endl;

	std::cout << " --- set(...) ---" << std::endl;
	double element[] = { 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0, 1.0, 0.0, 1.0, 2.0, 3.0};
	rmcp::dhmat a;
	a.set(element);

	std::cout << ">>> double element[] = { 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0, 1.0, 0.0, 1.0, 2.0, 3.0};" << std::endl;
	std::cout << ">>> rmcp::dhmat a;" << std::endl;
	std::cout << ">>> a.set(element);" << std::endl;
	std::cout << "a = " << a << std::endl;

	rmcp::dhmat c;
	c.set(a);

	std::cout << ">>> rmcp::dhmat c; c.set(a);" << std::endl;
	std::cout << "c = " << c << std::endl;

	rmcp::dmat3 R(element);
	R.inv();
	rmcp::dvec3 v(element+9);
	v *=2;

	std::cout << "R = " << R << std::endl;
	std::cout << "v = " << v << std::endl;
	rmcp::dhmat b;
	b.set(R, v);
	std::cout << ">>> rmcp::dhmat b; b.set(R, v);" << std::endl;
	std::cout << "b = " << b << std::endl;

	std::cout << "c = " << c << std::endl;
	c.set(R);
	std::cout << ">>> c.set(R);" << std::endl;
	std::cout << "c = " << c << std::endl;

	c.set(v);
	std::cout << ">>> c.set(v);" << std::endl;
	std::cout << "c = " << c << std::endl;

	util_pause();

	std::cout << " --- inv() : inverse matrix ---" << std::endl;

	std::cout << "a = " << a << std::endl;
	a.inv();
	std::cout << ">>> a.inv();" << std::endl;
	std::cout << "a = " << a << std::endl;
}

void dhmat_friend_test(void)
{
	std::cout << " ========== dhmat friend function ==========" << std::endl;
	std::cout << " --- inv() : inverse matrix ---" << std::endl;

	double element[] = { 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0, 1.0, 0.0, 1.0, 2.0, 3.0};
	rmcp::dhmat a(element);
	std::cout << "a = " << a << std::endl;
	rmcp::dhmat b = inv(a);
	std::cout << ">>> rmcp::dhmat b = inv(a);" << std::endl;
	std::cout << "b = " << b << std::endl;

}

void dvec3_test(void)
{
	dvec3_constructor_test();
	util_pause();
	dvec3_operator_test();
	util_pause();
	dvec3_friend_test();
	util_pause();
}

void dmat3_test(void)
{
	dmat3_constructor_test();
	util_pause();
	dmat3_operator_test();
	util_pause();
	dmat3_member_test();
	util_pause();
	dmat3_friend_test();
	util_pause();
}

void
eulerangle(void)
{
	rmcp::dvec3 angle(107.799* rmcp::PI / 180.0, 23.186 * rmcp::PI / 180.0, -136.837 * rmcp::PI / 180.0);
	rmcp::dvec3 v(0.0, 183.404, 0.0);
	rmcp::dmat3 R(1.0, 1.0, 1.0);

	rmcp::dhmat m(R, v);
	rmcp::dmat3 ro = rmcp::rotate(angle, rmcp::ROT_RPY);

	rmcp::dvec3 p(25.914, 43.31, 1444.618);

	rmcp::dhmat a(ro, p);

	rmcp::dvec3 t = a * v;

	util_cout("a = ", a);
	util_cout("target ", t);



}
void dhmat_test(void)
{
	/*
	dhmat_constructor_test();
	util_pause();
	dhmat_operator_test();
	util_pause();
	dhmat_member_test();
	util_pause();
	dhmat_friend_test();
	util_pause();
	*/

	eulerangle();
	util_pause();
}
