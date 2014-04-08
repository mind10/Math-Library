
#include "rmcp/dvector.h"
#include "rmcp/dmatrix.h"
#include "utility.h"

void dvector_constructor_test(void)
{
	std::cout << " === dvector constructor === " << std::endl;

	double element[] = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7 };

	/* the elements of dvector 'a' are garbage */
	rmcp::dvector a(4);
	util_cout("a = ", a);

	rmcp::dvector b(4, 1.0);
	util_cout("b = ", b);

	rmcp::dvector c(3, element);
	util_cout("c = ", c);

	rmcp::dvector d(c);	/* d = c */
	util_cout("d = ", d);

	/* if the row of a matrix is 1 or the column of a matrix 1,
	 * the matrix can be converted to a vector 
	 */

	rmcp::dmatrix m1(4, 1, element);
	rmcp::dvector e(m1);
	util_cout("e = ", e);

	rmcp::dmatrix m2(4, 1, element);
	rmcp::dvector f(m2);
	util_cout("f = ", f);

}

void dvector_operator_test(void)
{
	std::cout << " === dvector operator === " << std::endl;

	double element[] = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7 };

	std::cout << " --- access the part of a vector --- " << std::endl;

	rmcp::dvector a(4, 1.0);
	util_cout("a = ", a);

	a(1) = 3.4;
	util_cout("a = ", a);

	double a1 = a(1);
	std::cout << "a1 = " << a1 << std::endl << std::endl;

	a[2] = 2.1;
	util_cout("a = ", a);
	double a2 = a[2];
	std::cout << "a2 = " << a2 << std::endl << std::endl;
	std::cout << " On vector, () operator is the same as [] operator " << std::endl;

	util_pause();

	std::cout << " --- '+', '-' unary operator --- " << std::endl;
	util_cout("a = ", a);
	util_cout("+a = ", +a);
	util_cout("-a = ", -a);

	util_pause();

	std::cout << " --- '=' operator --- " << std::endl;
	util_cout("a = ", a);

	rmcp::dvector k = a;
	util_cout("k = ", k);

	double data[] = { 1.1, 1.2, 1.3, 1.4 };
	k = data;
	util_cout("k = ", k);

	util_pause();

	std::cout << " ---- '+', '-', '*', '/', '&' operator ---" << std::endl;
	a.resize(3);
	a = element;
	rmcp::dvector b(a);
	util_cout("a = ", a);
	util_cout("b = ", b);
	util_cout("a + b = ", a+b);
	util_cout("a - b = ", a-b);
	util_cout("a * 2.0 = ", a * 2.0);
	util_cout("2.0 * a = ", 2.0 * a);
	util_cout("a / 2.0 = ", a / 2.0);
	std::cout << "a * b (i.e. inner product) = " << a * b << std::endl;
	util_cout("a & b (i.e. cross product) = " , a & b);
	util_cout("a ^ b (i.e. matrix = column vector x row vector) = ", a ^ b);
	rmcp::dmatrix c(3, 3, 1.0);
	util_cout("c = ", c);
	util_cout("a * c (i.e. vector = row vector x matrix) = " , a * c );
	util_cout("c * a (i.e. vector = matrix x column vector) = " , c * a );

	util_pause();

	std::cout << " --- '+=', '-=', '*=', '/=' operator --- " << std::endl;
	util_cout("a = ", a);
	util_cout("b = ", b);
	a += b;
	util_cout("a = ", a);
	a -= b;
	util_cout("a = ", a);
	a *= 2.0;
	util_cout("a = ", a);
	a /= 2.0;
	util_cout("a = ", a);

	util_pause();
	std::cout << " --- '==' operator --- " << std::endl;
	std::cout << "a == b ? " << (a==b) << std::endl;
}

void dvector_memeber_test(void)
{
	std::cout << " === dvector member === " << std::endl;

	double element[] = { 1.3, 4.5, 6.4, 3.4, -5.6, -7.8, 4.5, -6.7, 5.4, 3.4, -3.4, 5.6, 7.1, 8.7, 6.7, 2.9 };

	rmcp::dvector v1(10, element);
	util_cout("v1 = ", v1);
	std::cout << v1.size() << std::endl << std::endl;

	double *ptr = v1.element();
	std::cout << ptr[0] << std::endl << std::endl;

	v1.set(1.0);
	util_cout("v1 = ", v1);

	v1.set(element);
	util_cout("v1 = ", v1);

	double contents[10];
	v1.get(contents);

	std::cout << contents[4] << std::endl << std::endl;

	v1.resize(8); /* after resize, the vector v1 is fulled with garbage */
	util_cout("v1 = ", v1);

	v1.resize(10, 2.0); /* after resize, the vector v1 is fulled with 2.0 */
	util_cout("v1 = ", v1);

	std::cout << " Is v1 zero vector? " << v1.zero() << std::endl << std::endl;

	rmcp::dvector v2(10, element);
	util_cout("v2 = ", v2);

	std::cout << "norm of v2 = " << v2.norm() << std::endl;
	std::cout << "norm of v2 = " << v2.normalize() << std::endl;
	util_cout("normalized v2 = ", v2);
	std::cout << "norm of normalized v1 = " << v2.norm() << std::endl;

	std::cout << "absolute max of v2 is " << v2.amax() << " where " << v2.amax_index() <<"th." << std::endl;
	std::cout << "absolute min of v2 is " << v2.amin() << " where " << v2.amin_index() <<"th." << std::endl;

	rmcp::dmatrix m(3, 4, element);
	rmcp::dvector g(2);
	g.diagonal(m);

	util_cout("m = ", m);
	util_cout("g (diagonal vector of m) = ", g);

	std::cout << "How to print dvector" << std::endl;
	std::cout << v2 << std::endl;


	std::cout << "How to sort the dvector" << std::endl;
	
	rmcp::dvector temp(10, element);
	util_cout("data = ", temp);
	
	rmcp::dvector masrt(temp);
	rmcp::dvector misrt(temp);
	rmcp::dvector seq1(10,0.0);
	rmcp::dvector seq2(10, 0.0);
	
//	masrt.max_sort(seq1);
//	misrt.min_sort(seq2);
	
	util_cout("sequence in decreasing order = ", seq1);
	util_cout("sort in decreasing order = ", masrt );
	util_cout("sequence in increasing order = ", seq2);
	util_cout("sort in increasing order = ", misrt );

}

void dvector_friend_test(void)
{
	std::cout << " === dvector friend === " << std::endl;

	double element[] = { 1.3, 4.5, 6.4, 3.4, -5.6, -7.8, 4.5, 6.7, 5.4, 3.4, 3.4, 5.6, 7.1, 8.7, 6.7, 2.9 };

	rmcp::dvector a(10, element);
	util_cout("a = ", a);

	std::cout << "norm of v2 = " << rmcp::norm(a) << std::endl;

	rmcp::dvector b(10, element);
	util_cout("a = ", a);

	std::cout << "inner product of a and b = " << rmcp::inner(a, b) << std::endl;

	b = rmcp::zeros(4); /* to make 'b' 4-dimensional zero vector */
	util_cout("b = ", b);
	b = rmcp::ones(4); /* to make 'b' 4-dimensional one vector */
	util_cout("b = ", b);
	b = rmcp::rand(4); 
	/* to make 'b' 4-dimensional random vector, -1 <= b <= 1 */
	util_cout("b = ", b);

	rmcp::dvector t(3, element);
	rmcp::dvector u(3, element+4);
	util_cout("t = ", t);
	util_cout("u = ", u);

	rmcp::dvector p(3);
	p = rmcp::cross(t, u);

	util_cout("p (cross(t, u)) = ", p);

	rmcp::dmatrix m(3, 3, element);
	util_cout("m = ", m);

	std::cout << "t * m * u = " << rmcp::quadratic(t, m, u) << std::endl;

}

void dvector_mean_test(void)
{
	double v0[] = { 0, 8, 12, 20};
	rmcp::dvector v(4, v0);

	util_cout("v = ", v);
	std::cout << "mean of v = " << v.mean() << std::endl;
}

void dvector_var_std_test(void)
{
	double v0[] = { 0, 8, 12, 20};
	rmcp::dvector v(4, v0);

	util_cout("v = ", v);
	std::cout << "variance (normal) of v = " << v.var() << std::endl;
	std::cout << "standard deviation (normal) of v = " << v.std() << std::endl;

	std::cout << "variance (divide by n) of v = " << v.var(1) << std::endl;
	std::cout << "standard deviation (divide by n) of v = " << v.std(1) << std::endl;

}

void dvector_cov_test(void)
{
	double element0[] = { 9, 15, 25, 14, 10, 18, 0, 16, 5, 19, 16, 20}; // 12 ea
	double element1[] = { 39, 56, 93, 61, 50, 75, 32, 85, 42, 70, 66, 80}; // 12 ea

	rmcp::dvector v0(12, element0);
	rmcp::dvector v1(12, element1);

	util_cout("v0 = ", v0);
	std::cout << "mean of v0 = " << v0.mean() << std::endl;

	util_cout("v1 = ", v1);
	std::cout << "mean of v1 = " << v1.mean() << std::endl;


	std::cout << "cov (normal) of v0, v1 = " << rmcp::cov(v0, v1) << std::endl; 
	std::cout << "cov (divide by n) of v0, v1 = " << rmcp::cov(v0, v1, 1) << std::endl; 

}

void dvector_sort_test(void)
{
	double element0[] = { 9, -15, 25, 14, -10, 18, 0, 16, 9, 20}; // 11 ea

	rmcp::dvector v(10, element0);

	util_cout("v = ", v);
	v.sort(0);
	util_cout("v (sort numbers in increasing order) = ", v);

	v.set(element0);
	util_cout("v = ", v);
	v.sort(1);
	util_cout("v (sort numbers in decreasing order) = ", v);

	unsigned int key[10];

	v.set(element0);
	util_cout("v = ", v);
	v.sort(key, 0);
	util_cout("v (sort numbers in increasing order) = ", v);

	for (int i = 0; i < 10; i ++) {
		std::cout << key[i] << " ";
	}
	std::cout << std::endl;

	v.set(element0);
	util_cout("v = ", v);
	v.sort(key, 1);
	util_cout("v (sort numbers in decreasing order) = ", v);

	for (int i = 0; i < 10; i ++) {
		std::cout << key[i] << " ";
	}
	std::cout << std::endl;
}

void dvector_test(void)
{
//	dvector_constructor_test();
//	util_pause();
//	dvector_operator_test();
//	util_pause();
//	dvector_memeber_test();
//	util_pause();
//	dvector_friend_test();
//	util_pause();

	dvector_mean_test();
	util_pause();

	dvector_var_std_test();
	util_pause();

	dvector_cov_test();
	util_pause();

	dvector_sort_test();
	util_pause();
}


