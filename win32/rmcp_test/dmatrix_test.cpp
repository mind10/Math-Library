#include "rmcp/dmatrix.h"
#include "rmcp/dvector.h"
#include "rmcp/dgqmatrix.h"
#include "utility.h"

void dmatrix_constructor_test(void)
{
	double element[] = { 1.3, 4.5, 6.4, 3.4, 5.6, 7.8, 4.5, 6.7,
		5.4, 3.4, 3.4, 5.6, 7.1, 8.7, 6.7, 2.9 };

	/* the elements of dvector 'a' are garbage */
	rmcp::dmatrix a(4, 4);
	util_cout("a = ", a);

	rmcp::dmatrix b(3, 4, 1.0);
	util_cout("b = ", b);

	rmcp::dmatrix c(4, 3, element);
	util_cout("c = ", c);

	rmcp::dmatrix d(c);
	util_cout("d = ", d);

//	d.inv();
	
	util_cout("d inv =", d);
	
	rmcp::dmatrix f(3,3, 1.0);

	util_cout("f = ", f) ;

	rmcp::dmatrix e[3];

}

void dmatrix_operator_test(void)
{
	std::cout << " ========== dmatrix operator ==========" << std::endl;

	std::cout << " --- access the part of a matrix --- " << std::endl;
	double element[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3};
	rmcp::dmatrix a(3, 3, element);
	util_cout("a = ", a);
	a(0, 0) = 5.0;
	std::cout << ">>> a(0, 0) = 5.0;  /* first row, first column */" << std::endl;
	util_cout("a = ", a);

	double a12 = a(1, 2);
	std::cout << ">>> double a12 = a(1, 2);  /* second row, third column */" << std::endl;
	std::cout << "a12 = " << a12 << std::endl << std::endl;
	a[2] = 7.0;
	std::cout << ">>> a[2] = 7.0; /* third element, in column-by-column order */" << std::endl;
	util_cout("a = ", a);

	double a6 = a[6];
	std::cout << ">>> double a6 = a[6];  /* seventh element, in column-by-column order */" << std::endl;
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
	rmcp::dmatrix k = a;
	util_cout("k = ", k);
	double data[] = { 1.1, 2.1, 3.1, 1.2, 2.2, 3.2, 1.3, 2.3, 3.3};
	k = data;
	util_cout("k = ", k);

	util_pause();

	std::cout << " ---- '+', '-', '*', '/', operator ---" << std::endl;
	rmcp::dmatrix b(a);
	util_cout("a = ", a);
	util_cout("b = ", b);
	util_cout("a + b = ", a + b);
	util_cout("a - b = ", a - b);
	util_cout("a * 2.0 = ", a * 2.0);
	util_cout("2.0 * a = ", 2.0 * a);
	util_cout("a / 2.0 = ", a / 2.0);
	util_cout("a * b = ", a * b);

	rmcp::dvector c(3);
	c[0] = 1.0;
	c[1] = 2.0;
	c[2] = 3.0;
	util_cout("c = ", c);
	util_cout("a * c = ", a * c);

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
	a *= b;
	util_cout("a = ", a);

	util_pause();

	std::cout << " --- '==' operator --- " << std::endl;
	std::cout << "a == b ? " << (a==b) << std::endl;

}

void dmatrix_member_test(void)
{
	std::cout << " === dmatrix member === " << std::endl;

	double element[] = { -1.3, 4.5, 1.4, 3.4, 5.6, 7.8, 4.5, 6.7,
		5.4, 7.4, 3.0, 105.6, 1.1, 8.7, 6.7, -3052.9 };

	rmcp::dmatrix m(3, 4, element);
	util_cout("m = ", m);
	std::cout << "the size of m = " << m.size() << std::endl << std::endl;
	std::cout << "the row of m = " << m.row() << std::endl << std::endl;
	std::cout << "the column of m = " << m.column() << std::endl << std::endl;

	double *ptr = m.element();
	std::cout << ptr[0] << std::endl << std::endl;

	m.set(1.0);
	util_cout("m = ", m);

	m.set(element);
	util_cout("m = ", m);

	double contents[12];
	m.get(contents);

	std::cout << contents[4] << std::endl << std::endl;

	rmcp::dvector a(1);
	a = m.get_column(1); /* return 2nd column vector */
	util_cout("a = ", a);

	m.set_column(a, 2); /* insert 'a' vector to third column of 'm' matrix */
	util_cout("m = ", m);

	a = m.get_row(1); /* return 2nd row vector */
	util_cout("a = ", a);

	m.set_row(a, 2); /* insert 'a' vector to third row of 'm' matrix */
	util_cout("m = ", m);

	rmcp::dmatrix p(1,1);

	p = m.get_partial(1, 2);  /* 'p' matrix has the elements of 'm' matrix from second row to last row and from third column to last column */
	util_cout("p = ", p);

	p = m.get_partial(0, 1, 2, 2); /* 'p' matrix has the element of 'm' matrix from fist row to third row and from second column to third column */ 
	util_cout("p = ", p);

	rmcp::dmatrix f(2, 2);
	f[0] = -1.0;
	f[1] = -2.0;
	f[2] = -3.0;
	f[3] = -4.0;

	util_cout("f = ", f);
	m.set_partial(f, 1, 2); /* insert a matrix 'f' to a matrix m on (1, 2) */
	util_cout("m = ", m);

	rmcp::dvector v(100.1, 100.2, 100.3);
	m.set_partial(v, 0, 3);
	util_cout("m = ", m);

	//m.resize(3, 3); /* after resize, the matrix m is fulled with garbage */
	//util_cout("m = ", m);

	//m.resize(3, 3, 2.0); /* after resize, the matrix m is fulled with 2.0 */
	//util_cout("m = ", m);

	//rmcp::dmatrix n(4, 3, element);
	//util_cout("n = ", n);
	//n.swap_column(0, 2);
	//util_cout("n = ", n);
	//n.swap_row(0, 1);
	//util_cout("n = ", n);

	//std::cout << " Is matrix 'n' zero matrix? " << n.zero() << std::endl << std::endl;
	//
	util_pause();

	rmcp::dmatrix d(4, 4, element);
	util_cout("d = ", d);

	std::cout << "the determinant of matrix 'd' = " << d.det() << std::endl;

	std::cout << "the norm of matrix 'd' = " << d.norm() << std::endl;

	d.normalize();
	util_cout("d = ", d);
	std::cout << "the trace of matrix 'd' = " << d.trace() << std::endl;


	d.transpose();
	util_cout("transposed matrix 'd' = " , d);

	util_cout(" state  1 =", d);
	d.det();
	util_cout(" state 2 = ", d);
	d.inv();
	util_cout("inversed matrix 'd' = " , d);

	
	//std::cout << "How to print dmatrix" << std::endl;
	//std::cout << d << std::endl;

	rmcp::dgqmatrix A(2);

	A(0,0)= 1.0;
	A(0,1)= 2.0;
	A(1,0)= 3.0;
	A(1,1)= 4.0;



//	double element_3[] = { 1,7.2,4.0,1,5,-3.3,13,5,71,11,5,20.2,3,9,3,4,3,-0.22,15,9};

//	double element_4[] = { 2.5, 0.5, 2.2, 1.9, 3.1, 2.3, 2, 1, 1.5, 1.1, 2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9};

//	rmcp::dvector avr(3);
//	avr = temp.mean();
//	util_cout(" mean = ", avr);

//	rmcp::dmatrix aaa(4,5, element_3);
//	util_cout(" aa= ", aaa);
//	rmcp::dmatrix ccc(4,5, 0.0);
//	rmcp::dmatrix ddd(4,5, 0.0);
//	rmcp::dvector eee(5, 0.0);
//	ccc = aaa.cov(eee, ddd);				// ccc = covariance, ddd= data_adjust, aaa= data

//	util_cout(" Data_Adjust = ", ddd);
//	util_cout(" cov = ", ccc);

//	std::cout << " === Principal Components Analysis === " << std::endl;

//	util_cout(" Data= ", aaa);

//	rmcp::dvector avr(5,0.0);
//	rmcp::dmatrix d_ad(4,5, 0.0);
//	rmcp::dmatrix pcs(6,5, 0.0);
//	rmcp::dvector per(5, 0.0);
//	rmcp::dvector eig_val(5,0.0);

//	aaa.pca(avr, d_ad,eig_val, pcs,per);

//	util_cout(" mean = ", avr);
//	util_cout(" Data_Adjust = ", d_ad);
//	util_cout(" EigenValue = ", eig_val);
//	util_cout(" PCS = ", pcs);
//	util_cout(" Percent = ", per);

//	rmcp::dmatrix final_data(4,5,0.0);

//	final_data =d_ad * pcs;

//	util_cout("final data = ", final_data);


//	rmcp::dmatrix original_data(5,4,0.0);
//	original_data = (pcs*final_data)+avr; 
/*
	rmcp::dmatrix aaa(10, 2, element_4);

	rmcp::dvector avr(2, 0.0);
	rmcp::dmatrix pcs(2, 2);
	rmcp::dmatrix data_adjusted(10, 2, 0.0);
	rmcp::dvector eigen_value(2, 0.0);
	rmcp::dvector per(2, 0.0);
	aaa.pca(avr, data_adjusted, eigen_value, pcs, per);

	std::cout << aaa << std::endl;
	std::cout << pcs << std::endl;
*/
}

void dmatrix_friend_test(void)
{
	std::cout << " === dmatrix friend === " << std::endl;

	double element[] = { 1.3, 4.5, 6.4, 3.4, 5.6, 7.8, 4.5, 6.7,
		5.4, 3.4, 3.4, 5.6, 7.1, 8.7, 6.7, 2.9 };

	
	rmcp::dmatrix a(3, 4, element);
	rmcp::dmatrix b(3, 4, element);
	rmcp::dmatrix c(4, 4);



	rmcp::multiply_AtB(c, a, b); /* c = transpose(a) * b */

	util_cout("a = ", a);
	util_cout("b = ", b);
	util_cout("c = ", c);

	rmcp::multiply_ABt(c, a, b); /* c = a * transpose(b) */
	util_cout("a = ", a);
	util_cout("b = ", b);
	util_cout("c = ", c);

	rmcp::dmatrix d(4, 3, element);
	rmcp::multiply_AB(c, a, d); /* c = a * d */
	util_cout("a = ", a);
	util_cout("d = ", d);
	util_cout("c = ", c);

	std::cout << "the norm of matrix 'd' = " << rmcp::norm(d) << std::endl;

	util_pause();

	d = rmcp::zeros(5, 5);
	util_cout("d = ", d);

	d = rmcp::eye(5, 5);
	util_cout("d = ", d);

	d = rmcp::ones(3, 3);
	util_cout("d = ", d);
	std::cout << "the trace of matrix 'd' = " << rmcp::trace(d) << std::endl;

	d = rmcp::rand(4, 5);
	util_cout("d = ", d);

	util_cout("transpose matrix of d = ", rmcp::transpose(d));

	util_cout("diagonal vector of d = ", rmcp::diagonal(d));

}

void dmatrix_eig_test()
{
	double element2[] = { 2, 11, 9, 5, 7, 9, 7, 7, 8, 30, 34, 4, 90, 12, 3, 70};

	rmcp::dmatrix m8(4, 4, element2);
	util_cout(" m8 = ", m8);
	rmcp::dmatrix eigenvector(4, 4);
	rmcp::dmatrix eigenvalue(4, 4);

	m8.eig(eigenvector, eigenvalue);
	util_cout(" eigenvector of m8 = ", eigenvector);
	util_cout(" eigenvalue of m8 = ", eigenvalue);

}
void dmatrix_mean_test()
{
	double element0[] = { 9, 15, 25, 14, 10, 18, 0, 16, 5, 19, 16, 20};

	rmcp::dmatrix m0(3, 4, element0);
	util_cout(" m0 = ", m0);
	util_cout(" mean of column = ", m0.mean());
	util_cout(" mean of row = ", m0.mean(1));
}

void dmatrix_var_test()
{
	double element0[] = { 9, 15, 25, 14, 10, 18, 0, 16, 5, 19, 16, 20}; // 12 ea
	double element1[] = { 39, 56, 93, 61, 50, 75, 32, 85, 42, 70, 66, 80}; // 12 ea

	rmcp::dvector v0(12, element0);
	rmcp::dvector v1(12, element1);

	rmcp::dmatrix m(12, 2);
	m.set_column(v0, 0);
	m.set_column(v1, 1);

	util_cout("m = ", m);
	util_cout("variance (normal, each column is a variable) of m = ", m.var());
	util_cout("variance (divide by n, each column is a variable) of m = ", m.var(0, 1));
	util_cout("variance (normal, each column is a variable) of m = ", m.var(1));
	util_cout("variance (divide by n, each column is a variable) of m = ", m.var(1, 1));
}

void dmatrix_cov_test()
{
	double element0[] = { 9, 15, 25, 14, 10, 18, 0, 16, 5, 19, 16, 20}; // 12 ea
	double element1[] = { 39, 56, 93, 61, 50, 75, 32, 85, 42, 70, 66, 80}; // 12 ea

	rmcp::dvector v0(12, element0);
	rmcp::dvector v1(12, element1);

	rmcp::dmatrix m(12, 2);
	m.set_column(v0, 0);
	m.set_column(v1, 1);

	util_cout("m = ", m);
	util_cout("covariance (normal, each column is a variable) of m = ", m.cov());
	util_cout("covariance (divide by n, each column is a variable) of m = ", m.cov(1));
}

void dmatrix_pca_test()
{

	rmcp::dmatrix pca_pcs[2];

	rmcp::dmatrix data(10, 2);
	data[0] = 2.5;
	data[1] =0.5;
	data[2] = 2.2;
	data[3] = 1.9;
	data[4] = 3.1;
	data[5] = 2.3;
	data[6] = 2;
	data[7] = 1;
	data[8] = 1.5;
	data[9] = 1.1;
	data[10] = 2.4;
	data[11] = 0.7;
	data[12] = 2.9;
	data[13] = 2.2;
	data[14] = 3.0;
	data[15] = 2.7;
	data[16] = 1.6;
	data[17] = 1.1;
	data[18] = 1.6;
	data[19] = 0.9;
	rmcp::dmatrix cov_(2,2);
	cov_ = data.cov(0);
	rmcp::dmatrix eig_vec(2,2);
	rmcp::dmatrix eig_val(2,2);
	rmcp::dvector eigval(2);
	cov_.eig(eig_vec, eig_val);
	eigval.diagonal(eig_val);
	
	util_cout(" data mean = ", data.mean(0));
	util_cout(" data cov = ", cov_);
	util_cout(" data eig_vec = ", eig_vec);
	util_cout(" data eig_val = ", eig_val);
	util_cout(" data eigval = " , eigval);

	double dd[] = {1, 3.4, 5, 4, 7, 1.5, 2.7, 6.7, 3.5, 2.8, 4, 7, 1.4, 8.4, 9.4, 9.2, 1.7, 6.3, 7.3, 8.3, 8.6, 3.5, 7.2, 8.2, 4};
	rmcp::dmatrix test_data(5,5,dd);

	rmcp::dmatrix test_eigval(5,5);
	rmcp::dmatrix test_eigvec(5,5);
	rmcp::dvector aaa(5);
	test_data.eig(test_eigvec, test_eigval);
	aaa.diagonal(test_eigval);
	unsigned int *key = new unsigned int[6];
	aaa.sort(key, 1);		// sort of eigenvalue
	rmcp::dvector key_(6);
	for(int i=0; i<6; i++) key_[i] = key[i];
	util_cout("test_eigvec = ", test_eigvec);
	util_cout("test_eigval = ", test_eigval);
	util_cout("test_sort_eigval =", key_);

#if 0	
	double data0[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	double data1[4] = { 1,2,0,3};

	rmcp::dvector v(4, data1);
	rmcp::dmatrix m(4,4,data0);
	rmcp::dmatrix m_(4,4, 0.0);
	unsigned int *key;
	key = new unsigned int[4];

	util_cout(" v = ", v);
	util_cout(" m = ", m);

	v.sort(key, 1);

	util_cout(" sorted v = ", v);

	rmcp::dvector temp;
	for (int i = 0; i < 4; i++) {
		temp = m.get_column(key[i]);
		//		PC.set_column(temp, key[i]);
		m_.set_column(temp, i);
	}

	util_cout("m_ = ", m_);

	delete key;
#endif
}
void dmatrix_test()
{
//	dmatrix_constructor_test();
//	util_pause();
//	dmatrix_operator_test();
//	util_pause();
//	dmatrix_member_test();
//	util_pause();
//	dmatrix_eig_test();
//	util_pause();
//	dmatrix_mean_test();
//	util_pause();
//	dmatrix_var_test();
//	util_pause();
//	dmatrix_cov_test();
//	util_pause();
	dmatrix_pca_test();
	util_pause();
//	dmatrix_friend_test();
//	util_pause();
}



