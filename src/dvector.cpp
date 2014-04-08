
#include "rmcp/dvector.h"
#include "rmcp/dmatrix.h"
#include "rmcp/dse3.h"
#include "rmcp/qsort_wkey.h"

#include <cmath>
#include <cassert>

using namespace rmcp;

dvector::dvector()	/* 3-dimensional vector with zero values */
{
	size_ = 3;
	element_ = new double[size_];
	element_[0] = 0.0;
	element_[1] = 0.0;
	element_[2] = 0.0;
}

dvector::dvector(const int size)
{
	assert (size > 0 &&
		   	"dvector::dvector(const int)");

	size_ = size;
	element_ = new double[size_];
}

dvector::dvector(const int size, const double element)
{
	assert (size > 0 &&
		   	"dvector::dvector(const int, const double)");

	size_ = size;
	element_ = new double[size_];

	for (int i = 0; i < size; i++) {
		element_[i] = element;
	}
}

dvector::dvector(const int size, const double element[])
{ 
	assert (size > 0 &&
		   	"dvector::dvector(const int, const double[])");

	size_ = size;
	element_ = new double[size_];
	cblas_dcopy(size_, element, 1, element_, 1);
}

dvector::dvector(const double element0, const double element1, const double element2)
{
	size_ = 3;
	element_ = new double[size_];
	element_[0] = element0;
	element_[1] = element1;
	element_[2] = element2;
}

dvector::dvector(const dvector& v)
{
	size_ = v.size_;
	element_ = new double[size_];
	cblas_dcopy(size_, v.element_, 1, element_, 1);
}

/*
dvector::dvector(const rmcp::dmatrix &m)
{
	assert( ((m.row_ == 1)||(m.column_ == 1)) &&
		"dvector::dvector(const rmcp::dmatrix &)" );

	size_ = m.size_;
	element_ = new double[size_];
	cblas_dcopy(size_, m.element_, 1, element_, 1);
}
*/

dvector::~dvector() 
{
	assert (element_ &&
			"dvector::~dvector()");
	
	delete[] element_;
}

double &
dvector::operator [] (const int i)
{
	assert (i >= 0 &&
			i < size_ &&
		   	"double &dvector::operator [] (const int)");

   	return element_[i];
}

const double &
dvector::operator [] (int i) const
{
	assert (i >= 0 &&
			i < size_ &&
		   	"const double &dvector::operator [] (int) const");

   	return element_[i] ;
}

double & 
dvector::operator () (const int i)
{
	assert (i >= 0 &&
		i < size_ &&
		"double &dvector::operator(const int)");

	return element_[i];
}

const double &
dvector::operator () (int i) const
{
	assert (i >= 0 &&
		i < size_ &&
		"const double &dvector::operator(int) const");

	return element_[i];
}

const dvector &
dvector::operator + (void) const
{
   	return *this;
}

dvector 
dvector::operator - (void) const	
{
	dvector vr(size_);
	
	for (int i = 0; i < vr.size_; i++) {
		vr.element_[i] = - element_[i];
	}

   	return vr;
} 

dvector &
dvector::operator = (const dvector &v)
{
	assert (element_ &&
			v.element_ &&
			"const dvector &dvector::operator = (const dvector &)");

	if (size_ != v.size_) {
		size_ = v.size_;
		delete element_;
		element_ = new double[size_];
	}

	cblas_dcopy(size_, v.element_, 1, element_, 1);
   	return *this;
}

dvector &
dvector::operator = (const double element[])
{
	cblas_dcopy(size_, element, 1, element_, 1);
   	return *this;
}

dvector &
dvector::operator += (const dvector& v)
{
	assert (size_ == v.size_ &&
			"const dvector &dvector::operator += (const dvector&)");

	for (int i = 0; i < size_; i++) {
		element_[i] += v.element_[i];
	}

	return *this;
}

dvector &
dvector::operator -= (const dvector& v)
{
	assert (size_ == v.size_ &&
			"const dvector &dvector::operator -= (const dvector&)");

	for (int i = 0; i < size_; i++) {
		element_[i] -= v.element_[i];
	}

   	return *this;
}

dvector & 
dvector::operator *= (double s)
{
	cblas_dscal(size_, s, element_, 1);
   	return *this;
}

dvector & 
dvector::operator /= (double s)
{
	assert ( s != 0.0 &&
			"const dvector &dvector::operator /= (double)");

	cblas_dscal(size_, 1.0 / s, element_, 1);
   	return *this;
}

dvector
dvector::operator * (double s) const
{
	dvector vr(*this);
	cblas_dscal(vr.size_, s, vr.element_, 1);
   	return vr;
}

dvector
dvector::operator / (double s) const
{
	assert ( s != 0.0 &&
			"dvector dvector::operator / (double)");

	dvector vr(*this);
	cblas_dscal(vr.size_, 1.0 / s, vr.element_, 1);
   	return vr;
}

dvector
dvector::operator + (const dvector& v) const
{
	assert (size_ == v.size_ &&
			"dvector dvector::operator + (const dvector &) const");

	dvector vr(*this);

	for (int i = 0; i < vr.size_; i++) {
		vr.element_[i] += v.element_[i];
	}

   	return vr;
}

dvector
dvector::operator - (const dvector& v) const
{
	assert (size_ == v.size_ &&
			"dvector dvector::operator - (const dvector&) const");

	dvector vr(*this);
	for (int i = 0; i < vr.size_; i++) {
		vr.element_[i] -= v.element_[i];
	}

   	return vr;
}

double
dvector::operator * (const dvector& v) const
{
	assert (size_ == v.size_ &&
			"double dvector::operator * (const dvector&) const");

	return cblas_ddot(size_, element_, 1, v.element_, 1);
}

dvector
dvector::operator * (const dmatrix &m) const
{
	assert (size_ == m.row_ &&
		"dvector dvector::operator * (const dmatrix &) const");

	dvector vr(m.column_);
	cblas_dgemv(CblasColMajor, CblasTrans, m.row_, m.column_, 1.0, m.element_, m.row_, element_, 1, 0.0, vr.element_, 1);
	return vr;
}

dmatrix
dvector::operator ^ (const dvector &v) const  // matrix = col vec x row vec
{
	assert (size_ == v.size_ &&
		"dmatrix dvector::operator ^ (const dvector &) const");

	dmatrix mr(v.size_, v.size_);

	for (int i = 0; i < size_; i++) {
		for (int j = 0; j < size_; j++) {
			mr.element_[i + size_ * j] = element_[i] * v.element_[j];
		}
	}

	return mr;
}

dvector
dvector::operator & (const dvector &v) const // cross product
{
	assert (size_ == 3 &&
		v.size_ == 3 &&
		"dvector dvector::operator ^ (const dvector &) const");

	dvector vr(3);
	vr.element_[0] = element_[1] * v.element_[2] - element_[2] * v.element_[1];
	vr.element_[1] = element_[2] * v.element_[0] - element_[0] * v.element_[2];
	vr.element_[2] = element_[0] * v.element_[1] - element_[1] * v.element_[0];
	return vr;
}

bool
dvector::operator == (const dvector& v) const
{
	assert (size_ == v.size_ &&
			"bool dvector::operator == (const dvector&) const");

	static bool re;
	re = true;
	for (int i = 0; i < size_; i++) {
		re = re && (element_[i] == v.element_[i]);
	}

	return re;
}

void
dvector::set(const double element[])
{
	cblas_dcopy(size_, element, 1, element_, 1);
	return;
}

void
dvector::set(const double element)
{
	for (int i = 0; i < size_; i++) {
		element_[i] = element;
	}

	return;
}

void
dvector::set(const double element0, const double element1, const double element2) /* for 3-dimensional vector */
{
	element_[0] = element0;
	element_[1] = element1;
	element_[2] = element2;
}

void
dvector::get(double element[])
{
	cblas_dcopy(size_, element_, 1, element, 1);
	return;
}

void
dvector::resize(const int size)
{
	assert (size > 0 &&
		   	"dvector::resize(const int)");

	if (size_ != size) {
		delete[] element_;
		size_ = size;
		element_ = new double[size_];
	}
}

void
dvector::resize(const int size, const double element)
{
	assert (size > 0 &&
		   	"dvector::resize(const int, const double)");

	if (size_ != size) {
		delete[] element_;
		size_ = size;
		element_ = new double[size_];
	}

	for (int i = 0; i < size; i++) {
		element_[i] = element;
	}
}

void
dvector::resize(const int size, double *element)
{
	assert (size > 0 &&
		"dvector::resize(const int, double *)");

	size_ = size;
	delete element_;
	element_ = element;
}

bool
dvector::zero(void) const
{
	static bool re;
	re = true;
	for (int i = 0; i < size_; i++) {
		re = re && (element_[i] == 0.0);
	}
	return re;
}

double
dvector::norm(void) const
{
	return cblas_dnrm2(size_, element_, 1);
}

double
dvector::square_sum(void) const
{
	double square_sum = 0.0;

	for (int i = 0; i < size_; i++) {
		square_sum += (element_[i] * element_[i]);
	}

	return square_sum;
}

/* if it's norm is less than VECTOR_NORM_EPS, error! */
double
dvector::normalize(void)
{
	double norm = this->norm();

	assert (std::fabs(norm) > rmcp::VECTOR_NORM_EPS &&
		"double dvector::normalize(void)");

	cblas_dscal(size_, 1.0 / norm, element_, 1);
	return norm;
}

double
dvector::amax(void) const
{
	return element_[cblas_idamax(size_, element_, 1)];
}

double
dvector::amin(void) const
{
	return element_[cblas_idamin(size_, element_, 1)];
}

int
dvector::amax_index(void) const
{
	return (int)cblas_idamax(size_, element_, 1);
}

int
dvector::amin_index(void) const
{
	return (int)cblas_idamin(size_, element_, 1);
}

void
dvector::diagonal(rmcp::dmatrix &m)
{
	if (m.row_ <= m.column_) {
		this->resize(m.row_);
	}
	else {
		this->resize(m.column_);
	}

	int j = 0;
	for (int i = 0; i < m.size_; i += (m.row_ + 1)) {
		element_[j] = m.element_[i];
		j++;
	}
}

double
dvector::mean(void) 
{
	assert (size_ > 0 && "double dvector::mean(void)");

	double sum = 0.0;
	for (int i = 0; i < size_; i++) {
		sum += element_[i];
	}

	return sum / size_;
}

double
dvector::var(int n) 
{
	assert (size_ > 0 && "double dvector::var(int)");

	double mean = this->mean();

	double var = 0.0;
	for (int i = 0; i < size_; i++) {
		var += (element_[i] - mean)*(element_[i] - mean);
	}

	switch (n) {
		case 0 :
			if (size_ > 1) var /= (size_-1);
			break;
		case 1 : 
			var /= size_;
			break;
	}

	return var;
}

double
dvector::std(int n) 
{
	assert (size_ > 0 && "double dvector::std(int)");
	
	return sqrt(this->var(n));
}

void 
dvector::sort(unsigned int key[], int dir)
{
	for (unsigned int i = 0; i < (unsigned int)size_; i++) {
		key[i] = i;
	}

	switch (dir) {
		case 0 :
			qsort_wkey(element_, size_, sizeof(double), key, dqsort_cmp_decrease);
			break;
		default :
			qsort_wkey(element_, size_, sizeof(double), key, dqsort_cmp_increase);
			break;
	}
}

//void
//dvector::sort(unsigned int *key, int dir)
//{
//	for (unsigned int i = 0; i < (unsigned int)size_; i++) {
//		key[i] = i;
//	}
//
//	switch (dir) {
//		case 0 :
//			qsort_wkey(element_, size_, sizeof(double), key, dqsort_cmp_decrease);
//			break;
//		default :
//			qsort_wkey(element_, size_, sizeof(double), key, dqsort_cmp_increase);
//			break;
//	}
//}

void 
dvector::sort(int dir)
{
	char id;

	switch (dir) {
		case 0 :
			id = 'I';
			break;
		default :
			id = 'D';
			break;
	}

	int info;

	dlasrt(&id, &size_, element_, &info);
}

/////////////////////////////////////////////////////////////
// friend functions
////////////////////////////////////////////////////////////

std::ostream &
rmcp::operator << (std::ostream& os, const dvector& v)
{
	os.setf(std::ios::fixed);
	os << "vector (dimension : " << v.size_ << ")" << std::endl; 
	os << "[" << std::endl;
	for (int i = 0; i < v.size_; i++) {
		if (v.element_[i] >= 0.0) {
				os << " ";
		}
		os << v.element_[i] << std::endl;
	}
	os << "];" << std::endl;
	return os;
}

dvector
rmcp::operator * (double s, const dvector &v)
{
	dvector vr(v);
	cblas_dscal(vr.size_, s, vr.element_, 1);
	return vr;
}

double
rmcp::norm(const dvector &v)
{
	return v.norm();
}

double
rmcp::inner(const dvector &x, const dvector &y)
{
	assert (x.size_ == y.size_ &&
		"double rmcp::inner(const dvector &, const dvector &)");

	return cblas_ddot(x.size_, x.element_, 1, y.element_, 1);
}

dvector
rmcp::zeros(const int size)
{
	dvector vr(size, 0.0);
	return vr; 
}

dvector
rmcp::rand(const int size)
{
	dvector vr(size);

	for (int i = 0; i < vr.size_; i++) {
		vr.element_[i] = 2.0 * ((double)std::rand() / (double)RAND_MAX - 0.5);
	}

	return vr;	
}

dvector
rmcp::ones(const int size)
{
	dvector vr(size, 1.0);

	return vr; 
}

dvector
rmcp::cross(const dvector &x, const dvector &y)
{
	assert (x.size_ == 3 &&
		y.size_ == 3 &&
		"dvector dvector::cross(const dvector &x, const dvector &y)");

	dvector vr(3);

	vr.element_[0] = x.element_[1] * y.element_[2] - x.element_[2] * y.element_[1];
	vr.element_[1] = x.element_[2] * y.element_[0] - x.element_[0] * y.element_[2];
	vr.element_[2] = x.element_[0] * y.element_[1] - x.element_[1] * y.element_[0];

	return vr;

}

double 
rmcp::cov(dvector &v1, dvector &v2, int n) // the sizes are equal.
{
	assert ( v1.size_ == v2.size_ && "double rmcp::cov(dvector, dvector)");

	double mean1 = v1.mean();
	double mean2 = v2.mean();
	
	double cov = 0.0;
	for (int i = 0; i < v1.size_; i++) {
		cov += (v1.element_[i] - mean1) * (v2.element_[i] - mean2);
	}

	switch (n) {
		default:
		case 0:
			cov /= (v1.size_ - 1);
			break;
		case 1:
			cov /= v1.size_;
			break;
	}
	return cov;
}

/* Lie Group Class(dso3, dSO3, dse3, dSE3) */

dvector &
dvector::operator %= (const dSE3 &T)
{
	static double tmp[3];
	tmp[0] = element_[0];
	tmp[1] = element_[1];
	tmp[2] = element_[2];

	element_[0] = 
		T.element_[0] * (tmp[0] - T.element_[12]) +
		T.element_[1] * (tmp[1] - T.element_[13]) +
		T.element_[2] * (tmp[2] - T.element_[14]);
	element_[1] = 
		T.element_[4] * (tmp[0] - T.element_[12]) + 
		T.element_[5] * (tmp[1] - T.element_[13]) +
		T.element_[6] * (tmp[2] - T.element_[14]);
	element_[2] = 
		T.element_[8] * (tmp[0] - T.element_[12]) + 
		T.element_[9] * (tmp[1] - T.element_[13]) +
		T.element_[10] * (tmp[2] - T.element_[14]);

	return *this;
}

dvector &
dvector::operator *= (const dSE3 &T)
{
	static double tmp[3];
	tmp[0] = element_[0];
	tmp[1] = element_[1];
	tmp[2] = element_[2];

	element_[0] = 
		T.element_[12] +
		T.element_[0] * tmp[0] +
		T.element_[4] * tmp[1] +
		T.element_[8] * tmp[2];
	element_[1] =
		T.element_[13] +
		T.element_[1] * tmp[0] +
		T.element_[5] * tmp[1] +
		T.element_[9] * tmp[2];
	element_[2] = 
		T.element_[14] +
		T.element_[2] * tmp[0] +
		T.element_[6] * tmp[1] +
		T.element_[10] * tmp[2];

	return  *this;
}

void
dvector::rotate(const rmcp::dSE3 &T, const dvector &v)
{
	element_[0] = 
		T.element_[0] * v.element_[0] +
		T.element_[1] * v.element_[1] +
		T.element_[2] * v.element_[2];
	element_[1] = 
		T.element_[4] * v.element_[0] +
		T.element_[5] * v.element_[1] +
		T.element_[6] * v.element_[2];
	element_[2] = 
		T.element_[8] * v.element_[0] +
		T.element_[9] * v.element_[1] +
		T.element_[10] * v.element_[2];

	return;
}

void
dvector::inv_rotate(const rmcp::dSE3 &T, const dvector &v)
{
	element_[0] =
		T.element_[0] * v.element_[0] +
		T.element_[4] * v.element_[1] +
		T.element_[8] * v.element_[2];
	element_[1] = 
		T.element_[1] * v.element_[0] +
		T.element_[5] * v.element_[1] +
		T.element_[9] * v.element_[2];
	element_[2] = 
		T.element_[2] * v.element_[0] +
		T.element_[6] * v.element_[1] +
		T.element_[10] * v.element_[2];

	return;
}

void
dvector::rotate(const rmcp::dSE3 &T)
{
	static double tmp[3];
	tmp[0] = element_[0];
	tmp[1] = element_[1];
	tmp[2] = element_[2];

	element_[0] = T.element_[0] * tmp[0] + T.element_[4] * tmp[1] + T.element_[8] * tmp[2];
	element_[1] = T.element_[1] * tmp[0] + T.element_[5] * tmp[1] + T.element_[9] * tmp[2];
	element_[2] = T.element_[2] * tmp[0] + T.element_[6] * tmp[1] + T.element_[10] * tmp[2];

}

void
dvector::inv_rotate(const rmcp::dSE3 &T)
{
	static double tmp[3];
	tmp[0] = element_[0];
	tmp[1] = element_[1];
	tmp[2] = element_[2];

	element_[0] = T.element_[0] * tmp[0] + T.element_[1] * tmp[1] + T.element_[2] * tmp[2];
	element_[1] = T.element_[4] * tmp[0] + T.element_[5] * tmp[1] + T.element_[6] * tmp[2];
	element_[2] = T.element_[8] * tmp[0] + T.element_[9] * tmp[1] + T.element_[10] * tmp[2];

}

dvector
rmcp::rotate(const rmcp::dSE3 &T, const dvector &v)
{
	rmcp::dvector vr;
	vr.rotate(T, v);
	return vr;
}

dvector
rmcp::inv_rotate(const rmcp::dSE3 &T, const dvector &v)
{
	rmcp::dvector vr;
	vr.inv_rotate(T, v);
	return vr;
}
