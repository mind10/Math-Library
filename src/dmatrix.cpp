
#include "rmcp/dmatrix.h"
#include "rmcp/dvector.h"

#include "rmcp/dse3.h"
#include "rmcp/dinertia.h"

#include <cmath>
#include <cassert>

using namespace rmcp;

dmatrix::dmatrix()
{
	row_ = 1;
	column_ = 1;
	size_ = 1;
	element_ = new double[size_];
}

dmatrix::dmatrix(const int row, const int column)
{
	assert (row > 0 &&
			column > 0 &&
			"dmatrix::dmatrix(const int, const int)");

	row_ = row;
	column_ = column;
	size_ = row * column;
	element_ = new double[size_];
}

dmatrix::dmatrix(const int row, const int column, const double element)
{
	assert (row > 0 &&
			column > 0 &&
			"dmatrix::dmatrix(const int, const int, const double)");

	row_ = row;
	column_ = column;
	size_ = row * column;
	element_ = new double[size_];

	for (int i = 0; i < size_; i++) {
		element_[i] = element;
	}
}

dmatrix::dmatrix(const int row, const int column, const double element[])
{
	assert (row > 0 &&
			column > 0 &&
			"dmatrix::dmatrix(const int, const int, const double[])");

	row_ = row;
	column_ = column;
	size_ = row * column;
	element_ = new double[size_];

	cblas_dcopy(size_, element, 1, element_, 1);
}

dmatrix::dmatrix(const dmatrix &m)
{
	assert (m.element_ &&
			"dmatrix::dmatrix(const dmatrix)");

	row_ = m.row_;
	column_ = m.column_;
	size_ = m.size_;
	element_ = new double[size_];

	cblas_dcopy(size_, m.element_, 1, element_, 1);
}

/*
dmatrix::dmatrix(const dvector &v)
{
	row_ = v.size_;
	column_ = 1;
	size_ = v.size_;
	element_ = new double[size_];

	cblas_dcopy(size_, v.element_, 1, element_, 1);
}
*/
dmatrix::~dmatrix()
{
	assert (element_ &&
			"dmatrix::~dmatrix()");

	delete[] element_;
}

double &
dmatrix::operator [] (int i)
{
	assert (i > -1 &&
			i < size_ &&
		   	"double &dmatrix::operator [] (const int)");
	
	return element_[i];
}

const double &
dmatrix::operator [] (int i) const
{
	assert (i > -1 &&
			i < size_ &&
		   	"const double &dmatrix::operator [] (int) const");

	return element_[i];
}

double & 
dmatrix::operator () (const int i, const int j)
{
	assert (i > -1 &&
		   	i < row_ &&
		    j > -1 &&
			j < column_ && 
			"double &dmatrix::operator(const int, const int)");

	return element_[i + j * row_];
}

const double &
dmatrix::operator () (int i, int j) const
{
	assert (i > -1 &&
		   	i < row_ &&
		    j > -1 &&
			j < column_ && 
			"const double &dmatrix::operator(int, int) const");

	return element_[i + j * row_];
}

const dmatrix &
dmatrix::operator + (void) const
{
   	return *this;
}

dmatrix 
dmatrix::operator - (void) const
{
	dmatrix mr(row_, column_);

	for (int i = 0; i < size_; i++) {
		mr.element_[i] = - element_[i];
	}

	return mr; 
}

dmatrix 
dmatrix::operator ~ (void) const
{ 
	dmatrix mr(column_, row_);

	for (int i = 0; i < column_; i++) {
		for (int j = 0; j < row_; j++) {
			mr.element_[i + j * column_] = element_[j + i * row_];
		}
	}

	return mr;
}

dmatrix &
dmatrix::operator = (const dmatrix &m)
{
	assert (element_ &&
			m.element_ &&
			"const dmatrix &dmatrix::operator = (const dmatrix &)");

	if (size_ != m.size_) {
		size_ = m.size_;
		delete element_;
		element_ = new double[size_];
	}

	row_ = m.row_;
	column_ = m.column_;
	cblas_dcopy(size_, m.element_, 1, element_, 1);

	return *this;
}

dmatrix &
dmatrix::operator = (const double element[])
{
	cblas_dcopy(size_, element, 1, element_, 1);
	return *this;
}

dmatrix &
dmatrix::operator += (const dmatrix &m)
{
	assert (m.row_ == row_ &&
		   	m.column_ == column_ &&
		   "const dmatrix &dmatrix::operator += (const dmatrix &)");

	for (int i = 0; i < size_; i++) {
		element_[i] += m.element_[i];
	}

	return *this;
}

dmatrix &
dmatrix::operator -= (const dmatrix &m)
{
	assert (m.row_ == row_ &&
		   	m.column_ == column_ &&
		   "const dmatrix &dmatrix::operator -= (const dmatrix &)");

	for (int i = 0; i < size_; i++) {
		element_[i] -= m.element_[i];
	}

	return *this;
}

dmatrix &
dmatrix::operator *= (double s)
{
	cblas_dscal(size_, s, element_, 1);
	return *this;
}

dmatrix &
dmatrix::operator *= (const dmatrix &m)
{
	assert (column_ == m.row_ &&
		"dmatrix &dmatrix::operator *= (const dmatrix &)");

	static double* mat;
	mat = new double[row_ * m.column_];

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, row_, m.column_, column_, 1.0, element_, row_, m.element_, m.row_, 0.0, mat, row_);

	delete element_;
	element_ = mat;
	column_ = m.column_;

	return *this;
}

dmatrix &
dmatrix::operator /= (double s)
{
	assert (s != 0.0 &&
			"const dmatrix &dmatrix::operator /= (double)");

	cblas_dscal(size_, 1.0 / s, element_, 1);
	return *this;
}

dmatrix 
dmatrix::operator + (const dmatrix &m) const
{
	assert (m.row_ == row_ &&
		   	m.column_ == column_ && 
			"dmatrix dmatrix::operator + (const dmatrix &) const");

	dmatrix mr(*this);

	for (int i = 0; i < size_; i++) {
		mr.element_[i] += m.element_[i];
	}

	return mr;
}

dmatrix 
dmatrix::operator - (const dmatrix &m) const
{
	assert (m.row_ == row_ &&
		   	m.column_ == column_ && 
			"dmatrix dmatrix::operator - (const dmatrix &)");

	dmatrix mr(*this);

	for (int i = 0; i < size_; i++) {
		mr.element_[i] -= m.element_[i];
	}

	return mr;
}

dmatrix 
dmatrix::operator * (const double s) const
{
	dmatrix mr(*this);
	cblas_dscal(mr.size_, s, mr.element_, 1);
	return mr;
}

dmatrix 
dmatrix::operator / (const double s) const
{
	assert (s != 0.0 &&
			"dmatrix dmatrix::operator / (double s) const");

	dmatrix mr(*this);
	cblas_dscal(size_, 1.0 / s, mr.element_, 1);
	return mr;
}

dmatrix 
dmatrix::operator * (const dmatrix &m) const
{
	assert (column_ == m.row_ &&
			"dmatrix dmatrix::operator * (const dmatrix &) const)");

	dmatrix mr(row_, m.column_);
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, mr.row_, mr.column_, column_, 1.0, element_, row_, m.element_, m.row_, 0.0, mr.element_, mr.row_);
	return mr;
}

dvector
dmatrix::operator * (const dvector &v) const
{
	assert (column_ == v.size_ &&
		"dvector dmatrix::operator * (const dvector &) const");

	dvector vr(row_);
	cblas_dgemv(CblasColMajor, CblasNoTrans, row_, column_, 1.0, element_, row_, v.element_, 1, 0.0, vr.element_, 1);
	return vr;
}

/*
dmatrix 
dmatrix::operator | (const dmatrix &m) const
{
	dmatrix re;
	multiply_ABt(re, *this, M);
	return re;	
}

dmatrix 
dmatrix::operator ^ (const dmatrix& M) const
{
	dmatrix re;
	multiply_AtB(re, *this, M);
	return re;
}



dmatrix 
dmatrix::operator % (const dmatrix &M) const
{
	dmatrix x;
	solve_AxeB(*this, x, M);
	return x;
}

dmatrix
dmatrix::operator & (const dmatrix &M) const
{
	dmatrix x;
	solve_AtxeB(*this, x, M);
	return x;
}
*/
bool
dmatrix::operator == (const dmatrix &m) const
{
	assert (row_ == m.row_ &&
			column_ == m.column_ &&
			"bool dmatrix::operator == (const dmatrix &m) const");

	bool re = true;

	for (int i = 0; i < size_; i++) {
		re = re && (m.element_[i] == element_[i]);
	}

	return re;  
}

void
dmatrix::set(const double element)
{
	for (int i = 0; i < size_; i++) {
		element_[i] = element;
	}
	return;
}

void
dmatrix::set(const double element[])
{
	cblas_dcopy(size_, element, 1, element_, 1);
	return;
}

void
dmatrix::set_column(const rmcp::dvector &v, const int idx)
{
	assert (idx >= 0 &&
		idx < column_ &&
		v.size_ == row_ &&
		"void dmatrix::set_column(const rmcp::dvector &, const int)");

	cblas_dcopy(row_, v.element_, 1, element_ + row_ * idx, 1);

}

void
dmatrix::set_row(const rmcp::dvector &v, const int idx)
{
	assert (idx >= 0 &&
		idx < row_ &&
		v.size_ == column_ &&
		"void dmatrix::set_row(const rmcp::dvector &, const int)");

	cblas_dcopy(column_, v.element_, 1, element_ + idx, row_);

}

void
dmatrix::set_partial(const rmcp::dmatrix &m, const int top, const int left)
{
	assert(top >= 0 &&
		left >= 0 &&
		m.row_ <= row_ - top &&
		m.column_ <= column_ - left &&
		"rmcp::dmatrix dmatrix::set_partial(const rmcp::dmatrix &m, const int top, const int left)");
	
	for (int i = 0; i < m.row_; i++) {
		for (int j = 0; j < m.column_; j++) {
			element_[top + i + (j + left) * row_] = m.element_[i + j * m.row_];
		}
	}
}

void
dmatrix::set_partial(const rmcp::dvector &v, const int top, const int left)
{
	assert(top >= 0 &&
		left >= 0 &&
		v.size_ <= row_ - top &&
		1 <= column_ - left &&
		"rmcp::dmatrix dmatrix::set_partial(const rmcp::dvector &v, const int top, const int left)");
	
	for (int i = 0; i < v.size_; i++) {
		element_[top + i + left * row_] = v.element_[i];
	}
}

void
dmatrix::get(double element[])
{
	cblas_dcopy(size_, element_, 1, element, 1);
	return;
}

rmcp::dvector
dmatrix::get_column(const int idx)
{
	assert (idx >= 0 &&
		idx < column_ &&
		"rmcp::dvector dmatrix::get_column(const int)");

	rmcp::dvector vr(row_);
	cblas_dcopy(row_, element_ + row_ * idx, 1, vr.element_, 1);
	return vr;
}

rmcp::dvector
dmatrix::get_row(const int idx)
{
	assert (idx >= 0 &&
		idx < row_ &&
		"rmcp::dvector dmatrix::get_row(const int)");

	rmcp::dvector vr(column_);
	cblas_dcopy(column_, element_ + idx, row_, vr.element_, 1);
	return vr;
}

rmcp::dmatrix
dmatrix::get_partial(const int top, const int left)
{
	assert(top >= 0 &&
		left >= 0 &&
		top < row_ &&
		left < column_ &&
		"rmcp::dmatrix dmatrix::get_partial(const int top, const int left)");

	return this->get_partial(top, left, row_ - 1, column_ - 1);

}

rmcp::dmatrix
dmatrix::get_partial(const int top, const int left, const int bottom, const int right)
{
	assert(top >= 0 &&
		left >= 0 &&
		bottom < row_ &&
		right < column_ &&
		top <= bottom &&
		left <= right &&
		"rmcp::dmatrix dmatrix::get_partial(const int top, const int left, const int bottom, const int right)");

	rmcp::dmatrix mr(bottom - top + 1, right - left + 1);

	for (int i = 0; i < mr.row_; i++) {
		for (int j = 0; j < mr.column_; j++) {
			mr.element_[i + j * mr.row_] = element_[top + i + (left + j) * row_];
		}
	}

	return mr;
}

void
dmatrix::resize(const int row, const int column)
{
	if (row * column != size_) {
		delete[] element_;
		size_ = row * column;
		element_ = new double[size_];
	}
	row_ = row;
	column_ = column;
}

void
dmatrix::resize(const int row, const int column, const double element)
{
	if (row * column != size_) {
		delete[] element_;
		size_ = row * column;
		element_ = new double[size_];
	}
	row_ = row;
	column_ = column;

	for (int i = 0; i < size_; i++) {
		element_[i] = element;
	}
}

void 
dmatrix::swap_column(const int u, const int v)
{
	assert (u >= 0 &&
			u < column_ &&
			v >= 0 &&
			v < column_ &&
			"void dmatrix::swap_column(const int, const int)");
	static double* tmp;
	tmp	= new double[row_];

	/*
	 * 1. u-th column -(copy)-> tmp
	 * 2. v-th column -(copy)-> u-th column
	 * 3. tmp -(copy)-> v-th column
	 */
	cblas_dcopy(row_, element_ + row_ * u, 1, tmp, 1);
	cblas_dcopy(row_, element_ + row_ * v, 1, element_ + row_* u, 1);
	cblas_dcopy(row_, tmp, 1, element_ + row_ * v, 1);
	
	delete tmp;
}

void
dmatrix::swap_row(const int u, const int v)
{
	assert (u >= 0 &&
			u < row_ &&
			v >= 0 &&
			v < row_ &&
			"void dmatrix::swap_row(const int, const int)");
	static double* tmp;
	tmp	= new double[column_];


	/*
	 * 1. u-th row -(copy)-> tmp
	 * 2. v-th row -(copy)-> u-th row
	 * 3. tmp -(copy)-> v-th row
	 */
	cblas_dcopy(column_, element_ + u, row_, tmp, 1);
	cblas_dcopy(column_, element_ + v, row_, element_ + u, row_);
	cblas_dcopy(column_, tmp, 1, element_ + v, row_);

	delete tmp;
}

bool
dmatrix::zero(void) const
{
	static bool re;
	re = true;
	for (int i = 0; i < size_; i++) {
		re = re && (element_[i] == 0.0);
	}
	return re;
}

rmcp::dvector
dmatrix::export_dvector(void)
{
	rmcp::dvector v(size_);

	cblas_dcopy(size_, element_, 1, v.element_, 1);

	return v;
}

double
dmatrix::norm(void) const
{
	return cblas_dnrm2(size_, element_, 1);
}

/* if the matrix's norm is less than MATRIX_NORM_EPS, error! */
double
dmatrix::normalize(void)
{
	static double norm;

	norm = this->norm();

	assert (std::fabs(norm) > rmcp::MATRIX_NORM_EPS &&
		"double dmatrix::normalize(void)");

	cblas_dscal(size_, 1.0 / norm, element_, 1);
	return norm;
}

double
dmatrix::trace(void) const
{
	assert (row_ == column_ &&
		"double dmatrix::trace(void) const");

	static double tr;
	tr = 0.0;
	for (int i = 0; i < size_; i += (row_ + 1)) {
		tr += element_[i];
	}
	return tr;
}

double
dmatrix::det(void)
{
	assert (row_ == column_ &&
		"void dmatrix::det(void)");

	static int *ipiv;
	static int info;
	static double *element;	
	double det = 1.0;

	ipiv = new int[row_];
	element = new double[size_];

	cblas_dcopy(size_, element_, 1, element, 1);
	dgetrf(&row_, &column_, element, &row_, ipiv, &info);

	for (int i = 0; i < row_; i++) {
		if ((i+1) != ipiv[i]) det = -det;
		det *= element[i + i * row_];
	}

	delete[] ipiv;
	delete[] element;

	return det;
}


void
dmatrix::transpose(void)
{
	static double* element;
	element = new double[size_];

	cblas_dcopy(size_, element_, 1, element, 1);

	static int tmp;
	tmp = row_;
	row_ = column_;
	column_ = tmp;

	for (int i = 0; i < row_; i++) {
		for (int j = 0; j < column_; j++) {
			element_[i + j * row_] = element[j + i * column_];
		}
	}

	delete element;
}

void
dmatrix::inv(void) 
{
	assert (row_ == column_ &&
			"void dmatrix::inv(void)");

	static int *ipiv;
	static int info;
	static double *work;

	ipiv = new int[row_];
	work = new double[row_];

	dgetrf(&row_, &column_, element_, &row_, ipiv, &info);
	dgetri(&row_, element_, &row_, ipiv, work, &row_, &info);

 	delete[] ipiv;
	delete[] work;

	return;
}

void
dmatrix::eig(dmatrix &v, dmatrix &d)
{
	static double* element;
	element = new double[size_];

	cblas_dcopy(size_, element_, 1, element, 1);

	int info;
	double work1;		// workspace 
	double *work2;		// workspace
	double *wr;			// real parts of eigenvalue
	double *wi;			// imaginary parts of eigenvalue
	double *vr;			// right eigenvector 
	double *vl;			// left eigenvector 

	wr = new double[row_];
	wi = new double[row_];
	vr = new double[size_];
	vl = new double[size_];

	char JOBVL = 'N';		// left eigenvector on='V',off='N'
	char JOBVR = 'V';		// right eigenvector on='V', off='N'
	int n = row_;			
	int lda = column_;
	int ldvr = row_;
	int ldvl = row_;
	int lwork = -1; // if lwork=-1, then a workspace query is assumed; the routine only calculates the optimal size of the work
					// array, returns this value as the first entry of the work array.

	dgeev(&JOBVL, &JOBVR, &n, element, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &work1, &lwork, &info);	// calculates the optimal size of the work array
	lwork = (int)work1;	// input the optimal value of size of the work array
		
	work2 = new double[lwork];	// renew the optimal size of the work array
	dgeev(&JOBVL, &JOBVR, &n, element, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work2, &lwork, &info);
	v.resize(row_, column_);
	v.set(vr);

	d.resize(row_, column_, 0.0);
	for (int i = 0; i < row_; i++) {
		d.element_[i + row_ * i] = wr[i];
	}

	delete[] work2;
	delete[] wr;
	delete[] wi;
	delete[] vr;
	delete[] vl;
}

rmcp::dvector
dmatrix::mean(int dim) 
{
	if (dim == 0) {
		rmcp::dvector avr(column_, 0.0);
		for (int i = 0; i < column_; i++) {	
			for (int j = 0; j < row_; j++) {
				avr.element_[i] += element_[j + row_ * i];
			}
			avr.element_[i] /= row_;
		}
		return avr;
	}
	else {
		dvector avr(row_, 0.0);
		for (int i = 0; i < row_; i++) {	
			for (int j = 0; j < column_; j++) {
				avr.element_[i] += element_[i + row_ * j];
			}
			avr.element_[i] /= column_;
		}
		return avr;
	}
}

rmcp::dvector
dmatrix::var(int dim, int n)
// dim = 0; column's variance, dim = 1; row's variance
// n = 0 : divided by size-1, n = 1: divided by size
{
	if (dim == 0) {
		rmcp::dvector var(column_, 0.0);
		rmcp::dvector mean(column_, 0.0);
		mean = this->mean(0);

		for (int i = 0; i < column_; i++) {	
			for (int j = 0; j < row_; j++) {
				var.element_[i] += (element_[j + row_ * i] -mean.element_[i])*(element_[j + row_ * i] - mean.element_[i]);
			}
		}

		switch (n) {
			default :
			case 0 :
				for (int i = 0; i < column_; i++) {
					var.element_[i] /= (row_-1);
				}
				break;
			case 1 :
				for (int i = 0; i < column_; i++) {
					var.element_[i] /= row_;
				}
				break;
		}
		return var;
	}
	else {
		dvector var(row_, 0.0);
		rmcp::dvector mean(row_, 0.0);
		mean = this->mean(1);
		for (int i = 0; i < row_; i++) {	
			for (int j = 0; j < column_; j++) {
				var.element_[i] += (element_[i + row_ * j] - mean.element_[i])*(element_[i + row_ * j] - mean.element_[i]);
			}
		}

		switch (n) {
			default :
			case 0 :
				for (int i = 0; i < row_; i++) {
					var.element_[i] /= (column_-1);
				}
				break;
			case 1 :
				for (int i = 0; i < row_; i++) {
					var.element_[i] /= column_;
				}
				break;
		}
		return var;
	}
}

rmcp::dmatrix
dmatrix::cov(int n)	
{
	dvector mean(column_,0.0);
	dmatrix cov(column_, column_);

	static double* element;
	element = new double[size_];
	
	mean = this->mean(0); // column's mean;

	for (int i = 0; i < column_; i++) {
		for (int j = 0; j < row_; j++) {
			element[j+i*row_] = element_[j+i*row_] - mean.element_[i];
		}
	}

	double sum;

	for (int i = 0; i < column_; i++) {
		for (int j = 0; j < i + 1; j++) {
			sum = 0.0;
			for (int k = 0; k < row_; k++) {
				sum += element[k+i*row_] * element[k+j*row_];
			}
			cov.element_[i+j*column_] = sum;
		}
	}

	switch (n) {
		default :
		case 0 :
			for (int i = 0; i < column_; i++) {
				for (int j = 0; j < i + 1; j++) {
					cov.element_[i+j*column_] /= (row_ -1);
				}
			}
			break;
		case 1 :
			for (int i = 0; i < column_; i++) {
				for (int j = 0; j < i + 1; j++) {
					cov.element_[i+j*column_] /= row_;
				}
			}
			break;
	}

	for (int i = 0; i < column_; i++) {
		for (int j = i+1; j < column_; j++) {
			cov.element_[i + j * column_] = cov.element_[j + i * column_];
		}
	}

	delete element;
	return cov;
}

void
dmatrix::pca(dmatrix &PC, rmcp::dvector &mean, rmcp::dvector &eigenvalue)
{

	dmatrix cov(column_,column_);	// covariance
	dmatrix eig_vector(column_, column_);  // eigenvector
	dmatrix eig_value(column_, column_);  // eigenvalue
	dvector temp(column_, 0.0);
	dvector d(column_);
	dvector rank(column_,0.0);

	double eig_sum = 0.0;

	mean = this->mean(0);
	cov = this->cov(0);

	cov.eig(eig_vector, eig_value);

	eigenvalue.diagonal(eig_value);

	unsigned int *key = new unsigned int[d.size_];
	eigenvalue.sort(key, 1);		// sort of eigenvalue
	
	for (int i = 0; i < column_; i++) {
		temp = eig_vector.get_column(key[i]);
		PC.set_column(temp, i);
	}

	delete key;
}


/////////////////////////////////////////////////////////////////////
// friend functions
/////////////////////////////////////////////////////////////////////

std::ostream &
rmcp::operator << (std::ostream &os, const dmatrix &m)
{
	os.setf(std::ios::fixed);
	os << "matrix (row x column : " << m.row_ << " x " << m.column_ <<")"  << std::endl; 
	os << "[" << std::endl;
	for (int i = 0; i < m.row_; i++)	{
		for (int j = 0; j < m.column_; j++) {
			if (m.element_[i + j * m.row_] >= 0.0) {
				os << " ";
			}
			os << m.element_[i + j * m.row_];
			if (j == m.column_ - 1) {
				os << " ;" << std::endl;
			}
			else os << "  ";
		}
	}
	os << "];" << std::endl;
	return os;
}

dmatrix
rmcp::operator * (double s, const dmatrix &m)
{
	dmatrix mr(m);
	cblas_dscal(mr.size_, s, mr.element_, 1);
	return mr;
}

double
rmcp::norm(const dmatrix &m)
{
	return m.norm();
}

void
rmcp::multiply_AB(dmatrix &re, const dmatrix &A, const dmatrix &B)
{
	assert (A.column_ == B.row_ &&
		"void rmcp::multiply_AB(dmatrix &, const dmatrix &, const dmatrix &");

	re.resize(A.row_, B.column_);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, re.row_, re.column_, A.column_, 1.0, A.element_, A.row_, B.element_, B.row_, 0.0, re.element_, re.row_);
}


void
rmcp::multiply_ABt(dmatrix &re, const dmatrix &A, const dmatrix &B)
{
	assert (A.column_ == B.column_ &&
		"void rmcp::multiply_ABt(dmatrix &, const dmatrix &, const dmatrix &");

	re.resize(A.row_, B.row_);

	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, re.row_, re.column_, A.column_, 1.0, A.element_, A.row_, B.element_, B.row_, 0.0, re.element_, re.row_);
}

void
rmcp::multiply_AtB(dmatrix &re, const dmatrix &A, const dmatrix &B)
{
	assert (A.row_ == B.row_ &&
		"void rmcp::multiply_ABt(dmatrix &, const dmatrix &, const dmatrix &");

	re.resize(A.column_, B.column_);

	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, re.row_, re.column_, A.row_, 1.0, A.element_, A.row_, B.element_, B.row_, 0.0, re.element_, re.row_);
}


double
rmcp::quadratic(const dvector &x, const dmatrix &A, const dvector &y)
{
	assert (x.size_ == A.row_ &&
		A.column_ == y.size_ &&
		"double rmcp::quadratic(const dvector &, const dmatrix &, const dvector &)");

	dvector vr(A.row_);

	// vr = A * y
	cblas_dgemv(CblasColMajor, CblasNoTrans, A.row_, A.column_, 1.0, A.element_, A.row_, y.element_, 1, 0.0, vr.element_, 1);
	
	// x * vr
	return cblas_ddot(x.size_, x.element_, 1, vr.element_, 1);
}

dvector
rmcp::diagonal(dmatrix &m)
{
	rmcp::dvector vr(1);
	vr.diagonal(m);
	return vr;
}

double
rmcp::trace(const dmatrix &m)
{
	return m.trace();
}

dmatrix 
rmcp::transpose(const dmatrix &m)
{
	dmatrix mr(m);
	mr.transpose();
	return mr;
}

dmatrix
rmcp::zeros(const int row, const int column)
{
	return rmcp::dmatrix(row, column, 0.0);
}

dmatrix
rmcp::ones(const int row, const int column)
{
	return rmcp::dmatrix(row, column, 1.0);
}

dmatrix
rmcp::rand(const int row, const int column)
{
	rmcp::dmatrix mr(row, column);
	for (int i = 0; i < mr.size_; i++) {
		mr.element_[i] = 2.0 * ( (double)std::rand() / (double)RAND_MAX - 0.5 );
	}
	return mr;	
}

dmatrix
rmcp::eye(const int row, const int column)
{
	rmcp::dmatrix mr(row, column, 0.0);

	int n = (mr.row_ < mr.column_) ? mr.row_ : mr.column_; // min(r,c)
	for (int i = 0; i < mr.size_; i += (mr.row_ + 1)) {
		mr.element_[i] = 1.0;
	}
	return mr;
}

/*
bool
rmcp::solve_AXeB(const dmatrix &A, dmatrix &X, const dmatrix &B)
{
	assert(A.row_ == B.row_ && "solve_AxeB(..)");

	if (A.row_ == A.column_) {
		dmatrix bfA(A);
		x = B;

		int lda = bfA.row_;
		int n = bfA.column_;
		int* ipvt = new int[n];
		int job = 0;	// solve a * x = b	
		int info;

		dgefa_(bfA.element_, &lda, &n, ipvt, &info);
		dgesl_(bfA.element_, &lda, &n, ipvt, x.element_, &job);

		delete[] ipvt;

		return true;
	}
	else {
			if ( A.row > A.col ) return SolveAxEqualB(A ^ A, x, A ^ B);
		bool flag = SolveAxEqualB(A | A, x, B);
		x = A ^ x;
		return flag;
		//return QRSolveAxEqualB(A, x, B);
		
		return true;
	}

}
*/


//// Derived

dmatrix::dmatrix(const dse3 &S)
{
	row_ = 4;
	column_ = 4;
	size_ = 16;
	element_ = new double[size_];

	element_[0] = 0.0;
	element_[1] = S.element_[2];
	element_[2] = - S.element_[1];
	element_[3] = 0.0;
	element_[4] = - S.element_[2];
	element_[5] = 0.0;
	element_[6] = S.element_[0];
	element_[7] = 0.0;
	element_[8] = S.element_[1];
	element_[9] = - S.element_[0];
	element_[10] = 0.0;
	element_[11] = 0.0;
	element_[12] = S.element_[3];
	element_[13] = S.element_[4];
	element_[14] = S.element_[5];
	element_[15] = 0.0;
}

dmatrix::dmatrix(const dinertia &I)
{
	row_ = 6;
	column_ = 6;
	size_ = 36;
	element_ = new double[size_];

	element_[0] = I.I_[0];
	element_[6] = I.I_[3];
	element_[12] = I.I_[4];
	element_[18] = 0.0;
	element_[24] = - I.r_[2];
	element_[30] = I.r_[1];

	element_[1] = I.I_[3];
	element_[7] = I.I_[1];
	element_[13] = I.I_[5];
	element_[19] = I.r_[2];
	element_[25] = 0.0;
	element_[31] = - I.r_[0];

	element_[2] = I.I_[4];
	element_[8] = I.I_[5];
	element_[14] = I.I_[2];
	element_[20] = - I.r_[1];
	element_[26] =  I.r_[0];
	element_[32] = 0.0;
	element_[3] = 0.0;
	element_[9] = I.r_[2];
	element_[15] = - I.r_[1];
	element_[21] = I.m_;
	element_[27] = 0.0;
	element_[33] = 0.0;
	element_[4] = - I.r_[2];
	element_[10] = 0.0;
	element_[16] = I.r_[0];
	element_[22] = 0.0;
	element_[28] = I.m_;
	element_[34] = 0.0;
	element_[5] = I.r_[1];
	element_[11] = - I.r_[0];
	element_[17] = 0.0;
	element_[23] = 0.0;
	element_[29] = 0.0;
	element_[35] = I.m_;
}

