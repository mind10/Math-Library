#ifndef _RMCP_DMATRIX_H_
#define _RMCP_DMATRIX_H_

#include <rmcp/rmcp.h>
#include <iostream>

namespace rmcp {

	class dmatrix : public rmcp::_dfullmatrix {

		friend class rmcp::dvector;
		friend class rmcp::dse3;
		friend class rmcp::dinertia;

		public:
			dmatrix();
			dmatrix(const int row, const int column);
			dmatrix(const int row, const int column, const double element);
			dmatrix(const int row, const int column, const double element[]);
			dmatrix(const dmatrix &m);
//			dmatrix(const dvector &v);
			virtual ~dmatrix();

			double &operator [] (const int i);
			const double &operator [] (int i) const;
			double &operator () (const int i, const int j);
			const double &operator () (int i, int j) const;
			const dmatrix &operator + (void) const;
			dmatrix operator - (void) const;
			dmatrix operator ~ (void) const;	// transpose
			dmatrix &operator = (const dmatrix &m);
			dmatrix &operator = (const double element[]);
			dmatrix &operator += (const dmatrix &m);
			dmatrix &operator -= (const dmatrix &m);
			dmatrix &operator *= (double s);
			dmatrix &operator *= (const dmatrix &m);
			dmatrix &operator /= (double s);
//			const dmatrix &operator /= (const dmatrix &m);
			dmatrix operator + (const dmatrix &m) const;
			dmatrix operator - (const dmatrix &m) const;
			dmatrix operator * (const double s) const;
			dmatrix operator / (const double s) const;
			dmatrix operator * (const dmatrix &m) const;
			dvector operator * (const dvector &v) const;
	/*		dmatrix operator / (const dmatrix &m) const; // A * inv(B)
			dmatrix operator % (const dmatrix &m) const; // inv(A) * B
			dmatrix operator | (const dmatrix &m) const; // A * ~B
			dmatrix operator ^ (const dmatrix &m) const; // ~A * B
			dmatrix operator & (const dmatrix &m) const; // inv(~A) * B
			dmatrix operator @ (const dmatrix &m) const; // A * inv(~B)
	*/		bool operator == (const dmatrix &m) const;

			int row(void) const { return row_; }
			int column(void) const { return column_; }
			int size(void) const { return size_; }
			double *element(void) const { return element_; }
		
			virtual void set(const double element);
			virtual void set(const double element[]);
			virtual void set_column(const rmcp::dvector &v, const int idx);
			virtual void set_row(const rmcp::dvector &v, const int idx);
			virtual void set_partial(const rmcp::dmatrix &m, const int top, const int left);
			virtual void set_partial(const rmcp::dvector &v, const int top, const int left);
			virtual void get(double element[]);
			virtual rmcp::dvector get_column(const int idx);
			virtual rmcp::dvector get_row(const int idx);
			rmcp::dmatrix get_partial(const int top, const int left);
			rmcp::dmatrix get_partial(const int top, const int left, const int bottom, const int right);

			rmcp::dvector export_dvector(void);

			void resize(const int row, const int column);
			void resize(const int row, const int column, const double element);
			void swap_column(const int u, const int v);
			void swap_row(const int u, const int v);
	
			bool zero(void) const;

			virtual double norm(void) const;
			virtual double normalize(void);
			double trace(void) const;
			virtual double det(void);
			virtual void transpose(void);
			virtual void inv(void);
			virtual void eig(dmatrix &v, dmatrix &d);
			rmcp::dvector mean(int dim = 0); // dim = 0; column's mean, dim = 1; row's mean
			rmcp::dvector var(int dim = 0, int n = 0); 
			// dim = 0; column's variance, dim = 1; row's variance
			// n = 0 : divided by size-1, n = 1: divided by size
			rmcp::dmatrix cov(int n = 0);		// return = symmetric matrix
			// n = 0 : divided by row_-1, n = 1: divided by row_ , row_ is the number of observations
			// each row is a observation, and each column a variable
			void pca(dmatrix &PC, rmcp::dvector &mean, rmcp::dvector &eigenvalue);
			
			friend std::ostream & operator << (std::ostream & os, const dmatrix &m);
			friend dmatrix operator * (double s, const dmatrix &m);
			friend void multiply_AB(dmatrix &re, const dmatrix &A, const dmatrix &B);
			friend void multiply_ABt(dmatrix &re, const dmatrix &A, const dmatrix &B);
			friend void multiply_AtB(dmatrix &re, const dmatrix &A, const dmatrix &B);
			friend double quadratic(const dvector &x, const dmatrix &A, const dvector &y);
		//	friend double MaxVec(const dmatrix &m, int* idx);
		//	friend double MinVec(const dmatrix &m, int* idx);
			friend dmatrix rand(const int row, const int column);
			friend dmatrix eye(const int row, const int column);
			friend double det(dmatrix &M);
		//	friend dmatrix inv(const dmatrix &m);
		//	friend bool	solve_AXeB(const dmatrix &A, dmatrix &X, const dmatrix &B);
		//	friend bool	solve_AtxeB(const dmatrix &A, dmatrix & x, const dmatrix & B);
		//	friend bool	solve_poAxeB(const dmatrix &A, dmatrix & x, const dmatrix & B);
		//	friend int GaussElimination(dmatrix &A, Idmatrix & row__pivot, Idmatrix & col_umn_pivot);
		//	friend int GaussElimination(dmatrix &A, Idmatrix & row__pivot, Idmatrix & col_umn_pivot, double eps);
		//	friend dmatrix Eig(dmatrix m);
		//	friend void	Eig(dmatrix & re, dmatrix &m);
		//	friend void	Eig(dmatrix &m, dmatrix &v, dmatrix & d);
		//	friend bool	QRSolveAxEqualB(const dmatrix &A, dmatrix & x, const dmatrix & B);
		//	friend bool	QRSolveAtxEqualB(const dmatrix & A, dmatrix & x, const dmatrix & B);
		//	friend dmatrix SVD(const dmatrix & M);
		//	friend void	SVD(const dmatrix & M, dmatrix & U, dmatrix & S, dmatrix & V);
		//	friend bool	SVDSolveAxEqualB(const dmatrix & A, dmatrix & x, const dmatrix & B);
		//	friend int Rank(const dmatrix & m, double singular_criteria);
		//	friend int Rank(const dmatrix & m);
		//	friend bool	SolveLCP(const dmatrix & A, dmatrix & f, const dmatrix & b);
	

			// dse3, dinteria
			dmatrix(const dse3 &S);
			dmatrix(const dinertia &I);

	};

	double norm(const dmatrix &m);
	dmatrix zeros(const int row, const int column);
	dmatrix ones(const int row, const int column);
	double trace(const dmatrix &m);
	dmatrix transpose(const dmatrix &m);
	dvector diagonal(dmatrix &m);

	void multiply_AB(dmatrix &re, const dmatrix &A, const dmatrix &B);
	void multiply_ABt(dmatrix &re, const dmatrix &A, const dmatrix &B);
	void multiply_AtB(dmatrix &re, const dmatrix &A, const dmatrix &B);
	double quadratic(const dvector &x, const dmatrix &A, const dvector &y);
	dmatrix rand(const int row, const int column);
	dmatrix eye(const int row, const int column);
	double det(dmatrix &M);
	
} // namespace rmcp

#endif // _RMCP_DMATRIX_H_

