#ifndef _RMCP_DSQOBJECTFN_H_
#define _RMCP_DSQOBJECTFN_H_

#include <rmcp/rmcp.h>

namespace rmcp {

	class dsqobjectfn {
		public:
			virtual double evaluate(double x[]) = 0;
			virtual void gradient(double x[], double grad[]) = 0;
			virtual void nonlinear_constraint(double x[], double con[], int err[]) = 0;
			virtual void nonlinear_constraint(double x[], int i, double &con, int &err) = 0;
			virtual void nonlinear_gradient(double x[], int i, double grad[]) = 0;
			virtual bool newx(double x[], double u[], int itstep, double **accinf) { return true; }
	};

}

/* 
	err : 0 or 1, 0 == FALSE, 1 == TRUE
	newx(...) is called after completion of an iterative steps
	and allows the user to control the process of optimization.

	x is the current solution.
	u is the vector of Lagrange multipliers.
	itstep is the step number of the last step.
	accinf is the array of accumulated information defined from [1][] to [itstep[]
	return to be set true for continuation and false for disruption of the optimization by the user. 
*/

#endif /* _RMCP_DSQOBJECTFN_H_ */
