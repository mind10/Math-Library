
#ifndef _RMCP_DSQOPTIMIZE_H_
#define _RMCP_DSQOPTIMIZE_H_

#include <rmcp/dsqobjectfn.h>
#include <rmcp/dvector.h>
#include <rmcp/dmatrix.h>

#include <string>

namespace rmcp {
	
	class dsqoptimize {

		private:
			rmcp::dsqobjectfn *fn_;	///< object function
			bool analytic_;				///< has analytic gradient function or not?
			bool silent_;				///< printing the outputs to file or console, or not
			int n_;					///< problem dimension n = dim(x)
			int nlin_;				///< number of linear constraints
			int nonlin_;				///< number of nonlinear constraints
			int iteration_max_;		///< maximum number of steps allowed
			int step_;
			double tau0_;
			double del0_;
			int nreset_;			///< nreset
			double epsdif_;
			double epsfcn_;
			double taubnd_;
			bool te0_;
			bool te1_;
			bool te2_;
			bool te3_;
			rmcp::dvector upper_;				///< upper limit of x
			rmcp::dvector lower_;				///< lower limit of x
			rmcp::dvector initial_x_;			///< starting value on optimization
			rmcp::dvector optimal_x_;			///< optimal solution x
			double optimal_value_;		///< optimal value of object function
			std::string	log_name_;			///< name of log file
			rmcp::dmatrix linear_constraint_;	///< coefficient matrix of linear constraints
			double elapsed_time_;		///< time during optimization
			int difftype_;					///< the method of getting numerical gradient
			/* ordinary forward difference = 1
			 symmetric difference = 2
			 sixth order approximation = 3 */

			int optite_;

		public:
			dsqoptimize(rmcp::dsqobjectfn *fn);
			virtual ~dsqoptimize();

			int difftype(void) { return difftype_; }
			double optimal_value(void) const { return optimal_value_; }
			double elapsed_time(void) const { return elapsed_time_; }
			rmcp::dvector optimal_solution(void) { return optimal_x_; }

			void set_difftype(int difftype) { difftype_ = difftype; }
			void set_objectfn(rmcp::dsqobjectfn *fn) { fn_ = fn; }
			void set_analytic(bool analytic) { analytic_ = analytic; }
			void set_silent(bool silent) { silent_= silent; }
			void set_upper(int i, double value);
			void set_lower(int i, double value);
			void set_lincon_upper(int i, double value);
			void set_lincon_lower(int i, double value);
			void set_nonlincon_upper(int i, double value);
			void set_nonlincon_lower(int i, double value);
			void set_linear_constraint(rmcp::dmatrix &constraint);
			void set_logfile(std::string name) { log_name_ = name; }
			void set_tau0(double tau0) { tau0_ = tau0; }
			void set_del0(double del0) { del0_ = del0; }
			void set_nreset(int nreset)	{ nreset_ = nreset; }
			void set_epsdif(double epsdif) { epsdif_ = epsdif; }
			void set_epsfcn(double epsfcn) { epsfcn_ = epsfcn; }
			void set_taubnd(double taubnd) { taubnd_ = taubnd; }
			void set_te0(bool te0) { te0_ = te0; }
			void set_te1(bool te1) { te1_ = te1; }
			void set_te2(bool te2) { te2_ = te2; }
			void set_te3(bool te3) { te3_ = te3; }
			int get_optite(void) { return optite_; }
			void init(int dim, int linear, int nonlinear, int iteration = 100, int step = 20);
			bool solve(rmcp::dvector initial_x);

	};

} // namespace rmcp 


#endif // _RMCP_DSQOPTIMIZE_H_


		/* functions from donlp2 */
		/*
		void user_init_size(void);
		void donlp2(void);
		void user_eval(DOUBLE xvar[],INTEGER mode);
		void newx ( DOUBLE x[] , DOUBLE u[] , INTEGER itstep , DOUBLE **accinf, LOGICAL *cont );
		void user_init(void);
		void setup(void);
		void solchk(void);
		void ef(DOUBLE x[],DOUBLE *fx);
		void egradf(DOUBLE x[],DOUBLE gradf[]);
		void econ(INTEGER type ,INTEGER liste[], DOUBLE x[],DOUBLE con[], LOGICAL err[]);
		void econgrad(INTEGER liste[] ,INTEGER shift ,  DOUBLE x[], DOUBLE **grad);
		void eval_extern(INTEGER mode);
		*/

