
#include "rmcp/dsqobjectfn.h"
#include "rmcp/o8para.h"

/* **************************************************************************** */
/*                              donlp2-intv size initialization                 */
/* **************************************************************************** */
void user_init_size(void) {

    #define  X extern
    #include "rmcp/o8comm.h"
    #include "rmcp/o8fint.h"
    #include "rmcp/o8cons.h"
    #undef   X
	
	extern int dsqoptimize_n;
	extern int dsqoptimize_nlin;
	extern int dsqoptimize_nonlin;
	extern int dsqoptimize_iterma;
	extern int dsqoptimize_nstep;

    /* problem dimension n = dim(x), nlin=number of linear constraints
       nonlin = number of nonlinear constraints  */
    
    n      = dsqoptimize_n;
    nlin   = dsqoptimize_nlin;
    nonlin = dsqoptimize_nonlin;
    iterma = dsqoptimize_iterma;
    nstep  = dsqoptimize_nstep;

}

/* **************************************************************************** */
/*                              donlp2-intv standard setup                      */
/* **************************************************************************** */
void user_init(void) {

	#define  X extern
	#include "rmcp/o8comm.h"
	#include "rmcp/o8fint.h"
	#include "rmcp/o8cons.h"
	#undef   X

	extern double *dsqoptimize_x;
	extern char dsqoptimize_name[40];
	extern bool dsqoptimize_analyt;
	extern bool dsqoptimize_silent;
	extern int dsqoptimize_difftype;
	extern double dsqoptimize_tau0;
	extern double dsqoptimize_del0;
	extern double dsqoptimize_epsdif;
	extern double dsqoptimize_epsfcn;
	extern double dsqoptimize_taubnd;
	extern double dsqoptimize_big;

	extern double *dsqoptimize_up;
	extern double *dsqoptimize_low;
	extern double *dsqoptimize_gres;

	static INTEGER  i, j;

	/* name is ident of the example/user and can be set at users will       */
	/* the first static character must be alphabetic. 40 characters maximum */

	strcpy(name, dsqoptimize_name);

	/* x is initial guess and also holds the current solution */
	/* problem dimension n = dim(x), nlin=number of linear constraints
	nonlin = number of nonlinear constraints  */

	analyt = (dsqoptimize_analyt) ? TRUE : FALSE;
	silent = (dsqoptimize_silent) ? TRUE : FALSE;
	difftype = dsqoptimize_difftype;

	epsdif = dsqoptimize_epsdif;   /* gradients exact to machine precision */
	/* if you want numerical differentiation being done by donlp2 then:*/
	epsfcn   = dsqoptimize_epsfcn; /* function values exact to machine precision */
	taubnd   = dsqoptimize_taubnd;

	del0 = dsqoptimize_del0;
	tau0 = dsqoptimize_tau0;

	big = dsqoptimize_big;

	/* initial value of x */
	for (i = 1 ; i <= n ; i++) {
		x[i] = dsqoptimize_x[i-1];
	}

	/*  set lower and upper bounds */
	for (i = 1; i <= n + nlin + nonlin; i++) {
		up[i] = dsqoptimize_up[i - 1];
		low[i] = dsqoptimize_low[i - 1];
	}

	/* store coefficients of linear constraints directly in gres */
	for ( i = 1; i <= nlin; i++) {
		for ( j = 1; j <= n; j++) {
			gres[j][i] = dsqoptimize_gres[i - 1 + (j - 1) * nlin];
		}
	}
	
	return;
}

/* **************************************************************************** */
/*                                 special setup                                */
/* **************************************************************************** */
void setup(void) {
	#define  X extern
	#include "rmcp/o8comm.h"
	#undef   X

	extern int dsqoptimize_nreset;
	extern bool dsqoptimize_te0;
	extern bool dsqoptimize_te1;
	extern bool dsqoptimize_te2;
	extern bool dsqoptimize_te3;

	nreset = dsqoptimize_nreset;

	te0 = (dsqoptimize_te0) ? TRUE : FALSE;
	te1 = (dsqoptimize_te1) ? TRUE : FALSE;
	te2 = (dsqoptimize_te2) ? TRUE : FALSE;
	te3 = (dsqoptimize_te3) ? TRUE : FALSE;
	return;
}

/* **************************************************************************** */
/*  the user may add additional computations using the computed solution here   */
/* **************************************************************************** */
void solchk(void) {
	#define  X extern
	#include "rmcp/o8comm.h"
	#undef   X
	#include "rmcp/o8cons.h"

	extern double *dsqoptimize_optimal_x;
	extern double dsqoptimize_fx;
	extern double dsqoptimize_runtim;

	static INTEGER i;

	for ( i = 1; i <= n; i++) {
		dsqoptimize_optimal_x[i - 1] = x[i];
	}

	dsqoptimize_fx = fx;
	dsqoptimize_runtim = runtim;

	return;
}

/* **************************************************************************** */
/*                               objective function                             */
/* **************************************************************************** */

void ef(DOUBLE x[],DOUBLE *fx) {
	#define  X extern
	#include "rmcp/o8fuco.h"
	#undef   X

	extern rmcp::dsqobjectfn *dsqoptimize_fn;

	*fx = dsqoptimize_fn->evaluate(x + 1);

	return;
}

/* **************************************************************************** */
/*                          gradient of objective function                      */
/* **************************************************************************** */
void egradf(DOUBLE x[],DOUBLE gradf[]) {
	#define  X extern
	#include "rmcp/o8fuco.h"
	#undef   X

	extern rmcp::dsqobjectfn *dsqoptimize_fn;
		
	dsqoptimize_fn->gradient(x + 1, gradf + 1);

	return;
}

/* **************************************************************************** */
/*                  compute the i-th nonlinear constraints                      */
/* **************************************************************************** */
void econ(INTEGER type ,INTEGER liste[], DOUBLE x[],DOUBLE con[], LOGICAL err[]) {
	#define  X extern
	#include "rmcp/o8fuco.h"
	#undef   X

	extern rmcp::dsqobjectfn *dsqoptimize_fn;

	static int i;

	switch (type) {
		case 1 :
			dsqoptimize_fn->nonlinear_constraint(x + 1, con + 1, err + 1);
			break;
		case 2:
			for (i = 1; i <= liste[0]; i++) {
				dsqoptimize_fn->nonlinear_constraint(x + 1, liste[i] - 1, con[liste[i]], err[liste[i]]);
			}
			break;
	}

	return;
}

/* **************************************************************************** */
/*          compute the gradient of the i-th equality constraint                */
/* **************************************************************************** */
void econgrad(INTEGER liste[] ,INTEGER shift, DOUBLE x[], DOUBLE **grad) {
	#define  X extern
	#include "rmcp/o8fuco.h"
	#undef   X

	extern rmcp::dsqobjectfn *dsqoptimize_fn;
	extern double *dsqoptimize_congrad;
	static int i, j;

	for (i = 1; i <= liste[0]; i++) {
		dsqoptimize_fn->nonlinear_gradient(x + 1, liste[i] - 1, dsqoptimize_congrad);
		for (j = 1; j <= n; j++) {
			grad[j][liste[i] + shift] = dsqoptimize_congrad[j - 1];
		}
	}
	return;
}


/* **************************************************************************** */
/*                        user functions (if bloc == TRUE)                      */
/* **************************************************************************** */
void eval_extern(INTEGER mode) {
	#define  X extern
	#include "rmcp/o8comm.h"
	#include "rmcp/o8fint.h"
	#undef   X
	#include "rmcp/o8cons.h"

	return;
}
