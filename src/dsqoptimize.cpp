
#include <rmcp/dsqoptimize.h>

#include <cassert>

using namespace rmcp;

dsqoptimize::dsqoptimize(rmcp::dsqobjectfn *fn)
:	n_(1),
	nlin_(1),
	nonlin_(0),
	upper_(n_ + nlin_ + nonlin_),
	lower_(n_ + nlin_ + nonlin_),
	initial_x_(n_),
	optimal_x_(n_),
	linear_constraint_(nlin_, n_) 
{
	fn_ = fn;
	analytic_ = false;
	silent_ = false;
	iteration_max_ = 0;
	step_ = 0;
	elapsed_time_ = 0.0;
	optimal_value_ = 0.0;
	nreset_ = 0;
	difftype_ = 3;
	tau0_ = 1.0e1;
	del0_ = 0.2e0;
	epsdif_ = 1.e-8;
	epsfcn_ = 1.e-16;
	taubnd_ = 1.0;
	te0_ = false;
	te1_ = false;
	te2_ = false;
	te3_ = false;

}

dsqoptimize::~dsqoptimize()
{

}

void
dsqoptimize::init(int dim, int linear, int nonlinear, int iteration, int step)
{
	assert( dim > 0 &&
			linear > -1 &&
			nonlinear > -1 &&
			iteration > 0 &&
			dim + linear + nonlinear > 0 &&
			step > 0 &&
			"void dsqoptimize::init(int, int, int, int, int)");

	n_ = dim;
	nlin_ = linear;
	nonlin_ = nonlinear;
	iteration_max_ = iteration;
	step_ = step;
	nreset_ = n_;

	upper_.resize(n_ + nlin_ + nonlin_, rmcp::SQOPTIMIZE_INF);
	lower_ .resize(n_ + nlin_ + nonlin_, -rmcp::SQOPTIMIZE_INF);
	initial_x_.resize(n_);
	optimal_x_.resize(n_);
	if (nlin_ > 0) linear_constraint_.resize(nlin_, n_);
	log_name_ = std::string("optimize");
}

void
dsqoptimize::set_upper(int i, double value)
{
	assert( i >= 0 &&
			i < n_ &&
		   "void dsqoptimize::set_upper(int i, double value)");

	upper_[i] = value;
	return;
}

void
dsqoptimize::set_lower(int i, double value)
{
	assert( i >= 0 &&
			i < n_ &&
			"void dsqoptimize::set_lower(int i, double value)");

	lower_[i] = value;
	return;
}

void
dsqoptimize::set_lincon_upper(int i, double value)
{
	assert( i >= 0 &&
			i < nlin_ &&
		   "void dsqoptimize::set_lincon_upper(int i, double value)");

	upper_[n_+ i] = value;

	return;
}

void
dsqoptimize::set_lincon_lower(int i, double value)
{
	assert( i >= 0 &&
			i < nlin_ &&
		   "void dsqoptimize::set_lincon_lower(int i, double value)");

	lower_[n_+ i] = value;

	return;
}

void
dsqoptimize::set_nonlincon_upper(int i, double value)
{
	assert( i >= 0 &&
			i < nonlin_ &&
			"void dsqoptimize::set_nonlincon_upper(int i, double value)");

	upper_[n_ + nlin_ + i] = value;

	return;
}

void
dsqoptimize::set_nonlincon_lower(int i, double value)
{
	assert( i >= 0 &&
			i < nonlin_ &&
			"void dsqoptimize::set_nonlincon_lower(int i, double value)");

	lower_[n_ + nlin_ + i] = value;

	return;
}

void
dsqoptimize::set_linear_constraint(rmcp::dmatrix &constraint)
{
	assert(nlin_ == constraint.row() &&
			n_ == constraint.column() &&
			"void dsqoptimize::set_linear_constraint(rmcp::dmatrix)");

	linear_constraint_ = constraint;
}


int dsqoptimize_n;
int dsqoptimize_nlin;
int dsqoptimize_nonlin;
int dsqoptimize_iterma;
int dsqoptimize_nstep;

double dsqoptimize_tau0;
double dsqoptimize_del0;
int dsqoptimize_nreset;
double dsqoptimize_epsdif;
double dsqoptimize_epsfcn;
double dsqoptimize_taubnd;
double dsqoptimize_big;
bool dsqoptimize_te0;
bool dsqoptimize_te1;
bool dsqoptimize_te2;
bool dsqoptimize_te3;

int dsqoptimize_optite;

double *dsqoptimize_x;
char dsqoptimize_name[40];

bool dsqoptimize_analyt;
bool dsqoptimize_silent;
int dsqoptimize_difftype;

double *dsqoptimize_up;
double *dsqoptimize_low;
double *dsqoptimize_gres;

double *dsqoptimize_optimal_x;
double dsqoptimize_fx;
double dsqoptimize_runtim;
rmcp::dsqobjectfn *dsqoptimize_fn;
double *dsqoptimize_congrad;

void donlp2(void);

bool
dsqoptimize::solve(rmcp::dvector initial_x)
{
	assert (fn_ &&
		n_ == initial_x.size() &&
		"bool dsqoptimize::solve(rmcp::dvector)");

	dsqoptimize_fn = fn_;

	dsqoptimize_n = n_;
	dsqoptimize_nlin = nlin_;
	dsqoptimize_nonlin = nonlin_;
	dsqoptimize_iterma = iteration_max_;
	dsqoptimize_nstep = step_;
	dsqoptimize_tau0 = tau0_;
	dsqoptimize_del0 = del0_;
	dsqoptimize_nreset = nreset_;
	dsqoptimize_epsdif = epsdif_;
	dsqoptimize_epsfcn = epsfcn_;
	dsqoptimize_taubnd = taubnd_;
	dsqoptimize_big = rmcp::SQOPTIMIZE_INF;
	dsqoptimize_te0 = te0_;
	dsqoptimize_te1 = te1_;
	dsqoptimize_te2 = te2_;
	dsqoptimize_te3 = te3_;

	dsqoptimize_x = initial_x.element();
	strcpy(dsqoptimize_name, log_name_.c_str());

	dsqoptimize_up = upper_.element();
	dsqoptimize_low = lower_.element();
	dsqoptimize_gres = linear_constraint_.element();
	dsqoptimize_optimal_x = optimal_x_.element();

	dsqoptimize_analyt = analytic_;
	dsqoptimize_silent = silent_;
	dsqoptimize_difftype = difftype_;	

	if (nonlin_ > 0) dsqoptimize_congrad = new double[n_];
	donlp2();
	if (nonlin_ > 0) delete dsqoptimize_congrad;

	optimal_value_ = dsqoptimize_fx;
	elapsed_time_ = dsqoptimize_runtim;
	optite_ = dsqoptimize_optite;

	return true;
}

