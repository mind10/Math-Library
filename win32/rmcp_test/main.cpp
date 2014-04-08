
#pragma comment(lib, "mkl_intel_c")
#pragma comment(lib, "mkl_intel_thread")
#pragma comment(lib, "mkl_core")
//#pragma comment(lib, "libiomp5md")

void dvector_test(void);
void dmatrix_test(void);
void dgqmatrix_test(void);
void dvec3_test(void);
void dmat3_test(void);
void dhmat_test(void);
void dsqoptimize_test(void);
void dso3_test(void);
void dSO3_test(void);
void dse3_test(void);
void dSE3_test(void);
//void dcspline_test(void);
void dquaternion_test(void);
void duquaternion_test(void);
void dvec3_test(void);
void dmat3_test(void);
void dhmat_test(void);

#include <cmath>
#include <iostream>

int main(void)
{
//	dvector_test();
//	dmatrix_test();
//	dgqmatrix_test();
//	dsqoptimize_test();
//	dvec3_test();
//	dmat3_test();
//	dhmat_test();
//	dso3_test();
//	dSO3_test();
	dse3_test();
//	dSE3_test();

//	dcspline_test();
//	dquaternion_test();
//	duquaternion_test();

//	dvec3_test();
//	dmat3_test();
//	dhmat_test();
//	dSO3_test();
	return 0;
}
