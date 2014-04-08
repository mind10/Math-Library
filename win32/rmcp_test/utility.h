#ifndef _RMCP_UTILITY_
#define _RMCP_UTILITY_

//#include "rmcp/dse3.h"
//#include "rmcp/dso3.h"

void util_pause(void);
void util_codeout(std::string str);
void util_cout(std::string str, rmcp::dvector v);
void util_cout(std::string str, rmcp::dmatrix m);

// test
//rmcp::dse3 aaa_;
//void test_dse3_func1(rmcp::dse3& aa) 
//{ 
//	aaa_ = aa; 
//	util_cout("use &", aaa_);
//}
//void test_dse3_func2(rmcp::dse3 aa) 
//{ 
//	aaa_ = aa;
//	util_cout("use ", aaa_);
//}

#endif