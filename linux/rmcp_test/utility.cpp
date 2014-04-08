
#include "rmcp/dvector.h"
#include "rmcp/dmatrix.h"

#include <iostream>
#include <string>

void util_pause(void)
{
	std::cout << "press enter key to continue to test";
	std::cin.get();
	std::cout << std::endl;
}

void util_codeout(std::string str)
{
	std::cout << ">>> " << str << std::endl;
}

void util_cout(std::string str, rmcp::dvector v)
{
	std::cout << str << v << std::endl;
}

void util_cout(std::string str, rmcp::dmatrix m)
{
	std::cout << str << m << std::endl;
}

