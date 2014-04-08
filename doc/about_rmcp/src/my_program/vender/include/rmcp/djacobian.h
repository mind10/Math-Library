#ifndef _RMCP_DJACOBIAN_H_
#define _RMCP_DJACOBIAN_H_

#include "rmcp/dmatrix.h"
#include "rmcp/dse3.h"

namespace rmcp {

	class djacobian : public rmcp::dmatrix {

		friend class rmcp::dse3;
		
		public:
			djacobian(int dof);
			djacobian(const djacobian &J);
			
			int dof(void)	{ return column_; }

			rmcp::dse3 get(const int i);
			void set(const int i, const rmcp::dse3 &S);
	};

}


#endif /* _RMCP_DJACOBIAN_H_ */
