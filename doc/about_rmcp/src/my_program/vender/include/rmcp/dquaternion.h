#ifndef _RMCP_DQUATERNION_H_
#define _RMCP_DQUATERNION_H_

#include "rmcp/rmcp.h"
#include <rmcp/dvector.h>

namespace rmcp {
	
	class dquaternion : public rmcp::dvector {

		public:
			dquaternion(void);
			dquaternion(double s, double x, double y, double z);
					
			void set(const double s, const double x, const double y, const double z);
			void get(double &s, double &x, double &y, double &z);
	};
} 

#endif /* _RMCP_DQUATERNION_H_ */
