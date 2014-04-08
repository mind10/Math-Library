#ifndef _RMCP_DORMATRIX_H_
#define _RMCP_DORMATRIX_H_

#include "rmcp/rmcp.h"
#include "rmcp/dgqmatrix.h"

namespace rmcp {

	class dormatrix : public rmcp::dgqmatrix {
		public:
			dormatrix(const int row);
			dormatrix(const int row, const double element);
			dormatrix(const int row, const double element[]);
			dormatrix(const dormatrix &m);

	};

} // namespace rmcp

#endif /* _RMCP_DORMATRIX_H_ */
