#ifndef _RMCP_DGQMATRIX_H_
#define _RMCP_DGQMATRIX_H_

#include "rmcp/dmatrix.h"

namespace rmcp {

	class dgqmatrix : public rmcp::dmatrix {

		public:
			dgqmatrix(const int row);
			dgqmatrix(const int row, const double element);
			dgqmatrix(const int row, const double element[]);
			dgqmatrix(const dgqmatrix &m);

			virtual bool identity(void);
			virtual void set_identity(void);

			virtual void inv(void);
	};

} // namespace rmcp

#endif /* _RMCP_DGQMATRIX_H_ */
