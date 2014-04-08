#ifndef _RMCP_DONMATRIX_H_
#define _RMCP_DONMATRIX_H_

#include "rmcp/dormatrix.h"

namespace rmcp {

	class donmatrix : public rmcp::dormatrix {

		public:
			donmatrix(const int row);
			donmatrix(const int row, const double element);
			donmatrix(const int row, const double element[]);
			donmatrix(const donmatrix &m);

			virtual void inv(void);

			friend donmatrix inv(const donmatrix &m);
	};

	donmatrix inv(const donmatrix &m);
	
} // namespace rmcp

#endif /* _RMCP_DONMATRIX_H_ */
