#ifndef _RMCP_DSSMATRIX_H_
#define _RMCP_DSSMATRIX_H_

#include "rmcp/dgqmatrix.h"

namespace rmcp {

	class dssmatrix : public rmcp::dgqmatrix {

		public:
			dssmatrix(const int row);
			dssmatrix(const int row, const double element);
			dssmatrix(const int row, const double element[]);
			dssmatrix(const dssmatrix &m);
	};

} // namespace rmcp

#endif /* _RMCP_DSSMATRIX_H_ */
