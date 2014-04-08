#ifndef _RMCP_DCSPLINE_H_
#define _RMCP_DCSPLINE_H_

#include <rmcp/dvector.h>

/* matlab help csape, csapi */

namespace rmcp {

	class dcspline {
	public:
		enum TYPE {
			CLAMPED = 0,
			NATURAL,
			NOTAKNOT
		};

	private :
		rmcp::dvector coeff_;

	public :
		dcspline();
		~dcspline() { };

		bool init(rmcp::dcspline::TYPE type, dvector &x, dvector &y, double yp1 = 0.0, double ypn = 0.0);

		bool interp(dvector &x, double new_x, double &new_y, double &new_dy, double &new_ddy);
		bool interp(dvector &x, dvector &new_x, dvector &new_y, dvector &new_dy, dvector &new_ddy);
	};

} // namespace rmcp 

#endif // _RMCP_CSPLINE_H_


