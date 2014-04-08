#ifndef _RMCP_DUQUATERNION_H_
#define _RMCP_DUQUATERNION_H_

#include <rmcp/rmcp.h>
#include <rmcp/dquaternion.h>

namespace rmcp {
	
	class duquaternion : public rmcp::dquaternion {

		protected :
		public :
			duquaternion();
			duquaternion(double s, double x, double y, double z);

			/* axis : 3 dimensional, angle = radian = [0 2pi] */

			void quaternion2euler(rmcp::dvector &euler_rpy);
			void quaternion2euler(rmcp::dvector axis, double angle);
			void euler2quaternion(rmcp::dvector euler_rpy);
			void euler2quaternion(rmcp::dvector axis, double angle);
			void get_rotation(rmcp::dvector axis, double angle);
			void convert(dSO3 &R);


	};
} 

#endif /* _RMCP_DUQUATERNION_H_ */
