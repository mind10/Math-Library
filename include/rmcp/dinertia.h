#ifndef _RMCP_DINERTIA_H_
#define _RMCP_DINERTIA_H_

#include "rmcp/rmcp.h"
#include "rmcp/dmatrix.h"

namespace rmcp {
	
	class dinertia {

		friend class rmcp::dvector;
		friend class rmcp::dmatrix;
		friend class rmcp::dainertia;
		friend class rmcp::dse3;

		protected:
			double I_[6], r_[3], m_;

		public:
			dinertia(void);
			dinertia(double m);
			dinertia(double mass, double Ixx, double Iyy, double Izz);
			dinertia(double mass, double Ixx, double Iyy, double Izz, double Ixy,
					double Ixz, double Iyz, double r0, double r1, double r2);
			dinertia(const dinertia &I);

			dse3 operator * (const dse3 &acc) const;
			dse3 operator * (const rmcp::dvector &acc) const;
			dainertia operator + (const dainertia &J) const;
			dainertia operator - (const dainertia &J) const;
			dinertia operator + (const dinertia &J) const;
			dinertia& operator = (const dinertia &I);
			dinertia& operator *= (double x);
			bool operator == (const dinertia &I) const;

			double mass(void) const;
			rmcp::dvector com(void) const;
			void get_inertia(double I[]);

			void set_mass(double mass);
			void set_com(const rmcp::dvector &v);
			void set_inertia(double Ixx, double Iyy, double Izz, double Ixy,
					double Ixz, double Iyz);

			dinertia transform(const dSE3 &T) const;
			void to_array(double I[]) const;

			friend dinertia operator * (double x, const dinertia &I);
			friend dainertia inv(const dinertia &I);
					
	};


	dainertia inv(const dinertia &I);

	class dainertia {
		friend class dinertia;

		public:
			dainertia(void);
			dainertia(double d);
			dainertia(const dinertia& I);
			dainertia(double a0, double a3, double a4, double a6, double a7, double a8,
					double b0, double b1, double b2, double b3, double b4, double b5,
					double b6, double b7, double b8, double c0, double c3, double c4,
					double c6, double c7, double c8);
			dainertia(const dainertia& J);

			const dainertia& operator + (void) const;
			dainertia operator - (void) const;
			dse3 operator * (const dse3& a) const;
			dainertia operator + (const dainertia& J) const;
			dainertia operator + (const dinertia& J) const;
			dainertia operator - (const dainertia& J) const;
			dainertia operator - (const dinertia& J) const;
			dainertia& operator += (const dainertia& J);
			dainertia& operator += (const dinertia& J);
			dainertia& operator -= (const dainertia& J);
			dainertia& operator -= (const dinertia& J);
			bool operator == (const dainertia& J) const;
			dse3 operator % (const dse3& t) const;

			void subtract_kronecker_product(const dse3& x, const dse3& y);
			void add_transform(const dainertia& J, const dSE3& T);
			dainertia transform(const dSE3& T) const;
			void to_array(double I[]) const;
			rmcp::dgematrix to_dgematrix(void) const;
			
			friend dainertia kronecker_product(const dse3& x, const dse3& y);
			friend dainertia kronecker_product(const dse3& x, const dse3& y, double d);
			friend dainertia inv(const dinertia& I);
			
		private:
			double J_[36];
	};

	dainertia kronecker_product(const dse3& x, const dse3& y);
	dainertia kronecker_product(const dse3& x, const dse3& y, double d);
	dainertia inv(const dinertia& I);

} 

#endif /* _RMCP_DINERTIA_H_ */
