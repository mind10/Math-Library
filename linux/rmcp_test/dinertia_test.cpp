#include "Ram/Ram.h"
#include <ctime>

using namespace std;
using namespace Ram;

int main(void)
{
	clock_t time;

	/*
	 *	Ram::Matrix
	 */

	// Constructor
	Ram::Matrix M1(4, 4);
	cout << "Ram::Matrix Constructor - case 1" << endl << M1 << endl;

	Ram::Matrix M2(4, 3, 8.623);
	cout << "Ram::Matrix constructor - case 2" << endl << M2 << endl;

	Ram::Matrix M3(M2);
	cout << "Ram::Matrix constructor - case 3" << endl << M3 << endl;

	double M4_data[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
	Ram::Matrix M4(3, 4, M4_data);
	cout << "Ram::Matrix constructor - case 4" << endl << M4 << endl;


	/* 
	 * Ram::SE3
	/*
	Vec3 w1(1.0);
	Vec3 w2(1,2,3);
	Vec3 q1(0.4);
	Vec3 q2(1.0);

	Vec3 v1 = - cross(w1, q1);
	Vec3 v2 = - cross(w2, q2);

	cout << "v : " << v2 << endl;
	cout << "v, norm : " << v2.norm() << endl;

	se3 A1(w1,v1);
	se3 A2(w2,v2);

	cout << "se3 : " << A1 << endl;

	dse3 F;

	cout << "dse3 : " << F << endl;

	SE3 M(Vec3(10, 3, 5));

	double theta1 = 1.0;
	double theta2 = 1.0;
	SE3 T = exp(A1, theta1) * exp(A2, theta2) * M;
	cout << "SE3 : " << T << endl;
	T.exp(A1, theta1);

	cout << "SE3 : " << T << endl;
	//cin >> theta1;

	SO3 R = SO3(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0);
	Vec3 w = log(R);
	cout << "SO3 : " << R << endl;
	cout << w << endl;
	
	Matrix RM = M.to_Matrix();

	RM = rand(4, 4);
	double d[] = { 1, 4.2, 7, 2.1, 5, 8, 3, 6.456, 10};
	Matrix DR(3, 3, d);

	cout << M << endl;
	cout << "DR = " << DR << endl;
	cout << "inv(DR) = " << inv(DR) << endl;
	double deter = det(DR);
	cout << "det : " << deter << endl;

	double r[] = { 4.6, 8.1, 2.0};

	Matrix x;
	Matrix DB(3,1, r);
	solve_AtxeB(DR, x, DB);

	cout << " x: " << x << endl;

	double a = 4.5;
	
//	cout << "d * RM : " << 4.5 / RM << endl;

	RM = w.to_Matrix();

	cout << RM << endl;

	cout << R.to_Matrix() << endl;

//	AInertia AIner(1, 4, 5, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 10, 1, 4, 5,7, 7, 8);
	Inertia iner(3.45, 5.67, 34.3, 445.4, 54.4, 54, 54, 61, 344, 0.45);
	cout << iner.to_Matrix() << endl;

	cout << inv(iner).to_Matrix() << endl;

	Matrix ABC(3, 3, 4);
	Matrix ABCD(ABC);
	Matrix ABCDE(3, 3, 5);
	Matrix ABCDEF(3, 2, 4);

	cout << (ABC == ABCD) << '\t' << (ABC == ABCDE) << '\t' << (ABC == ABCDEF) << endl;
	


		int hold;
	cin >> hold;
*/

	/* 
	 * Ram::Inertia
	 */

	// Constructor
	Inertia I(1.3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.1, 0.2, 0.3);
	std::cout << I.to_Matrix();
	std::cout << std::endl;

	double inertia[6];

	// Member Functions : mass, center_of_mass, get_inertia
	std::cout << I.mass() << std::endl;
	std::cout << I.center_of_mass();
	I.get_inertia(inertia);
	std::cout << inertia[0] << '\t' << inertia[1] << '\t' << inertia[2] << '\t'
		<< inertia[3] << '\t' << inertia[4] << '\t' << inertia[5] << std::endl;
	std::cout << std::endl;

	// Member Functions : set_mass
	I.set_mass(0.8);

	std::cout << I.to_Matrix();
	std::cout << I.mass() << std::endl;
	std::cout << I.center_of_mass();
	I.get_inertia(inertia);
	std::cout << inertia[0] << '\t' << inertia[1] << '\t' << inertia[2] << '\t'
		<< inertia[3] << '\t' << inertia[4] << '\t' << inertia[5] << std::endl;
	std::cout << std::endl;

	// Member Functions : set_center_of_mass
	I.set_center_of_mass(Vec3(3.3, 4.4, 5.5));

	std::cout << I.to_Matrix();
	std::cout << I.mass() << std::endl;
	std::cout << I.center_of_mass();
	I.get_inertia(inertia);
	std::cout << inertia[0] << '\t' << inertia[1] << '\t' << inertia[2] << '\t'
		<< inertia[3] << '\t' << inertia[4] << '\t' << inertia[5] << std::endl;
	std::cout << std::endl;

	// Member Functions : set_inertia
	I.set_inertia(9.0, 8.0, 7.0, 6.0, 5.0, 4.0);

	std::cout << I.to_Matrix();
	std::cout << I.mass() << std::endl;
	std::cout << I.center_of_mass();
	I.get_inertia(inertia);
	std::cout << inertia[0] << '\t' << inertia[1] << '\t' << inertia[2] << '\t'
		<< inertia[3] << '\t' << inertia[4] << '\t' << inertia[5] << std::endl;
	std::cout << std::endl;
	
	/* 
	 * Friend Functions
	 */

	// det

	double sample[] = { 1.3, 4.5, 6.4, 3.4, 5.6, 7.8, 4.5, 6.7, 5.4, 3.4, 3.4, 5.6, 7.1, 8.7, 6.7, 2.9 };
	Matrix detM(4, 4, sample);
	double determinent;

	time = clock();
	for (int i = 0; i < 1234567; i++ ) {
		determinent = det(detM);
	}
	time = clock() - time;

	cout << "detM = " << detM << endl;
	cout << "determinant of detM : " << determinent << "  elpased time : " << time << endl;

	return 0;
}
