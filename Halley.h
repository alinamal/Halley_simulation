#pragma once
#include <vector>
#include <fstream>

// enum class BC{
// 	STIFF,
// 	LOOSE
// };

class Halley{
public:
	Halley();

	// performs single step of the explicit Euler scheme; 
	// the variables x,y,vx,vy are updated 
	void stepEuler( double dt );
	// solves the problem using Euler scheme in time (0, tmax)   
	void solveEuler( std::ofstream &f , double dt );

	// performs single step of the RK4 scheme; 
	// the variables x,y,vx,vy are updated 
	void stepRK4( double dt );
	// performs single step of the RK4 scheme; 
	// the variables x2,y2,vx2,vy2 are updated
	void stepRK4_2( double dt );
	// solves the problem using RK4 scheme in time (0, tmax)   
	void solveRK4( std::ofstream &f , double dt); 

	// solves the problem using RK4 scheme in time (0, tmax)   
	// using the automatic time step adjustment
	void solveRK4automatic( std::ofstream &f , double dt );
	
private:
	const double G; // gravitational constant
	const double M; // mass o Sun
	const double au; // astronomic unit
	double x;	// x coordinate of the comet position
	double y;	// y coordinate of the comet position
	double vx;	// x coordinate of the comet velocity
	double vy;	// y coordinate of the comet velocity

	double x2;	// same as x, used for the automatic time step adjustment
	double y2;	// as above
	double vx2;	// as above
	double vy2;	// as above

	// calculates the acceleration. The determinant is symmetric for x,y
	// swapped, the x or y coordinate is calculated depending on the order
	// of x and y as arguments coord1 and coord2. 
	double acceleration(double coord1, double coord2, double dt);

	// insert the initial values for x,y,vx,vy
	void fillInitialConditions();
	// insert the initial values for x2,y2,vx2,vy2
	void fillInitialConditions_2();
	// calculates the values of every k1, k2, k3, and k4 for RK4 scheme
	void f_rk4(const double u[], double k[], double dt);
};