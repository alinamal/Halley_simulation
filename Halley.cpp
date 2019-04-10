#include "Halley.h"
#include <iostream>
#include <cmath>

Halley::Halley():G(6.6741e-11), M(1.989e30), au(149597870700){
}

// moje błędy: pomnożyłam tutaj przez dt, niepotrzebnie jeśli w stepEuler / stepRK4 mnozymy
double Halley::acceleration(double coord1, double coord2, double dt){
	return -G * M * coord1 / pow( coord1*coord1 + coord2*coord2, 1.5 );
}

// moje błędy: nie pomnożyłam przez dt wszystkich dodawanych wartości
void Halley::stepEuler(double dt){
	double x_old = x;
	double y_old = y;
	x = x + dt*vx;
	y = y + dt*vy;

	vx = vx + acceleration(x_old, y_old,dt) * dt;
	vy = vy + acceleration(y_old, x_old,dt) * dt;
}

void Halley::solveEuler( std::ofstream &f, double dt ){
	fillInitialConditions();
	
	double t = 0;
	double tmax = 2000000;
	std::cout << "Wyliczamy z dt = " << dt << ", tmax = " << tmax << std::endl;
	int nr_it=0;

	for(t=0;t<tmax+dt/2/3600;t+=dt/3600){
		if(nr_it%10==0){
				f << t << "\t" << x/au << "\t" << y/au << "\t" << vx << "\t" << vy << std::endl;
		}
		// i teraz obliczenia po kolei :)
		stepEuler(dt);
		nr_it++;
	}
}


void Halley::f_rk4(const double u[], double k[], double dt){
	k[0] = u[2];
	k[1] = u[3];
	k[2] = acceleration( u[0], u[1], dt );
	k[3] = acceleration( u[1], u[0], dt );
}

void Halley::stepRK4(double dt){
	const int rksize = 4;
	double u[rksize] = {x, y, vx, vy};
	double u_rk[rksize]{};

	double k1[rksize]{};
	double k2[rksize]{};
	double k3[rksize]{};
	double k4[rksize]{};

	f_rk4( u, k1, dt );
	for (int i = 0; i < rksize; ++i) u_rk[i] = u[i] + dt/2*k1[i];
	f_rk4( u_rk, k2, dt);
	for (int i = 0; i < rksize; ++i) u_rk[i] = u[i] + dt/2*k2[i];
	f_rk4( u_rk, k3, dt);
	for (int i = 0; i < rksize; ++i) u_rk[i] = u[i] + dt*k3[i];
	f_rk4( u_rk, k4, dt);

	for (int i = 0; i < rksize; ++i) u[i] = u[i] + dt/6*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
	x = u[0];
	y = u[1];
	vx = u[2];
	vy = u[3];
}
void Halley::stepRK4_2(double dt){
	const int rksize = 4;
	double u[rksize] = {x2, y2, vx2, vy2};
	double u_rk[rksize]{};

	double k1[rksize]{};
	double k2[rksize]{};
	double k3[rksize]{};
	double k4[rksize]{};

	f_rk4( u, k1, dt );
	for (int i = 0; i < rksize; ++i) u_rk[i] = u[i] + dt/2*k1[i];
	f_rk4( u_rk, k2, dt);
	for (int i = 0; i < rksize; ++i) u_rk[i] = u[i] + dt/2*k2[i];
	f_rk4( u_rk, k3, dt);
	for (int i = 0; i < rksize; ++i) u_rk[i] = u[i] + dt*k3[i];
	f_rk4( u_rk, k4, dt);

	for (int i = 0; i < rksize; ++i) u[i] = u[i] + dt/6*( k1[i] + 2*k2[i] + 2*k3[i] + k4[i] );
	x2 = u[0];
	y2 = u[1];
	vx2 = u[2];
	vy2 = u[3];
}

void Halley::solveRK4( std::ofstream &f, double dt ){
	fillInitialConditions();
	
	double t = 0;
	double tmax = 2000000; // w godzinach, nie sekundach
	std::cout << "RK4 z dt = " << dt << ", tmax = " << tmax << std::endl;
	int nr_it=0;

	for(t=0;t<tmax+dt/2/3600;t+=dt/3600){
		if(nr_it%10==0){
			f << t << "\t" << x/au << "\t" << y/au << "\t" << vx << "\t" << vy << std::endl;
		}
		// i teraz obliczenia po kolei :)
		stepRK4(dt);
		nr_it++;
	}
}

// moje błędy: zapomniałam, że x,y,vx,vy aktualizujemy tylko jeśli krok
// czasowy został zaakceptowany
void Halley::solveRK4automatic( std::ofstream &f, double dt ){
	fillInitialConditions();
	fillInitialConditions_2();

	double t = 0;
	double tol = 100000;
	// double tmax = 60 * 60 * 80000;
	double tmax = 2500000;
	std::cout << "RK4 automatic = " << dt << ", tmax = " << tmax << std::endl;
	int nr_it=0;
	double epsilon_x;
	double epsilon_y;
	double epsilon;

	while(t<tmax+dt/2/3600){
		double xt = x, yt = y, vxt = vx, vyt= vy;
		double xt2 = x2, yt2 = y2, vxt2 = vx2, vyt2= vy2;
		// i teraz obliczenia po kolei :)
		stepRK4(dt);

		stepRK4_2(dt/2);
		stepRK4_2(dt/2);

		epsilon_x = fabs(x2-x)/(pow(2,4)-1);
		epsilon_y = fabs(y2-y)/(pow(2,4)-1);

		epsilon = std::max(epsilon_y, epsilon_x);

		if( epsilon<tol ){
			t+=dt/3600;
			x=x2;
			y=y2;
			vx=vx2;
			vy=vy2;
			f << t << "\t" << x/au << "\t" << y/au << "\t" << vx << "\t" << vy << "\t" << dt << std::endl;
		}
		else{ // if error not accepted, we set back the values prior to stepRK4.
			x=xt; y=yt; vx= vxt; vy=vyt;
			x2=xt2; y2=yt2; vx2= vxt2; vy2=vyt2;
		}
		dt = dt * 0.9 * pow(tol/epsilon, 1./5.);

		nr_it++;
	}
}

void Halley::fillInitialConditions(){
	x = 0;
	y = 0.586*au;
	vx = 54600;///au*8901554;
	vy = 0;
}

void Halley::fillInitialConditions_2(){
	x2 = 0;
	y2 = 0.586*au;
	vx2 = 54600;///au*8901554;
	vy2 = 0;
}