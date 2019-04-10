#include "Halley.h"
#include <iostream>
#include <fstream>
#include <cmath>

int main(){

	Halley zad1;
	std::ofstream file;
	// file.open("rozw1a.txt", std::ofstream::out);
	double dt = 60 * 15;
	// zad1.solveEuler(file, dt);
	// file.close();


	file.open("rozw2.txt", std::ofstream::out);
	dt = 60 * 6000;
	zad1.solveRK4(file, dt);
	file.close();


	file.open("rozw2auto.txt", std::ofstream::out);
	dt = 60 * 6000;
	zad1.solveRK4automatic(file, dt);
	file.close();

}