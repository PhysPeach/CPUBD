#ifndef System_HPP
#define System_HPP

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "Parameters.hpp"
#include "gaussian_rand.hpp"
#include "Particle.hpp"

//from MT.h
double genrand_real3();

class System{
private:
	int id;
	double dt;
	double t; //Time
	double T; //Temparature
	double thermalFuctor; //sqrt(2 * ZT * T /dt)
	double L; //Length of Box

	Particle* p;

	int** list; //for list
    double* positionMemory[D];

	ostringstream positionFileName;
	ofstream positionFile;

	inline void setdt_T(double setdt, double setT) { dt = setdt; T = setT; thermalFuctor = sqrt(2 * setT / setdt); return; }
	inline double getV(){double V = 1.0; for(char d = 0; d < D; d++){V *= L;} return V;}
	inline double getDnsty() { return N / getV(); }
	inline double getvg(char); //velocity of the System (must be 0)
	double getU();
	inline double getK();

	inline short getM(){ return (int)(L / (6 * a0));}  //for cell

	inline void periodic(int);
	void culc_Interaction(); 
	inline void tDvlpBD();
	void culc_harmonicInteraction();
	void tHarmonicDvlp();

	void updateCell2D();
	void judgeUpdateCell();
	
	void equilibrateSys(double);

	void makeInitPosition();

	void recordSys();

public:
	System(int ID);
	~System();

	void initSys();

	void procedure();
	void benchmark(unsigned int loop);
};

#endif
