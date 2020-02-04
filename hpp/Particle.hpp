#ifndef Particle_HPP
#define Particle_HPP
#include <string>
#include <sstream>

using namespace std;

#include "Parameters.hpp"
#include "gaussian_rand.hpp"

class Particle{
private:
	double diam; //diameter
protected:
	double xv[2 * D]; //Position, Velosity

public:
	double force[D]; // Force from outside

	void vEvoLD(double, double);
	void halfvEvoMD(double);
	void xEvo(double);


	inline void resetx(char d, double x){ xv[d] = x; return; }
    inline void resetv(char d, double v){ xv[d+D] = v; return; }
	inline double getx(char i){ return xv[i]; }
	inline double getv(char i){ return xv[D + i]; }
	inline double getf(char i){ return force[i]; }
	inline double getDiam(){ return diam; }
	double getvv(){
		double vv = 0;
		for (char i = 0; i < D; i++){
			vv += xv[D + i] * xv[D + i];
		}
		return vv;
	}
	double getxx(){
		double xx = 0;
		for (char i = 0; i < D; i++){
			xx += xv[i] * xv[i];
		}
		return xx;
	}

	void init(double*, double); //Initialisation
};

#endif

