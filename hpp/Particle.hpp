#ifndef Particle_HPP
#define Particle_HPP
#include <string>
#include <sstream>

using namespace std;

#include "Parameters.hpp"

//from gaussian_rand.hpp
double gaussian_rand();

class Particle{
private:
	double diam; //diameter
protected:
	double xv[2 * D]; //Position, Velosity

public:
	double force[D]; // Force from outside

	void vDvlpBD(double, double);
	void vDvlpOD(double, double);
	void xDvlp(double);


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
