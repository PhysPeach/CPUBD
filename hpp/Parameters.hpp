#ifndef Parameters_HPP
#define Parameters_HPP

#include <cmath>
using namespace std;

//Mathematical Constant
const double pi = 4 * atan(1.0);

//Time Constants
const double dt_init = 0.01;
const double dt_BD = 0.003;
const double dt_MD = 0.001;
extern double tmax;

//Physical Constants

//Dimentions
const char D = 2; //�����2

//Particles Diameter
const double a0 = 1;
const double a1 = a0 * 1;
const double a2 = a0 * 1.4;

//Number of Particle

extern int N;

//Numbers of Initial Condition
extern int IDs;
extern int IDe;

//Density of Particles (Volume of Box = N / Density)
const double dnsty = 0.8;

//Initial Lattice Constant
const double linit = 0.5*(a1 + a2);


//Friction Parameter: (ZT/m)sqrt(m * a^2/ep)
//ZT = 1.0;

//Temparature Parameter: k_b * T / ep
extern double Tfin;

#endif

//Checked
