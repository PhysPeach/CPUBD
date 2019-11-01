#include <cmath>
using namespace std;

//funcs in MT.h

double genrand_real1(void);

//gaussian_rand
double gaussian_rand() {

	static bool iset = true;
	static double gset;
	double fac, rsq, v1, v2;

	if (iset) {
		do {
			//get 2 uniform random numbers v = [-1, 1]
			v1 = 2 * genrand_real1() - 1; //genrand_real1 = [0,1]
			v2 = 2 * genrand_real1() - 1;

			//length of random number
			rsq = v1 * v1 + v2 * v2;

		} while (rsq >= 1.0 || rsq == 0.0);

		fac = sqrt(-2.0 * log(rsq) / rsq);

		gset = fac * v1;
		iset = false;
		return fac * v2;
	}
	else {
		iset = true;
		return gset;
	}
}