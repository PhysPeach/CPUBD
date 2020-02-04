#include <iostream>
#include <ctime>
using namespace std;

#include "../hpp/MT.hpp"
#include "../hpp/gaussian_rand.hpp"
#include "../hpp/Parameters.hpp"
#include "../hpp/Particle.hpp"
#include "../hpp/System.hpp"

template <typename T>
T pow_tmp(T x, int y) {
	T xx = (T)1;
	if (y > 0) {
		for (int i = 1; i <= y; i++) {
			xx *= x;
		}
	}
	else {
		for (int i = -1; i >= y; i--) {
			xx /= x;
		}
	}

	return xx;
}

//extern
double Tfin;
int N;
double tmax;
int IDs;
int IDe;

int main(int argc, char** argv){
    
	IDs = atoi(argv[1]);
	IDe = atoi(argv[2]);
    N = atoi(argv[3]);
	Tfin = atof(argv[4]);
	char tscale = atoi(argv[5]);

	tmax = pow_tmp(2., tscale);
	
	cout << "--Parameters--" << endl;
	cout << "IDstart: " << IDs << endl;
	cout << "IDend  : " << IDe << endl;
	cout << "N = " << N << endl;
	cout << "Tfin = " << Tfin << endl;
	cout << "tmax = " << tmax << endl;
	cout << "--------------" << endl << endl;

	cout << "Set Random Number" << endl;
	time_t seed;
	time(&seed);
	init_genrand((unsigned long)seed);

	cout << "see logfile" << endl;
	
	for (int id = IDs; id <= IDe; id++) {
		System s(id);
		s.initSys();
        s.getDataLD();
        s.connectLDtoMD();
		s.getDataMD();
	}

	return 0;
}
