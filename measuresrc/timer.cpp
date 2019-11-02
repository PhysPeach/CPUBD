#include <iostream>
#include <ctime>
using namespace std;

#include "../hpp/MT.hpp"
#include "../hpp/timer.hpp"
#include "../hpp/gaussian_rand.hpp"
#include "../hpp/Parameters.hpp"
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
    
	IDs = 0;
	IDe = 0;
    N = atoi(argv[1]);
	Tfin = 1.;
	tmax = 1.;
	
	cout << "--MeasureTime--" << endl;
	cout << "N = " << N << endl;
	cout << "---------------" << endl << endl;

	cout << "Set Random Number" << endl;
	time_t seed;
	time(&seed);
    init_genrand((unsigned long)seed);

    System s(0);
    s.initSys();
    const unsigned int loop = 1000;
    double endTime;
    cout <<"starting benchmark" << endl;
    measureTime();
    s.benchmark(loop);
    endTime = measureTime();

    cout << "---Results---" << endl;
    cout << "Time: " << endTime << "ms" << endl;
    cout << "Loop: " << loop << "steps" << endl;
    cout << "   -> " << endTime/(double)loop << "ms/step" << endl;
    cout << "-------------" << endl;

	return 0;
}