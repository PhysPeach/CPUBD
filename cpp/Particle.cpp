#include "../hpp/Particle.hpp"

void Particle::init(double* xvinit, double diameter){

	diam = diameter;

	for (char i = 0; i < 2 * D; i++){
		xv[i] = xvinit[i];
	}

	return;
}

void Particle::vDvlpBD(double dt, double thermalFuctor){

	for (char i = 0; i < D; i++){
		//v += dt*F : F = -ZT * v + force_outside + fluctuation
		//ZT is ignored.
		//thermalFuctor = sqrt(2 * ZT * T /dt)
		xv[D + i] += dt * (-xv[D + i] + force[i] + thermalFuctor * gaussian_rand()); //()=F(t)
	}

	return;
}
void Particle::vDvlpOD(double dt, double thermalFuctor) {
	for (char i = 0; i < D; i++) {
		//x += F*dt/ZT : F = force_outside + fluctuation
		//ZT is ignored.
		//thermalFuctor = sqrt(2 * ZT * T /dt)
		xv[D + i] = force[i] + thermalFuctor * gaussian_rand(); //()=F(t)
	}
}
void Particle::xDvlp(double dt) {

	for (char i = 0; i < D; i++) {
		xv[i] += dt * xv[D + i];
	}

	return;
}

//OK
