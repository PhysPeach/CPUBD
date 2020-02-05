#include "../hpp/Particle.hpp"

void Particle::init(double* xvinit, double diameter){

	diam = diameter;

	for (char i = 0; i < 2 * D; i++){
		xv[i] = xvinit[i];
	}

	return;
}

void Particle::vEvoLD(double dt, double thermalFuctor){

	for (char i = 0; i < D; i++){
		//v += dt*F : F = -ZT * v + force_outside + fluctuation
		//ZT is ignored.
		//thermalFuctor = sqrt(2 * ZT * T /dt)
		xv[D + i] += dt * (-xv[D + i] + force[i] + thermalFuctor * gaussian_rand()); //()=F(t)
	}

	return;
}
void Particle::xEvoLD(double dt) {

	for (char i = 0; i < D; i++) {
		xv[i] += dt * xv[D + i];
	}

	return;
}
void Particle::halfvEvoMD(double dt){

	for (char i = 0; i < D; i++){
		//v += dt*F : F = force_outside
		//ZT is ignored.
		//thermalFuctor = sqrt(2 * ZT * T /dt)
		xv[D + i] += 0.5 * dt * force[i]; //()=F(t)
	}

	return;
}
void Particle::xEvoMD(double dt) {

	for (char i = 0; i < D; i++) {
		xv[i] += dt * (xv[D + i] + 0.5 * dt * force[i]);
	}

	return;
}

//OK
