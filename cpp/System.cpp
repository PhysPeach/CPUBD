#include "../hpp/System.hpp"

template <typename T>
T pow_tmp(T x, int y){
	T xx = (T)1;
	if (y > 0){
		for (int i = 1; i <= y; i++) {
			xx *= x;
		}
	}
	else{
		for (int i = -1; i >= y; i--) {
			xx /= x;
		}
	}

	return xx;
}

int f(int i, int M){
	if (i < 0){
		return i + M;
	}
	else if (i >= M){
		return i - M;
	}
	else{
		return i;
	}
}

//protected
inline double System::getvg(char d) {
	double vg = 0;
	for (int n = 0; n < N; n++) {
		vg += p[n].getv(d) / N;
	}

	return vg;
}
double System::getU() {

	double U = 0;

	double xik[D];
	double rik;
	double aik;
	double ar;

	int k;

	for (int i = 0; i < N; i++) {
		for (int j = 1; j <= list[i][0]; j++) {

			k = list[i][j];

			//culc position difference
			rik = 0;
			for (char d = 0; d < D; d++) {
				xik[d] = p[i].getx(d) - p[k].getx(d);

				//Periodic between 2 particles;
				if (xik[d] > (0.5 * L)) {
					xik[d] -= L;
				}
				if (xik[d] < -(0.5 * L)) {
					xik[d] += L;
				}
				rik += xik[d] * xik[d];
			}
			rik= sqrt(rik);

			//culc diam difference
			aik = (p[i].getDiam() + p[k].getDiam()) / 2.0;
			ar = aik / rik;
			if (rik < 3 * aik) {
				      //(aik/rik)^12
				U += (ar*ar*ar*ar*ar*ar*ar*ar*ar*ar*ar*ar - 1/(3*3*3*3*3*3*3*3*3*3*3*3)) / N; //including Cutoff
			}
		}
	}

	return U;
}
inline double System::getK() {

	double K = 0;

	//m = 1
	for (int i = 0; i < N; i++) {
		K += p[i].getvv() / (2. * N);
	}

	return K;
}

inline void System::periodic(int i){
	//If pi is outside of the box, this func fix the pos of pi

	for (char d = 0; d < D; d++){
		if (p[i].getx(d) > L){
			p[i].resetx(d, p[i].getx(d) - L);
		}
		else if (p[i].getx(d) < 0){
			p[i].resetx(d, p[i].getx(d) + L);
		}
	}

	return;
}
void System::culc_Interaction(){
	//culc the every interactions between pi and pk

	double xik[D];
	//double xik, yik;
	double rik2; //rik^2
	double aik2; //aik^2
	double ar2; //(aik/rik)^2
	double f_rik; // f/rik
	int k;

	//Init force
	for (int n = 0; n < N; n++){
		for (int d = 0; d < D; d++){
			p[n].force[d] = 0.0;
		}
	}

	//culc_Interaction between i and j
	int i, j;
	for (i = 0; i < N; i++){
		for (j = 1; j <= list[i][0]; j++){

			k = list[i][j];

			//Position Difference i - j
			rik2 = 0;
			for (char d = 0; d < D; d++){
				xik[d] = p[k].getx(d) - p[i].getx(d);


				//Periodic between 2 particles;
				if (2 * xik[d] >  L)
					xik[d] -= L;
				if (2 * xik[d] < -L)
					xik[d] += L;

				rik2 += xik[d] * xik[d];
			}

			//culc Diameter Average^2
			aik2 = (p[i].getDiam() + p[k].getDiam()) / 2.;
			aik2 *= aik2;

			//culc_force
			ar2 = aik2 / rik2;

			//if rik < 3 * aik, continue culcration
			if (1 < 9 * ar2){
				            //(aik/rik)^12
				f_rik = -12 * (ar2*ar2*ar2*ar2*ar2*ar2)/rik2;
				// f/rik = -12 * (aik/rik)^12/rik^2

				for (int d = 0; d < D; d++) {
					//f_xik = f * xik/rik = f * xik * rik / (rik^2)

					p[i].force[d] += f_rik * xik[d];
					p[k].force[d] -= f_rik * xik[d];
				}
			}
		}
	}

	return;
}
inline void System::tEvoLD() {
	//temporal development of the system

	culc_Interaction();

	for (int i = 0; i < N; i++) {
		p[i].vDvlpBD(dt, thermalFuctor);
	}

	for (char d = 0; d < D; d++) {
		double vtmp = getvg(d);
		for (int n = 0; n < N; n++) {
			p[n].resetv(d, p[n].getv(d) - vtmp);
		}
	}
	for (int i = 0; i < N; i++) {
		p[i].xDvlp(dt);
		periodic(i);
	}
	judgeUpdateCell();
	return;
}
void System::culc_harmonicInteraction() {
	//culc the every interactions between pi and pk by using harmonic potential

	double xik[D];
	double rik2; //rik^2
	double aik;
	double aik2; //aik^2
	double f_rik; // f/rik
	int k;

	//Init force
	for (int n = 0; n < N; n++) {
		for (int d = 0; d < D; d++) {
			p[n].force[d] = 0.0;
		}
	}

	//culc_Interaction between i and j
	for (int i = 0; i < N; i++) {
		for (int j = 1; j <= list[i][0]; j++) {

			k = list[i][j];

			//Position Difference i - j
			rik2 = 0;
			for (char d = 0; d < D; d++) {
				xik[d] = p[k].getx(d) - p[i].getx(d);

				//Periodic between 2 particles;
				if (xik[d] > (0.5 * L)) {
					xik[d] -= L;
				}
				if (xik[d] < -(0.5 * L)) {
					xik[d] += L;
				}

				rik2 += xik[d] * xik[d];
			}

			//culc Diameter Average
			aik = (p[i].getDiam() + p[k].getDiam()) / 2.;
			aik2 = aik * aik;
			
			//culc_force
			if (rik2 < aik2) {
				// f/rik = -(1 - rik/aik)/(aik*rik)
				f_rik = 50 * (1/aik2 - 1 /(aik*sqrt(rik2)));
				

				for (int d = 0; d < D; d++) {
					p[i].force[d] += f_rik * xik[d];
					p[k].force[d] -= f_rik * xik[d];
				}
			}
		}
	}

	return;
}
inline void System::tHarmonicEvo() {
	//temporal development of the system by using harmonic potential

	culc_harmonicInteraction();

	for (int i = 0; i < N; i++) {
		p[i].vDvlpBD(dt, 0);
	}
	double vtmp;
	for (char d = 0; d < D; d++) {
		vtmp = getvg(d);
		for (int n = 0; n < N; n++) {
			p[n].resetv(d, p[n].getv(d) - vtmp);
		}
	}
	for (int i = 0; i < N; i++) {
		p[i].xDvlp(dt);
		periodic(i);
	}
	judgeUpdateCell();
	return;
}

void System::updateCell2D(){
	//determine positionMemory[D][N]
	//determine list[N][Num of pj around pi]

	int nx[D];
	double xij[D];
	double rrij;
	short Mtmp = getM();
	int M2 = Mtmp * Mtmp; //M^2

	//map[Num of cell][Num of cell list]
	int** map = new int*[M2];
	for (int i = 0; i < M2; i++){
		map[i] = new int[(int)(1.3 * 9 * N / M2)];
	}

	for (int i = 0; i < M2; i++){
		map[i][0] = 0;
	}

	for (int n = 0; n < N; n++){
		for (char d = 0; d < D; d++){
			nx[d] = f((int)(p[n].getx(d) * Mtmp / L), Mtmp);
		}
		for (int j = nx[1] - 1; j <= nx[1] + 1; j++) {
			for (int i = nx[0] - 1; i <= nx[0] + 1; i++) {
				map[f(i, Mtmp) + Mtmp*f(j, Mtmp)][0]++;
				map[f(i,Mtmp) + Mtmp*f(j, Mtmp)][map[f(i, Mtmp) + Mtmp*f(j, Mtmp)][0]] = n;
			}
		}
	}
	int j;
	for (int i = 0; i < N; i++){
		list[i][0] = 0;
		for (char d = 0; d < D; d++){
			nx[d] = f((int)(p[i].getx(d) * Mtmp / L), Mtmp);
		}

		for (int k = 1; k <= (map[nx[0] + Mtmp*nx[1]][0]); k++){
			j = map[nx[0] + Mtmp*nx[1]][k];

			if (j > i){
				rrij = 0;
				for (char d = 0; d < D; d++){
					xij[d] = p[i].getx(d) - p[j].getx(d);
					if (xij[d] < -(L / 2.0)){
						xij[d] += L;
					}
					else if (xij[d] > (L / 2.0)){
						xij[d] -= L;
					}
					rrij += xij[d] * xij[d];
				}
				if (rrij < (L*L) / M2){
					list[i][0]++;
					list[i][list[i][0]] = j;
				}
			}
		}
	}

	for (int i = 0; i < M2; i++){
		delete[] map[i];
	}
	delete[] map;
	for (int n = 0; n < N; n++){
		for (char d = 0; d < D; d++){
			positionMemory[d][n] = p[n].getx(d);
		}
	}

	return;
}
void System::judgeUpdateCell(){
	//if pi pos is far from pi0, activate updateCell()

	double dr2; //dr^2
	double dx[D];
    double delta2 = ((L/getM())-3*a0) / 2;
	delta2 *= delta2; //delta^2
	for (int n = 1; n < N; n++){
		dr2 = 0;
		
		for (char d = 0; d < D; d++){
			dx[d] = positionMemory[d][n] - p[n].getx(d);
			//Periodic between 2 particles;
			if (dx[d] >(0.5*L)){
				dx[d] -= L;
			}
			if (dx[d] < -(0.5*L)){
				dx[d] += L;
			}
			dr2 += dx[d] * dx[d];
		}
		if (dr2 >= delta2){
			updateCell2D();
			return;
		}
	}
	return;
}
inline void System::equilibrateSys(double teq) {
	cout << "Equilibrate the System: ID = " << id << endl;
	cout << "Time = " << teq << endl;
	unsigned int s = teq/dt;
	for (unsigned int nt = 0; nt < s; nt++) {
		tEvoLD();
	}
	cout << " -> Edone: ID = "<< id << endl;
	return;
}

void System::makeInitPosition() {

	cout << "Make InitPosition: ID = " << id << endl;

	double xvtmp[2 * D];

	//set velosity = 0
	for (char d = 0; d < D; d++) {
		xvtmp[d + D] = 0;
	}

	//avoiding super overraps
	double Ltmp = L - (a1 + a2) / 2;

	//set position
	for (int n = 0; n < N; n += 2) {
		for (char d = 0; d < D; d++) {
			xvtmp[d] = Ltmp * genrand_real3();
		}
		p[n].init(xvtmp, a1);
	}
	for (int n = 1; n < N; n += 2) {
		for (char d = 0; d < D; d++) {
			xvtmp[d] = Ltmp * genrand_real3();
		}
		p[n].init(xvtmp, a2);
	}

	struct stat st;
    std::ostringstream diamName;
    diamName << "../traj/N" << N << "/diam.data";
    if(stat(diamName.str().c_str(), &st) != 0){
        std::ofstream diamFile;
        diamFile.open(diamName.str().c_str());
        std::cout << "created " << diamName.str() << std::endl;
        for(unsigned int n = 0; n < N; n++){
            diamFile << p[n].getDiam() << std::endl;
	    }
        diamFile.close();
    }

	//set posMem and list
	updateCell2D();

	//remove overraps
	for (int nt = 0; nt < (int)(20 / dt); nt++) {
		tHarmonicEvo();
	}
	judgeUpdateCell();

	cout << " -> MIdone: ID = "<< id << endl;
	return;
}

void System::recPos(std::ofstream *of) {
	for (int n = 0; n < N; n++) {
		for (char d = 0; d < D; d++) {
			*of << p[n].getx(d) << " ";
		}
	}
	*of << endl;

	return;
}

//public
System::System(int ID) {

	//default setting
	p = new Particle[N];

	id = ID;
	dt = dt_MD;
	t = 0;
	T = Tfin;
	Eav = 0;
	thermalFuctor = sqrt(2 * T / dt);
	L = sqrt(N / dnsty);

	LDDir = "/LD";
    MDDir = "/MD";
    EDir = "/E";
    posDir = "/pos";
	velDir = "/vel";

	//for list
	for (char d = 0; d < D; d++) {
		positionMemory[d] = new double[N];
	}

	int Mtmp = (int)(sqrt(N / dnsty) / (6 * a0)); // M ~= L/(2*a_cut)
	cout << "M = " << Mtmp << endl;
	//list[Num of Particle][Num of list]
	list = new int* [N];
	for (int i = 0; i < N; i++) {
		list[i] = new int[(int)(1.5 * N * pi / (Mtmp * Mtmp))];
	}

	//Make dir tree
    std::ostringstream trajName;
    struct stat st;
    trajName << "../traj/N" << N;
    if(stat(trajName.str().c_str(), &st) != 0){
        mkdir(trajName.str().c_str(), 0755);
        std::cout << "created " << trajName.str() << std::endl;
    }
    trajName << "/T" << Tfin;
    NTDir = trajName.str();
    if(stat(NTDir.c_str(), &st) != 0){
        mkdir(NTDir.c_str(), 0755);
        std::cout << "created " << NTDir << std::endl;

        mkdir((NTDir + LDDir).c_str(), 0755);
        mkdir((NTDir + LDDir + EDir).c_str(), 0755);
        mkdir((NTDir + LDDir + posDir).c_str(), 0755);
        mkdir((NTDir + LDDir + velDir).c_str(), 0755);
        
        mkdir((NTDir + MDDir).c_str(), 0755);
        mkdir((NTDir + MDDir + EDir).c_str(), 0755);
        mkdir((NTDir + MDDir + posDir).c_str(), 0755);
        mkdir((NTDir + MDDir + velDir).c_str(), 0755);
    }

	cout << "Created System: ID = " << id << endl;

	return;
}
System::~System() {

	for (int i = 0; i < N; i++) {
		delete[] list[i];
	}
	delete[] list;

	for (char d = 0; d < D; d++) {
		delete[] positionMemory[d];
	}

	delete[] p;
	cout << "Finish! ID = " << id << endl << endl;
}

void System::initSys() {

	cout << "Start Initialisation: ID = " << id << endl;

	setdt_T(dt_init, Tfin);

	makeInitPosition();

	setdt_T(dt_BD, Tfin);

	equilibrateSys(tmax);
	setdt_T(dt_BD, Tfin);
	cout << "End Initialisation: ID = " << id << endl;

	return;
}

void System::getDataLD() {

	std::cout << "Starting LD time loop: ID = " << id << std::endl;
    unsigned int Nt;
    unsigned int ntAtOutput;

	std::ofstream tFile;
	std::ofstream eFile;
	std::ofstream posFile;

	if(id == 1){
		std::cout << "getting liniarPlot datas in 5 secs" << std::endl;

		std::string tLinpltName = "/tliniar.data";
        tFile.open((NTDir + LDDir + tLinpltName).c_str());

        std::ostringstream eLinpltName;
        eLinpltName << NTDir + LDDir + EDir << "/liniar.data";
        eFile.open(eLinpltName.str().c_str());

        std::ostringstream posLinpltName;
        posLinpltName << NTDir + LDDir + posDir << "/liniar.data";
        posFile.open(posLinpltName.str().c_str());

        Nt = 5./dt;
        ntAtOutput = 0;
        for(unsigned int nt = 0; nt < Nt; nt++){
            tEvoLD();
            if(nt >= ntAtOutput){
                if(id == 1){
                    tFile << nt * dt << std::endl;
                }
                eFile << getK() << " " << getU() << " " << std::endl;
                recPos(&posFile);
                ntAtOutput += 0.1/dt;
            }
        }
        posFile.close();
        eFile.close();
        tFile.close();
        std::string tLogpltName = "/tlog.data";
        tFile.open((NTDir + LDDir + tLogpltName).c_str());
	}

	std::cout << "getting logPlot datas" << std::endl;
    std::ostringstream eLogpltName;
    eLogpltName << NTDir + LDDir + EDir << "/id" << id << ".data";
    eFile.open(eLogpltName.str().c_str());

    std::ostringstream posLogpltName;
    posLogpltName << NTDir + LDDir + posDir << "/id" << id << ".data";
    posFile.open(posLogpltName.str().c_str());

    Nt = tmax/dt;
    ntAtOutput = 10;
    unsigned int ntAtTakingAverage = 0;
    unsigned int NextNtAtTakingAverage = Nt>>7;
    if(NextNtAtTakingAverage == 0){
        NextNtAtTakingAverage = 1;
    }
    unsigned int numOfEnsemble = 0;
    Eav = 0;
    for(unsigned int nt = 0; nt <= Nt; nt++){
        tEvoLD();
        if(nt >= ntAtOutput){
            if(id == 1){
                tFile << nt * dt << std::endl;
            }
            eFile << getK() << " " << getU() << std::endl;
            recPos(&posFile);
            ntAtOutput *= 1.3;
        }
        if(nt >= ntAtTakingAverage){
            Eav += getK() + getU();
            numOfEnsemble++;
            ntAtTakingAverage += NextNtAtTakingAverage;
        }
    }
    Eav /= (float)numOfEnsemble;
    std::cout <<"Ensemble: " <<numOfEnsemble << ", Eav = " << Eav << std::endl;
    if(id == 1){
        tFile.close();
    }
    eFile.close();
    posFile.close();
    std::cout << "Every LD steps have been done: ID = " << id << std::endl;
    return;
}
void System::benchmark(unsigned int loop){
	for(unsigned int nt = 0; nt < loop; nt++) {
		tEvoLD();
	}
}