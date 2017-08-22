#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

const int nz = 1000; //! 解析領域の分割数
const int nstep = 2000; //! 計算のステップ数

const std::string filename = "fdtd_data.dat";

const double c = 2.9979246e8; //! 光速
const double eps0 = 8.8541878e-12; //! 真空の誘電率
const double mu0 = 1.256637e-6; //! 真空の透磁率
const double z0 = 376.73031; //! 真空のインピーダンス

std::vector<double> ex(nz+1,0.0), hy(nz+1,0.0); //! 電界, 磁界を記憶する配列
std::vector<double> ae(nz + 1, 0.0), be(nz + 1, 0.0); //! 電界磁界更新式の係数
double bm; //! 電界磁界更新式の係数
double z, dz, dt, t; //! 位置, セルサイズ, 時間ステップ, 時刻
double exlold, exrold; //! k=1, K=nz-1の1ステップ前の電界
double czl, czr; //! Murの係数
double zp, a; //! 励振パルスのパラメータ
double d, epsr, eps, sgm, a0, b, v;
std::vector<double> epsd(nz + 1, 0.0), sgmd(nz + 1, 0.0);
int nd, kd;

void setup(void);
double pulse(double z, double tm);
void e_cal(void);
void mur(void);
void h_cal(void);

int main() {

	std::ofstream outputfile(filename);

	setup();

	t = dt;

	for (int n = 1; n <= nstep; n++) {
		e_cal();
		mur();
		t += 0.5*dt;
		h_cal();
		t += 0.5*dt;

		if (n % 50 == 0) {
			for (int k = 0; k <= nz; k++) {
				if (k % 10 == 0) {
					z = k*dz;
					outputfile << z << " " << ex[k] << std::endl;
				}
			}
			outputfile << std::endl << std::endl;
			std::cout << n << "/" << nstep << std::endl;
		}
		
	}

	outputfile.close();

	return 0;
}

void setup(void) {
	d = 0.1;
	nd = 50;
	dz = d / nd;
	epsr = 3.0;
	kd = nz / 2;

	for (int k = 0; k <= nz; k++) {
		if ((k < kd) || (k >= kd + nd)) {
			epsd[k] = 1.0;
			sgmd[k] = 0.0;
		}
		else {
			epsd[k] = epsr;
			sgmd[k] = 0.0;
		}
	}

	dt = dz / c;

	for (int k = 1; k <= nz-1; k++) {
		eps = 0.5*(epsd[k]+epsd[k-1])*eps0;
		sgm = 0.5*(sgmd[k] + sgmd[k - 1]);
		b = dt / eps;
		a0 = 0.5*sgm*b;
		ae[k] = (1.0 - a0) / (1.0 + a0);
		be[k] = b / (1.0 + a0) / dz;
	}

	bm = dt / mu0 / dz;

	zp = 100.0*dz;
	a = 20.0*dz;

	for (int k = 0; k <= nz; k++) {
		z = k * dz;
		ex[k] = pulse(z, 0.0);
	}

	for (int k = 0; k <= nz; k++) {
		z = (k+0.5) * dz;
		hy[k] = pulse(z, 0.5*dt)/z0;
	}

	exlold = ex[1];
	exrold = ex[nz - 1];
	v = c / sqrt(epsd[0]);
	czl = (v*dt - dz) / (v*dt+dz);
	v = c / sqrt(epsd[nz-1]);
	czr= (v*dt - dz) / (v*dt + dz);
}

double pulse(double z, double tm) {
	return exp(-1.0*pow((z-zp-c*tm)/a,2));
}

void e_cal(void) {
	for (int k = 1; k <= nz-1; k++) {
		ex[k] = ae[k] * ex[k] - be[k] * (hy[k]-hy[k-1]);
	}
}

void mur(void) {
	ex[0] = exlold + czl*(ex[1]-ex[0]);
	ex[nz] = exrold + czr*(ex[nz-1] - ex[nz]);
	exlold = ex[1];
	exrold = ex[nz-1];
}

void h_cal(void) {
	for (int k = 1; k <= nz - 1; k++) {
		hy[k] = hy[k] - bm * (ex[k+1] - ex[k]);
	}
}