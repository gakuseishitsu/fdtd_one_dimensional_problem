#include <iostream>
#include <vector>

int nz, nstep; //! ðÍÌæÌª, vZÌXebv
std::vector<double> ex, hy; //! dE, ¥EðL¯·ézñ
std::vector<double> ae, be; //! dE¥EXV®ÌW
double bm; //! dE¥EXV®ÌW
double c, dz, dt, t; //! õ¬, ZTCY, ÔXebv, 
double exload, exrold; //! k=1, K=nz-1Ì1XebvOÌdE
double czl, czr; //! MurÌW
double zp, a; //! ãUpXÌp[^

int main() {

	std::cout << "hello world" << std::endl;

	c = 2.9979246e8;

	nz = 1000;
	nstep = 2000;

	ex.resize(nz);
	hy.resize(nz);

	return 0;
}