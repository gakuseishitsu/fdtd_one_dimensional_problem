#include <iostream>
#include <vector>

int nz, nstep; //! 解析領域の分割数, 計算のステップ数
std::vector<double> ex, hy; //! 電界, 磁界を記憶する配列
std::vector<double> ae, be; //! 電界磁界更新式の係数
double bm; //! 電界磁界更新式の係数
double c, dz, dt, t; //! 光速, セルサイズ, 時間ステップ, 時刻
double exload, exrold; //! k=1, K=nz-1の1ステップ前の電界
double czl, czr; //! Murの係数
double zp, a; //! 励振パルスのパラメータ

int main() {

	std::cout << "hello world" << std::endl;

	c = 2.9979246e8;

	nz = 1000;
	nstep = 2000;

	ex.resize(nz);
	hy.resize(nz);

	return 0;
}