#include <iostream>
#include <vector>

int nz, nstep; //! ��͗̈�̕�����, �v�Z�̃X�e�b�v��
std::vector<double> ex, hy; //! �d�E, ���E���L������z��
std::vector<double> ae, be; //! �d�E���E�X�V���̌W��
double bm; //! �d�E���E�X�V���̌W��
double c, dz, dt, t; //! ����, �Z���T�C�Y, ���ԃX�e�b�v, ����
double exload, exrold; //! k=1, K=nz-1��1�X�e�b�v�O�̓d�E
double czl, czr; //! Mur�̌W��
double zp, a; //! ��U�p���X�̃p�����[�^

int main() {

	std::cout << "hello world" << std::endl;

	c = 2.9979246e8;

	nz = 1000;
	nstep = 2000;

	ex.resize(nz);
	hy.resize(nz);

	return 0;
}