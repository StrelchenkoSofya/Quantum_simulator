#pragma once
#include<math.h>
#include <stdlib.h>
#include <time.h>
#include<iostream>
#include<vector>
#include<complex>
#include<omp.h>

using namespace std;

class Quantum {
	vector<complex<double>> q;
	size_t size;
	size_t rl_sz;
public:
	Quantum();
	Quantum(size_t n);

	size_t getSize();
	size_t getRlSize();
	static void setNumThreads(int n);

	void CNOT(size_t ind1, size_t ind2);
	void CP(size_t ind1, size_t ind2, double r);
	void I();
	void H(size_t ind);
	void H_all();
	void X(size_t ind);
	void Y(size_t ind);
	void Z(size_t ind);
	void P(double r, size_t ind);
	void Rx(double r, size_t ind);
	void Ry(double r, size_t ind);
	void Rz(double r, size_t ind);
	void controlled_swap(size_t control, size_t target1, size_t target2, size_t ancilla);

	vector<size_t> measure(size_t count);

	friend std::ostream& operator<<(std::ostream& ostr, const Quantum& other);
};

void QFT_once(Quantum& q, size_t start, size_t fin);
void QFT(Quantum& q);
void QFT(Quantum& q, size_t st, size_t fin);
void IQFT_once(Quantum& q, size_t start, size_t fin);
void IQFT(Quantum& q);
void IQFT(Quantum& q, size_t st, size_t fin);
void add(Quantum& q, size_t st_a, size_t fin_a, size_t st_b, size_t fin_b);
void add(Quantum& q, int a, size_t st, size_t fin);
void controlled_add(Quantum& q, int a, size_t control, size_t st, size_t fin);
void add_mod(Quantum& q, int a, int N, size_t st, size_t fin, size_t ancilla);
void c2_add_mod(Quantum& q, int a, int N, size_t control1, size_t control2, size_t st, size_t fin, size_t ancilla);
void cmult(Quantum& q, int a, int N, size_t control, size_t st_a, size_t fin_a, size_t st_b, size_t fin_b, size_t ancilla);
int mod_inverse(int a, int N);
void c_ua(Quantum& q, int a, int N, size_t control, size_t st_x, size_t fin_x, size_t st_ancilla, size_t fin_ancilla, size_t temp_ancilla, size_t work_ancilla);