#include"Quantum.h"

Quantum::Quantum() : size(1), rl_sz(2), q(vector<complex<double>>(2)) {
	q[0] = complex<double>(1, 0);
#pragma omp parallel for
	for (int i = 1; i < rl_sz; ++i) {
		q[i] = complex<double>(0, 0);
	}
}
Quantum::Quantum(size_t n) : size(n), rl_sz(1ULL << n), q(vector<complex<double>>(1ULL << n)) {
	q[0] = complex<double>(1, 0);
#pragma omp parallel for
	for (int i = 1; i < rl_sz; i++) {
		q[i] = complex<double>(0, 0);
	}
}

size_t Quantum::getSize() {
	return size;
}

size_t Quantum::getRlSize() {
	return rl_sz;
}

void Quantum::setNumThreads(int n) {
	if (n > 0) {
		omp_set_num_threads(n);
	}
}

void Quantum::CNOT(size_t ind1, size_t ind2) {
	if (ind1 >= size)throw - 1;
	if (ind2 >= size)throw - 1;
	int mask1 = 1ULL << ind1;
	int mask2 = 1ULL << ind2;
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		complex<double> tmp0;
		if ((i & mask1) && (!(i & mask2))) {
			tmp0 = q[i];
			int j = i;
			while (j != (i | mask2)) j++;
			q[i] = q[j];
			q[j] = tmp0;
		}
	}
}
void Quantum::CP(size_t ind1, size_t ind2, double r) {
	if (ind1 >= size)throw - 1;
	if (ind2 >= size)throw - 1;
	int mask1 = 1ULL << ind1;
	int mask2 = 1ULL << ind2;
	complex<double> e(cos(r), sin(r));
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if ((i & mask1) && (i & mask2)) {
			q[i] = q[i] * e;
		}
	}
}
void Quantum::I() {
	return;
}
void Quantum::H(size_t ind) {
	if (ind >= size)throw - 1;
	vector<complex<double>> copy = q;
	int m = 1ULL << ind;
	double s2 = sqrt(2);
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if (i & m) {
			q[i] = (copy[i - m] - copy[i]) / s2;
			q[i - m] = (copy[i - m] + copy[i]) / s2;
		}
	}
	return;
}
void Quantum::H_all() {
#pragma omp parallel for
	for (int i = 0; i < size; ++i) {
		H(i);
	}
}
void Quantum::X(size_t ind) {
	if (ind >= size)throw - 1;
	vector<complex<double>> copy = q;
	int m = 1ULL << ind;
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if (i & m) {
			q[i] = copy[i - m];
			q[i - m] = copy[i];
		}
	}
	return;
}

void Quantum::Y(size_t ind) {
	if (ind >= size)throw - 1;
	vector<complex<double>> copy = q;
	int m = 1ULL << ind;
	complex<double> cmp = complex<double>(0, 1.0);
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if (i & m) {
			q[i] = copy[i - m] * cmp * (-1.0);
			q[i - m] = copy[i] * cmp;
		}
	}
	return;
}

void Quantum::Z(size_t ind) {
	if (ind >= size)throw - 1;
	vector<complex<double>> copy = q;
	int m = 1ULL << ind;
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if (i & m)q[i] = q[i] * (-1.0);
	}
	return;
}

void Quantum::P(double r, size_t ind) {
	if (ind >= size)throw - 1;
	complex<double> e(cos(r), sin(r));
	int m = 1ULL << ind;
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if (i & m) {
			q[i] = q[i] * e;
		}
	}
	return;
}

void Quantum::Rx(double r, size_t ind) {
	if (ind >= size)throw - 1;
	vector<complex<double>> copy = q;
	int m = 1ULL << ind;
	double c = cos(r / 2);
	double s = sin(r / 2);
	complex<double> cmp = complex<double>(0, 1.0);
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if (i & m) {
			q[i] = copy[i] * c - copy[i - m] * s * cmp;
			q[i - m] = copy[i - m] * c - copy[i] * s * cmp;
		}
	}
	return;
}

void Quantum::Ry(double r, size_t ind) {
	if (ind >= size)throw - 1;
	vector<complex<double>> copy = q;
	int m = 1ULL << ind;
	double c = cos(r / 2);
	double s = sin(r / 2);
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if (i & m) {
			q[i] = copy[i] * c + copy[i - m] * s;
			q[i - m] = copy[i - m] * c - copy[i] * s;
		}
	}
	return;
}

void Quantum::Rz(double r, size_t ind) {
	if (ind >= size)throw - 1;
	if (ind >= size)throw - 1;
	int m = 1ULL << ind;
#pragma omp parallel for
	for (int i = 0; i < rl_sz; ++i) {
		if (i & m) q[i] = q[i] * complex<double>(cos(r / 2), sin(r / 2));
		else q[i] = q[i] * complex<double>(cos(- r / 2), sin(- r / 2));
	}
	return;
}
void Quantum::controlled_swap(size_t control, size_t target1, size_t target2, size_t ancilla) {
	CNOT(target2, ancilla);
	CNOT(target1, target2);
	CP(control, target1, 4.0 * atan(1));
	CP(control, target2, 4.0 * atan(1));
	CNOT(target1, target2);
	CNOT(target2, ancilla);
}

vector<size_t> Quantum::measure(size_t count) {
	vector<size_t> res(rl_sz);
	vector<double> brd(rl_sz);
	double tmp = 0;
	for (size_t i = 0; i < rl_sz; ++i) {
		res[i] = 0;
		double m = abs(q[i]);
		tmp += (m * m);
		brd[i] = tmp;
	}
	brd[rl_sz - 1] = 1;
#pragma omp parallel for
	for (int i = 0; i < count; ++i) {
		double rnd = (double)(rand()) / RAND_MAX;
		int k = 0;
		while (rnd > brd[k])k++;
		res[k]++;
	}
	/*q.clear();
	size = 0;
	rl_sz = 0;*/
	return res;
}

std::ostream& operator<<(std::ostream& ostr, const Quantum& other) {
	for (int i = 0; i < other.rl_sz; i++) {
		ostr << "(" << real(other.q[i]) << "," << imag(other.q[i]) << "); ";
	}
	return ostr;
}

void QFT_once(Quantum& q, size_t start, size_t fin) {
	size_t size = q.getSize();
	if (start >= size)throw - 1;
	if (fin >= size)throw - 1;
	if (start > fin)throw - 1;
	q.H(start);
	int st = 2;
	double pi = atan(1) * 4;
	for (size_t i = start + 1; i <= fin; ++i) {
		q.CP(start, i, pi / st);
		st = st << 1;
	}
}
void QFT(Quantum& q) {
	size_t size = q.getSize();
	for (size_t i = 0; i < size; ++i) {
		QFT_once(q, i, size - 1);
	}
}
void QFT(Quantum& q, size_t st, size_t fin) {
	for (size_t i = st; i <= fin; ++i) {
		QFT_once(q, i, fin);
	}
}
void IQFT_once(Quantum& q, size_t start, size_t fin) {
	size_t size = q.getSize();
	if (start >= size)throw - 1;
	if (fin >= size)throw - 1;
	if (start > fin)throw - 1;
	int st = 1 << (fin - start);
	double pi = atan(1) * 4 * (-1);
	for (int i = fin; i > start; --i) {
		q.CP(start, i, pi / st);
		st = st >> 1;
	}
	q.H(start);
}
void IQFT(Quantum& q) {
	size_t size = q.getSize();
	for (int i = (int)size - 1; i >= 0; --i) {
		IQFT_once(q, i, size - 1);
	}
}
void IQFT(Quantum& q, size_t st, size_t fin) {
	for (int i = fin; i >= (int)st; --i) {
		IQFT_once(q, i, fin);
	}
}


void add(Quantum& q, size_t st_a, size_t fin_a, size_t st_b, size_t fin_b) {
	if ((st_a > fin_a || st_b > fin_b) || (st_a >= st_b && st_a <= fin_b) || (st_b >= st_a && st_b <= fin_a))throw - 1;
	size_t sz = q.getSize();
	if (fin_a >= sz || fin_b >= sz)throw - 1;
	int n = 1ULL << (fin_b - st_b + 1);
	double two_pi = 8.0 * atan(1);
	for (size_t i = st_b; i <= fin_b; ++i) {
		for (size_t j = st_a; j <= fin_a; ++j) {
			q.CP(j, i, two_pi * (1ULL << (j - st_a)) / n);
		}
	}
}
void add(Quantum& q, int a, size_t st, size_t fin) {
	if (st > fin || fin >= q.getSize()) throw - 1;
	size_t n = fin - st + 1;
	unsigned long long zn = 1ULL << n;
	double two_pi = 8.0 * atan(1);
	for (size_t i = st; i <= fin; ++i) {
		double theta = two_pi * a * (1ULL << (i - st)) / zn;
		q.P(theta, i);
	}
}
void controlled_add(Quantum& q, int a, size_t control, size_t st, size_t fin) {
	if (st > fin || fin >= q.getSize() || control >= q.getSize()) throw - 1;
	unsigned long long zn = 1ULL << (fin - st + 1);
	double two_pi = 8.0 * atan(1);
	for (size_t i = st; i <= fin; ++i) {
		double theta = two_pi * a * (1ULL << (i - st)) / zn;
		q.CP(control, i, theta);
	}
}
void add_mod(Quantum& q, int a, int N, size_t st, size_t fin, size_t ancilla) {
	if (st > fin || fin >= q.getSize() || ancilla >= q.getSize()) throw - 1;
	size_t msb = fin;
	add(q, a, st, fin);
	add(q, -N, st, fin);
	IQFT(q, st, fin);
	q.CNOT(msb, ancilla);
	QFT(q, st, fin);
	controlled_add(q, N, ancilla, st, fin);
	IQFT(q, st, fin);
	q.CNOT(msb, ancilla);
	QFT(q, st, fin);
}
void c2_add_mod(Quantum& q, int a, int N, size_t control1, size_t control2, size_t st, size_t fin, size_t ancilla) {
	if (st > fin || fin >= q.getSize() || ancilla >= q.getSize() || control1 >= q.getSize() || control2 >= q.getSize())throw - 1;
	size_t msb = fin;
	double two_pi = 8.0 * atan(1);
	for (size_t i = st; i <= fin; ++i) {
		double t = two_pi * a * (1ULL << (i - st)) / (1ULL << (fin - st + 1));
		q.CP(control2, i, t);
		q.CNOT(control1, control2);
		q.CP(control2, i, -t);
		q.CNOT(control1, control2);
		q.CP(control1, i, t);
	}
	add(q, -N, st, fin);
	IQFT(q, st, fin);
	q.CNOT(msb, ancilla);
	QFT(q, st, fin);
	controlled_add(q, N, ancilla, st, fin);
	IQFT(q, st, fin);
	q.CNOT(msb, ancilla);
	QFT(q, st, fin);
}
void cmult(Quantum& q, int a, int N, size_t control, size_t st_a, size_t fin_a, size_t st_b, size_t fin_b, size_t ancilla) {
	size_t n = fin_a - st_a + 1;
	for (size_t i = 0; i < n; ++i) {
		c2_add_mod(q, (a * (1ULL << i)) % N, N, control, st_a + i, st_b, fin_b, ancilla);
	}
}

int mod_inverse(int a, int N) {
	int t = 0, t1 = 1;
	int r = N, r1 = a;
	while (r1 != 0) {
		int quot = r / r1;
		int temp_t = t - quot * t1;
		t = t1;
		t1 = temp_t;
		int temp_r = r - quot * r1;
		r = r1;
		r1 = temp_r;
	}
	if (r > 1) return -1;
	if (t < 0) t += N;
	return t;
}

void c_ua(Quantum& q, int a, int N,
	size_t control,
	size_t st_x, size_t fin_x,
	size_t st_ancilla, size_t fin_ancilla,
	size_t temp_ancilla,
	size_t work_ancilla) {
	size_t n = fin_x - st_x + 1;
	cmult(q, a, N, control, st_x, fin_x, st_ancilla, fin_ancilla, work_ancilla);
	for (size_t i = 0; i < n; ++i) {
		q.controlled_swap(control, st_x + i, st_ancilla + i, temp_ancilla);
	}
	int a_inv = mod_inverse(a, N);
	if (a_inv == -1) throw - 1;
	cmult(q, a_inv, N, control, st_x, fin_x, st_ancilla, fin_ancilla, work_ancilla);
}