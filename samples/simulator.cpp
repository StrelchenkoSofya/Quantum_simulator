#include <iostream>
#include<cmath>
#include <chrono>
#include "Quantum.h"

int main() {
	for (size_t n = 1; n <= 12; ++n) {
		Quantum::setNumThreads(n);
		int s = 2;
		int i = 2;
		cout << "Proccessors " << n << endl;
		for (size_t s = 2; s <= 25; s += i) {
			Quantum q(s);
			auto start = std::chrono::steady_clock::now();
			QFT(q);
			add(q, 500000, 0, s - 1);
			IQFT(q);
			vector<size_t> m = q.measure(1000);
			auto end = std::chrono::steady_clock::now();
			chrono::duration<double> elapsed = end - start;
			cout << s << ": " << elapsed.count() << endl;
			if (s == 14)i = 1;
		}
	}

}