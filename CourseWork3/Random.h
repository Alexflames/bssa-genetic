#pragma once
#include <vector>
#include <ctime>
#include <algorithm>

using std::vector;
using std::exception;
using std::max;
using std::min;

class Random {
private:
	static const int a = 10'000;
public:
	Random() {
		srand(time(NULL));
	}

	int choice(vector<int> v) {
		int i = rand() % int(v.size());
		return v[i];
	}

	static double getRandomLowerOne() {
		double x = rand() % a;
		return x / double(a);
	}

	template<typename T>
	T choice(vector<T> v, vector<double> d) {
		if (v.size() != d.size()) {
			throw exception("The size of vectors must be equal");
		}
		double x = getRandomLowerOne();
		for (int i = 0; i < v.size(); i++) {
			x -= d[i];
			if (x <= 0) {
				return v[i];
			}
		}
		return v[v.size() - 1];
	}

	static int randint(int a, int b) {
		int maxV = max(a, b), minV = min(a, b);
		int x = rand() % (maxV - minV + 1);
		return minV + x;
	}
	
};

