#pragma once
#include <vector>
#include"Random.h"

using std::vector;
using std::pair;
using std::make_pair;

double getProb(int n) {
	double c = 1.3;
	//return c * log(n) / double(n);
	return 0.32;
}

vector<pair<int, int>> GnpRandomGraph(int n) {
	vector<pair<int, int>> e;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			double x = getProb(n);
			if (Random::getRandomLowerOne() < x) {
				e.push_back(make_pair(i, j));
			}
		}
	}
	return e;
}