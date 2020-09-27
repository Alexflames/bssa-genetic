#pragma once
#include <set>
#include <vector>

using std::set;
using std::vector;

template <typename T>
set<T> operator-(const set<T>& a, const set<T>& b) {
	set<T> result;
	for (auto it = a.begin(); it != a.end(); it++) {
		if (b.find(*it) == b.end()) {
			result.insert(*it);
		}
	}
	return result;
}

template <typename T>
vector<T> makeVector(const set<T>& a) {
	vector<T> result;
	for (auto it = a.begin(); it != a.end(); it++) {
		result.push_back(*it);
	}
	return result;
}

template <typename T>
set<T> makeSet(const vector<T>& a) {
	set<T> result;
	for (int i = 0; i < a.size(); i++) {
		result.insert(a[i]);
	}
	return result;
}


