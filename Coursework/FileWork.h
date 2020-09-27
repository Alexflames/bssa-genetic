#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>

using std::vector;
using std::pair;
using std::make_pair;
using std::string;
using std::ifstream;
using std::exception;
using std::getline;
using std::stringstream;
using std::max;

vector<pair<int, int>> readFromFileWhereEdges(string file, int& N, int& M) {
	ifstream in(file);
	if (!in) {
		throw exception("Проблема с файлом");
	}
	vector<pair<int, int>> e;
	string inputLine;
	int n = 0, m = 0;
	while (getline(in, inputLine)) {
		if (inputLine[0] == '%') {
			continue;
		}
		m++;
		stringstream s(inputLine);
		int x, y;
		s >> x >> y;
		x--, y--;
		n = max(n, max(x, y));
		e.push_back(make_pair(x, y));
	}
	N = n;
	M = m;
	return e;
}

vector<pair<int, int>> readFromFileNM(string file, int& N, int& M) {
	ifstream in(file);
	if (!in) {
		throw exception("Проблема с файлом");
	}
	vector<pair<int, int>> e;
	int n, m;
	in >> n >> m;
	for (int i = 0; i < m; i++) {
		int x, y;
		in >> x >> y;
		e.push_back(make_pair(x, y));
	}
	N = n;
	M = m;
	return e;
}