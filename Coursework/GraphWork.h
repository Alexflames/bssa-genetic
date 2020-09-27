#pragma once
#include <vector>
#include <queue>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <string>

using std::vector;
using std::queue;
using std::pair;
using std::exception;
using std::string;
using std::to_string;
using std::cin;
using std::cout;
using std::endl;

using std::max;

class Graph {
private:
	int N, M;
	vector<vector<int>> adjacency;
	vector<vector<int>> distance;
	vector<vector<vector<int>>> path;
	vector<bool> calculatedEccentricity;
	vector<int> eccentricity;
	vector<bool> subGraph;
	bool subGraphMode = false;
	void initializeMatrix() {
		adjacency.resize(this->N, vector<int>());
		distance.resize(this->N, vector<int>(this->N, INT_MAX));
		path.resize(this->N, vector<vector<int>>(this->N));
		calculatedEccentricity.resize(this->N, false);
		eccentricity.resize(this->N, INT_MAX);
		
		for (int i = 0; i < this->N; i++) {
			this->distance[i][i] = 0;
			this->path[i][i] = vector<int>(1, i);
		}


	}

public:
	Graph(int n, int m, vector<pair<int, int>> adj) {
		this->N = n, this->M = m;
		initializeMatrix();
		for (int i = 0; i < m; i++) {
			int x = adj[i].first, y = adj[i].second;
			this->adjacency[x].push_back(y);
			this->adjacency[y].push_back(x);
		}
	}

	void bfsFromVertex(int x) {
		if (x < 0 || x >= this->N) {
			throw exception("Vertex doesn't exist");
		}
		if (this->subGraphMode && !this->subGraph[x]) {
			throw exception("Graph is subgraph mode and vertex isn't in subgraph");
		}
		queue<int> q;
		q.push(x);
		vector<int> d(this->N, INT_MAX);
		d[x] = 0;
		while (!q.empty()) {
			int u = q.front();
			q.pop();
			vector<int> currentPath = this->path[x][u];
			for (unsigned int i = 0; i < this->adjacency[u].size(); i++) {
				int v = this->adjacency[u][i];
				if (this->subGraphMode && !this->subGraph[v]) {
					continue;
				}
				if (d[v] > INT_MAX / 2) {
					currentPath.push_back(v);
					this->path[x][v] = currentPath;
					this->path[v][x] = currentPath;
					currentPath.pop_back();
					d[v] = d[u] + 1;
					q.push(v);
				}
			}
		}
		for (unsigned int i = 0; i < d.size(); i++) {
			this->distance[x][i] = d[i];
			this->distance[i][x] = d[i];
		}
	}
	
	void printDistance() {
		cout << '\t';
		for (int i = 0; i < this->N; i++) {
			cout << i << '\t';
		}
		cout << std::endl;
		for (int i = 0; i < this->N; i++) {
			cout << i << '\t';
			for (int j = 0; j < this->N; j++) {
				if (this->distance[i][j] > INT_MAX / 2) {
					cout << "INF";
				}
				else {
					cout << this->distance[i][j];
				}
				cout << '\t';
			}
			cout << std::endl;
		}
	}

	void printPath(int x, int y) {
		if (this->path[x][y].empty()) {
			return;
		}
		cout << x << " ~~~> " << y << ": ";
		for (unsigned int i = 0; i < this->path[x][y].size(); i++) {
			cout << path[x][y][i] << " ";
		}
		cout << std::endl;
	}

	vector<int> getPath(int x, int y) {
		if (this->path[x][y].empty()) {
			this->bfsFromVertex(x);
		}
		return this->path[x][y];
	}

	int getEccentricity(int x) {
		if (x >= this->N) {
			throw exception("Vertex doesn't exist in graph");
		}

		if (!this->calculatedEccentricity[x]) {
			this->bfsFromVertex(x);
			this->calculatedEccentricity[x] = true;
			int result = -INT_MAX;
			for (int i = 0; i < this->N; i++) {
				result = max(this->distance[x][i], result);
			}
			this->eccentricity[x] = result;
		}
		return this->eccentricity[x];
	}

	vector<int> getNeighbour(int x) {
		return this->adjacency[x];
	}

	int Size() {
		return this->N;
	}

	void setSubGraph(vector<int> v)  {
		this->subGraph.resize(this->N);
		for (unsigned int i = 0; i < v.size(); i++) {
			this->subGraph[v[i]] = true;
		}
		this->subGraphMode = true;
	}
	void resetSubGraph() {
		this->subGraph.clear();
		this->subGraphMode = false;
	}
	int getDistance(int u, int v) {
		return this->distance[u][v];
	}

	// This can make some error
	void setDistance(int u, int v, int d) {
		this->distance[u][v] = d;
		this->distance[v][u] = d;
	}

	void printDistanceMatrix() {
		for (int i = 0; i < this->distance.size(); i++) {
			for (int j = 0; j < this->distance[i].size(); j++) {
				cout << this->distance[i][j] << " ";
			}
			cout << endl;
		}
	}
	vector<vector<int>> getDistanceMatrix() {
		return this->distance;
	}

};



