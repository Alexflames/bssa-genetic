#pragma once

#include <vector>
#include <algorithm>
#include <set>
#include <time.h>

#include "GraphWork.h"
#include "Random.h"
#include "Strassen.h"


using std::begin;
using std::end;
using std::max_element;
using std::min_element;
using std::pair;
using std::vector;
using std::set;

//#define PRINT_GA
#define PRINT_TIME

vector<int> getCentralVertex(vector<int> e) {
	int minV = *min_element(e.begin(), e.end());
	#ifdef PRINT_GA
		cout << "R = " << minV << " " << "Vertex = ";
	#endif
	vector<int> res;
	for (int i = 0; i < e.size(); i++) {
		if (e[i] == minV) {
			res.push_back(i);
		}
	}
	return res;
}

void printMatrix(const vector<vector<int>>& X) {
	for (int i = 0; i < X.size(); i++) {
		for (int j = 0; j < X[i].size(); j++) {
			cout << X[i][j] << " ";
		}
		cout << endl;
	}
	cout << "---------------\n";
}

void printVector(const vector<int>& v) {
	//for (int i = 0; i < v.size(); i++) {
	//	cout << v[i] << " ";
	//}
	cout << "Central vertices counted: " << v.size() << endl;
}

class GeneticAlgorithm {
private:
	Graph* G;
	int populationSize;
	double mutationP, crossP;
	vector<int> population;
	Random rd;

	void startAlgorithm(double& time) {
		int stepN = 20;
		double start = clock();
		for (int i = 0; i < stepN; i++) {
			this->evolutionStep();
		}
		double finish = clock();
		string message = "Time of genetic algorithm: ";
		time = (finish - start) / CLOCKS_PER_SEC;
		#ifdef PRINT_TIME
			cout << "" << (finish - start) / CLOCKS_PER_SEC << endl;
		#endif
	}


public:
	GeneticAlgorithm(Graph* g, int popSize, double pC, double pM) {
		this->G = g;
		this->populationSize = popSize;
		this->crossP = pC;
		this->mutationP = pM;
		for (int i = 0; i < popSize; i++) {
			population.push_back(rand() % g->Size());
		}
	}

	void makeSelection() {
		vector<int> e;
		for (int i = 0; i < population.size(); i++) {
			e.push_back(G->getEccentricity(population[i]));
		}

		vector<double> probability;
		double maxV = *max_element(e.begin(), e.end()), leftSum = 0;
		for (int i = 0; i < e.size(); i++) {
			leftSum += maxV / double(e[i]);
		}
		double x = 1.0 / leftSum;
		for (int i = 0; i < e.size(); i++) {
			probability.push_back(maxV / double(e[i]) * x);
		}
		vector<int> nextPopulation;
		for (int i = 0; i < this->populationSize; i++) {
			nextPopulation.push_back(rd.choice(population, probability));
		}
		this->population = nextPopulation;
	}

	void crossing() {
		vector<int> crossed;
		for (int i = 0; i < populationSize; i++) {
			if (rd.getRandomLowerOne() < crossP) {
				int ind1 = rd.randint(0, populationSize - 1);
				int ind2 = rd.randint(0, populationSize - 1);
				int u = population[ind1], v = population[ind2];
				vector<int> path = G->getPath(u, v);
				crossed.push_back(path[path.size() / 2]);
			}
		}
		this->population.insert(population.end(), crossed.begin(), crossed.end());
	}

	void mutation() {
		for (unsigned int i = 0; i < this->population.size(); i++) {
			if (rd.getRandomLowerOne() < mutationP) {
				vector<int> n = G->getNeighbour(population[i]);
				population[i] = rd.choice(n);	
			}
		}
	}

	void printPopulation() {
		for (int i = 0; i < this->population.size(); i++) {
			cout << population[i] << " ";
		}
		cout << endl;
	}

	void evolutionStep() {
		this->crossing();
		this->makeSelection();
		this->mutation();
		#ifdef PRINT_GA
			this->printPopulation();
		#endif
	}

	int getBestResult(double& time) {
		this->startAlgorithm(time);
		vector<int> e;
		for (int i = 0; i < population.size(); i++) {
			int x = G->getEccentricity(population[i]);
			e.push_back(x);
		}
	 	vector<int> ind = getCentralVertex(e);
		#ifdef PRINT_GA
			for (int i = 0; i < ind.size(); i++) {
				cout << population[ind[i]] << " ";
			}
			cout << endl;
		#endif
		return this->G->getEccentricity(population[ind[0]]);
	}

	vector<int> getPopulation() {
		return this->population;
	}
	
};

class ExactAlgorithm {
private:
	int n, m;
	vector<pair<int, int>> adj;
	vector<vector<int>> adjTable;

	vector<vector<int>> APSP(const vector<vector<int>>& A) {
		bool stopecurtion = true;
		cout << "Step into" << endl;
		for (unsigned int i = 0; i < A.size(); i++) {
			for (unsigned int j = 0; j < A.size(); j++) {
				if (i != j && !A[i][j]) {
					stopecurtion = false;
					break;
				}
			}
			if (!stopecurtion) {
				break;
			}
		}
		if (stopecurtion) {
			return A;
		}
		vector<vector<int>> Z = A * A;
		vector<vector<int>> B(A.size(), vector<int>(A.size(), 0));
		for (int i = 0; i < B.size(); i++) {
			for (int j = 0; j < B.size(); j++) {
				if (i != j && (A[i][j] == 1 || Z[i][j] > 0)) {
					B[i][j] = 1;
				}
			}
		}
		vector<vector<int>> T = APSP(B);
		vector<vector<int>> X = T * A;
		vector<int> degree(A.size());
		for (int i = 0; i < A.size(); i++) {
			int sum = 0;
			for (int j = 0; j < A.size(); j++) {
				sum += A[i][j];
			}
			degree[i] = sum;
		}
		vector<vector<int>> D(A.size(), vector<int>(A.size()));
		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < A.size(); j++) {
				if (X[i][j] >= T[i][j] * degree[j]) {
					D[i][j] = 2 * T[i][j];
				}
				else {
					D[i][j] = 2 * T[i][j] - 1;
				}
			}
		}
		return D;
		cout << "Step out" << endl;
	}

	vector<vector<int>> goodAPSP(const vector<vector<int>>& A) {
		cout << "Step" << endl;
		int n = A.size();
		vector<vector<int>> Z = A * A;
		vector<int> inner(n);
		vector<vector<int>> B(n, inner);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j && (A[i][j] == 1 || Z[i][j] > 0)) {
					B[i][j] = 1;
				}
				else {
					B[i][j] = 0;
				}
			}
		}

		bool goodMatrix = true;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i != j && !B[i][j]) {
					goodMatrix = false;
				}
			}
		}
		vector<vector<int>> D(n, inner);
		if (goodMatrix) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					D[i][j] = 2 * B[i][j] - A[i][j];
				}
			}
			return D;
		}
		vector<vector<int>> T = goodAPSP(B);
		vector<vector<int>> X = T * A;
		vector<int> degree(n);
		for (int i = 0; i < n; i++) {
			int sum = 0;
			for (int j = 0; j < n; j++) {
				sum += A[i][j];
			}
			degree[i] = sum;
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (X[i][j] > T[i][j] * degree[j]) {
					D[i][j] = 2 * T[i][j];
				}
				else {
					D[i][j] = 2 * T[i][j] - 1;
				}
			}
		}
		return D;
	}

	vector<int> makeDominatingSet(const vector<vector<int>>& A) {
		if (A.size() <= 0) {
			return vector<int>();
		}
		vector<bool> U(this->n);
		vector<int> C;
		vector<bool> usedS(A.size());
		while (true) {
			int maxV = 0, index = -1;
			for (int i = 0; i < A.size(); i++) {
				int current = 0;
				if (!usedS[i]) {
					for (int j = 0; j < A[i].size(); j++) {
						int x = A[i][j];
						if (!U[x]) {
							current++;
						}
					}
				}
				if (current > maxV) {
					maxV = current;
					index = i;
				}
			}
			if (index == -1) {
				break;
			}
			for (int i = 0; i < A[index].size(); i++) {
				U[A[index][i]] = true;
			}
			C.push_back(index);
		}

		set<int> Union;
		for (int i = 0; i < C.size(); i++) {
			Union.insert(A[C[i]].begin(), A[C[i]].end());
		}
		vector<int> result;
		for (auto i = Union.begin(); i != Union.end(); i++) {
			result.push_back(*i);
		}
		return result;
	}
	
public:
	ExactAlgorithm(int n, int m, vector<pair<int, int>> a) {
		this->n = n;
		this->m = m;
		this->adj = a;
		this->adjTable.resize(n);
		for (unsigned i = 0; i < a.size(); i++) {
			int u = a[i].first, v = a[i].second;
			this->adjTable[u].push_back(v);
			this->adjTable[v].push_back(u);
		}
	}

	double TrivialAlgorithm() {
		cout << "Trivial algorithm" << endl;
		Graph* g = new Graph(n, m, adj);
		vector<int> e;
		double start = clock();
		for (int i = 0; i < n; i++) {
			int x = g->getEccentricity(i);
			e.push_back(x);
		}
		double finish = clock();
		double diff = (finish - start) / CLOCKS_PER_SEC;
		cout << "Time of trivial algorithm: " << diff << endl;
		printVector(getCentralVertex(e));
		delete g;
		return diff;
	}
	
	void SeidelsAlgorithm() {
		cout << "Seidels Algorithm" << endl;
		vector<int> inner(this->n);
		vector<vector<int>> A(this->n, inner);
		for (int i = 0; i < adj.size(); i++) {
			int u = adj[i].first, v = adj[i].second;
			A[u][v] = 1;
			A[v][u] = 1;
		}

		vector<vector<int>> result = APSP(A);
		vector<int> e(this->n);
		for (int i = 0; i < result.size(); i++) {
			e[i] = *max_element(result[i].begin(), result[i].end());
		}
		printVector(getCentralVertex(e));
	}

	void AingworthAlgorithm() {
		cout << "Aingworth algorithm" << endl;
		double start = clock();
		Graph g(n, m, adj);
		int s = sqrt(n * log(n));
		vector<int> L;
		vector<int> H;
		for (int i = 0; i < n; i++) {
			int x = g.getNeighbour(i).size();
			if (x < s) {
				L.push_back(i);
			}
			else {
				H.push_back(i);
			}
		}
		cout << "S = " << s << endl;
		cout << "L = " << L.size() << endl;
		cout << "H = " << H.size() << endl;
		g.setSubGraph(L);
		for (int i = 0; i < L.size(); i++) {
			g.bfsFromVertex(L[i]);
		}
		g.resetSubGraph();
		vector<vector<int>> A;
		for (int i = 0; i < H.size(); i++) {
			vector<int> Ne = g.getNeighbour(H[i]);
			Ne.push_back(H[i]);
			A.push_back(Ne);
		}
		cout << "Making domain set" << endl;
		vector<int> D = makeDominatingSet(A);
		cout << "Domain set maked" << endl;
		cout << D.size() << endl; 
		for (int i = 0; i < D.size(); i++) {
			g.bfsFromVertex(D[i]);
		}
		vector<bool> Diff(this->n, true);
		for (int i = 0; i < D.size(); i++) {
			Diff[D[i]] = false;
		}

		cout << "Making Table" << endl;
		vector<vector<int>> d = g.getDistanceMatrix();
		for (int u = 0; u < n; u++) {
			if (!Diff[u]) {
				continue;
			}
			for (int v = u + 1; v < n; v++) {
				if (!Diff[v]) {
					continue;
				}
				int dUV = d[u][v];
				int minV = INT_MAX;
				for (int w = 0; w < D.size(); w++) {
					int d1 = d[D[w]][u];
					int d2 = d[D[w]][v];
					if (d1 + d2 < minV) {
						minV = d1 + d2;
					}
				}
				d[u][v] = min(minV, dUV);
				d[v][u] = min(minV, dUV);
			}
		}
		double finish = clock();
		cout << "Time of AingworthAlgorithm: " << (finish - start) / CLOCKS_PER_SEC << endl;
		vector<int> e;
		for (int i = 0; i < n; i++) {
			int x = *max_element(d[i].begin(), d[i].end());
			e.push_back(x);
		}
		printVector(getCentralVertex(e));
	}
};











