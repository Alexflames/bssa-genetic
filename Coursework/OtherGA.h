#pragma once

#include <vector>
#include <set>
#include <algorithm>
#include "GraphWork.h"
#include "Random.h"
#include "SetOperation.h"
#include "GeneticAlgorithm.h"


using std::set;
using std::max_element;

class SimpleGeneticAlgorithm {
private:
	Graph* G;
	int populationSize;
	int chromoLength;
	double mutationP, crossP;
	Random rd;
	vector<set<int>> population;
	
	void startAlgorithm(double& time) {
		int stepN = 20;
		double start = clock();
		for (int i = 0; i < stepN; i++) {
			this->evolutionStep();
		}
		double finish = clock();
		time = (finish - start) / CLOCKS_PER_SEC;
		string message = "Time of OTHER genetic algorithm: ";
#ifdef PRINT_TIME
		cout << "" << time << endl;
#endif
	}

	
	int function(vector<int> ch) {
		int sum = 0;
		for (auto it = ch.begin(); it != ch.end(); it++) {
		
			sum += this->G->getEccentricity(*it);
		}
		return sum;
	}

	int function(set<int> ch) {
		int sum = 0;
		for (auto it = ch.begin(); it != ch.end(); it++) {
			sum += this->G->getEccentricity(*it);
		}
		return sum;
	}


public:
	SimpleGeneticAlgorithm(Graph* g, int popSize, int chromoLength, double pC, double pM) {
		this->G = g;
		this->populationSize = popSize;
		this->crossP = pC;
		this->mutationP = pM;
		this->chromoLength = chromoLength;
		for (int i = 0; i < popSize; i++) {
			set<int> chromo;
			while (chromo.size() < chromoLength) {
				int v = rand() % g->Size();
				chromo.insert(v);
			}
			population.push_back(chromo);
		}
	}

	void evolutionStep() {
		this->makeSelection();
		this->crossing();
		this->N4N();
#ifdef PRINT_GA
		for (int i = 0; i < populationSize; i++) {
			cout << "{ ";
			for (auto it = population[i].begin(); it != population[i].end(); it++) {
				cout << *it << " ";
			}
			cout << "}" << endl;
		}
		cout << "----|----|----|----|----" << endl;
#endif
	}


	void N4N() {
		const int COUNT = 4;
		int percent10 = 0.1 * populationSize;
		vector<bool> used(populationSize, false);
		vector<int> indexChromo;
		for (int i = 0; i < percent10; i++) {
			while (true) {
				int c = rand() % populationSize;
				if (!used[c]) {
					used[c] = true;
					indexChromo.push_back(c);
					break;
				}
			}
		}
		for (int k = 0; k < indexChromo.size(); k++) {
			vector<int> best = makeVector(population[indexChromo[k]]);
			vector<int> X = best;
			int currentBest = function(X);
			for (int i = 0; i < X.size(); i++) {
				vector<int> neight = this->G->getNeighbour(X[i]);
				for (int j = 0, t = 0; j < neight.size() && t < COUNT; j++) {
					auto it = find(X.begin(), X.end(), neight[j]);
					if (it != X.end()) {
						continue;
					}
					t++;
					int tmp = X[i];
					X[i] = neight[j];
					int currentFunction = function(X);
					if (currentFunction < currentBest) {
						best = X;
						currentBest = currentFunction;
					}
					X[i] = tmp;
				}
			}
			population[indexChromo[k]] = makeSet(best);
		}
	}


	void mutation() {
		for (int i = 0; i < populationSize; i++) {
			set<int> chromo = population[i];
			while (true) {
				int v = rand() % this->G->Size();
				if (chromo.find(v) == chromo.end()) {
					chromo.erase(chromo.begin());
					chromo.insert(v);
					break;
				}
			}
		}
	}

	void makeSelection() {
		vector<int> fitnessFunction;

		for (int i = 0; i < this->populationSize; i++) {
			int sum = 0;
			set<int> chromo = this->population[i];
			for (auto it = chromo.begin(); it != chromo.end(); it++) {
				sum += G->getEccentricity(*it);
			}
			fitnessFunction.push_back(sum);
		}
		
		vector<double> probability;
		double maxV = *max_element(fitnessFunction.begin(), fitnessFunction.end()), leftSum = 0;
		for (int i = 0; i < fitnessFunction.size(); i++) {
			leftSum += maxV / double(fitnessFunction[i]);
		}
		double x = 1.0 / leftSum;
		for (int i = 0; i < fitnessFunction.size(); i++) {
			probability.push_back(maxV / double(fitnessFunction[i]) * x);
		}

		vector<set<int>> nextPopulation;
		for (int i = 0; i < this->populationSize; i++) {
			nextPopulation.push_back(rd.choice(population, probability));
		}
		this->population = nextPopulation;
	}





	void crossing() {
		for (int i = 0; i < this->populationSize / 2; i++) {
			if (rd.getRandomLowerOne() >= crossP) {
				continue;
			}
			int i1 = rand() % populationSize;
			int i2 = rand() % populationSize;
			set<int> vp1 = population[i1] - population[i2];
			set<int> vp2 = population[i2] - population[i1];
			if (vp1.size() == 0 || vp2.size() == 0) {
				continue;
			}
			if (vp1.size() != vp2.size()) {
				throw exception("«р€ ты так думал");
			}
			
			int r = rand() % vp1.size();
			vector<int> v1 = makeVector(population[i1]);
			vector<int> v2 = makeVector(population[i2]);
			vector<int> xch1 = makeVector(vp2);
			vector<int> xch2 = makeVector(vp1);

			for (int j = 0; j < r; j++) {
				v1[rand() % v1.size()] = xch1[j];
				v2[rand() % v2.size()] = xch2[j];
			}
			population[i1] = makeSet(v1);
			population[i2] = makeSet(v2);
		}
	}

	int getBestResult(double& time) {
		this->startAlgorithm(time);
		int index = 0, minSum = INT_MAX;
		for (int i = 0; i < this->populationSize; i++) {
			int sum = 0;
			for (auto it = population[i].begin(); it != population[i].end(); it++) {
				sum += this->G->getEccentricity(*it);
			}
			if (sum < minSum) {
				minSum = sum, index = i;
			}
		}
		vector<int> bestChromo = makeVector(population[index]);
		vector<int> e;
		for (int i = 0; i < bestChromo.size(); i++) {
			e.push_back(G->getEccentricity(bestChromo[i]));
		}
		vector<int> bestIndex = getCentralVertex(e);
#ifdef PRINT_GA
		for (int i = 0; i < bestIndex.size(); i++) {
			cout << bestChromo[bestIndex[i]] << " ";
		}
		cout << endl;
#endif
		return this->G->getEccentricity(bestChromo[bestIndex[0]]);
	}
};

