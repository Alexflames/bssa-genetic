#include <fstream>
#include "GraphWork.h"
#include "Random.h"
#include "GeneticAlgorithm.h"
#include "OtherGA.h"
#include "FileWork.h"

using std::make_pair;
using std::ifstream;


//#define GEOM_GRAPH
//#define BA_GRAPH
#define ERDOSH_GRAPH

#define SIMPLE_GA
//#define OUR_GA

const string BA_FILE[] = { "BA_Graph\\BarabasiAlbertGraph1_M2.txt", 
				    	   "BA_Graph\\BarabasiAlbertGraph2_M2.txt",
						   "BA_Graph\\BarabasiAlbertGraph3_M2.txt",
						   "BA_Graph\\BarabasiAlbertGraph4_M2.txt",
						   "BA_Graph\\BarabasiAlbertGraph5_M2.txt",
						   "BA_Graph\\BarabasiAlbertGraph6_M2.txt",
						   "BA_Graph\\BarabasiAlbertGraph7_M2.txt" };

const string GEOM_FILE[] = { "GEOM_Graph\\GeometricGraph1_R01.txt",
						     "GEOM_Graph\\GeometricGraph2_R01.txt",
						     "GEOM_Graph\\GeometricGraph3_R01.txt",
						     "GEOM_Graph\\GeometricGraph4_R01.txt",
						     "GEOM_Graph\\GeometricGraph5_R01.txt",
						     "GEOM_Graph\\GeometricGraph6_R01.txt",
							 "GEOM_Graph\\GeometricGraph7_R01.txt" };

const string ER_FILE[] = { "ERDOSRENYI_Graph\\ErdosRenyi1_P001.txt",
						   "ERDOSRENYI_Graph\\ErdosRenyi2_P001.txt",
						   "ERDOSRENYI_Graph\\ErdosRenyi3_P001.txt", 
						   "ERDOSRENYI_Graph\\ErdosRenyi4_P001.txt",
						   "ERDOSRENYI_Graph\\ErdosRenyi5_P001.txt",
						   "ERDOSRENYI_Graph\\ErdosRenyi6_P001.txt",
						   "ERDOSRENYI_Graph\\ErdosRenyi7_P001.txt" };


const int TEST = 5;
const int ITER = 10;

const int RADIUS_BA[]   = { 4, 4, 5, 4, 5, 5, 5 };
const int RADIUS_GEOM[] = { 9, 9, 8, 8, 8, 8, 8 };
const int RADIUS_ER[]   = { 5, 4, 4, 4, 3, 3, 3 };

int main() {
	srand(time(NULL) % INT_MAX);
	int n, m;
#ifdef BA_GRAPH
	vector<pair<int, int>> e = readFromFileNM(BA_FILE[TEST], n, m);
#endif
#ifdef GEOM_GRAPH
	vector<pair<int, int>> e = readFromFileNM(GEOM_FILE[TEST], n, m);
#endif

#ifdef ERDOSH_GRAPH
	vector<pair<int, int>> e = readFromFileNM(ER_FILE[TEST], n, m);
#endif

	cout << "N = " << n << " M = " << m << endl;
	ExactAlgorithm alg(n, m, e);
 	double t1 = alg.TrivialAlgorithm();
	cout << "Aingworth " << badTimeFunction(t1, n, m) << endl;
	
	//alg.SeidelsAlgorithm();
	//alg.AingworthAlgorithm();
	#ifdef BA_GRAPH
		const int REAL_R = RADIUS_BA[TEST];
	#endif

	#ifdef GEOM_GRAPH
		const int REAL_R = RADIUS_GEOM[TEST];
	#endif
	#ifdef ERDOSH_GRAPH
			const int REAL_R = RADIUS_ER[TEST];
	#endif

#ifdef OUR_GA
		{
			int error = 0;
			double avg_time = 0.0;
			for (int i = 0; i < ITER; i++) {
				Graph* g1 = new Graph(n, m, e);
				cout << "Start genetic algorithm" << endl;
				GeneticAlgorithm genAlg(g1, 50, 0.7, 0.1);
				double time = 0.0;
				int R = genAlg.getBestResult(time);
				avg_time += time;
				if (R != REAL_R)
					error++;
			}
			cout << "AVG Time = " << avg_time / double(ITER) << endl;
			cout << "Error = " << double(error) / double(ITER) << endl;
		}
#endif
	
#ifdef SIMPLE_GA
		{
			int error = 0;
			double avg_time = 0.0;
			for (int i = 0; i < ITER; i++) {
				Graph* g1 = new Graph(n, m, e);
				cout << "Start OTHER genetic algorithm" << endl;
				SimpleGeneticAlgorithm sgen(g1, 50, 10, 0.7, 0.1);
				double time = 0.0;
 				int R = sgen.getBestResult(time);
				avg_time += time;
				if (R != REAL_R)
					error++;
			}
			cout << "AVG Time = " << avg_time / double(ITER) << endl;
			cout << "Error = " << double(error) / double(ITER) << endl;
		}
#endif
}


