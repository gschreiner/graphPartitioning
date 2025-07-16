#ifndef MINILOUVAINCG_H
#define MINILOUVAINCG_H

#include <unordered_set>
#include <list>
#include <vector>
#include "../graph/largegraph.h"

#include <ilcplex/ilocplex.h>
#include <ilconcert/iloenv.h>

struct MLVCGNode{
    unsigned internalEdges = 0;
    unsigned degree = 0;
    unordered_set<unsigned> nodes;
   
    MLVCGNode(unsigned node, unsigned edges, unsigned deg):
        internalEdges(edges), degree(deg){
        nodes.insert(node);
    }

};

class SignedMiniLouvainCG
{
public:
    LargeGraph * graph;
    vector<MLVCGNode> coarsedNodes;
    float lambda = 0.5;
    IloNum **Wij;
	
    float currentValue = 0.0;
    vector<unsigned> currentNodes; //binary array
    list<unsigned> currentLNodes; //list of nodes that belongs to the cluster

    float sumLambdaS = 0.0;

    vector<unsigned> a;//solution which is constructed by 'constructSolutionArray' method

    float maxValue = 0.0;
    vector<unsigned> maxNodes;

	//problem parameters:
	IloNum W = 0;
	float PAlpha;
	int   K;



    SignedMiniLouvainCG(LargeGraph * graph, IloNum W, float PAlpha,  int   K);

    void execute(						IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5,
										float strength=0.0);
    float gainIfIEnter(unsigned node, 
										IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5);
    float gainIfIExit(unsigned node, 
										IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5);
   
    float calculateAll(					IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5);//verify correctness
    void constructSolutionArray();
    void shuffleCurrentSolution(		IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5,
										float strength);
    void copySolutionArray(vector<unsigned> &target);
};

#include "minilouvaincg.cpp"

#endif // MINILOUVAINCG_H
