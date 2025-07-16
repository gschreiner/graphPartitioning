#include "minilouvaincg.h"

#include <algorithm>

using namespace std;

SignedMiniLouvainCG::SignedMiniLouvainCG(LargeGraph *graph, IloNum W, float PAlpha,  int   K)
{
    this->graph = graph;
    this->W = W;
    this->PAlpha = PAlpha;
    this->K = K;

    
    for (unsigned i=0;i<graph->numberOfNodes;i++){
        this->maxNodes.push_back(0);
        this->currentNodes.push_back(0);
        a.push_back(0);
    }
}




float SignedMiniLouvainCG::gainIfIEnter(unsigned node, 
										IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5){
											
											
    //first part
    int  m=graph->numberOfEdges;		
    this->currentNodes[node] = 1;
    
    //verifying the balance limit 
    float limit = -(1.0 +PAlpha)* W/K;
    for (int i=0; i < graph->numberOfNodes; i++) {
		limit += this->currentNodes[i]*Wv[i];
	}
	if (limit > 0){//do not add
		this->currentNodes[node] = 0;
		return 1000;
	}
			
		
		
	float value = 0.0;
	for (int e=0;e<m;e++){
		unsigned a =  graph->edges[e].v1;
		unsigned b =  graph->edges[e].v2;
		if (this->currentNodes[a] + this->currentNodes[b] == 1){
			value += Ce[e];
		}
	}
	for (unsigned i=0;i<graph->numberOfNodes;i++){
		value -= this->currentNodes[i]*lambda1[i];
	}
	value -= lambda4[0];
	value += lambda5[0];
	
	
    this->currentNodes[node] = 0;

    return /*max*/ value-this->currentValue ;
}

float SignedMiniLouvainCG::gainIfIExit(	unsigned node, 
										IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5
									  ){
										  
	int  m=graph->numberOfEdges;										  
    //first part
    this->currentNodes[node] = 0;
    
    float value = 0.0;
    for (int e=0;e<m;e++){
		unsigned a =  graph->edges[e].v1;
		unsigned b =  graph->edges[e].v2;
		if (this->currentNodes[a] + this->currentNodes[b] == 1){
			value += Ce[e];
		}
	}
	for (unsigned i=0;i<graph->numberOfNodes;i++){
		value -= this->currentNodes[i]*lambda1[i];
	}
	value -= lambda4[0];
	value += lambda5[0];
    this->currentNodes[node] = 1;
    

    return /*max*/ value-this->currentValue ;
}


float SignedMiniLouvainCG::calculateAll(
										IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5){

 int  m=graph->numberOfEdges;		
    
    float value = 0.0;
    for (int e=0;e<m;e++){
		unsigned a =  graph->edges[e].v1;
		unsigned b =  graph->edges[e].v2;
		if (this->currentNodes[a] + this->currentNodes[b] == 1){
			value += Ce[e];
		}
	}
	for (unsigned i=0;i<graph->numberOfNodes;i++){
		value -= this->currentNodes[i]*lambda1[i];
	}
	value -= lambda4[0];
	value += lambda5[0];
    

    return /*max*/ value ;

}

void SignedMiniLouvainCG::constructSolutionArray(){
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            this->a[*ita] = this->currentNodes[i];
            ita++;
        }
    }
}

void SignedMiniLouvainCG::shuffleCurrentSolution(
										IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5,
										float strength){
    unsigned nchanges = this->graph->numberOfNodes * strength;
    float value;
    for (unsigned ic=0;ic<nchanges;ic++){
        unsigned node = rand()%this->coarsedNodes.size();
        if (this->currentNodes[node] == 0){
            value = this->gainIfIEnter(node, Wv, Ce, lambda1, lambda4, lambda5);
            this->currentNodes[node] = 1;
            //an improvement happens
            this->currentValue += value;
            //each subnode belongs to the solution
            unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
            while (itnews != this->coarsedNodes[node].nodes.end()){
                this->currentLNodes.push_back(*itnews);
                itnews++;
            }

        }else{
            value = this->gainIfIExit(node, Wv, Ce, lambda1, lambda4, lambda5);
            this->currentNodes[node] = 0;
            //an improvement happens
            this->currentValue += value;
            //each subnode belongs to the solution
            unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
            while (itnews != this->coarsedNodes[node].nodes.end()){
                this->currentLNodes.remove(*itnews);
                itnews++;
            }

        }

    }

}

void SignedMiniLouvainCG::copySolutionArray(vector<unsigned> &target){
    for (unsigned i=0;i<this->coarsedNodes.size();i++){
        unordered_set<unsigned>::iterator ita = this->coarsedNodes[i].nodes.begin();
        while (ita != this->coarsedNodes[i].nodes.end()){
            target[*ita] = this->currentNodes[i];
            ita++;
        }
    }
}
