#include "minifastgreedybli.h"

#include <algorithm>

using namespace std;

SignedMiniFastGreedyBLI::SignedMiniFastGreedyBLI(LargeGraph *graph,float alpha, IloNum W, float PAlpha,  int   K) 
    : SignedMiniLouvainCG(graph, W, PAlpha, K)
{
    this->alpha = alpha;
}


void SignedMiniFastGreedyBLI::execute(IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5){


    //choosing the alpha
    this->coarsedNodes.clear();
    this->currentLNodes.clear();
    ////starting a new search
    this->currentValue = 0.0;
    this->maxValue = 0.0;

    //... all nodes are not in the new cluster (new column)
    for (unsigned i=0;i<graph->numberOfNodes;i++){
        this->maxNodes[i]=0;
        this->currentNodes[i]=0;

    }

    //... start the coarsened nodes that are not coarsened yeat
    this->coarsedNodes.clear();
    for (unsigned i=0;i<graph->numberOfNodes;i++){
        unsigned deg = this->graph->degreeOfNode[i];
        this->coarsedNodes.push_back(MLVCGNode(i, 0, deg));
    }

    this->sumLambdaS = 0.0;
    float lastValue=0.0;


    ////mini louvain method
    vector<unsigned> randomOrder(this->coarsedNodes.size());

    unsigned notImproved = 0;

    while( notImproved < graph->numberOfNodes){
		//start random solution
		shuffleCurrentSolution(Wv, Ce, lambda1, lambda4, lambda5, 0.5);

        ////level phase
        int notImprovedLvL = 0;
        while (notImprovedLvL < 2){//graph->numberOfNodes){
            for (int i=0 ; i<this->coarsedNodes.size() ; i++){
                randomOrder[i]=i;
            }
            for (int i=0 ; i<this->coarsedNodes.size()-1 ; i++) {
                int randpos = rand()%(this->coarsedNodes.size()-i)+i;
                swap(randomOrder[i], randomOrder[randpos]);
            }

            for(unsigned inode = 0; inode < this->coarsedNodes.size(); inode++){
                unsigned node = randomOrder[inode];
                float value;
                if (this->currentNodes[node] == 0){
                    value = this->gainIfIEnter(node, Wv, Ce, lambda1, lambda4, lambda5);


                    if(value < -0.0000001){
						
                        this->currentNodes[node] = 1;
                        //an improvement happens
                        notImprovedLvL = 0;

                        this->currentValue += value;
                        //each subnode belongs to the solution
                        unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
                        while (itnews != this->coarsedNodes[node].nodes.end()){
                            this->currentLNodes.push_back(*itnews);
                            itnews++;
                        }

                    }
                }else{
                    value = this->gainIfIExit(node, Wv, Ce, lambda1, lambda4, lambda5);
                    if(value < -0.0000001){

                        this->currentNodes[node] = 0;
                        //an improvement happens
                        notImprovedLvL = 0;

                        this->currentValue += value;
                        //each subnode belongs to the solution
                        unordered_set<unsigned>::iterator itnews = this->coarsedNodes[node].nodes.begin();
                        while (itnews != this->coarsedNodes[node].nodes.end()){
                            this->currentLNodes.remove(*itnews);
                            itnews++;
                        }
                    }
                }
                
                node++;
            }
            notImprovedLvL++;

        }
        if (lastValue <= this->currentValue){
            notImproved++;
        }
        lastValue = this->currentValue;

        if (this->currentValue < this->maxValue){
            this->maxValue = this->currentValue;
            this->copySolutionArray(this->maxNodes);
        }
        //shuffleCurrentSolution(lambda, this->alpha);

    }

}


