#ifndef MINIFASTGREEDYBLI_H
#define MINIFASTGREEDYBLI_H

#include "minilouvaincg.h"



class SignedMiniFastGreedyBLI: public SignedMiniLouvainCG
{
public:
    SignedMiniFastGreedyBLI(LargeGraph *graph, float alpha, IloNum W, float PAlpha,  int   K);
    void execute(IloNumArray &Wv,
										IloNumArray &Ce,
										IloNumArray &lambda1,
										//IloNumVar   &lambda3 <=0,
										IloNumArray &lambda4,
										IloNumArray &lambda5);

    float alpha = 0.0;
};


#include "minifastgreedybli.cpp"
#endif // MINIFASTGREEDYBLI_H
