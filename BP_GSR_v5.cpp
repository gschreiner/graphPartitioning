#define GROUND_TRUTH

#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <cstdlib>
#include <time.h>
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloenv.h>
#include <chrono>
#include <list>





#include "./graph/largegraph.h"
#include "./graph/solution.h"
#include "./utils/modularitylg.h"


using namespace std;
using namespace std::chrono;
using namespace boost::heap;

ILOSTLBEGIN

#define RC_EPS 1.0e-9


//number of columns
long long int ncInitial = 0;
long long int ncAPHeuristic = 0;
long long int ncAPExact = 0;

long long int totalTimeInitial      = 0;
long long int totalTimeLP           = 0;
long long int totalTimeAPExact      = 0;
long long int totalTimeAPHeuristic  = 0;
long double solveRMP = 0.0;
long double solveAux = 0.0;
string experimento = "0";

IloEnv env;

bool PAR_VERBOSE_ITER = false;


IloNumArray Wv(env);
IloNumArray Ce(env);
string parameter = "";
IloNum BigM = 0;

chrono::system_clock::time_point antes;

#define EPSILON 1E-6
bool DBL_EQL(double a, double b) {
    return fabs(a - b) < EPSILON;
}




static IloNum Wt(LargeGraph &graph,IloNumArray &A){
	IloNum wt= 0;
	for (IloInt i=0; i<A.getSize(); i++){
		wt += A[i]*Wv[i];
	}
	return wt;
}



//"prepareCe" prepares each edge weight
void prepareCe(LargeGraph &graph, IloNumArray &Ce){	
    ////memory allocation
    IloNum maxi = 0;   
    for (IloInt i=0; i < graph.numberOfEdges; i++) {
		unsigned u = graph.edges[i].v1;
		unsigned v = graph.edges[i].v2;     		
        Ce.add(graph.getAdj(u, v));
        maxi =  max(maxi, Ce[i]);

    }
    BigM = max(BigM, maxi*graph.numberOfEdges);

}

//"prepareWv" prepares each edge weight
void prepareWv(LargeGraph &graph, IloNumArray &Wv){
    //memory allocation
    //Wv = new IloNumArray (graph.numberOfNodes);
    IloNum maxi = 0;
    for (IloInt i=0; i < graph.numberOfNodes; i++) {
        Wv.add(atof(graph.labelOfNode[i].c_str()));
        maxi =  max(maxi, Wv[i]);
    }
    BigM = max(BigM, maxi*graph.numberOfNodes);  

}





using namespace std;

struct TSolution{
    pair<IloInt, IloInt> fixed;
    long long int parentId;
    long long int id;
    IloNum value;
    bool operator<(TSolution const & rhs) const{
        return value > rhs.value; ///inverti do density pq agora temos uma minimização
    }
};


/*
(    graph,
                               PAlpha,
                               heap,
                               data,
                               bpnodes,
                               lowerbound,
                               bestZs,
                               A,
                               Z,
                               masterOpt,
                               masterCplex,
                               mauxOpt,
                               AA,
                               consLambda1,
							   consLambda2,
							   consLambda3,
							   consLambda4,
							   consLambda5,
                               masterObj
*/


void print(IloNumArray &xZ, IloArray<IloNumArray> & A ){
		string comms = "";
    for (IloInt i=0;i<xZ.getSize(); i++){
        string comm = "";
        if ((xZ)[i] > 0.0){
			
            for (IloInt v=0;v<A[i].getSize(); v++){
                if (A[i][v] > 0.9){
                    if (comm == ""){
                        comm = to_string(v);
                    }else{
                        comm += ","+to_string(v);
                    }
                }
            }
            if (comm != ""){
                if (comms  == ""){
                    comms+= "["+comm+"]";
                }else{
                    comms+= ",["+comm+"]";
                }
                comms += "("+to_string(xZ[i])+")";
                
            }
        }
    }
    cout<<"Fim: "<<comms<<endl;

}




IloNumArray * solver(
            LargeGraph                  &graph,
            float                       &PAlpha,
            int 						&K,
			IloNum						&W,
            fibonacci_heap<TSolution>   &heap,
            TSolution                   &bpnode,
            vector<TSolution *>         &bpnodes,
            IloNum                      &lowerbound,
            IloNumArray                 *&bestZs,
            IloArray<IloNumArray>       &A,
            IloNumVarArray              &Z,
            IloNumVarArray              &Y,
            IloModel                    &masterOpt,
            IloCplex                    &masterCplex,
            IloNumVarArray              &AA,
            IloRangeArray               &consLambda1,
            IloRangeArray               &consLambda2,
            IloRangeArray               &consLambda3,
            IloRangeArray               &consLambda4,
            IloObjective                &masterObj,
            IloNumArray 				&disableds
        ){

	masterOpt.remove(consLambda2);
	IloAdd(masterOpt ,consLambda2);
	masterOpt.remove(consLambda3);
	IloAdd(masterOpt ,consLambda3);

    //node and solution
    IloConstraintArray nodeConstraints(env);
cout<<"\n\n*******************\n";
    long long int id = bpnode.id;
    cout<<id<<" <- ";
    while (bpnodes[id]->parentId != -1){
//cout<<"# 182"<<endl;                           
cout<<"::"<<bpnodes[id]->fixed.first<<endl;
        IloConstraint cons = Z[bpnodes[id]->fixed.first] == bpnodes[id]->fixed.second;
        masterOpt.add(cons);
        nodeConstraints.add(cons);
        id = bpnodes[id]->parentId;
        cout<<id<<" <- ";
    }
cout<<"\n*******************\n\n";
    IloInt it, i, t, tt, j, u, v, e;
    IloInt  m=graph.numberOfEdges, n=graph.numberOfNodes;

    //parameters and results (price and new column)
    IloNumArray lambda1(env, graph.numberOfNodes); //dual variables from model
    //IloNumArray lambda2(env, graph.numberOfEdges*A.getSize()); //dual variables from model
    IloNumVarArray   lambda2(env, m, 0, +IloInfinity, ILOFLOAT);
    //IloNumVarArray   R(env, m, 0, 1, ILOINT);
    IloNumVarArray   Yl(env, m, 0, 1, ILOINT);
    //IloNumVarArray   S(env, m, 0, 1, ILOINT);
    //IloNumArray lambda3(env,  1); //dual variables from model
    
    IloNumVar   lambda3(env, -IloInfinity, 0.0, ILOFLOAT);
    IloNumArray lambda4(env,  1); //dual variables from model
    IloNumArray newCommunity(env, graph.numberOfNodes);

    IloNum minValue = IloInfinity;
    IloNumArray * Zs;
    //                 _===    _     
    //////////////////||===//|| ||//||>>///////////////////////////////
	//////////////////||/////||_||//||\\///////////////////////////////
    for(it=1;;it++){
		/*for (i=0;i<A.getSize();i++){
            masterOpt.add(Z[i] >= 0);
        }*/
        
        IloCplex masterCplex(masterOpt);
        masterCplex.setParam(IloCplex::Param::Threads,1);
        masterCplex.setParam(IloCplex::TiLim,36000);
        
		cout<<"Starting iteration: "<<it<<""<<endl;	
		masterCplex.exportModel("antes.lp");
		cout<<"\t sOlution printed"<<endl;	
        //solve the reduced master problem
        masterCplex.solve();
		cout<<"Master value: "<< masterCplex.getObjValue()<<endl;
IloNumArray xZ(env, Z.getSize());
masterCplex.getValues(xZ, Z);
cout<<"Sol.: "<< xZ <<endl;

for (int i=0; i< xZ.getSize(); i++){
	if (xZ[i] < 0){
		cout<<"Deu caca"<<endl;
		
		
		exit(0);
	}
}
		//exit(0);
		
        if (minValue > masterCplex.getObjValue() ){
            minValue = masterCplex.getObjValue();
        }
        
        //print the data of each iteration
        if (PAR_VERBOSE_ITER){
            //instance,method,dvalue,iteration
            string file = "evolutionByIter.csv";
            string texto = parameter+",BP,"+to_string(minValue)+","+to_string(it)+","+experimento+"\n";
            ofstream ofs (file, std::fstream::out | std::fstream::app);
            ofs << texto;
            ofs.close();
        }
		IloModel mauxOpt(env);
        //storing the dual variables
        for (i=0;i<n;i++){
            lambda1[i] = masterCplex.getDual(consLambda1[i]);
        }
        lambda4[0] = masterCplex.getDual(consLambda4[0]);

/*
cout<<"Lambda_1: "<<lambda1<<endl;	
cout<<"Lambda_2: "<<lambda2<<endl;	
cout<<"Lambda_3: "<<lambda3<<endl;	
cout<<"Lambda_4: "<<lambda4<<endl;	*/
		
		
		
		
		//Composing the model for the auxiliary problem
		//Objective function
		IloExpr mauxOFExp(env);
		//the main part of the objective function
		for (IloInt e=0; e<m; e++){
			mauxOFExp +=   Ce[e] * Yl[e];
		}
		for (IloNum v=0;v<n;v++){
			mauxOFExp -= AA[v]*lambda1[v];
		}
		for (e=0;e<m;e++){
			mauxOFExp -=  -Y[e]*lambda2[e];
		}
		for (IloInt i=0; i < n; i++) {
			mauxOFExp -= Wv[i]*AA[i]*lambda3;
		}
		mauxOFExp -= lambda4[0];
		IloObjective obj = IloAdd(mauxOpt, IloMinimize(env, mauxOFExp  ));
		
		//executar o corte diversas vezes até exaurir o valor de K.
		int V = 0; //vertices ja cobertos
		int k=0;//quantidade de clusters
		/*bool aux_opt = false;
		while (k<K){	*/
			//constraints
			for (e=0;e<m;e++){
				IloInt a =  graph.edges[e].v1;
				IloInt b =  graph.edges[e].v2;
				IloInt u = min(a,b);
				IloInt v = max(a,b);
				mauxOpt.add(Yl[e] >= AA[u] - AA[v]);
				mauxOpt.add(Yl[e] >= AA[v] - AA[u]);
				mauxOpt.add(Yl[e] <= AA[v] + AA[u]);
				mauxOpt.add(Yl[e] <= 2 -AA[v] - AA[u]);
			}
			mauxOFExp.end();
			
			//solving the column generator
			IloCplex mauxCplex(mauxOpt);
			mauxCplex.setParam(IloCplex::Param::Threads,1);
			mauxCplex.setParam(IloCplex::TiLim,36000);
			mauxCplex.exportModel("aux.lp");
			mauxCplex.solve();
			cout<<"Auxiliary value: "<< mauxCplex.getObjValue()<<endl;
			
			
			
			
			   
			//getting the new column
			//copy the new column
			//mauxCplex.getValues(newCommunity, AA);
			for (IloInt i=0; i<AA.getSize(); i++){
				try{
					newCommunity[i] = mauxCplex.getValue(AA[i]);
				}catch (IloAlgorithm::NotExtractedException& e){ 
					newCommunity[i] = 0;
				}
				//cout<<"\t - "<<i<<" = "<<newCommunity[i]<<endl;
			}
			cout<<" - New column = "<<newCommunity<<endl;
				


			//if solution is optimal
			//if (k==0 && mauxCplex.getObjValue() >= 0.0){//RC_EPS){
			if (mauxCplex.getObjValue() >= 0.0){//RC_EPS){
				Zs =  new IloNumArray(env, Z.getSize());
				masterCplex.getValues((*Zs), Z);
	cout<<"Zs: "<<*Zs<<endl;
				print(*Zs, A );
				mauxOpt.remove(obj);
				obj.end();
				mauxCplex.end();
				//aux_opt = true;//
				break;
			}

			//pass AA to master model (inserting a new community in model.)
			//
			
			Z.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
			string name = "Z"+to_string(Z.getSize());
			Z[Z.getSize()-1].setName(name.c_str());
			
			//gerar todos os novos Zs
			
			A.add(IloNumArray(env,graph.numberOfNodes));
			ncAPExact++;
			t=A.getSize()-1;
			for (j=0;j<n;j++){
				A[t][j] = newCommunity[j];
			}
	/*cout<<"****************************"<<endl;        
	cout<<"* Updating constraints"<<endl;
	cout<<"****************************"<<endl;        */
			//updating constraints:
			//(l1)
			for (i=0;i<n;i++){
				consLambda1[i].setExpr(
					consLambda1[i].getExpr()+
					A[t][i]*Z[t]
				);
			}
			
			//(l2)
			for (e=0;e<m;e++){
					IloInt a =  graph.edges[e].v1;
					IloInt b =  graph.edges[e].v2;
					IloInt u = min(a,b);
					IloInt v = max(a,b);
					IloExpr exp(env);
					//exp += Y[e] + A[t][u]*A[t][v]*Z[t];
					exp += Y[e] -  A[t][u]*(1-A[t][v])*Z[t] - A[t][v]*(1-A[t][u])*Z[t];
					IloRange rng (env, 0.0, exp, +IloInfinity);
					exp.end();
					consLambda2.add(rng);
			}
			masterOpt.remove(consLambda2);
			IloAdd(masterOpt ,consLambda2); 		
			
			//(l3)
			IloExpr exp3(env);
			exp3 += Wt(graph, A[t])*Z[t]-(1.0 +PAlpha)* W / K;
			IloRange rng3 (env, -IloInfinity, exp3, 0.0);
			exp3.end();
			consLambda3.add(rng3);
			masterOpt.remove(consLambda3);
			IloAdd(masterOpt ,consLambda3);
			
			consLambda4[0].setExpr(consLambda4[0].getExpr()+ Z[t]);
					
			obj.end();
			mauxCplex.end();
			
			/*k++;
			for (int i=0; i<n; i++ ){
				if (newCommunity[i] > 0.9){
					V++;
					mauxOpt.add(AA[i] == 0);
				}
			}*/
			/*if (k==K && V < n){
				break;
			}*/
			
		//}//k<K
		//if(aux_opt == true){break;}
//cout<<"PASS # 369 -- it("<<it<<")"<<endl;
//if(it > 8) {cout<<"EXIT # 335"<<endl;exit(0);  }
	}

//cout<<"EXIT # 371"<<endl;exit(0);            
	//removing constraints
	for (IloInt x = 0; x < nodeConstraints.getSize(); x++){
		 masterOpt.remove( nodeConstraints[x] );
	}

    //BRANCH
    bool isIntegral = true;
    IloNum bestFrac = 0.0;
    IloInt idFrac   = -1;
    for (IloInt c=0; c<Zs->getSize(); c++){
        if (!DBL_EQL((*Zs)[c], 0.0) && !DBL_EQL((*Zs)[c], 1.0)){
            isIntegral = false;
            IloInt intPart = (*Zs)[c];
            if (bestFrac < (*Zs)[c] - intPart){
                bestFrac = (*Zs)[c];
                idFrac = c;
            }
        }
    }
//cout<<"# 375: "<<idFrac<<endl;
    
    if (isIntegral && minValue < lowerbound){
        cout<<"INTEGRAL!"<<endl;
        lowerbound = minValue;
        if (bestZs != NULL){
            delete bestZs;
        }
        bestZs = Zs;
    }else{
		cout<<"BRANCHING: "<<minValue<<endl;
        //branching
        TSolution *left = new TSolution, *right = new TSolution;
        left->value = minValue;// - coef *(*Zs)[idFrac] + 0.0 *coef;
        right->value = minValue;// - coef *(*Zs)[idFrac] + 1.0 *coef;
        right->parentId = left->parentId = bpnode.id;
        left->id = bpnodes.size();
        right->id = bpnodes.size()+1;
        left->fixed.first = idFrac;
        right->fixed.first = idFrac;
        left->fixed.second = 0;
        right->fixed.second = 1;
        bpnodes.push_back(left);
        bpnodes.push_back(right);
        heap.push(*left);
        heap.push(*right);
    }



}






int main(int argc, char **argv){
    srand (time(NULL));

	//parameters
    if (argc < 3){
			cout<<"Error: Three parameters are required: instance file and alpha value."<<endl;
			return 0;
    }
	parameter    = argv[1];
    float PAlpha = atof(argv[2]);
    int   K      = atof(argv[3]);
    
	//instance
    string filepath = "./instances/"+parameter;
    cout<<"Reading the graph at file '"<<filepath<<"'."<<endl;
    LargeGraph graph(filepath);
	
	cout<<"Graph G=(V, E) loaded: "<<endl;
	cout<<" - |V|: "<<graph.numberOfNodes<<endl;
	cout<<" - |E|: "<<graph.numberOfEdges<<endl;


	//best solution
	IloNumArray * bestZs = NULL;

    //constructing starting solution
    IloNum minValue = +IloInfinity;
	
	//starting the chronometer
	chrono::system_clock::time_point antes  = chrono::system_clock::now();
	
    env.setOut(env.getNullStream() );
    try {
        IloInt  i, j, t, tt, m=graph.numberOfEdges, n=graph.numberOfNodes;
        IloNum W;
		//preparing weights
		prepareCe (graph, Ce);
        prepareWv(graph, Wv);


        /// MASTER PROBLEM ///
        IloModel masterOpt (env);

        //creating the constraint coeficients (and two clusters)
        IloArray<IloNumArray> A(env);
        unsigned commStCount = K;
        int csize = n / K ;
        for (t=0; t<K-1; t++){
			A.add(IloNumArray(env,graph.numberOfNodes ));        
			for (unsigned v = 0; v< n; v++){
				if (v >= t*csize && v < (t+1)*csize)
					A[t][v] = 1;
				else
					A[t][v] = 0;
			}
		}
		A.add(IloNumArray(env,graph.numberOfNodes ));  
		for (unsigned v = 0; v< n; v++){
			A[t][v] = 0;
		}
		for (unsigned v = t*csize; v< n; v++){
			A[t][v] = 1;
		}

        //creating vars from model
        //... Z
        IloNumVarArray Z(env, A.getSize(), 0.0, IloInfinity, ILOFLOAT);//deveria ser IloInt, mas por causa do Dual deve-se relaxar essa restrição
        for (IloInt i=0; i<A.getSize(); i++){
			string name = "Z_"+to_string(i+1);
			Z[i].setName(name.c_str());
		}
		//... Ye
		IloNumVarArray Y(env, m, 0.0, IloInfinity, ILOFLOAT);//deveria ser IloInt, mas por causa do Dual deve-se relaxar essa restrição
		for (IloInt i=0; i<m; i++){
			string name ="Y_"+to_string(i+1); 
			Y[i].setName(name.c_str());
		}
	

        //define objective
        IloExpr masterOFExp(env);
        //the main part of the objective function
        for (IloInt e=0; e<m; e++){
            masterOFExp +=   Ce[e] * Y[e];
        }
        IloObjective masterObj = IloAdd(masterOpt, IloMinimize(env, masterOFExp) );
        masterOFExp.end();

        //define master constraints
        //(1) Each node belongs to a single cluster
        IloRangeArray  consLambda1(env);
        for (i=0;i<n;i++){
            IloExpr  masterConsExp(env);
            for (t=0;t<A.getSize();t++){
                masterConsExp += A[t][i]*Z[t] ;

            }
            IloRange rng (env, 1.0, masterConsExp ,1.0);
            masterConsExp.end();
            consLambda1.add(rng);
        }///... add new sum in each constraint when adding a new cluster
        IloAdd(masterOpt ,consLambda1);

        //(2) Is edge Ye intercluster?
        IloRangeArray  consLambda2(env);
        for (t=0;t<A.getSize();t++){        
			for (IloInt e=0; e<m; e++){
				IloInt a =  graph.edges[e].v1;
				IloInt b =  graph.edges[e].v2;
				IloInt u = min(a,b);
				IloInt v = max(a,b);
				//masterOpt.add(Y[e] + A[t][u]*A[t][v]*Z[t] >= 1);
				IloExpr exp(env);
				//exp += Y[e] + A[t][u]*A[t][v]*Z[t];
				exp += Y[e] -  A[t][u]*(1-A[t][v])*Z[t] - A[t][v]*(1-A[t][u])*Z[t];
				IloRange rng (env, 0.0, exp, +IloInfinity);
				exp.end();
				consLambda2.add(rng);
			}
		}
		IloAdd(masterOpt, consLambda2);
		
        //(3) Ensuring the balance suggested by Geomar
        IloRangeArray  consLambda3(env);
        W = 0;
        for (IloInt u=0; u<n; u++){
			W += Wv[u];
		}
        for (t=0;t<A.getSize();t++){
			IloExpr exp(env);
			exp += Wt(graph, A[t])*Z[t] - (1.0 +PAlpha)* W/K;
			IloRange rng (env, -IloInfinity, exp, 0.0);
			exp.end();
			consLambda3.add(rng);
		}
		IloAdd(masterOpt ,consLambda3);
		
		//(4) Ensuring the number of clusters
        IloRangeArray  consLambda4(env);
        IloExpr exp5(env);
        for (t=0;t<A.getSize();t++){
			exp5 += Z[t];
		}
		IloRange rng (env, K, exp5, K);
		exp5.end();
		consLambda4.add(rng);		
		IloAdd(masterOpt, consLambda4);
		
		        
        IloCplex masterCplex(masterOpt);
        masterCplex.setParam(IloCplex::Param::Threads,1);
        masterCplex.setParam(IloCplex::TiLim,36000);


        ///auxiliar problem	
        //define AA_v \in {0,1}^|V|: our a_v of the auxiliar problem
        IloNumVarArray AA(env, graph.numberOfNodes, 0,1, ILOINT);//-IloInfinity,IloInfinity, ILOFLOAT);

        
        //defining the starting node of the BP
        TSolution solInicial;
        solInicial.value = +IloInfinity;
        solInicial.parentId = -1;

        vector<TSolution *> bpnodes;
        solInicial.id = bpnodes.size();
        bpnodes.push_back(&solInicial);

        //the queue
        fibonacci_heap<TSolution> heap;
        heap.push(solInicial);

        IloNum lowerbound = +IloInfinity;
        
        //clusters que não geram descendentes
        IloNumArray disableds(env);

        bool terminate = false;
        long maxIT = 200;
        bestZs =  new IloNumArray(env, A.getSize());
        for (int i=0; i<bestZs->getSize();i++){
            (*bestZs)[i]=1;
        }
        IloInt it = 0;
        while (!terminate && maxIT > 0){
			it ++;
			maxIT --;
            TSolution data = heap.top();
            heap.pop();

            if (data.value > lowerbound){
                break;
            }
			

            solver(            graph,
                               PAlpha,
                               K,
                               W,
                               heap,
                               data,
                               bpnodes,
                               lowerbound,
                               bestZs,
                               A,
                               Z,
                               Y,
                               masterOpt,
                               masterCplex,
                               AA,
                               consLambda1,
							   consLambda2,
							   consLambda3,
							   consLambda4,
                               masterObj,
                               disableds
                           );
cout<<"# 632"<<endl;                           
            if (heap.empty()){
                break;
            }

        }

cout<<(*bestZs)<<endl;

cout<<"# 642"<<endl;                           



    string comms = "";
    for (IloInt i=0;i<bestZs->getSize(); i++){
        string comm = "";
        if ((*bestZs)[i] > 0.9){
            for (IloInt v=0;v<A[i].getSize(); v++){
                if (A[i][v] > 0.9){
                    if (comm == ""){
                        comm = to_string(v);
                    }else{
                        comm += ","+to_string(v);
                    }
                }
            }
            if (comm != ""){
                if (comms  == ""){
                    comms+= "["+comm+"]";
                }else{
                    comms+= ",["+comm+"]";
                }
            }
        }
    }
    
    cout<<"Grupos: "<<comms<<endl;

    /*Solution *sole = new Solution (&graph, comms);
    cout<<"\n*** True: "<<mlg.calculateDensity(sole, Plambda)<<endl;*/
    


//grafo,tipo,todaheuristica,todoprocesso,maxdensity,status,tempo1vez,valor1vez,tempo2vez,valor2vez,tempoRMP,tempoAux,it

    //masterCplex.exportModel("depois.lp");
    //cout<<"MaxDensity: "<<maxDensity;cout.flush();
    string file = "GeomarExato.csv";
    string texto = parameter+";NP_GSR;";
    texto += to_string(PAlpha)+";";
    texto += to_string(K)+";";
    texto += to_string(chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now()-antes).count()/1000.0)+";";
    texto += to_string(lowerbound)+";";
    texto += to_string(masterCplex.getStatus())+";";
    texto += to_string(it)+";";
    texto += comms+"\n";

    ofstream ofs (file, std::fstream::out | std::fstream::app);
    ofs << texto;
    ofs.close();




    env.end();
    //cout<<"\n\n";
    return 0;



}
catch (IloAlgorithm::CannotExtractException& e)
{ std::cerr <<
"CannoExtractException: " << e << std::endl; IloExtractableArray failed = e.getExtractables();

for (IloInt i = 0; i < failed.getSize(); ++i) std::cerr <<
"\t" << failed[i] << std::endl;
// Handle exception ...
}





catch (IloException& ex) {
   cerr << "Error: " << ex << endl;
}
catch (...) {
   cerr << "Error" << endl;
}


    return 0;
}
