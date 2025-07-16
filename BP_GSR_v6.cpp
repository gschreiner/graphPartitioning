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
#include <set>
#include <unordered_set>





#include "./graph/largegraph.h"
#include "./graph/solution.h"
#include "./utils/modularitylg.h"
#include "./heuristics/minifastgreedybli.h"


using namespace std;
using namespace std::chrono;
using namespace boost::heap;

ILOSTLBEGIN

#define RC_EPS 1.0e-9

//number of threads
int N_THREADS = 78; ///2;

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

bool P_HEURISTIC_CG = true; //heuristic column generation active?
bool PAR_VERBOSE_ITER = false;

IloNum W = 0;
float PAlpha;
int   K;
IloNumArray Wv(env);
IloNumArray Ce(env);
string parameter = "";
IloNum BigM = 0;

chrono::system_clock::time_point antes;

#define EPSILON 1E-6
bool DBL_EQL(double a, double b) {
    return fabs(a - b) < EPSILON;
}

/*
 * Estrutura para apoiar a verificação se cluster já existem
 */ 
//conjunto de clusters existente
struct TCluster{
	string nodes;
	
	TCluster (IloNumArray &A){
		nodes="";
		for (int i=0; i< A.getSize(); i++){
			int value = A[i];
			nodes+= to_string(value);
		}
	}
	
	bool operator<(TCluster const & rhs) const{
		return nodes < rhs.nodes; 
	}
	bool operator==(const TCluster& other) const {
		return nodes == other.nodes;
    }

};
namespace std {
  template <>
  struct hash<TCluster>
  {
    std::size_t operator()(const TCluster& k) const
    {
      return (hash<string>()(k.nodes));
    }
  };
}

unordered_set < TCluster > clusters;
bool clusterExists (unordered_set < TCluster > &clusters, IloNumArray &A){
	TCluster temp(A);
	return clusters.find(A) != clusters.end();
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
        Wv.add(graph.valueOfNode[i]);
        maxi =  max(maxi, Wv[i]);
    }
    BigM = max(BigM, maxi*graph.numberOfNodes);  

}


static IloNum coefOF( LargeGraph &graph, IloNumArray &A){


    set<int> nodes;

    for (IloInt v=0; v<A.getSize(); v++){
        if(A[v] >= 0.9){
           nodes.insert(v);
        }

    }
    IloNum coef = 0;
	for (IloInt e=0; e<graph.numberOfEdges; e++){
		//cout<<"\t -- Ce["<<e<<"]"<<Ce[e]<<endl;
		int a =  graph.edges[e].v1;
		int b =  graph.edges[e].v2;
		if(A[a] >= 0.9){
           if (nodes.find(b) == nodes.end()){
				coef+= Ce[e];
				//cout<<"\tCe["<<e<<"]"<<Ce[e]<<endl;
		   }
        }else if(A[b] >= 0.9){
				   if (nodes.find(a) == nodes.end()){
						coef+= Ce[e];
				   }
				   //cout<<"\tCe["<<e<<"]"<<Ce[e]<<endl;
				}
        
	}
	return coef;

}





using namespace std;

struct TSolution{
    pair<IloInt, IloInt> fixed;
    long long int parentId;
    long long int id;
    IloNum value;
    IloNum rank;
    bool operator<(TSolution const & rhs) const{
        return rank > rhs.rank; ///inverti do density pq agora temos uma minimização
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




//IloNumArray * solver(
void  solver(
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
            IloNumVar              		&H,
            IloModel                    &masterOpt,
            IloCplex                    &masterCplex,
            IloNumVarArray              &AA,
            IloRangeArray               &consLambda1,
            IloRangeArray               &consLambda2,
            IloRangeArray               &consLambda3,
            IloRangeArray               &consLambda4,
            IloRangeArray               &consLambda5,
            IloObjective                &masterObj,
            IloNumArray 				&disableds
        ){

	//masterOpt.remove(consLambda2);
	IloAdd(masterOpt ,consLambda2);
	//masterOpt.remove(consLambda3);
	IloAdd(masterOpt ,consLambda3);

    //node and solution
    IloConstraintArray nodeConstraints(env);
cout<<"\n\n*******************\n";
    long long int id = bpnode.id;
    cout<<" Value: "<<bpnode.value;
    cout<<id<<" <- ";
    string sconstraints= "";
    while (bpnodes[id]->parentId != -1){
//cout<<"# 182"<<endl;                           
cout<<"::"<<bpnodes[id]->fixed.first<<endl;
        IloConstraint cons = Z[bpnodes[id]->fixed.first] == bpnodes[id]->fixed.second;
        sconstraints += "Z["+to_string(bpnodes[id]->fixed.first)+"] == "+to_string(bpnodes[id]->fixed.second)+"\n";
        masterOpt.add(cons);
        nodeConstraints.add(cons);
        id = bpnodes[id]->parentId;
        cout<<id<<" <- ";
    }
    cout<<"\n\n"<<sconstraints;
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
    IloNumArray lambda5(env,  1); //dual variables from model
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
        cout<<"Starting iteration: "<<it<<""<<endl;	
        
        cout<<"Preparing master model..."<<endl;	
        IloCplex masterCplex(masterOpt);
        masterCplex.setParam(IloCplex::Param::Threads,N_THREADS);
        masterCplex.setParam(IloCplex::TiLim,36000);
        cout<<"Master model prepared."<<endl;
		
		
		//masterCplex.exportModel("antes.lp");
		
        //solve the reduced master problem
        cout<<"Solving master model..."<<endl;	
        masterCplex.solve();
        cout<<"Master model solved."<<endl;	
        
        cout<<"Master status: "<< masterCplex.getStatus()<<endl;
        if(masterCplex.getStatus() == IloAlgorithm::Infeasible){
			//removing constraints
			for (IloInt x = 0; x < nodeConstraints.getSize(); x++){
				 masterOpt.remove( nodeConstraints[x] );
			}
			return;
		}
		cout<<"Master value: "<< masterCplex.getObjValue()<<endl;
		
IloNumArray xZ(env, Z.getSize());
masterCplex.getValues(xZ, Z);
cout<<"Sol.: "<< xZ <<endl;

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
        lambda5[0] = masterCplex.getDual(consLambda5[0]);


cout<<"Lambda_1: "<<lambda1<<endl;	
//cout<<"Lambda_2: "<<lambda2<<endl;	
//cout<<"Lambda_3: "<<lambda3<<endl;	
cout<<"Lambda_4: "<<lambda4<<endl;	
cout<<"Lambda_5: "<<lambda5<<endl;	
		//Heuristic Auxiliary solver
		//if heuristic CG is active
		bool heuristicColumnFound = false;
		if (P_HEURISTIC_CG == true){
			
			SignedMiniFastGreedyBLI heu_gc(&graph, 0.7, W, PAlpha, K);
			heu_gc.execute(Wv, Ce, lambda1, lambda4, lambda5);
			cout<<"Heuri: \t"<<heu_gc.maxValue<<" \t["<<endl;
			for (int v : heu_gc.maxNodes){
				cout<<v<<", ";
			} 
			cout<<endl;
			
			//verificando se cluster já existe:
			IloNumArray cluster(env,graph.numberOfNodes);
			for (j=0;j<graph.numberOfNodes;j++){
				cluster[j] = heu_gc.maxNodes[j];
			}
			if (clusterExists(clusters, cluster)){
				heu_gc.maxValue = 1;
				cout<<"\t\t######### existe igual #######"<<endl;
			}
		
		//if the heuristic has obtained a new promising column
			int checksum = 0;
			for (int v : heu_gc.maxNodes){
				checksum+=v;
			} 

			
			if ( (heu_gc.maxValue < 0) && (checksum > 0 ) && (checksum< graph.numberOfNodes)  ){
				heuristicColumnFound = true;
				
				//adiciona cluster
				clusters.insert(TCluster(cluster));
				
				
				//gerar todos os novos Zs
				cout<<"Preparing new column..."<<endl;
				cout<<"\t preparing new Z..."<<endl;
				Z.add(IloNumVar(env, 0.0, IloInfinity, ILOFLOAT));
				string name = "Z"+to_string(Z.getSize());
				Z[Z.getSize()-1].setName(name.c_str());
				
				
				cout<<"\t preparing new A..."<<endl;
				//A.add(IloNumArray(env,graph.numberOfNodes));
				A.add(cluster);
				ncAPHeuristic++;
				t=A.getSize()-1;
				
				//updating constraints:
				cout<<"\t updating constraing l1..."<<endl;
				//(l1)
				for (i=0;i<n;i++){
					consLambda1[i].setExpr(
						consLambda1[i].getExpr()+
						A[t][i]*Z[t]
					);
				}
				
				//(l2)
				cout<<"\t updating constraing l2..."<<endl;
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
				//masterOpt.remove(consLambda2);
				IloAdd(masterOpt ,consLambda2); 		
				
				//(l3)
				cout<<"\t updating constraing l3..."<<endl;
				IloExpr exp3(env);
				exp3 += Wt(graph, A[t])*Z[t]-(1.0 +PAlpha)* W / K;
				IloRange rng3 (env, -IloInfinity, exp3, 0.0);
				exp3.end();
				consLambda3.add(rng3);
				//masterOpt.remove(consLambda3);
				IloAdd(masterOpt ,consLambda3);
				
				cout<<"\t updating constraing l4 and l5..."<<endl;
				consLambda4[0].setExpr(consLambda4[0].getExpr()+ Z[t]);
				consLambda5[0].setExpr(consLambda5[0].getExpr()- Z[t]);
						
				cout<<"New column added."<<endl;

			}
		}
		
		if (P_HEURISTIC_CG ==false || heuristicColumnFound == false){
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
			mauxOFExp += lambda5[0];
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
				//Do not pass the balance limit
				IloExpr mauxBalLimit(env);
				for (IloInt i=0; i < graph.numberOfNodes; i++) {
					mauxBalLimit += AA[i]*Wv[i];
				}
				mauxBalLimit += - (1.0 +PAlpha)* W/K;
				mauxOpt.add(mauxBalLimit <= 0);

				mauxOFExp.end();
				
			
/*			IloExpr expZero(env);
			for (IloNum v=0;v<n;v++){
				expZero+= AA[v];
			}
			mauxOpt.add(expZero >= 1);
			mauxOpt.add(expZero <= n-1);
*/				
				//solving the column generator
				IloCplex mauxCplex(mauxOpt);
				mauxCplex.setParam(IloCplex::Param::Threads,N_THREADS);
				mauxCplex.setParam(IloCplex::TiLim,36000);
				//mauxCplex.exportModel("aux.lp");
				mauxCplex.solve();
				cout<<"Auxiliary value: "<< mauxCplex.getObjValue()<<endl;

//if (it == 2) exit(0);
				
				
				
				
				   
				//getting the new column
				//copy the new column
				//mauxCplex.getValues(newCommunity, AA);
				for (IloInt i=0; i<AA.getSize(); i++){
					//try{
						newCommunity[i] = mauxCplex.getValue(AA[i]);
					//}catch (IloAlgorithm::NotExtractedException& e){ 
					//	newCommunity[i] = 0;
					//}
					//cout<<"\t - "<<i<<" = "<<newCommunity[i]<<endl;
				}
				cout<<" - New column = "<<newCommunity<<endl;
					
				//cluster
				clusters.insert(TCluster(newCommunity));

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
				consLambda5[0].setExpr(consLambda5[0].getExpr()- Z[t]);
						
				obj.end();
				mauxCplex.end();
			}
			masterCplex.end();
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
    if (idFrac == -1 && isIntegral == false){
		isIntegral = true;
	}
//cout<<"# 375: "<<idFrac<<endl;
    
    if (isIntegral){
		if (minValue < lowerbound ){
			cout<<"INTEGRAL!("<<minValue<<")"<<endl;
			lowerbound = minValue;
			if (bestZs != NULL){
				delete bestZs;
			}
			bestZs = Zs;
		}
    }else{
		cout<<"BRANCHING: value: "<<minValue<<", id: "<<idFrac<<endl;
        //branching
        TSolution *left = new TSolution, *right = new TSolution;
        IloNum coef = coefOF(graph, A[idFrac]);
        /*left->value = -coef;//minValue;// - coef *(*Zs)[idFrac] + 0.0 *coef;
        right->value = coef; //minValue;// - coef *(*Zs)[idFrac] + 1.0 *coef;*/
        left->value = minValue - coef *(*Zs)[idFrac] + 0.0 *coef;
        right->value = minValue - coef *(*Zs)[idFrac] + 1.0 *coef;        
        /*if (right->value == 0 ){
			right->value = -BigM;
		}*/
		
		
		left->rank = (coef+1)*BigM*BigM;
		right->rank = coef;
		
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
        //cout<<"***: "<< left->value <<" / "<< right->value <<" : "<<coef<<endl;
        //exit(0)        ;
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
    PAlpha = atof(argv[2]);
    K      = atof(argv[3]);
    
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
        
		//preparing weights
		prepareCe (graph, Ce);
        prepareWv(graph, Wv);
        W = 0; //total weights
        for (IloInt u=0; u<n; u++){
			W += Wv[u];
		}
		


        /// MASTER PROBLEM ///
        IloModel masterOpt (env);

        //creating the constraint coeficients (and two clusters)
        IloArray<IloNumArray> A(env);
        unsigned commStCount = K;
        int csize = n / K ;
        
        /*for (t=0; t<K-1; t++){
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
		}*/
		
		//greedy initial columns creator
		
		list<int> *ini_c_sets= new list<int>[K];
		int * ini_c_sum = new int[K];
		for (int i=0; i<K; i++){
			ini_c_sum[i]=0;
		}
		fibonacci_heap<pair<float, int>> weights_and_nodes;
		for (IloInt u=0; u<n; u++){
			weights_and_nodes.push(pair<float, int>(Wv[u], u) );
		}
		int ini_c=0;
		while ( weights_and_nodes.empty() == false){
			pair<float, int> data = weights_and_nodes.top(); 
			weights_and_nodes.pop(); 
			while (ini_c_sum[ini_c] >= (1.0 +PAlpha)* W/K){
				ini_c = (ini_c + 1) %K;
			}
			ini_c_sets[ini_c].push_back(data.second);
			ini_c_sum[ini_c] += data.first;
			ini_c = (ini_c + 1) %K;
		}
		IloNumArray X(env,graph.numberOfNodes );
		for (t=0; t<K; t++){
			A.add(IloNumArray(env,graph.numberOfNodes ));        
			for (unsigned v = 0; v< n; v++){
					A[t][v] = 0;
			}
			for (int u : ini_c_sets[t]){
				A[t][u] = 1;
			}
			cout<<A[t]<<endl;
		}
		for (t=0;t<A.getSize();t++){
			cout<<"Soma dos W["<<t+1<<"]: "<<Wt(graph, A[t])<<endl;
		}
		cout<<"Soma dos Ws: "<<W<<endl;
		IloNum solValue = 0;
		for (t=0;t<A.getSize();t++){
			solValue += coefOF(graph, A[t]);
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
		IloNumVar H(env, 0.0, IloInfinity, ILOFLOAT);
		{	string name ="H_"+to_string(i+1); 
			H.setName(name.c_str());
		}

        //define objective
        IloExpr masterOFExp(env);
        //the main part of the objective function
        for (IloInt e=0; e<m; e++){
            masterOFExp +=   Ce[e] * Y[e];
        }
        masterOFExp += BigM * H; //Garantir equilíbrio no número de clusters
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

        for (t=0;t<A.getSize();t++){
			IloExpr exp(env);
			exp += Wt(graph, A[t])*Z[t] - (1.0 +PAlpha)* W/K;
			IloRange rng (env, -IloInfinity, exp, 0.0);
			exp.end();
			consLambda3.add(rng);
		}
		IloAdd(masterOpt ,consLambda3);
		
		//(4 and 5) Ensuring the number of clusters
        IloRangeArray  consLambda4(env);
        IloExpr exp5(env);
        exp5 += H - K;
        for (t=0;t<A.getSize();t++){
			exp5 += Z[t];
		}
		IloRange rng1 (env, 0, exp5, IloInfinity);
		exp5.end();
		consLambda4.add(rng1);		
		IloAdd(masterOpt, consLambda4);
		//(...5) Ensuring the number of clusters
        IloRangeArray  consLambda5(env);
        IloExpr exp6(env);
        exp6 += H + K;
        for (t=0;t<A.getSize();t++){
			exp6 -= Z[t];
		}
		IloRange rng2 (env, 0, exp6, IloInfinity);
		exp6.end();
		consLambda5.add(rng2);		
		IloAdd(masterOpt, consLambda5);
		        
        IloCplex masterCplex(masterOpt);
        masterCplex.setParam(IloCplex::Param::Threads,N_THREADS);
        masterCplex.setParam(IloCplex::TiLim,36000);
        
        


        ///auxiliar problem	
        //define AA_v \in {0,1}^|V|: our a_v of the auxiliar problem
        IloNumVarArray AA(env, graph.numberOfNodes, 0,1, ILOINT);//-IloInfinity,IloInfinity, ILOFLOAT);

        
        //defining the starting node of the BP
        TSolution solInicial;
        solInicial.value = solValue;//+IloInfinity;
        solInicial.parentId = -1;
        solInicial.rank = 0;

        vector<TSolution *> bpnodes;
        solInicial.id = bpnodes.size();
        bpnodes.push_back(&solInicial);

        //the queue
        fibonacci_heap<TSolution> heap;
        heap.push(solInicial);

        IloNum lowerbound = solValue;//+IloInfinity;
        
        
        //cout<<lowerbound<<endl;
        //exit(0);
        
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
cout<<"# 846"<<endl;  
			it ++;
			maxIT --;
            TSolution data = heap.top();
            heap.pop();

            if (data.value > lowerbound){
                break;
            }
cout<<"# 855"<<endl;  			

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
                               H,
                               masterOpt,
                               masterCplex,
                               AA,
                               consLambda1,
							   consLambda2,
							   consLambda3,
							   consLambda4,
							   consLambda5,
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
