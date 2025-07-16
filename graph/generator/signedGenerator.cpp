#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <unordered_set>
#include <chrono>
#include <cstdlib>

using namespace std;

struct Edge{
    unsigned v1;
    unsigned v2;
    int      value;
    Edge(){
        v1=0;
        v2=0;
        value=0;
    }
    Edge(unsigned v1, unsigned v2, int value):
        v1(v1), v2(v2), value(value){}
};

int main(int argc, char*argv[]){
    srand(time(NULL));

    //number of nodes
    int     nodes      = atoi(argv[1]);
    //density of edges inside commmunities
    float   density    = atof(argv[2]);
    //percentual of negative links inside communities
    float   percMinus  = atof(argv[3]);
    //percentual of random relinks
    float   percRelink = atof(argv[4]);

    //name
    string  name       = argv[5];

    int sizeCommunities = nodes/log2f(nodes);

    //communities
    vector< vector<unsigned> > communities;
    for (unsigned node=0; node<nodes;node++){
        int com = node / sizeCommunities;
        if (com == communities.size()){
            communities.push_back(vector<unsigned> ());
        }
        communities[com].push_back(node);
    }
cout<<"50\n";cout.flush();
    //creating edges
    vector< Edge > edgelist;
    unordered_set<string> edgesBelongingList;
    for (unsigned ic = 0; ic< communities.size(); ic++){

        if (communities [ic].size() == 1){
            continue;
        }

        int nodesInside = communities [ic].size();

        int nedges = nodesInside*density;
        for (int e=0; e<nedges;e++){
            unsigned iu = rand()%communities [ic].size();
            unsigned iv = iu;
            while (iv == iu) {
                iv = rand()%communities [ic].size();
            }
            unsigned v1 = communities [ic][iu];
            unsigned v2 = communities [ic][iv];
            int value = 1;

            //if it was not inserted (verification)
            string key = to_string(v1)+","+to_string(v2);
            if (edgesBelongingList.find(key) != edgesBelongingList.end()){
                e--;
                continue;
            }else{
                edgesBelongingList.insert(key);
                edgesBelongingList.insert(to_string(v2)+","+to_string(v1));
            }
            edgelist.push_back(Edge(v1, v2, value));

        }
    }
cout<<"85\n";cout.flush();
    //rewiring
    for (int e=0;e<edgelist.size();e++){
        float ran = rand()%100;
        if ( ran/100 >= percRelink){
            continue;
        }
        unsigned iu = rand()%nodes;
        unsigned iv = iu;
        while (iv == iu) {
            iv = rand()%nodes;
        }
        unsigned v1 = iu;
        unsigned v2 = iv;
        int value = 1;

        //if it was not inserted (verification)
        string key = to_string(v1)+","+to_string(v2);
        if (edgesBelongingList.find(key) != edgesBelongingList.end()){
            e--;
            continue;
        }else{
            edgesBelongingList.insert(key);
            edgesBelongingList.insert(to_string(v2)+","+to_string(v1));
        }
        edgelist[e].v1 = v1;
        edgelist[e].v2 = v2;

    }

    //distribution of negative values
    for (int e=0;e<edgelist.size();e++){
        float ran = rand()%100;
        if ( ran/100 >= percMinus){
            continue;
        }
        edgelist[e].value = -1;
    }

    //printing graph file
    ofstream file;
    file.open(name+".net");
    file<<"*Vertices "<<nodes<<"\n";
    for (int i=0;i<nodes; i++){
        file<<i+1<<"\n";
    }
    file<<"*Edges\n";
    for (int i=0;i<edgelist.size(); i++){
        file<<edgelist[i].v1+1<<" "<<edgelist[i].v2+1<<" "<<edgelist[i].value<<"\n";
    }
    file.close();


    //printing community file
    ofstream cfile;
    cfile.open(name+".cty");
    for (vector<unsigned> cmty:  communities){
        for(unsigned node: cmty){
            cfile<<node<<";";
        }
        cfile<<"\n";
    }
    file.close();

    return 0;
}
