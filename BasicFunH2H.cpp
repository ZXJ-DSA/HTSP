/*
 * BasicFunH2H.cpp
 *
 *  Created on: 31 Aug 2023
 *      Author: zhangmengxuan, Xinjie ZHOU
 */
#include "headH2H.h"

extern vector<int> NodeOrder_;//nodeID order

//// Graph RW
//function for reading graph
void Graph::ReadGraph(string filename){
    ifstream inGraph(filename);
    if(!inGraph){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Reading graph..."<<endl;
    Timer tt;
    tt.start();
    string line;
    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
    node_num=stoi(vs[0]); edge_num=0;
    int tempENum=stoi(vs[1]);
    getline(inGraph,line);
    //graph g initialize
    Neighbor.assign(node_num, vector<pair<vertex,int>>());
    set<int> vertices;

    int ID1,ID2, weight;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]); ID2=stoi(vs[1]); weight=stoi(vs[2]);
//        weight=1;
        if(ID1>=0&&ID1<node_num && ID2>=0&&ID2<node_num && weight>0){
            edge_num++;
            Neighbor[ID1].emplace_back(ID2,weight);
//            Neighbor[ID2].emplace_back(ID1,weight);
            if(vertices.find(ID1)==vertices.end()){//if not found
                vertices.insert(ID1);
            }
        }
        else{
            cout<<"Wrong data! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }

        if(inGraph.eof()) break;
        getline(inGraph,line);
    }
    inGraph.close();
    tt.stop();
    cout<<"Finish Reading! Vertex number: "<<node_num<<"; Edge number: "<<edge_num<<". Time: "<<tt.GetRuntime()<<" s."<< endl;
    if(vertices.size()!=node_num){
        cout<<"vertices is wrong! "<<vertices.size()<<" "<<node_num<<endl; exit(1);
    }
    set<int> mcc;
    DFS_CC(Neighbor,vertices,mcc,node_num);
}
//function of vertex order reading
void Graph::WriteOrder(string filename){
    //Write order file to disk
    ofstream OF(filename);
    if(!OF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Writing vertex order..."<<endl;
    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        OF << i << "\t" << NodeOrder[i] << endl;//ID, order
    }

    OF.close();
    cout<<"Write done."<<endl;
}

void Graph::ReadOrder(string filename){
//    filename="/Users/zhouxj/Documents/1-Research/Datasets/NY/NY_NC_64/vertex_order";
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }
    int nodeNum;
    string line;
    getline(inFile,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    nodeNum=stoi(vs[0]);

    if(node_num==0){
        node_num=nodeNum;
    }else{
        if(nodeNum != node_num){
            cout<<"Wrong vertex number: "<<nodeNum<<" "<<node_num<<endl;
            exit(1);
        }
    }
    cout<<"Vertex number: "<<node_num<<endl;


    NodeOrder.assign(node_num,-1);
    vNodeOrder.assign(node_num,-1);

    getline(inFile,line);

    int ID, order, num=0;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs.size()!=2){
            cout<<"Wrong syntax! "<<line<<endl; exit(1);
        }

        ID=stoi(vs[0]); order=stoi(vs[1]);

        NodeOrder[ID] = order;
        vNodeOrder[order] = ID;
        num++;
        if(inFile.eof())
            break;
        getline(inFile,line);
    }
    if(num!=nodeNum){
        cout<<"Inconsistent! "<<num<< " "<<nodeNum<<endl;
    }

    NodeOrder_=NodeOrder;
}

//// Update RW
void Graph::ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldw;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    IF>>num;
    for(int i=0;i<num;i++){
        IF>>ID1>>ID2>>oldw;
        TestData.push_back(make_pair(make_pair(ID1, ID2), oldw));
    }
    IF.close();
}
void Graph::ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldW,newW;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    string line;
    getline(IF,line);

    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" "),boost::token_compress_on);
        ID1=stoi(vs[0]); ID2=stoi(vs[1]); oldW=stoi(vs[2]); newW=stoi(vs[3]);
        TestData.push_back(make_pair(make_pair(ID1, ID2), make_pair(oldW,newW)));
        if(IF.eof())
            break;
        getline(IF,line);
    }

    IF.close();
}

void Graph::ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData){
    TestData.clear();

    int num, ID1, ID2, oldW,newW1,newW2;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    string line;
    getline(IF,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    num=stoi(vs[0]);
    getline(IF,line);

    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs.size()==5){
            ID1=stoi(vs[0]); ID2=stoi(vs[1]); oldW=stoi(vs[2]); newW1=stoi(vs[3]); newW2=stoi(vs[4]);
            TestData.push_back(make_pair(make_pair(ID1, ID2), make_tuple(oldW,newW1, newW2)));
        }else{
            cout<<"Wrong input! vs.size: "<<vs.size()<<" "<<line<<endl;
        }

        if(IF.eof())
            break;
        getline(IF,line);
    }

    IF.close();
}

//// Dijkstra
//Dijkstra's algorithm
int Graph::Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor){
	benchmark::heap<2, int, int> pqueue(node_num);
	pqueue.update(ID1,0);

	vector<bool> closed(node_num, false);
	vector<int> distance(node_num, INF);
	vector<int> prece(node_num, 0);
	distance[ID1]=0;
	int topNodeID, topNodeDis;
	int NNodeID,NWeigh;

	int d=INF;//initialize d to infinite for the unreachable case

	while(!pqueue.empty()){
		pqueue.extract_min(topNodeID, topNodeDis);
		if(topNodeID==ID2){
			d=distance[ID2];
			break;
		}
		closed[topNodeID]=true;

		for(auto it=Neighbor[topNodeID].begin();it!=Neighbor[topNodeID].end();it++){
			NNodeID=(*it).first;
			NWeigh=(*it).second+topNodeDis;
			if(!closed[NNodeID]){
				if(distance[NNodeID]>NWeigh){
					distance[NNodeID]=NWeigh;
					pqueue.update(NNodeID, NWeigh);
					prece[NNodeID]=topNodeID;
				}
			}
		}
	}
    //retrieve path
//    RetrievePath(ID1, ID2, prece);
	return d;
}

int Graph::BiDijkstra(int ID1, int ID2, vector<vector<pair<vertex,int>>> &Neighbor) {
    if(ID1==ID2) return 0;
    //if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;//to avoid the incorrectness caused by the isolated vertex
    benchmark::heap<2, int, int> queueF(node_num), queueB(node_num);
    queueF.update(ID1,0);
    queueB.update(ID2,0);

    vector<bool> closedF(node_num, false), closedB(node_num, false);
    vector<int> distanceF(node_num, INF), distanceB(node_num, INF);

    distanceF[ID1]=0;
    distanceB[ID2]=0;
    int topNodeIDF, topNodeDisF, topNodeIDB, topNodeDisB;
    int NNodeIDF,NWeighF, NNodeIDB, NWeighB;

    int d=INF;//initialize d to infinite for the unreachable case

    while(!queueF.empty() || !queueB.empty()){
//        if(queueF.top_key()+queueB.top_key()>=d){//termination condition 1
//            return d;
//        }
        //forward
        queueF.extract_min(topNodeIDF, topNodeDisF);
        closedF[topNodeIDF]=true;
        if(closedB[topNodeIDF]){//termination condition 2
            return d;
        }
        for(auto itF=Neighbor[topNodeIDF].begin();itF!=Neighbor[topNodeIDF].end();itF++){
            NNodeIDF=(*itF).first;
            NWeighF=(*itF).second+topNodeDisF;
            if(closedB[NNodeIDF] && NWeighF+distanceB[NNodeIDF]<d){
                d=NWeighF+distanceB[NNodeIDF];
            }
            if(!closedF[NNodeIDF]){
                if(distanceF[NNodeIDF]>NWeighF){
                    distanceF[NNodeIDF]=NWeighF;
                    queueF.update(NNodeIDF, NWeighF);
                }
            }
        }
        //backward
        queueB.extract_min(topNodeIDB, topNodeDisB);
        closedB[topNodeIDB]=true;
        if(closedF[topNodeIDB]){//termination condition 2
            return d;
        }
        for(auto itB=Neighbor[topNodeIDB].begin();itB!=Neighbor[topNodeIDB].end();itB++){
            NNodeIDB=(*itB).first;
            NWeighB=(*itB).second+topNodeDisB;
            if(closedF[NNodeIDB] && NWeighB+distanceF[NNodeIDB]<d){
                d=NWeighB+distanceF[NNodeIDB];
            }
            if(!closedB[NNodeIDB]){
                if(distanceB[NNodeIDB]>NWeighB){
                    distanceB[NNodeIDB]=NWeighB;
                    queueB.update(NNodeIDB, NWeighB);
                }
            }
        }
    }
    return d;
}

//function of retrieving the shortest path
void Graph::RetrievePath(int ID1, int ID2, vector<int> & prece){

	//path retrieval
	vector<int> path;
	path.clear();
	path.push_back(ID2);
	int preID=prece[ID2];
	while(preID!=ID1){
		path.push_back(preID);
		preID=prece[preID];
	}
	path.push_back(ID1);

    pair<int,int> highestVertex(-1,0);//ID, order
	cout<<"path from "<<ID1<<" to "<<ID2<<": "<<endl;
	for(int i=path.size()-1;i>-1;i--){
		cout<<" "<<path[i]<<"("<<NodeOrder[path[i]]<<") ";//<<endl;
        if(NodeOrder[path[i]] > highestVertex.second){
            highestVertex.second = NodeOrder[path[i]];
            highestVertex.first = path[i];
        }
        if(i>0){
            for(int j=0;j<Neighbor[path[i]].size();++j){
                if(Neighbor[path[i]][j].first == path[i-1]){
                    cout<<Neighbor[path[i]][j].second<<endl;
                    break;
                }
            }
        }
	}
	cout<<endl;
    cout<<"Highest-order vertex: "<<highestVertex.first<<" ("<<highestVertex.second<<")"<<endl;
}


//// Other preprocessing
//function of generating update edges
void Graph::UpdateGene(int num, string filename){
    vector<pair<pair<int,int>, pair<int,int>>> UpdateData;

    set<pair<int,int>> Edges;
    vector<pair<pair<int,int>,int>> ENodeID;
    int ID1,ID2,wei;
    for(int i=0;i<Neighbor.size();i++){
        ID1=i;
        for(int j=0;j<Neighbor[i].size();j++){
            ID2=Neighbor[i][j].first;
            wei=Neighbor[i][j].second;
            if(ID1<ID2 && Edges.find(make_pair(ID1,ID2))==Edges.end()){
                Edges.insert(make_pair(ID1,ID2));
                ENodeID.push_back(make_pair(make_pair(ID1,ID2),wei));
            }
            else if(ID2<ID1 && Edges.find(make_pair(ID2,ID1))==Edges.end()){
                Edges.insert(make_pair(ID2,ID1));
                ENodeID.push_back(make_pair(make_pair(ID2,ID1),wei));
            }
        }
    }

    ofstream OF(filename);
    OF<<num<<endl;
    set<int> eid;
    for(int k=0;k<num;k++){
        int edgeid=rand()%ENodeID.size();
        if(eid.find(edgeid)==eid.end()){
            OF<<ENodeID[edgeid].first.first<<" "<<ENodeID[edgeid].first.second<<" "<<ENodeID[edgeid].second<<endl;
            eid.insert(edgeid);
        }else{
            k--;
        }
    }
    OF.close();
}
//function of generating OD pairs
void Graph::ODGene(int num, string filename){
    set<pair<int,int>> ODpair;
    vector<pair<int,int>> ODpairVec;

    srand (time(NULL));
    int s, t;
    for(int i=0;i<num;i++){
        s=rand()%node_num;
        t=rand()%node_num;
        if(ODpair.find(make_pair(s,t))==ODpair.end()){
            ODpairVec.push_back(make_pair(s,t));
            ODpair.insert(make_pair(s,t));
            ODpair.insert(make_pair(t,s));
        }else{
            i--;
        }
    }
    cout<<"generated OD pair number "<<ODpairVec.size()<<endl;

    ofstream OF(filename);
    OF<<ODpairVec.size()<<endl;
    for(int k=0;k<ODpairVec.size();k++){
        OF<<ODpairVec[k].first<<" "<<ODpairVec[k].second<<endl;
    }
    OF.close();
}
//function of connectivity checking
void Graph::StainingMethod(int ID){
    queue<int> Q;

    vector<bool> Stained;
    Stained.assign(node_num, false);

    Q.push(ID);
    Stained[ID]=true;
    int frontid, neiid;
    while(!Q.empty()){
        frontid=Q.front();
        Q.pop();
        for(int k=0;k<Neighbor[frontid].size();k++){
            neiid=Neighbor[frontid][k].first;
            if(!Stained[neiid]){
                Q.push(neiid);
                Stained[neiid]=true;
            }
        }
    }

    int stainNum=0;
    for(int i=0;i<node_num;i++){
        if(Stained[i])
            stainNum+=1;
    }
    //cout<<"Stained Number "<<stainNum<<endl;
    if(stainNum != node_num){
        cout<<"Incorrect!!! stain number: "<<stainNum<<" ; node number: "<<node_num<<endl;
    }

    vector<int> VertexInverted;
    VertexInverted.assign(node_num, -1);
    int j=0;
    for(int i=0;i<node_num;i++){
        if(Stained[i]){
            VertexInverted[i]=j;
            j+=1;
        }
    }
    //cout<<"Check j= "<<j<<", stainNum= "<<stainNum<<endl;

    int Orinode_num=node_num;
    node_num=stainNum;
    vector<vector<pair<vertex,int>>> Neighbor1=Neighbor;
    Neighbor.clear();
    Neighbor.assign(node_num, vector<pair<vertex,int>>());
    int InvertedID, nei, Invertednei, wei;
    for(int ID=0;ID<Orinode_num;ID++){
        if(VertexInverted[ID]!=-1){
            InvertedID=VertexInverted[ID];
            for(int k=0;k<Neighbor1[ID].size();k++){
                nei=Neighbor1[ID][k].first;
                wei=Neighbor1[ID][k].second;
                if(VertexInverted[nei]!=-1){
                    Invertednei=VertexInverted[nei];
                    Neighbor[InvertedID].push_back(make_pair(Invertednei,wei));
                }
            }
        }
    }

}
//function of checking the connectivity, set_A: the vertex set
vector<int> Graph::DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_LCC, int node_num) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    unordered_set<int> set_B;//nodes visited for current component
    set_B.clear();
    int item_id,temp_id;
    vector<bool> flag_visited(node_num,false);
    bool flag_finish = false;
    int temp_num = 0;
    int component_i = 0;
    pair<unordered_set<int>,int> LCC;
    vector<int> CCs;//the vertex size of each connected component

//    for(int i=0;i<node_num;++i){
//        set_A.insert(i);
//    }
    int seed = *set_A.begin();
    stack_A.push(seed);
    set_A.erase(seed);
    set_B.insert(seed);
    flag_visited[seed] = true;
    //Get the connected components by DFS
    while(!set_A.empty()) {//if not finish
        temp_num = 0;
        while (!stack_A.empty()) {
            item_id = stack_A.top();
            stack_A.pop();
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                temp_num += 1;
                if (!flag_visited[temp_id]) {//if not visited
                    stack_A.push(temp_id);
                    set_A.erase(temp_id);
                    set_B.insert(temp_id);
                    flag_visited[temp_id] = true;
                }
            }
        }
        if (set_B.size() > LCC.first.size()) {
            LCC.first.clear();
            LCC.first = set_B;
            LCC.second = temp_num;// /2
        }
        assert(!set_B.empty());
        CCs.push_back(set_B.size());
//        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
//        }
        ++component_i;
        set_B.clear();
        if (!set_A.empty()) {
            stack_A.push(*set_A.begin());
            set_B.insert(*set_A.begin());
            flag_visited[*set_A.begin()] = true;
            set_A.erase(*set_A.begin());
        } else {
            break;
        }
    }
    if(component_i==1){
        cout<<"This graph has only one connected component. ";
        cout<<"Nodes size of graph: "<< LCC.first.size() << " ; ";
        cout<<"Edges size of graph: "<< LCC.second << endl;
    }else{
        cout<<"!!! This graph has "<< component_i <<" connected component!"<<endl;
        cout<<"Nodes size of the largest connected component is: "<<LCC.first.size()<<endl;
        cout<<"Edges size of the largest connected component is: "<<LCC.second<<endl;
    }
    for(auto it=LCC.first.begin();it!=LCC.first.end();++it){
        set_LCC.insert(*it);
    }
    std::sort(CCs.begin(), CCs.end());
    return CCs;
//    return component_i;
}

//function of checking the connectivity, set_A: the vertex set
vector<int> Graph::DFS_CC(vector<vector<pair<int,int>>> & Edges, set<int> set_A, set<int> & set_LCC, int node_num) {
    /// DFS for connected component
    stack<int> stack_A;
//    set<int> set_A;//nodes waiting for visited
    unordered_set<int> set_B;//nodes visited for current component
    set_B.clear();
    int item_id,temp_id;
    vector<bool> flag_visited(node_num,false);
    bool flag_finish = false;
    int temp_num = 0;
    int component_i = 0;
    pair<unordered_set<int>,int> LCC;
    vector<int> CCs;//the vertex size of each connected component

//    for(int i=0;i<node_num;++i){
//        set_A.insert(i);
//    }
    int seed = *set_A.begin();
    stack_A.push(seed);
    set_A.erase(seed);
    set_B.insert(seed);
    flag_visited[seed] = true;
    //Get the connected components by DFS
    while(!set_A.empty()) {//if not finish
        temp_num = 0;
        while (!stack_A.empty()) {
            item_id = stack_A.top();
            stack_A.pop();
            for (auto it = Edges[item_id].begin(); it != Edges[item_id].end(); ++it) {
                temp_id = it->first;
                temp_num += 1;
                if (!flag_visited[temp_id]) {//if not visited
                    stack_A.push(temp_id);
                    set_A.erase(temp_id);
                    set_B.insert(temp_id);
                    flag_visited[temp_id] = true;
                }
            }
        }
        if (set_B.size() > LCC.first.size()) {
            LCC.first.clear();
            LCC.first = set_B;
            LCC.second = temp_num;// /2
        }
        assert(!set_B.empty());
        CCs.push_back(set_B.size());
//        if(!set_B.empty() && set_B.size() < mcc.first.size()){
//            cout<<"Smaller connected component with vertex size "<<set_B.size()<<": ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<*it<<" ";
//            }
//            cout<<"; degree: ";
//            for(auto it=set_B.begin();it!=set_B.end();++it){
//                cout<<Edges[*it].size()<<" ";
//            }
//            cout<<endl;
//        }
        ++component_i;
        set_B.clear();
        if (!set_A.empty()) {
            stack_A.push(*set_A.begin());
            set_B.insert(*set_A.begin());
            flag_visited[*set_A.begin()] = true;
            set_A.erase(*set_A.begin());
        } else {
            break;
        }
    }
    if(component_i==1){
        cout<<"This graph has only one connected component. ";
        cout<<"Nodes size of graph: "<< LCC.first.size() << " ; ";
        cout<<"Edges size of graph: "<< LCC.second << endl;
    }else{
        cout<<"!!! This graph has "<< component_i <<" connected component!"<<endl;
        cout<<"Nodes size of the largest connected component is: "<<LCC.first.size()<<endl;
        cout<<"Edges size of the largest connected component is: "<<LCC.second<<endl;
    }
    for(auto it=LCC.first.begin();it!=LCC.first.end();++it){
        set_LCC.insert(*it);
    }
    std::sort(CCs.begin(), CCs.end());
    return CCs;
//    return component_i;
}

