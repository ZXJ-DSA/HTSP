/*
 * BasicFun.cpp
 *
 *  Created on: 13 Oct 2022
 *      Author: zhangmengxuan
 */
#include "head.h"

extern vector<int> NodeOrder_;//nodeID order

//// Graph RW
//function for reading graph
void Graph::ReadGraph(string filename){
    ifstream inGraph(filename);
    if(!inGraph){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    cout<<"Reading graph... "<<filename<<endl;
    Timer tt;
    tt.start();
    string line;
    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    node_num=stoi(vs[0]); edge_num=stoull(vs[1]);
    int tempENum=stoi(vs[1]);
    //graph g initialize
    Neighbor.assign(node_num, vector<pair<vertex,int>>());
    vector<map<int,int>> NeighborMapT;
    NeighborMapT.assign(node_num,map<int,int>());
    set<int> vertices;
    unsigned long long edgeNum=0;
    int ID1,ID2, weight;
    while(getline(inGraph,line)){
        if(line.empty()) continue;
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        ID1=stoi(vs[0]); ID2=stoi(vs[1]); weight=stoi(vs[2]);
//        weight=1;
        if(ID1>=0&&ID1<node_num && ID2>=0&&ID2<node_num && weight>0){
            if(NeighborMapT[ID1].find(ID2)==NeighborMapT[ID1].end()){//if not found
                NeighborMapT[ID1].insert({ID2,weight});
                NeighborMapT[ID2].insert({ID1,weight});
                edgeNum+=2;
            }else{
                if(NeighborMapT[ID1][ID2]!=weight){
                    cout<<"Inconsistent graph data!!! "<<ID1<<" "<<ID2<<" "<<NeighborMapT[ID1][ID2]<<" "<<weight<<endl; exit(1);
                }
            }
            if(vertices.find(ID1)==vertices.end()){//if not found
                vertices.insert(ID1);
            }
        }
        else{
            cout<<"Wrong data! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }

    }
    inGraph.close();
    if(edgeNum != edge_num){
        cout<<"Nominate edge number: "<<edge_num<<" ; real edge number: "<<edgeNum<<endl;
        edge_num=edgeNum;
        exit(1);
    }
    for(int i=0;i<NeighborMapT.size();++i){
        for(auto it=NeighborMapT[i].begin();it!=NeighborMapT[i].end();++it){
            Neighbor[i].emplace_back(it->first,it->second);
        }
    }
    tt.stop();
    cout<<"Finish Reading! Vertex number: "<<node_num<<"; Edge number: "<<edge_num<<". Time: "<<tt.GetRuntime()<<" s."<< endl;
    if(vertices.size()!=node_num){
        cout<<"vertices is wrong! "<<vertices.size()<<" "<<node_num<<endl; exit(1);
    }
    set<int> mcc;
    DFS_CC(Neighbor,vertices,mcc,node_num);

    /// Read the coordinates
//    ReadCoordinate(sourcePath+dataset+".co");
}

//function for reading coordinates
void Graph::ReadCoordinate(string filename) {
    ifstream inGraph(filename);
    if(!inGraph){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    cout<<"Reading coordinates..."<<endl;
    Timer tt;
    tt.start();
    string line;
    getline(inGraph,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    int nodeNum=stoi(vs[0]);
    if(node_num!=nodeNum){
        cout<<"Inconsistent vertex number! "<<nodeNum<<" "<<node_num<<endl; exit(1);
    }
    getline(inGraph,line);
    //graph g initialize
    GraphLocation.assign(node_num, pair<int,int>());
    set<int> vertices;

    int ID, lon, lat;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        ID=stoi(vs[0]); lon=stoi(vs[1]); lat=stoi(vs[2]);
        if(ID>=0&&ID<node_num){
            GraphLocation[ID].first=lon, GraphLocation[ID].second=lat;
            if(vertices.find(ID)==vertices.end()){//if not found
                vertices.insert(ID);
            }else{
                cout<<"Already exist! "<<ID<<endl; exit(1);
            }
        }
        else{
            cout<<"Wrong data! "<<ID<<" "<<lon<<" "<<lat<<endl; exit(1);
        }

        if(inGraph.eof()) break;
        getline(inGraph,line);
    }
    inGraph.close();
    tt.stop();
    cout<<"Finish Reading! Vertex number: "<<nodeNum<<". Time: "<<tt.GetRuntime()<<" s."<< endl;
    if(vertices.size()!=node_num){
        cout<<"vertices is wrong! "<<vertices.size()<<" "<<node_num<<endl; exit(1);
    }
}

//function of writing the whole graph into disk
void Graph::WriteGraph(string graphfile){
    cout<<"Writing graph..."<<endl;
    ofstream OF(graphfile);
    if(!OF){
        cout<<"Cannot open "<<graphfile<<endl;
        exit(1);
    }
    OF<<node_num<<" "<<edge_num<<endl;
    int tempE=0;
    for(int ID1=0;ID1<node_num;ID1++){
        for(int j=0;j<Neighbor[ID1].size();j++){
            int ID2=Neighbor[ID1][j].first;
            int wei=Neighbor[ID1][j].second;
            OF<<ID1<<" "<<ID2<<" "<<wei<<endl;
            tempE++;
        }
    }
    OF.close();
    cout<<"Done."<<endl;
}
//function of writing partition result of PostMHL
void Graph::WritePostMHLPartiResult(string filename) {
    //Write order file to disk
    ofstream OF(filename);
    if(!OF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Writing PostMHL partition result..."<<endl;
    int vNum=0;
    OF<<partiNum<<endl;
    for(int pid = 0; pid < partiNum; ++pid){
        OF << PartiVertex[pid].size()<<endl;
        vNum+=PartiVertex[pid].size();
        for(int i=0;i<PartiVertex[pid].size();++i){
            OF << PartiVertex[pid][i] << " ";//ID, order
        }
        OF << endl;
    }
    vNum+=OverlayVertex.size();
    OF<<OverlayVertex.size()<<endl;
    for(int i=0;i<OverlayVertex.size();++i){
        OF << OverlayVertex[i] << " ";//ID, order
    }
    OF << endl;
    OF.close();
    if(vNum != node_num){
        cout<<"Inconsistent vertex number! "<<vNum<<" "<<node_num<<endl; exit(1);
    }
    cout<<"Write done."<<endl;

}
void Graph::ReadPostMHLPartiResult(string filename){
//    filename="/Users/zhouxj/Documents/1-Research/Datasets/NY/NY_NC_64/vertex_order";
    ifstream inFile(filename, ios::in);
    if (!inFile) { // if not exist
        cout << "Fail to open file" << filename << endl;
        exit(1);
    }
//    cout<<"Partition file: "<<filename<<endl;
    int nodeNum;
    string line;
    getline(inFile,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    partiNum=stoi(vs[0]);
    cout<<"Partition number: "<<partiNum<<endl;
    PartiVertex.assign(partiNum,vector<int>());

    int ID, order, num=0;
    for(int pid=0;pid<partiNum;++pid){
        getline(inFile,line);
        assert(!line.empty());
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        int vNum=stoi(vs[0]);

        getline(inFile,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vNum!=vs.size()-1){
            cout<<"PID "<<pid<<" . Inconsistent partition vertex number! "<<vNum<<" "<<vs.size()<<endl;
            exit(1);
        }
        for(int i=0;i<vs.size();++i){
            if(!vs[i].empty()){
                ID=stoi(vs[i]);
                if(ID>=0 && ID<node_num){
                    PartiVertex[pid].push_back(ID);
                }else{
                    cout<<"Wrong vertex id. "<<ID<<endl; exit(1);
                }
            }
//            else{
//                cout<<i<<" "<<vs[i]<<endl;
//                continue;
//            }


        }

        if(inFile.eof())
            break;
    }
    getline(inFile,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    int vNum=stoi(vs[0]);

    getline(inFile,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    if(vNum!=vs.size()-1){
        cout<<"Overlay. Inconsistent vertex number! "<<vNum<<" "<<vs.size()<<endl;
        exit(1);
    }
    for(int i=0;i<vs.size();++i){
        if(!vs[i].empty()){
            ID=stoi(vs[i]);
            if(ID>=0 && ID<node_num){
                OverlayVertex.push_back(ID);
            }else{
                cout<<"Wrong vertex id. "<<ID<<endl; exit(1);
            }
        }


    }
    cout<<"Read done."<<endl;
    inFile.close();
}

//function of vertex order writing
void Graph::WriteOrder(string filename){
    //Write order file to disk
    ofstream OF(filename);
    if(!OF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Writing vertex order... "<< filename<<endl;
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
    cout<<"Order file: "<<filename<<endl;
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



//function for comparing the orders
void Graph::CompareOrder(string filename1, string filename2){
    cout<<"Comparing orders..."<<endl;
    ifstream IF1(filename1);
    if(!IF1){
        cout<<"Cannot open file "<<filename1<<endl;
        exit(1);
    }
    ifstream IF2(filename2);
    if(!IF2){
        cout<<"Cannot open file "<<filename2<<endl;
        exit(1);
    }
    vector<int> order1, order2;//Label1 is the ground truth
    order1.assign(node_num,-1);
    order2.assign(node_num,-1);
    string line;
    int ID,ord;

    //read label 1
    getline(IF1,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    int node_num=stoi(vs[0]);
    assert(node_num == node_num);
    getline(IF1,line);
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs.size()==2){
            ID=stoi(vs[0]), ord=stoi(vs[1]);
            order1[ID] = ord;
        }else{
            cout<<"Wrong syntax! vs.size(): "<<vs.size() <<" "<< line<<endl;
            exit(1);
        }


        if(IF1.eof())
            break;
        getline(IF1,line);
    }
    IF1.close();

    //read label 2
    getline(IF2,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    node_num=stoi(vs[0]);
    assert(node_num == node_num);
    getline(IF2,line);
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        if(vs.size()==2){
            ID=stoi(vs[0]), ord=stoi(vs[1]);
            order2[ID] = ord;
        }else{
            cout<<"Wrong syntax! vs.size(): "<<vs.size() <<" "<< line<<endl;
            exit(1);
        }
        if(IF2.eof())
            break;
        getline(IF2,line);
    }
    IF2.close();

    int ord1,ord2;

    for(int i=0;i<node_num;++i){
        if(order1[i]!=order2[i]){
            cout<<"Inconsistent! "<<i<<" "<<order1[i]<<" "<<order2[i]<<endl;
        }
    }
    cout<<"Done."<<endl;
}

void Graph::GraphPartitionRead(string filename){
    //read the partitioned graphs
//    NeighborMap.assign(node_num,unordered_map<vertex,int>());
    Neighbor.assign(node_num,vector<pair<vertex,int>>());
    NeighborsParti.assign(node_num, vector<pair<vertex,int>>());
    NeighborsOverlay.assign(node_num,unordered_map<vertex,int>());
//    NeighborsOverlay.assign(node_num,vector<pair<vertex,int>>());
    PartiTag.assign(node_num, make_pair(-1,false));
    IDMap.assign(node_num,-1);

    bool flag_minus = false;


    ifstream IF1(filename+"/subgraph_vertex");
    if(!IF1){
        cout<<"Cannot open file "<<"subgraph_vertex"<<endl;
        exit(1);
    }

    int pnum2;
    IF1>>pnum2;
    if(algoParti == "NC"){
//        flag_minus = true;
        partiNum = pnum2;
    }else if(algoParti == "SC" || algoParti == "MT"){
//        flag_minus = true;
//        pnum2 = pnum;
    }
    cout<<"Partition number: "<<pnum2<<endl;

    vector<set<OrderCompMax>> orderPartiV(partiNum,set<OrderCompMax>());

    for(int k=0;k<pnum2;k++){
        int vernum,ID;
        IF1>>vernum;
        for(int i=0;i<vernum;i++){
            IF1>>ID;
//            if(flag_minus){
//                ID = ID-1;
//            }
            assert(ID>=0);
            /*if(ID>=nodenum)
                cout<<"ID "<<ID<<" ,partition ID "<<k<<endl;*/

            if(ID>=0 && ID<node_num){
                if(PartiTag[ID].first==-1){
                    PartiTag[ID].first=k;
                    orderPartiV[k].insert(OrderCompMax(ID));
                }else{
                    cout<<"vertex already in one partition!"<<ID<<" "<<PartiTag[ID].first<<" "<<k<<endl;
                }
            }else{
                cout<<"Wrong vertex ID! "<<ID<<endl; exit(1);
            }

        }
    }
    //further check that each vertex is in one and only one partition
    for(int vid=0;vid<node_num;vid++){
        if(PartiTag[vid].first==-1){
            cout<<"vertex "<<vid<<" not within any partition"<<endl; exit(1);
        }
    }
    //record the vertex to PartiVertex in vertex order: from lower-rank vertex to higher-rank vertex
    PartiVertex.assign(partiNum,vector<vertex>());
    BoundVertex.assign(partiNum,vector<vertex>());
    for(int pid=0;pid<partiNum;++pid){
        int id=0;
        for(auto it=orderPartiV[pid].begin();it!=orderPartiV[pid].end();++it){
//            cout<<pid<<": "<<it->ID<<" "<<NodeOrder[it->ID]<<" "<<PartiTag[it->ID].second<<endl;
            PartiVertex[pid].emplace_back(it->ID);
            IDMap[it->ID] = id;
            ++id;
        }
    }
    orderPartiV.clear();

    ifstream IF(filename+"/subgraph_edge");
    if(!IF){
        cout<<"Cannot open file "<<"subgraph_edge"<<endl;
        exit(1);
    }

    int pnum1;
    IF>>pnum1;
//    if(algoName == "NC"){
//        pnum1 = pnum2;
//    }else if(algoName == "SC" || algoName == "MT"){
//        pnum1 = pnum;
//    }
    for(int k=0;k<pnum1;k++){
        int edgenum0,ID1,ID2,weight;
        IF>>edgenum0;
        for(int i=0;i<edgenum0;i++){
            IF>>ID1>>ID2>>weight;
//            if(flag_minus){
//                ID1 = ID1-1; ID2 = ID2-1;
//            }
            if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
                NeighborsParti[ID1].emplace_back(ID2,weight);
                Neighbor[ID1].emplace_back(ID2,weight);

//                if(NeighborMap[ID1].find(ID2)==NeighborMap[ID1].end()){//not found
//                    NeighborMap[ID1].insert(make_pair(ID2,weight));
//                }else{
//                    cout<<"Wrong for subgraph_edge! edge ("<<ID1<<", "<<ID2<<") already exist!"<<endl; exit(1);
//                }
            }else{
                cout<<"Wrong for subgraph_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
            }


        }
    }

    //read the cut edges
//    CutEdges.clear();
//    set<int> ss; ss.clear();
//    BoundVerSet.assign(pnum,ss);
    set<OrderCompMax> orderOverlayV; orderOverlayV.clear();
    ifstream IF2(filename+"/cut_edges");
    if(!IF2){
        cout<<"Cannot open file "<<"cut_edges"<<endl;
        exit(1);
    }

    int ednum,ID1,ID2,weight;
    int boundaryNum=0;
    IF2>>ednum;
    for(int i=0;i<ednum;i++){
        IF2>>ID1>>ID2>>weight;
//        if(flag_minus){
//            ID1 = ID1-1; ID2 = ID2-1;
//        }

        if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                cout<<"two end points of cut edge are in the same partition"<<endl; exit(1);
            }
            if(!PartiTag[ID1].second || !PartiTag[ID2].second){

            }
            if(!PartiTag[ID1].second){//if not boundary
                PartiTag[ID1].second=true;
                boundaryNum++;
                orderOverlayV.insert(OrderCompMax(ID1));
//                BoundVertex[PartiTag[ID1].first].emplace_back(ID1);
            }
            if(!PartiTag[ID2].second){
                PartiTag[ID2].second=true;
                boundaryNum++;
                orderOverlayV.insert(OrderCompMax(ID2));
//                BoundVertex[PartiTag[ID2].first].emplace_back(ID2);
            }

            NeighborsOverlay[ID1].insert({ID2,weight});
//            NeighborsOverlay[ID1].emplace_back(ID2,weight);
            Neighbor[ID1].emplace_back(ID2,weight);
//            if(NeighborMap[ID1].find(ID2)==NeighborMap[ID1].end()){//not found
//                NeighborMap[ID1].insert(make_pair(ID2,weight));
//            }else{
//                cout<<"Wrong for cut_edge! edge ("<<ID1<<", "<<ID2<<") already exist!"<<endl; exit(1);
//            }
        }else{
            cout<<"Wrong for cut_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }
    }

//    for(int pid=0;pid<partiNum;++pid){
//        for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){
//            if(pid==57||pid==43){
//                cout<<pid<<": "<<*it<<" "<<NodeOrder[*it]<<" "<<PartiTag[*it].second<<endl;
//            }
//        }
//    }

    for(auto it=orderOverlayV.begin();it!=orderOverlayV.end();++it){
//        cout<<PartiTag[it->ID].first<<": "<<it->ID<<" "<<NodeOrder[it->ID]<<" "<<NeighborsOverlay[it->ID].size()<<endl;
        OverlayVertex.emplace_back(it->ID);
        BoundVertex[PartiTag[it->ID].first].emplace_back(it->ID);
    }

    for(int i=0;i<boundaryNum;++i){
        ID1=vNodeOrder[node_num-i-1];
        if(!PartiTag[ID1].second){
            cout<<boundaryNum<<", "<<node_num-i<<": "<<ID1<<"("<<NodeOrder[ID1]<<") is not boundary vertex! "<<endl; exit(1);
        }
    }
    vector<int> bNums;
    for(int i=0;i<partiNum;++i){
        bNums.emplace_back(BoundVertex[i].size());
    }

    cout<<"Overall boundary vertex number: "<<boundaryNum<<"("<<OverlayVertex.size()<<") ; Average boundary vertex number for each partition: "<<boundaryNum/partiNum<<" ; Maximum boundary number: "<<*max_element(bNums.begin(),bNums.end())<<endl;

    //further check the edges
    DD.assign(node_num,0);
    unsigned long long calEdgeNum=0;
    for(int i=0;i<Neighbor.size();i++){
        DD[i]=NeighborsParti[i].size();
        calEdgeNum+=Neighbor[i].size();
        if(IDMap[i]==-1){
            cout<<"Wrong ID map for "<<i<<" "<<IDMap[i]<<endl; exit(1);
        }
    }
    if(edge_num == 0){
        edge_num = calEdgeNum;
    }else{
        if(edge_num!=calEdgeNum){
            cout<<"Inconsistent edge number! "<<edge_num<<" "<<calEdgeNum<<endl; exit(1);
        }
    }
    cout<<"calculated edge number "<<calEdgeNum<<", graph edge number "<<edge_num<<endl;
    cout<<"Partition data finish reading!"<<endl;
}

//// Index RW

//function of writing core label to disk
void Graph::WriteTreeIndexOverlay(string file){
    cout<<"Writing tree index into disk..."<<endl;
    ofstream OF(file+".tree");
    if(!OF){
        cout<<"Cannot open "<<file+".tree"<<endl;
        exit(1);
    }

    OF<<Tree.size()<<endl;
    int tempE=0;
    int ID2,wei;
    for(int ID1=0;ID1<Tree.size();ID1++){
        OF<<Tree[ID1].uniqueVertex;
//        OF<<" "<<Tree[ID1].treeroot;
        OF<<" "<<Tree[ID1].pa;
        OF<<" "<<Tree[ID1].height;
        OF<<" "<<Tree[ID1].hdepth;
        OF<<" "<<Tree[ID1].ch.size();
        for(int i=0;i<Tree[ID1].ch.size();++i){
            OF<<" "<<Tree[ID1].ch[i];
        }
        OF<<" "<<Tree[ID1].vert.size();//vert
        for(int i=0;i<Tree[ID1].vert.size();++i){
            OF<<" "<<Tree[ID1].vert[i].first<<" "<<Tree[ID1].vert[i].second.first<<" "<<Tree[ID1].vert[i].second.second<<" "<<Tree[ID1].pos[i];
        }
        OF<<" "<<Tree[ID1].vAncestor.size()-1;//dis
        for(int i=0;i<Tree[ID1].vAncestor.size()-1;++i){
            OF<<" "<<Tree[ID1].vAncestor[i]<<" "<<Tree[ID1].dis[i]<<" "<<Tree[ID1].cnt[i]<<" "<<Tree[ID1].FN[i];
        }
//        OF<<" "<<Tree[ID1].disInf.size();//disInf
//        for(auto it=Tree[ID1].disInf.begin();it!=Tree[ID1].disInf.end();++it){
//            OF<<" "<<it->first<<" "<<it->second<<" "<<Tree[ID1].FNInf[it->first];
//        }
        OF<<endl;
    }
    OF.close();
    /// Write SCconNodesMT
    ofstream OF2(file+".sc");
    if(!OF2){
        cout<<"Cannot open "<<file+".sc"<<endl;
        exit(1);
    }
    int cid;
//    OF2<<node_num<<endl;
    for(int ID1=0;ID1<node_num;ID1++){
        if(!SCconNodesMT[ID1].empty()){
            OF2<<ID1<<" "<<SCconNodesMT[ID1].size();
            for(auto it=SCconNodesMT[ID1].begin();it!=SCconNodesMT[ID1].end();++it){
                ID2=it->first;
                OF2<<" "<<ID2<<" "<<it->second.size();
                for(auto it2=it->second.begin();it2!=it->second.end();++it2){
                    cid = it2->first; wei = it2->second;
                    OF2<<" "<<cid<<" "<<wei;
                }
            }
        }

    }
    OF2.close();
    cout<<"Done."<<endl;
}
//function of reading core label to disk
void Graph::ReadTreeIndex(string file) {
    cout<<"Reading tree index..."<<endl;
    /// read label
    ifstream IF(file+".tree");
    if(!IF){//if the label file does not exist, construct it
        cout<<"Cannot open file "<<file+".tree"<<endl;
        exit(1);
    }
    Timer tt;
    tt.start();


    int tree_num;
    int tempE=0;
    int ID1,ID2,weight,sz;
    string line;

    getline(IF,line);
    vector<string> vs;
    boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
    tree_num = stoi(vs[0]);

    Tree.reserve(tree_num);

    for(int id=0;id<tree_num;++id){
        getline(IF,line);
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        Tree[id].uniqueVertex=stoi(vs[0]);
//        Tree[id].treeroot=stoi(vs[1]);
        Tree[id].pa=stoi(vs[2]);
        Tree[id].height=stoi(vs[3]);
        Tree[id].hdepth=stoi(vs[4]);
        int index=5;
        //ch
        int sz=stoi(vs[index++]);
        for(int i=0;i<sz;++i){
            Tree[id].ch.emplace_back(stoi(vs[index++]));
        }
        //vert
        int vertS=stoi(vs[index++]);
        for(int i=0;i<vertS;++i){
            pair<int,pair<int,int>> temp;
            temp.first=stoi(vs[index++]);
            temp.second.first=stoi(vs[index++]);
            temp.second.second=stoi(vs[index++]);
            Tree[id].vert.emplace_back(temp);
            Tree[id].pos.emplace_back(stoi(vs[index++]));
        }
        //dis
        int ancS=stoi(vs[index++]);
        Tree[id].pos.emplace_back(ancS);
        for(int i=0;i<ancS;++i){
            Tree[id].vAncestor.emplace_back(stoi(vs[index++]));
            Tree[id].dis.emplace_back(stoi(vs[index++]));
            Tree[id].cnt.emplace_back(stoi(vs[index++]));
            Tree[id].FN.push_back(stoi(vs[index++]));
        }
        Tree[id].vAncestor.emplace_back(Tree[id].uniqueVertex);
        Tree[id].dis.emplace_back(0);
        //disInf
        int intS=stoi(vs[index++]);

//        for(int i=0;i<ancS;++i){
//            ID2=stoi(vs[index++]);
//            Tree[id].disInf.insert({ID2,stoi(vs[index++])});
//            Tree[id].FNInf.insert({ID2,stoi(vs[index++])});
//        }
        Tree[id].vAncestor.emplace_back(Tree[id].uniqueVertex);
        Tree[id].dis.emplace_back(0);
    }

    IF.close();
    /// read SuppPartiID
    ifstream IF2(file+".sc");
    if(!IF2){
        cout<<"Cannot open file "<<file+".sc"<<endl;
        exit(1);
    }

    SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());
    getline(IF2,line);
    int sz1,sz2,cid,wei;
    while(!line.empty()){
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        ID1=stoi(vs[0]); sz1=stoi(vs[1]);
        int index=2;
        for(int i=0;i<sz1;++i){
            ID2= stoi(vs[index++]); sz2=stoi(vs[index++]);
            vector<pair<int,int>> tempV;
            for(int j=0;j<sz2;++j){
                cid=stoi(vs[index++]); wei= stoi(vs[index++]);
                tempV.emplace_back(cid,wei);
            }
            SCconNodesMT[ID1].insert({ID2,tempV});
        }

        if(IF2.eof())
            break;
        getline(IF2,line);
    }

    IF2.close();


    tt.stop();
    cout<<"Time for index reading: "<<tt.GetRuntime()<<" s."<<endl;
}
//function of writing core label to disk
void Graph::WriteTreeIndexParti(string file){
    cout<<"Writing tree index into disk..."<<endl;
    ofstream OF(file+".tree");
    if(!OF){
        cout<<"Cannot open "<<file+".tree"<<endl;
        exit(1);
    }

    OF<<Tree.size()<<endl;
    int tempE=0;
    int ID2,wei;
    for(int ID1=0;ID1<Tree.size();ID1++){
        OF<<Tree[ID1].uniqueVertex;
//        OF<<" "<<Tree[ID1].treeroot;
        OF<<" "<<Tree[ID1].pa;
        OF<<" "<<Tree[ID1].height;
        OF<<" "<<Tree[ID1].hdepth;
        OF<<" "<<Tree[ID1].ch.size();
        for(int i=0;i<Tree[ID1].ch.size();++i){
            OF<<" "<<Tree[ID1].ch[i];
        }
        OF<<" "<<Tree[ID1].vert.size();//vert
        for(int i=0;i<Tree[ID1].vert.size();++i){
            OF<<" "<<Tree[ID1].vert[i].first<<" "<<Tree[ID1].vert[i].second.first<<" "<<Tree[ID1].vert[i].second.second<<" "<<Tree[ID1].pos[i];
        }
        OF<<" "<<Tree[ID1].vAncestor.size()-1;//dis
        for(int i=0;i<Tree[ID1].vAncestor.size()-1;++i){
            OF<<" "<<Tree[ID1].vAncestor[i]<<" "<<Tree[ID1].dis[i]<<" "<<Tree[ID1].cnt[i]<<" "<<Tree[ID1].FN[i];
        }
//        OF<<" "<<Tree[ID1].disInf.size();//disInf
//        for(auto it=Tree[ID1].disInf.begin();it!=Tree[ID1].disInf.end();++it){
//            OF<<" "<<it->first<<" "<<it->second<<" "<<Tree[ID1].FNInf[it->first];
//        }
        OF<<endl;
    }
    OF.close();
    /// Write SCconNodesMT
    ofstream OF2(file+".sc");
    if(!OF2){
        cout<<"Cannot open "<<file+".sc"<<endl;
        exit(1);
    }
    int cid;
//    OF2<<node_num<<endl;
    for(int ID1=0;ID1<node_num;ID1++){
        if(!SCconNodesMT[ID1].empty()){
            OF2<<ID1<<" "<<SCconNodesMT[ID1].size();
            for(auto it=SCconNodesMT[ID1].begin();it!=SCconNodesMT[ID1].end();++it){
                ID2=it->first;
                OF2<<" "<<ID2<<" "<<it->second.size();
                for(auto it2=it->second.begin();it2!=it->second.end();++it2){
                    cid = it2->first; wei = it2->second;
                    OF2<<" "<<cid<<" "<<wei;
                }
            }
        }

    }
    OF2.close();
    cout<<"Done."<<endl;
}
//// Update RW
void Graph::ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData){
    TestData.clear();
    set<pair<int,int>> updateEdges;
    int num, ID1, ID2, oldw;
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    IF>>num;
    for(int i=0;i<num;i++){
        IF>>ID1>>ID2>>oldw;
        if(ID1<ID2){
            if(updateEdges.find(make_pair(ID1,ID2))==updateEdges.end()){//if not found
                TestData.push_back(make_pair(make_pair(ID1, ID2), oldw));
            }else{
                cout<<"Repeated edge update: "<<ID1<<" "<<ID2<<" "<<oldw<<endl;
            }
        }else{
            if(updateEdges.find(make_pair(ID2,ID1))==updateEdges.end()){//if not found
                TestData.push_back(make_pair(make_pair(ID2, ID1), oldw));
            }else{
                cout<<"Repeated edge update: "<<ID2<<" "<<ID1<<" "<<oldw<<endl;
            }
        }

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

int Graph::BiDijkstra(int ID1, int ID2, vector<vector<pair<vertex, int>>> &Neighbor) {
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

//A* algorithm
int Graph::Astar(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor){
    if(ID1==ID2) return 0;
    benchmark::heap<2, int, int> pqueue(node_num);
    int heurisDis=EuclideanDis(ID1,ID2);
    pqueue.update(ID1,heurisDis);

    vector<bool> closed(node_num, false);
    vector<int> distance(node_num, INF);//actual distance from source

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
            NWeigh=(*it).second+distance[topNodeID];
            heurisDis=EuclideanDis(NNodeID,ID2);
            if(!closed[NNodeID]){
                if(distance[NNodeID]>NWeigh){
                    distance[NNodeID]=NWeigh;
                    pqueue.update(NNodeID, NWeigh+heurisDis);
                }
            }
        }
    }
    return d;
}
//naive version
//int Graph::EuclideanDis(int s, int t){
//    int aveLat=(GraphLocation[s].second+GraphLocation[t].second)/2000000;
//    double coeLon=cos(aveLat*0.017453)*0.111319;//
////    int lon=(int)(abs(GraphLocation[s].first-GraphLocation[t].first)*0.083907);
////    cout<<"coeLon: "<<coeLon<<" "<<cos(60)<<" "<<cos(3.1415/3)<<endl;
//    int lon=(int)(abs(GraphLocation[s].first-GraphLocation[t].first)*coeLon);
//    int lat=(int)(abs(GraphLocation[s].second-GraphLocation[t].second)*0.111319);
//    int d=sqrt(lon*lon + lat*lat);
//    return d;
//}
//optimized version
int Graph::EuclideanDis(int s, int t){
    int lon=(int)(abs(GraphLocation[s].first-GraphLocation[t].first)*0.083907);//longitude
    int lat=(int)(abs(GraphLocation[s].second-GraphLocation[t].second)*0.111319);//latitude

    int min,max;
    min=(lat>lon)?lon:lat;
    max=(lat>lon)?lat:lon;
    int approx=max*1007+min*441;
    if(max<(min<<4))
        approx-=max*40;
    return (approx+512)>>10;
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
		cout<<" "<<path[i]<<"("<<PartiTag[path[i]].first<<","<<PartiTag[path[i]].second<<","<<NodeOrder[path[i]]<<") ";//<<endl;
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
//Dijkstra's search in core
int Graph::DijkstraCore(int ID1, int ID2){
    if(!PartiTag[ID1].second || !PartiTag[ID2].second){
        cout<<"Wrong! There is non-boundary vertex! "<<ID1<<"("<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].second<<")"<<endl;
        return -1;
    }
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

        for(auto it=NeighborsOverlay[topNodeID].begin();it!=NeighborsOverlay[topNodeID].end();it++){
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

    srand (0);
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
void Graph::QueryGenerationParti(bool ifSame){
    string partitionfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    string orderfile;

    if(algoChoice==3){
        orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
        int partiNumber=partiNum;
        ReadOrder(orderfile);
        GraphPartitionRead(partitionfile);//read partitions
//    ODGeneSameParti(100000,sourcePath+dataset+"partitions/"+dataset+"_"+algoParti+"_"+ to_string(partiNum)+"/SameParti.query");//same partition
//    ODGeneParti(100000,sourcePath+dataset+"partitions/"+dataset+"_"+algoParti+"_"+ to_string(partiNum)+"/Simulate.query",0.9);//real-world simulation
//    ODGeneCrossParti(100000,sourcePath+dataset+"partitions/"+dataset+"_"+algoParti+"_"+ to_string(partiNum)+"CrossParti.query");//cross-partition simulation
        for(int i=100;i>=50;i-=10){
            double portion=i;
            portion/=100;
            ODGeneParti(100000,sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+ to_string(partiNumber)+"/sameParti_"+to_string(i)+".query",portion);
        }
    }else if(algoChoice==5){
        orderfile= sourcePath + "tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".order";
        int partiNumber=partiNum;

        ifstream IF(orderfile);
        if(!IF){//if the label file does not exist, construct it
            ReadGraph(sourcePath+dataset);//
            TreeDecompositionPartitioning(partiNum,bRatioUpper,bRatioLower);
            WritePostMHLPartiResult(sourcePath + "tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".partiInfo");
        }else{
            ReadOrder(orderfile);
            ReadPostMHLPartiResult(sourcePath + "tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".partiInfo");//read partitions
            IF.close();
        }

        for(int i=100;i>=50;i-=10){
            ODGeneParti(100000, sourcePath + "tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".sameParti_"+to_string(i)+".query",i);
        }
    }

    exit(0);
}
//Function of generating realistic query
void Graph::ODGeneParti(int num, string filename, int portion){
    set<pair<int,int>> ODpair;
    vector<pair<int,int>> ODpairVec;

    srand (0);
    int s, t;
//    double portion=0.99;
    int numSameParti=portion*num/100;
    int numCrossParti=(100-portion)*num/100;
    vector<pair<int,int>> ODpairVecSP, ODpairVecCP;

    cout<<"Generating mixture queries ("<<portion<<"% same-partition queries)..."<<endl;

    for(int i=0;i<numSameParti;i++){
        int pid=rand()%partiNum;
        s=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        t=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        while(s==t){
            t=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        }

        if(ODpair.find(make_pair(s,t))==ODpair.end()){
            ODpairVecSP.emplace_back(s,t);
            ODpair.insert(make_pair(s,t));
            ODpair.insert(make_pair(t,s));
        }else{
            i--;
        }
    }
    cout<<"ODpairVecSP: "<<ODpairVecSP.size()<<endl;
    for(int i=0;i<numCrossParti;i++){
        int pid1=rand()%partiNum;
        int pid2=rand()%partiNum;
        while(pid1==pid2){
            pid2=rand()%partiNum;
        }
        s=PartiVertex[pid1][rand()%PartiVertex[pid1].size()];
        t=PartiVertex[pid2][rand()%PartiVertex[pid2].size()];

        if(ODpair.find(make_pair(s,t))==ODpair.end()){
            ODpairVecCP.emplace_back(s,t);
            ODpair.insert(make_pair(s,t));
            ODpair.insert(make_pair(t,s));
        }else{
            i--;
        }
    }
    cout<<"ODpairVecCP: "<<ODpairVecCP.size()<<endl;
    int base=100-portion;
    cout<<"Base: "<<base<<endl;
    int sp_i=0; int cp_i=0;
    for(int i=0;i<num;++i){
        if(i%100<base && i%100>=0){
            ODpairVec.emplace_back(ODpairVecCP[cp_i]);
            cp_i++;
        }else{
            ODpairVec.emplace_back(ODpairVecSP[sp_i]);
            sp_i++;
        }
    }
    cout<<"generated OD pair number "<<ODpairVec.size()<<" "<<sp_i<<" "<<cp_i<<endl;

    ofstream OF(filename);
    OF<<ODpairVec.size()<<endl;
    for(int k=0;k<ODpairVec.size();k++){
        OF<<ODpairVec[k].first<<" "<<ODpairVec[k].second<<endl;
    }
    OF.close();
    cout<<"Finish."<<endl;
}
//Function of generating realistic query
void Graph::ODGeneSameParti(int num, string filename){
    set<pair<int,int>> ODpair;
    vector<pair<int,int>> ODpairVec;

    srand (0);
    int s, t;
    cout<<"Generating same-partition queries..."<<endl;

    for(int i=0;i<num;i++){
        int pid=rand()%partiNum;
        s=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        t=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        while(s==t){
            t=PartiVertex[pid][rand()%PartiVertex[pid].size()];
        }

        if(ODpair.find(make_pair(s,t))==ODpair.end()){
            ODpairVec.emplace_back(s,t);
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
    cout<<"Finish."<<endl;
}
//Function of generating realistic query
void Graph::ODGeneCrossParti(int num, string filename){
    set<pair<int,int>> ODpair;
    vector<pair<int,int>> ODpairVec;

    srand (0);
    int s, t;
    int pid1,pid2;
    cout<<"Generating cross-partition queries..."<<endl;
    for(int i=0;i<num;i++){
        pid1=rand()%partiNum;
        pid2=rand()%partiNum;
        while(pid1==pid2){
            pid2=rand()%partiNum;
        }
        s=PartiVertex[pid1][rand()%PartiVertex[pid1].size()];
        t=PartiVertex[pid2][rand()%PartiVertex[pid2].size()];

        if(ODpair.find(make_pair(s,t))==ODpair.end()){
            ODpairVec.emplace_back(s,t);
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
    cout<<"Finish."<<endl;
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
