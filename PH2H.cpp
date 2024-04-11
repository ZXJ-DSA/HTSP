/*
 * Construction.cpp
 *
 *  Created on: 24 August 2023
 *      Author: Xinjie ZHOU
 */
#include "head.h"
#include "PH2H.hpp"

/// Index Construction
void Graph::IndexConstruction(){
    if(algoChoice==1){
        cout<<"System Index: CH + H2H."<<endl;
        HybridSPIndexConstruct();
    }else if(algoChoice==2){
        cout<<"System Index: PH2H."<<endl;
        HybridPSPIndexConstruct();
//        PH2HIndexConstruct();
    }else if(algoChoice==3){
        cout<<"System Index: PCH + PH2H."<<endl;
        PMHLIndexConstruct();
    }else if(algoChoice==4){
        cout<<"System Index: PCH + PH2H with Optimizations."<<endl;
        PMHLIndexConstructOpt();
    }else if(algoChoice==5){
        cout<<"System Index: PostMHL index."<<endl;
        cout<<"Bandwidth: "<<bandWidth<<endl;
        PostMHLIndexConstruct();
    }else if(algoChoice==0){
        cout<<"A* search."<<endl;
        ReadGraph(sourcePath+dataset);//
    }
}
//function of hybrid multi-stage SP index construction
void Graph::HybridSPIndexConstruct(){
    string orderfile=sourcePath+dataset+".order";
    orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";

    double runT1=0, runT2=0, runT3=0;
    Timer tt;

    ReadGraph(sourcePath+dataset);//

    tt.start();
    MDEContraction(orderfile);
    tt.stop();
    runT1=tt.GetRuntime();
    cout<<"Time for MDE contraction: "<<runT1<<" s."<<endl;

    tt.start();
    makeTree();
    tt.stop();
    runT2=tt.GetRuntime();
    cout<<"Time for Tree construction: "<<runT2<<" s."<<endl;

    tt.start();
    makeIndex();
    tt.stop();
    runT3=tt.GetRuntime();
    cout<<"Time for Index building: "<<runT3<<" s."<<endl;

    cout<<"Overall index construction time: "<<runT1+runT2+runT3<<" s."<<endl;

    IndexsizeH2H();
}

//function for MDE contraction
void Graph::MDEContraction(string orderfile){
    cout<<"MDE contraction..."<<endl;
    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }
    ifstream IF(orderfile);
//    if(true){
    if(!IF.is_open()){/// if no order file, use MDE to generate order
        cout<<"Cannot open vertex ordering file "<<orderfile<<endl;
        int Twidth=0;//tree width
        //initialize SCconNodesMT
        SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());

        //initialize E
        map<int,pair<int,int>> m;
        E.assign(node_num,m);
        for(int i=0;i<Neighbor.size();i++){
            for(int j=0;j<Neighbor[i].size();j++)
                E[i].insert(make_pair(Neighbor[i][j].first,make_pair(0,1)));
        }

        _DD_.assign(node_num,0);
        DD.assign(node_num,0);

        set<DegComp> Deg;
        int degree;
        for(int i=0;i<node_num;i++){
            degree=Neighbor[i].size();
            if(degree!=0){
                _DD_[i]=degree;
                DD[i]=degree;
                Deg.insert(DegComp(i));
            }
        }

        vector<bool> exist; exist.assign(node_num,true);
        vector<bool> change; change.assign(node_num,false);

        vector<pair<int,pair<int,int>>> vect;
        NeighborCon.assign(node_num,vect); //NeighborCon.clear();
        //SCconNodes.clear();

        //cout<<"Begin to contract"<<endl;
        int count=0;

        while(!Deg.empty()){
            if(count%10000==0)
                cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
            count+=1;
            int x=(*Deg.begin()).x;

            while(true){
                if(change[x]){
                    Deg.erase(DegComp(x));
                    _DD_[x]=DD[x];
                    Deg.insert(DegComp(x));
                    change[x]=false;
                    x=(*Deg.begin()).x;
                }else
                    break;
            }

            vNodeOrder.push_back(x);
            Deg.erase(Deg.begin());
            exist[x]=false;

            vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

            for(auto it=E[x].begin();it!=E[x].end();it++){
                if(exist[(*it).first]){
                    Neigh.push_back(*it);
                }
            }

            if(Neigh.size()>Twidth)
                Twidth=Neigh.size();

            NeighborCon[x].assign(Neigh.begin(),Neigh.end());

            //multi threads for n^2 combination
            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteEOrderGenerate(x,y);
                change[y]=true;
            }

            int stepf=Neigh.size()/threadnum;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*stepf;
                if(i==threadnum-1)
                    p.second=Neigh.size();
                else
                    p.second=(i+1)*stepf;
                threadf.add_thread(new boost::thread(&Graph::NeighborComOrderGenerate, this, boost::ref(Neigh), p));
            }
            threadf.join_all();
        }

        NodeOrder.assign(node_num,-1);
        for(int k=0;k<vNodeOrder.size();k++){
            NodeOrder[vNodeOrder[k]]=k;
        }
        ofstream ofile(orderfile);
        ofile << node_num << endl;
        for(int i = 0; i < NodeOrder.size(); i++)
            ofile << i << " " << NodeOrder[i] << endl;
        ofile.close();
        cout<<"Finish Contract"<<" , treewidth "<<Twidth<<endl;
        exit(0);
    }
    else{///if there is an order file
        cout<<"Reading vertex ordering..."<<endl;
        NodeOrder.assign(node_num, -1);
        vNodeOrder.assign(node_num, -1);
        int num, nodeID, nodeorder;
        string line;
        getline(IF,line);
        vector<string> vs;
        boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
        num=stoi(vs[0]);
        for(int i=0;i<num;i++){
            vs.clear();
            getline(IF,line);
            boost::split(vs,line,boost::is_any_of(" \t"),boost::token_compress_on);
            if(vs.size()!=2){
                cout<<"Wrong syntax for ordering. "<<line<<endl; exit(1);
            }
            nodeID=stoi(vs[0]);nodeorder=stoi(vs[1]);

            NodeOrder[nodeID]=nodeorder;
            if(nodeorder!=-1){
                vNodeOrder[nodeorder]=nodeID;
            }else{
                cout<<"Wrong order! "<<nodeID<<" "<<nodeorder<<endl; exit(1);
            }
        }
        IF.close();
        unordered_set<int> vertices; vertices.clear();
        for(int i=0;i<node_num;++i){
            if(vertices.find(vNodeOrder[i])==vertices.end()){//if not found
                vertices.insert(vNodeOrder[i]);
            }
        }
        if(vertices.size()!=node_num){
            cout<<"Order wrong! "<<vertices.size()<<" "<<node_num<<endl; exit(1);
        }
        vertices.clear();

        for(int i=0;i<2;++i){
            int id=vNodeOrder[node_num-1-i];
            cout<<"Order "<<node_num-1-i<<": "<<id<<" "<<Neighbor[id].size()<<endl;
        }
        for(int i=0;i<2;++i){
            int id=vNodeOrder[i];
            cout<<"Order "<<i<<": "<<id<<" "<<Neighbor[id].size()<<endl;
        }

        vector<pair<int,pair<int,int>>> vect;
        NeighborCon.assign(node_num,vect);

        SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());//record the supportive vertices of a shortcut, only record edge once by leveraging the ID positions of endpoints

        //initialize E
        map<int,pair<int,int>> m;
        E.assign(node_num,m);
        for(int i=0;i<Neighbor.size();i++){
            for(int j=0;j<Neighbor[i].size();j++)
                E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
        }

        vector<bool> exist; exist.assign(node_num,true);
        //vector<bool> change; change.assign(nodenum,false);

        //cout<<"Begin to contract"<<endl;
        for(int nodeorder=0;nodeorder<node_num;nodeorder++){//start from the most important vertex
            int x=vNodeOrder[nodeorder];
            if(x!=-1){//to identify and exclude the isolated vertices
                exist[x]=false;

                vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

                for(auto it=E[x].begin();it!=E[x].end();it++){
                    if(exist[(*it).first]){
                        Neigh.push_back(*it);
                    }
                }
                NeighborCon[x].assign(Neigh.begin(),Neigh.end());

                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteEorder(x,y);
                    //change[y]=true;
                }

                int ID1,ID2;
                if(Neigh.size()<=100){
//                if(true){
                    //single thread
                    for(int i=0;i<Neigh.size();i++){
                        ID1=Neigh[i].first;
                        for(int j=i+1;j<Neigh.size();j++){
                            ID2=Neigh[j].first;
                            int temp=Neigh[i].second.first+Neigh[j].second.first;
                            insertEorder(ID1,ID2,temp);
                            if(ID1<ID2){
                                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                                    SCconNodesMT[ID1].insert({ID2,vector<pair<int,int>>()});
                                }
                                SCconNodesMT[ID1][ID2].emplace_back(x,temp);//only record onece
                            }
                            else if(ID2<ID1){
                                if(SCconNodesMT[ID2].find(ID1)==SCconNodesMT[ID2].end()){//if not found
                                    SCconNodesMT[ID2].insert({ID1,vector<pair<int,int>>()});
                                }
                                SCconNodesMT[ID2][ID1].emplace_back(x,temp);
                            }

                        }
                    }
                }else{
                    if(Neigh.size()>threadnum){
                        int step=Neigh.size()/threadnum;
                        boost::thread_group thread;
                        for(int i=0;i<threadnum;i++){
                            pair<int,int> p;
                            p.first=i*step;
                            if(i==threadnum-1)
                                p.second=Neigh.size();
                            else
                                p.second=(i+1)*step;
                            thread.add_thread(new boost::thread(&Graph::NeighborComorderH2H, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }else{
                        boost::thread_group thread;
                        for(int i=0;i<Neigh.size();i++){
                            pair<int,int> p;
                            p.first=i; p.second=(i+1);
                            thread.add_thread(new boost::thread(&Graph::NeighborComorderH2H, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }

                }

            }
            else{
                cout<<"Wrong order! "<<x<<" "<<nodeorder<<endl; exit(1);
            }
        }
    }
    NodeOrder_ = NodeOrder;
}

void Graph::insertEMTorder(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }
}

void Graph::NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
//            vSm[ID1]->wait();
//            vSm[ID2]->wait();
            insertEMTOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
//            vSm[ID1]->notify();
//            vSm[ID2]->notify();
        }
    }
//    sm->notify();
}

void Graph::insertEMTOrderGenerate(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//        DD2[u]++;
    }
}

void Graph::deleteEorder(int u,int v){
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        //DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        //DD[v]--;
    }
}

void Graph::insertEorder(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        //DD[u]++;
        //DD2[u]++;
    }
    else{
        if(E[u][v].first>w)
            E[u][v]=make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        //DD[v]++;
        //DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}

//function for H2H tree construction
void Graph::makeTree(){
    cout<<"Building H2H tree..."<<endl;
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(node_num,vecemp);//record the tree node id whose unique vertex involves this vertex as neighbor

    rank.assign(node_num,0);
    //Tree.clear();
    int len=vNodeOrder.size()-1;
    heightMax=0;

    Node rootn;
    int x=vNodeOrder[len];//from the highest vertex
    //cout<<"len "<<len<<" , ID "<<x<<endl;
    while(x==-1){//to skip those vertices whose ID is -1
        len--;
        x=vNodeOrder[len];
        //cout<<"len "<<len<<" , ID "<<x<<endl;
    }
    rootn.vert=NeighborCon[x];
    rootn.uniqueVertex=x;
    rootn.pa=-1;
    rootn.height=1;
    rank[x]=0;
    Tree.push_back(rootn);
    len--;

    int nn;
    for(;len>=0;len--){
        int x=vNodeOrder[len];
        Node nod;
        nod.vert=NeighborCon[x];
        nod.uniqueVertex=x;
        int pa=match(x,NeighborCon[x]);
        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        nod.hdepth=Tree[pa].height+1;

        for(int i=0;i<NeighborCon[x].size();i++){
            nn=NeighborCon[x][i].first;
            VidtoTNid[nn].push_back(Tree.size());
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1)
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
        }
        if(nod.height>heightMax) {
            heightMax=nod.height;
        }
        rank[x]=Tree.size();

        Tree.push_back(nod);
        //cout<<"len "<<len<<" , ID "<<x<<endl;
    }
}
int Graph::match(int x,vector<pair<int,pair<int,int>>> &vert){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(rank[vert[i].first]>rank[nearest])
            nearest=vert[i].first;
    }
    int p=rank[nearest];
    return p;
}
//function of H2H index construction
void Graph::makeIndex(){
    cout<<"Building H2H index..."<<endl;
    makeRMQ(toRMQ, RMQIndex, Tree);//build LCA index

    //initialize
    vector<int> list; //list.clear();
    list.push_back(Tree[0].uniqueVertex);
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);

    for(int i=0;i<Tree[0].ch.size();i++){
        makeIndexDFS(Tree[0].ch[i],list);
    }

}

//function for computing the index size
void Graph::IndexsizeH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0;
    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*2*sizeof(int);//dis
        m3+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].cnt.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count

    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"Distance label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"H2H label size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"CH Update information size: "<<(double)m4/1024/1024<<" MB"<<endl;
    cout<<"H2H Update information size: "<<(double)(m2+m4)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::PMHLIndexConstruct() {
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
    orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    ReadOrder(orderfile);

    string partitionfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
    tt.start();
    /// Partition index and Overlay graph construction
//    Construct_PartiIndex(false);
    Construct_PartiIndex(true, true);
//    PMHLConstructPartiIndexCH(true);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition index construction time: "<<runT1<<" s"<<endl;

    /// Overlay graph construction
    tt.start();
    Construct_OverlayGraphNoAllPair(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
    Construct_OverlayIndex(true);
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


    /// Partition index repair
    if(algoUpdate>=PH2H_Post){
        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true, true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
//        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)

    if(algoUpdate==PH2H_Cross){
        tt.start();
//        ConstructExtensionLabelsNoAllPair();//correct version
        ConstructExtensionLabelsNoAllPairTopDown();
        tt.stop();
        runT5=tt.GetRuntime();
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
    }

    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}

void Graph::PMHLIndexConstructOpt() {
    double runT0, runT1, runT2, runT3, runT4, runT5;
    runT0=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
    orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    ReadOrder(orderfile);

    string partitionfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
//    tt.start();
//    PreConstructAllPairs(true);
//    tt.stop();
//    runT0 = tt.GetRuntime();
//    cout<<"All-pair boundary edges insertion time: "<<runT0<<" s"<<endl;

    tt.start();
    /// Partition index and Overlay graph construction
    Construct_PartiIndex(true, false);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition shortcuts construction time: "<<runT1<<" s"<<endl;

    /// Overlay graph construction
    tt.start();
    Construct_OverlayGraphNoAllPair(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
    Construct_OverlayIndex(true);
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay shortcuts construction time: "<<runT3<<" s."<< endl;
    algoQuery=PCH_No;

    /// Partition index repair
    if(algoUpdate>=PH2H_Post){
        /// partition shortcuts update and label construction
        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true, true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;
        algoQuery=PH2H_Post;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)

    if(algoUpdate==PH2H_Cross){
        tt.start();
//        ConstructExtensionLabelsNoAllPair();//correct version
        ConstructExtensionLabelsNoAllPairTopDown();
        tt.stop();
        runT5=tt.GetRuntime();
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
        algoQuery=PH2H_Cross;
    }

    cout<<"Overall Construction Time: "<<runT0+runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}

void Graph::PostMHLIndexConstruct() {
    double runT0, runT1, runT2, runT3, runT4, runT5;
    runT0=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
    Timer tt;
    ReadGraph(sourcePath+dataset);//

    vSm.reserve(node_num);
    for(int i = 0; i < node_num; i++)
    {
        Semaphore* s = new Semaphore(1);
        vSm.push_back(s);
    }
//    tt.start();
//    PreConstructAllPairs(true);
//    tt.stop();
//    runT0 = tt.GetRuntime();
//    cout<<"All-pair boundary edges insertion time: "<<runT0<<" s"<<endl;

    tt.start();
    /// Partition and Overlay graph construction
//    TreeDecompositionPartitioningNaive();
    TreeDecompositionPartitioning(partiNum,bRatioUpper,bRatioLower);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Tree-Decomposition based Partitioning time: "<<runT1<<" s"<<endl;
    // write the partition results to disk
    string partiF = sourcePath + "tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".partiInfo";
    ifstream IF(partiF);
    if(!IF){//if the label file does not exist, construct it
        WritePostMHLPartiResult(sourcePath + "tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".partiInfo");
    }
    IF.close();

    algoQuery=PCH_No;


    /// Overlay index construction
    tt.start();
    PostMHLIndexConstructOverlay();
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;
//    algoQuery=PCH_No;

    /// Partition index construction, single thread
    if(algoUpdate>=PH2H_Post){
        double tPost=0;
        tt.start();
        PostMHLIndexConstructPost(true,tPost);
//        PostMHLIndexConstructPost(false);
        tt.stop();
        runT4 = tt.GetRuntime();
        algoQuery=PH2H_Post;
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;
    }
    if(algoUpdate==PH2H_Cross){
        double tCross=0;
        tt.start();
        PostMHLIndexConstructExtend(true,tCross);
//        PostMHLIndexConstructExtend(false);
        tt.stop();
        runT5=tt.GetRuntime();
        algoQuery=PH2H_Cross;
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
    }


    cout<<"Overall Construction Time: "<<runT0+runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePostMHL();
}



// partitioned SP
void Graph::HybridPSPIndexConstruct(){
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
    orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    ReadOrder(orderfile);

    string partitionfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
    tt.start();
    /// Partition index and Overlay graph construction
    Construct_PartiIndex(true, true);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition index construction time: "<<runT1<<" s"<<endl;


    /// Overlay graph construction
    tt.start();
    Construct_OverlayGraph(true);
//    Construct_OverlayGraphNoAllPair(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
    Construct_OverlayIndex(true);
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


    /// Partition index repair
    if(algoUpdate>=PH2H_Post){
        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true, true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)

    if(algoUpdate==PH2H_Cross){
        tt.start();
        ConstructExtensionLabels();
        tt.stop();
        runT5=tt.GetRuntime();
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
    }

    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}

void Graph::PH2HIndexConstruct(){
    double runT1, runT2, runT3, runT4, runT5;
    runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;

    /// Read order and partitions
    string orderfile;
    orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum)+"/vertex_orderMDE2";
    ReadOrder(orderfile);

    string partitionfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
    GraphPartitionRead(partitionfile);//read partitions

//    vSm.reserve(node_num);
//    for(int i = 0; i < node_num; i++)
//    {
//        Semaphore* s = new Semaphore(1);
//        vSm.push_back(s);
//    }

    Timer tt;
    tt.start();
    /// Partition index and Overlay graph construction
    Construct_PartiIndex(true,true);
    tt.stop();
    runT1 = tt.GetRuntime();
    cout<<"Partition index construction time: "<<runT1<<" s"<<endl;


    /// Overlay graph construction
    tt.start();
    Construct_OverlayGraph(true);
    tt.stop();
    runT2 = tt.GetRuntime();
    cout<<"Overlay graph construction time: "<<runT2<<" s."<< endl;

    /// Overlay index construction
    tt.start();
    Construct_OverlayIndex(true);
    tt.stop();
    runT3 = tt.GetRuntime();
    cout<<"Overlay index construction time: "<<runT3<<" s."<< endl;


    /// Partition index repair
    if(algoUpdate>=PH2H_Post){
        tt.start();
        ConstructPartitionPost(true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition construction: "<<tt.GetRuntime()<<" s."<<endl;

        tt.start();
        ConstructPartitionPostIndex(true, true);
        tt.stop();
        runT4 += tt.GetRuntime();
        cout<<"Time for post partition index construction: "<<tt.GetRuntime()<<" s."<<endl;

//        tt.start();
//        Repair_PartiIndex(true);
////        Repair_PartiIndex(false);
//        tt.stop();
//        runT3 = tt.GetRuntime();
        cout<<"Partition index repair time: "<<runT4<<" s"<<endl;
    }
//    Compute_tree_label(ifParallel, ifOpt);//Construct periphery index (H2H label + interface label)

    if(algoUpdate==PH2H_Cross){
        tt.start();
        ConstructExtensionLabels();
        tt.stop();
        runT5=tt.GetRuntime();
        cout<<"Extended label construction time: "<<runT5<<" s"<<endl;
    }

    cout<<"Overall Construction Time: "<<runT1+runT2+runT3+runT4+runT5<<" s."<<endl;

    IndexSizePH2H();//index (labeling+pruning point) size computation
}


//function for computing the index size
void Graph::IndexSizePH2H(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;

    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*sizeof(int);//dis
        m1+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].dis.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count
    }
    //overlay
    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m2+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    if(algoUpdate>=PH2H_Post){
//        if(algoUpdate==PCH_Post){
//            for(int pid=0;pid<TreesPost.size();++pid){
//                for(int i=0;i<TreesPost[pid].size();i++){
////                    m3+=TreesPost[pid][i].dis.size()*sizeof(int);//dis
////                    m3+=TreesPost[pid][i].pos.size()*sizeof(int);//pos
////                    m4+=TreesPost[pid][i].dis.size()*sizeof(int);//cnt
//                    m4+=TreesPost[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
//                }
//            }
//        }else
        if(algoUpdate==PH2H_Post){
            for(int pid=0;pid<TreesPost.size();++pid){
                for(int i=0;i<TreesPost[pid].size();i++){
                    m3+=TreesPost[pid][i].dis.size()*sizeof(int);//dis
                    m3+=TreesPost[pid][i].pos.size()*sizeof(int);//pos
                    m4+=TreesPost[pid][i].dis.size()*sizeof(int);//cnt
                    m4+=TreesPost[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
                }
            }
        }
        if(algoUpdate==PH2H_Cross){
            for(int pid=0;pid<TreesPost.size();++pid){
                for(int i=0;i<TreesPost[pid].size();i++){
                    m5+=TreesPost[pid][i].dis.size()*sizeof(int);//dis
                    m5+=TreesPost[pid][i].pos.size()*sizeof(int);//pos
                    m5+=TreesPost[pid][i].vAncestor.size()*sizeof(int);//cnt
                }
            }
        }

        for(int i=0;i< SCconNodesMTPost.size();i++){
            for(auto it=SCconNodesMTPost[i].begin(); it!=SCconNodesMTPost[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }

    }
    else{
        if(algoUpdate==PH2H_No || algoUpdate>=PH2H_Post){
            //partitions
            for(int pid=0;pid<Trees.size();++pid){
                for(int i=0;i<Trees[pid].size();i++){
                    m3+=Trees[pid][i].dis.size()*sizeof(int);//dis
                    m3+=Trees[pid][i].pos.size()*sizeof(int);//pos
                    m4+=Trees[pid][i].dis.size()*sizeof(int);//cnt
                    m4+=Trees[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
                }
            }
        }else if(algoUpdate==PCH_No){
            for(int pid=0;pid<Trees.size();++pid){
                for(int i=0;i<Trees[pid].size();i++){
//                    m3+=Trees[pid][i].dis.size()*sizeof(int);//dis
//                    m3+=Trees[pid][i].pos.size()*sizeof(int);//pos
//                    m4+=Trees[pid][i].dis.size()*sizeof(int);//cnt
                    m4+=Trees[pid][i].vert.size()*3*sizeof(int);//neighID/weight/count
                }
            }
        }
        //partitions
        for(int i=0;i< SCconNodesMTP.size();i++){
            for(auto it=SCconNodesMTP[i].begin(); it!=SCconNodesMTP[i].end(); it++){
                m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
            }
        }
    }



    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4+m5;
    cout<<"Distance labeling size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"Overlay graph index size: "<<(double)(m1+m2)/1024/1024<<" MB"<<endl;
    cout<<"Partition graphs index size "<<(double)(m3+m4)/1024/1024<<" MB"<<endl;
    cout<<"Extended label size "<<(double)(m5)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

//function for computing the index size
void Graph::IndexSizePostMHL(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0,m5=0;

    //Overlay index
    for(int i=0;i<Tree.size();i++){
        m1+=Tree[i].dis.size()*2*sizeof(int);//dis
        m1+=Tree[i].pos.size()*sizeof(int);//pos
        m2+=Tree[i].cnt.size()*sizeof(int);//cnt
        m2+=Tree[i].vert.size()*3*sizeof(int);//neighID/weight/count
        m3+=Tree[i].disInf.size()*2*sizeof(int);
        m3+=Tree[i].disPost.size()*2*sizeof(int);
    }
    //overlay
    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m2+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }


    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4+m5;
    cout<<"Distance labeling size: "<<(double)(m1+m3)/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::PH2HVertexOrdering(int type){
    ReadGraph(sourcePath+dataset);
    int pNum=partiNum;
    switch (type) {
        case 0:{//MDE partition + distant MDE overlay
            cout<<"MDE ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuild();
            PartitionOrderingBuildMDE(true);
//            PartitionOrderingBuildMDE(false);
            OrderingAssemblyMDE(pNum);
            break;
        }
        case 1:{
            cout<<"Boundary-first ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuildBoundaryFirst(node_num, NeighborSketch);
            PartitionOrderingBuildMDE(true);
            OrderingAssemblyBoundaryFirst(pNum);
            break;
        }
        case 2:{
            cout<<"Boundary-first MDE ordering."<<endl;
            SketchGraphBuild();
            OverlayOrderingBuild();
            PartitionOrderingBuildMDE(true);
            OrderingAssemblyMDEBoundaryFirst(pNum);
            break;
        }
        default:{
            cout<<"Wrong ordering type! "<<type<<endl; exit(1);
        }
    }
    exit(0);
}
void Graph::OrderingAssemblyMDEBoundaryFirst(int pNum){
    string filename=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE2";

    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int PID,ID;
    NodeOrder.assign(node_num,-1);
    vNodeOrder.clear();
    int order_i=0;
    set<int> vertices;
    /// For non-boundary vertex
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            if(!PartiTag[ID].second){// if not boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }
        }
    }
    /// For boundary vertex
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            if(PartiTag[ID].second){// if boundary
                vNodeOrder.emplace_back(ID);
                vertices.insert(ID);
                ++order_i;
            }
        }
    }

//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num || vertices.size()!=node_num || vNodeOrder.size()!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<vertices.size()<<" "<<vNodeOrder.size()<<" "<<node_num<<endl; exit(1);
    }



    for(int i=0;i<node_num;++i){
        NodeOrder[vNodeOrder[i]]=i;
    }

    for(int i=0;i<nodeNumOverlay;++i){
        ID=vNodeOrder[node_num-i-1];
        if(!PartiTag[ID].second){
            cout<<nodeNumOverlay<<", "<<node_num-i<<": "<<ID<<"("<<NodeOrder[ID]<<") is not boundary vertex!"<<endl; exit(1);
        }
    }

    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}
//function of MDE ordering assemblying
void Graph::OrderingAssemblyMDE(int pNum){
    string filename=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(pNum)+"/vertex_orderMDE";

    ofstream OF(filename,ios::out);
    if(!OF.is_open()){
        cout<<"Cannot open file "<<filename<<endl;
        exit(1);
    }
    int PID,ID;
    NodeOrder.assign(node_num,-1);
    int order_i=0;
    set<int> vertices;
    for(int pid=0;pid<partiNum;++pid){
        PID=vNodeOrderOverlay[pid];
//        cout<<pid<<": "<<PID<<endl;
        for(auto it=vNodeOrderParti[PID].begin();it!=vNodeOrderParti[PID].end();++it){
//            cout<<PID<<": "<<it->first<<" "<<it->second<<endl;
//            ID=it->second;
            ID=*it;
            vertices.insert(ID);
            NodeOrder[ID]=order_i;
            ++order_i;
        }
    }
//    cout<<verticesV.size()<<" "<<vertices.size()<<endl;
    if(order_i!=node_num){
        cout<<"Wrong order number! "<<order_i<<" "<<node_num<<endl; exit(1);
    }
    OF<<node_num<<endl;
    for(int i = 0; i < NodeOrder.size(); i++){
        if(NodeOrder[i]==-1){
            cout<<"Wrong order! "<<i<<"("<<PartiTag[i].first<<") "<<NodeOrder[i]<<endl; exit(1);
        }
        OF << i << " " << NodeOrder[i] << endl;
    }
    OF.close();
    cout<<"Finish "<<endl;
}
//function of boundary-first assemblying
void Graph::OrderingAssemblyBoundaryFirst(int pNum){
    string orderfile=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(pNum)+"/vertex_order2";
    set<int> vcheck;//use to check the redundant ordered vertex
    vcheck.clear();
    vNodeOrder.clear();
    int pid,ID;
    //for vertex within partition
    for(int k=0;k<vNodeOrderOverlay.size();k++){
        pid=vNodeOrderOverlay[k];

        for(auto it=vNodeOrderParti[k].begin();it!=vNodeOrderParti[k].end();++it){
//            ID=it->second;
            ID=*it;
            if(!PartiTag[ID].second){//if not boundary vertex
                vNodeOrder.push_back(ID);
                if(vcheck.find(ID)!=vcheck.end())
                    cout<<"wrong: redundant vertex ordered"<<endl;
                vcheck.insert(ID);
            }
        }
    }

    //for boundary vertex
    for(int k=0;k<vNodeOrderOverlay.size();k++){
        pid=vNodeOrderOverlay[k];
        for(int i=0;i<BoundVertex[pid].size();i++){
            ID=BoundVertex[pid][i];
            vNodeOrder.push_back(ID);
            if(vcheck.find(ID)!=vcheck.end())
                cout<<"wrong: redundant vertex ordered"<<endl;
            vcheck.insert(ID);
        }
    }

    //cout<<"total number of ordered vertex "<<vNodeOrder.size()<<endl;
    if(vNodeOrder.size()!=node_num)
        cout<<"Something wrong happened: some vertices do not have the vertex order!"<<endl;

    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }

    ofstream OF(orderfile);
    if(!OF){
        cout<<"Cannot open file "<<orderfile<<endl;
        exit(1);
    }
    OF<<NodeOrder.size()<<endl;
    for(int i=0;i<NodeOrder.size();i++){
        OF<<i<<" "<<NodeOrder[i]<<endl;
    }
    OF.close();
    cout<<"Finished."<<endl;
}

void Graph::SketchOrder(vector<vector<pair<int,int>>> Neighbor, vector<int> &vNodeOrderSketch){
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(1,1)));
    }

    _DD_.assign(node_num,0);
    DD.assign(node_num,0);

    set<DegComp> Deg;
    int degree;
    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();
        if(degree!=0){
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
        }
    }

    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect); //NeighborCon.clear();
    //SCconNodes.clear();

    //cout<<"Begin to contract"<<endl;
    int count=0;

    while(!Deg.empty()){
        //if(count%10==0) cout<<"count "<<count<<endl;
        count+=1;
        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }

        vNodeOrderSketch.push_back(x);
        Deg.erase(Deg.begin());
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();i++){
            for(int j=i+1;j<Neigh.size();j++){
                insertEOrderGenerate(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                change[Neigh[i].first]=true;
                change[Neigh[j].first]=true;
            }
        }
    }

    /*NodeOrderSketch.assign(nodenum,-1);
    for(int k=0;k<vNodeOrderSketch.size();k++){
        NodeOrderSketch[vNodeOrderSketch[k]]=k;
        cout<<"order "<<k<<", vertex "<<vNodeOrderSketch[k]<<", degree "<<NeighborSketch[vNodeOrderSketch[k]].size()<<endl;
    }*/
    //cout<<"Finish Contract"<<endl;
}

void Graph::SketchGraphBuild(){
    NeighborsParti.assign(node_num, vector<pair<vertex,int>>());
    NeighborsOverlay.assign(node_num,unordered_map<vertex,int>());
//    NeighborsOverlay.assign(node_num,vector<pair<vertex,int>>());
    PartiTag.assign(node_num, make_pair(-1,false));

    bool flag_minus = false;

    string filename=sourcePath+"partitions/"+dataset+"_"+algoParti+"_"+to_string(partiNum);
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

    PartiVertex.assign(partiNum,vector<vertex>());
    for(int k=0;k<pnum2;k++){
        int vernum,ID;
        IF1>>vernum;
        for(int i=0;i<vernum;i++){
            IF1>>ID;
//            if(flag_minus){
//                ID = ID-1;
//            }

            if(ID>=0 && ID<node_num){
                if(PartiTag[ID].first==-1){
                    PartiTag[ID].first=k;
                    PartiVertex[k].emplace_back(ID);
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
    int nNum=0;
    for(int pid=0;pid<partiNum;++pid){
        nNum+=PartiVertex[pid].size();
    }
    if(nNum!=node_num){
        cout<<"Inconsistent node number! "<<nNum<<" "<<node_num<<endl; exit(1);
    }
    //record the vertex to PartiVertex in vertex order: from lower-rank vertex to higher-rank vertex

    ifstream IF(filename+"/subgraph_edge");
    if(!IF){
        cout<<"Cannot open file "<<"subgraph_edge"<<endl;
        exit(1);
    }

    int pnum1;
    IF>>pnum1;
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
            }else{
                cout<<"Wrong for subgraph_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
            }


        }
    }

    vector<int> vecint;
    vecint.clear();
    NeighborSketch.assign(partiNum, vecint);
//    vector<set<int>> NeighborSketchS;
    NeighborSketchS.assign(partiNum, set<int>());

    BoundVertex.assign(partiNum,vector<vertex>());
    //read the cut edges
    ifstream IF2(filename+"/cut_edges");
    if(!IF2){
        cout<<"Cannot open file "<<"cut_edges"<<endl;
        exit(1);
    }

    int ednum,ID1,ID2,weight;
    int boundaryNum=0;
    int PID1, PID2;

    IF2>>ednum;
    for(int i=0;i<ednum;i++){
        IF2>>ID1>>ID2>>weight;
//        if(flag_minus){
//            ID1 = ID1-1; ID2 = ID2-1;
//        }

        if(ID1>=0 && ID1 <node_num && ID2>=0 && ID2 <node_num && weight>0){
            PID1=PartiTag[ID1].first, PID2=PartiTag[ID2].first;
            if(PartiTag[ID1].first==PartiTag[ID2].first){
                cout<<"two end points of cut edge are in the same partition"<<endl; exit(1);
            }
            if(!PartiTag[ID1].second){
                PartiTag[ID1].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID1].first].emplace_back(ID1);
            }
            if(!PartiTag[ID2].second){
                PartiTag[ID2].second=true;
                boundaryNum++;
                BoundVertex[PartiTag[ID2].first].emplace_back(ID2);
            }

            NeighborsOverlay[ID1].insert({ID2,weight});
//            NeighborsOverlay[ID1].emplace_back(ID2,weight);

            if(NeighborSketchS[PID1].find(PID2)==NeighborSketchS[PID1].end()){//if not found, i.e., PID2 is not in the NeighborSketchS of PID1
                NeighborSketch[PID1].push_back(PID2);
                NeighborSketch[PID2].push_back(PID1);

                NeighborSketchS[PID1].insert(PID2);
                NeighborSketchS[PID2].insert(PID1);
            }
        }else{
            cout<<"Wrong for cut_edge! "<<ID1<<" "<<ID2<<" "<<weight<<endl; exit(1);
        }
    }

    nodeNumOverlay=boundaryNum;


    /*for(int k=0;k<pnum;k++){
        cout<<k<<" "<<NeighborSketch[k].size()<<endl;
    }*/
}

void Graph::OverlayOrderingBuild(){

    int lastParti = -1;
    map<int,pair<int,int>> m;
    E.assign(partiNum,m);
    for(int i=0;i<NeighborSketch.size();i++){
        for(auto it=NeighborSketch[i].begin();it!=NeighborSketch[i].end();++it){
            E[i].insert(make_pair(*it,make_pair(1,1)));
        }
    }
    _DD_.assign(partiNum,0);
    DD.assign(partiNum,0);
    _DD2_.assign(partiNum,0);

    set<DegComp> Deg;
    set<DegComp2> Deg2;
    int ID,degree;
    for(ID=0;ID<partiNum;ID++){
        degree=NeighborSketch[ID].size();
        if(degree!=0){
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
        }else{
            cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
        }
    }

    vector<bool> exist; exist.assign(partiNum,true);
    vector<bool> change; change.assign(partiNum,false);

    vector<set<int>> NeighborCon(partiNum,set<int>());
    vector<int> neix;
    neix.clear();

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    int x;

    set<int> pSet;
    while(!Deg.empty()){

        x=(*Deg.begin()).x;
        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }
        Deg.erase(Deg.begin());



        if(lastParti!=-1 && NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if x is lastParti's neighbor
            neix.clear();
            while(NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if found
                _DD2_[x]++;
                neix.emplace_back(x);

                if(Deg.empty()){
                    break;
                }else{
                    x=Deg.begin()->x;
                    Deg.erase(Deg.begin());
                }
            }

            if(NeighborCon[lastParti].find(x) != NeighborCon[lastParti].end()){//if found, i.e., x is the neighbor of lastParti
                if(neix.size()>1){//if neix has more than one element
                    if(Deg.empty()){
                        Deg2.clear();
                        for(int i=0;i<neix.size();++i){
                            Deg2.insert(DegComp2(neix[i]));
                        }
                        x=Deg2.begin()->x;
                        Deg2.erase(Deg2.begin());
                        if(!Deg2.empty()){
                            for(auto it=Deg2.begin();it!=Deg2.end();++it){
                                ID=it->x;
                                Deg.insert(DegComp(ID));
                            }
                            Deg2.clear();
                        }
                    }else{
                        cout<<"Wrong! "<<endl; exit(1);
                    }
                }
            }//if not the neighbor
            else{
                if(!neix.empty()){
                    for(int i=0;i<neix.size();++i){
                        Deg.insert(DegComp(neix[i]));
                    }
                }
            }
        }

//        cout<<x<<" "<<Deg.size()<<endl;
        vNodeOrderOverlay.emplace_back(x);
        pSet.insert(x);
        lastParti = x;
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.emplace_back(*it);
                NeighborCon[x].insert(it->first);
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderOverlay.size() != partiNum || pSet.size()!=partiNum){
        cout<<"Inconsistent size for sketch graph! "<<vNodeOrderOverlay.size()<<" "<<pSet.size() <<" "<< partiNum<<endl; exit(1);
    }

//    exit(0);
}

//based on SketchOrder function, to guarantee the order of neighboring vertices are not contiguous
void Graph::OverlayOrderingBuildBoundaryFirst(int nodenum, vector<vector<int>> Neighbor){
    map<int,pair<int,int>> m;
    E.assign(partiNum,m);
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j],make_pair(1,1)));
    }

    _DD_.assign(partiNum,0);
    DD.assign(partiNum,0);

    set<DegComp> Deg;
    int degree;
    for(int i=0;i<Neighbor.size();i++){
        degree=Neighbor[i].size();
        if(degree!=0){
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
        }
    }

    vector<bool> exist; exist.assign(partiNum,true);
    vector<bool> change; change.assign(partiNum,false);

    vector<vector<pair<int,pair<int,int>>>> NeighborCon;
    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(partiNum,vect); //NeighborCon.clear();
    //SCconNodes.clear();

    //cout<<"Begin to contract"<<endl;
    int count=0;
    int x;
    int lastx;
    vector<int> neix;
    neix.clear();
    while(!Deg.empty()){
        //if(count%10==0) cout<<"count "<<count<<endl;
        count+=1;

        while(true){
            if(Deg.empty()){//in case that all the remaining vertices are the neighbors of current vertex x
                x=neix[0];
                for(int j=1;j<neix.size();j++){
                    Deg.insert(DegComp(neix[j]));
//					cout<<"insert back/// "<<neix[j]<<endl;
                }
                neix.clear();
                break;///
            }
            else
                x=(*Deg.begin()).x;

            while(true){
                if(change[x]){
                    Deg.erase(DegComp(x));
                    _DD_[x]=DD[x];
                    Deg.insert(DegComp(x));
                    change[x]=false;
                    x=(*Deg.begin()).x;
                }else
                    break;
            }

            if(count==1){
                lastx=x;
                break;
            }else if(NeighborSketchS[lastx].find(x)==NeighborSketchS[lastx].end()){//if not found
                lastx=x;
                break;
            }else{
                Deg.erase(DegComp(x));
                neix.push_back(x);
//				cout<<"erase "<<x<<endl;
            }

//            if(Deg.empty())////
//                break;
        }

        if(neix.size()!=0){
            for(int j=0;j<neix.size();j++){
                Deg.insert(DegComp(neix[j]));
//				cout<<"insert back "<<neix[j]<<endl;
            }
        }
        neix.clear();

        vNodeOrderOverlay.push_back(x);
        Deg.erase(DegComp(x));
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.push_back(*it);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();i++){
            for(int j=i+1;j<Neigh.size();j++){
                insertEOrderGenerate(Neigh[i].first,Neigh[j].first,Neigh[i].second.first+Neigh[j].second.first);
                change[Neigh[i].first]=true;
                change[Neigh[j].first]=true;
            }
        }
    }

    /*NodeOrderSketch.assign(nodenum,-1);
    for(int k=0;k<vNodeOrderSketch.size();k++){
        NodeOrderSketch[vNodeOrderSketch[k]]=k;
        cout<<"order "<<k<<", vertex "<<vNodeOrderSketch[k]<<", degree "<<NeighborSketch[vNodeOrderSketch[k]].size()<<endl;
    }*/
    //cout<<"Finish Contract"<<endl;
}
//MDE-based vertex ordering for partition
void Graph::PartitionOrderingBuildMDE(bool ifParallel){
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsParti.size();i++){
        for(auto it=NeighborsParti[i].begin();it!=NeighborsParti[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    _DD_.assign(node_num,0);
    DD.assign(node_num,0);
//    vNodeOrderParti.assign(partiNum,map<int,int>());
    vNodeOrderParti.assign(partiNum,vector<int>());

    if(ifParallel){
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrderingV, this, boost::ref(processID[j]) ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::PartitionOrdering, this, j));
            }
            thread.join_all();
        }
    }
    else{
        for(int i=0;i<partiNum;++i){
            PartitionOrdering(i);
        }
    }

}

void Graph::PartitionOrderingV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        PartitionOrdering(p[i]);
    }
}

void Graph::PartitionOrdering(int pid){

    set<DegComp> Deg;
    int ID,degree;
    for(int i=0;i<PartiVertex[pid].size();i++){
        ID = PartiVertex[pid][i];
        degree=NeighborsParti[ID].size();
        if(degree!=0){
            _DD_[ID]=degree;
            DD[ID]=degree;
            Deg.insert(DegComp(ID));
        }else{
            cout<<"Wrong! degree is zero. "<<ID<<" "<<degree<<endl; exit(1);
        }
    }

    vector<bool> exist; exist.assign(node_num,true);
    vector<bool> change; change.assign(node_num,false);

    int count=0;
    int Twidth=0;
    int order_i=0;
    int ID1, ID2;
    while(!Deg.empty()){
//        if(count%10000==0)
//            cout<<"count "<<count<<" , treewidth "<<Twidth<<endl;
        count+=1;
        int x=(*Deg.begin()).x;

        while(true){
            if(change[x]){
                Deg.erase(DegComp(x));
                _DD_[x]=DD[x];
                Deg.insert(DegComp(x));
                change[x]=false;
                x=(*Deg.begin()).x;
            }else
                break;
        }

//        vNodeOrderParti[pid].insert({order_i,x});
        vNodeOrderParti[pid].push_back(x);
        order_i++;
        Deg.erase(Deg.begin());
        exist[x]=false;

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();

        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(exist[(*it).first]){
                Neigh.emplace_back(*it);
            }
        }

        if(Neigh.size()>Twidth)
            Twidth=Neigh.size();

        //multi threads for n^2 combination
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteEOrderGenerate(x,y);//delete x from y's neighbor
            change[y]=true;
        }

        for(int i=0;i<Neigh.size();++i){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();++j){
                ID2=Neigh[j].first;
                insertEOrderGenerate(ID1, ID2, 1);//for order generation, only need unweighted graph
            }
        }
    }
    if(vNodeOrderParti[pid].size() != PartiVertex[pid].size()){
        cout<<"Inconsistent size! "<< pid <<" "<<vNodeOrderParti[pid].size() <<" "<< PartiVertex[pid].size()<<endl; exit(1);
    }
}



/// Functions for MDE contraction
void Graph::deleteEOrderGenerate(int u,int v){//delete u from v's neighbor
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}

void Graph::insertEOrderGenerate(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//        DD2[u]++;
    }
    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        DD[v]++;
//        DD2[u]++;
    }
}


/// Query Processing
//function for correctness check
void Graph::CorrectnessCheck(int runtimes){
    Timer tt;
    double runT=0;
    srand (time(NULL));
    int s, t, d1, d2, d3;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ... AlgoQuery: "<<algoQuery<<".";
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
//        int pid=rand()%partiNum;
//        s=PartiVertex[pid][rand()%PartiVertex[pid].size()];
//        t=PartiVertex[pid][rand()%PartiVertex[pid].size()];
//        s=92291,t=22007;//PMHL
//        s=192997,t=163009;//PCH
//        s=116238,t=116089;//PostMHL
//        s=115541,t=116089;//PostMHL
//        s=115541,t=116033;//PostMHL//115999
//        s=174536,t=179173;//PostMHL
//        s=82961, t=91878;//PostMHL
//        s=67975, t=67432;//PostMHL
//        s=13076, t=13104;//PostMHL//13104, 23022
//        cout<<"Query "<<i<<": "<<s<<" "<<t<<endl;

//        if(runtimes == 1){
//            cout<<"s: "<<s<<" ; t: "<<t<<endl;
//        }

        if(algoChoice==1){
            tt.start();
            d2=QueryNP(s,t);
            tt.stop();
        }else if(algoChoice==2){
            tt.start();
            d2=Query(s,t);
            tt.stop();
        }else if(algoChoice==3){
            tt.start();
            d2=QueryPMHL(s,t);
            tt.stop();
        }
        else if(algoChoice==4){
            tt.start();
            d2=QueryPMHLOpt(s,t);
            tt.stop();
        }
        else if(algoChoice==5){
            tt.start();
            d2=QueryPostMHL(s,t);
            tt.stop();
        }else if(algoChoice==0){
            tt.start();
            d2= Astar(s,t,Neighbor);
            tt.stop();
        }
        runT+=tt.GetRuntime();

        d1=Dijkstra(s,t,Neighbor);
//        cout<<"Algorithm "<<algoQuery<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTag[s].first<<" "<<PartiTag[t].first<<"; Boundary Tag: "<<PartiTag[s].second<<" "<<PartiTag[t].second<<endl;
        if(d1!=d2){

            if(algoChoice==1){
                cout<<"Correct Test. InCorrect! Algorithm "<<algoQuery<<", "<<i<<": "<<s<<" "<<t<<": "<<d2<<" "<<d1<<endl;
            }else if(algoChoice==2 || algoChoice==3 || algoChoice==4){
                cout<<"Correct Test. InCorrect! Algorithm "<<algoQuery<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTag[s].first<<" "<<PartiTag[t].first<<"; Boundary Tag: "<<PartiTag[s].second<<" "<<PartiTag[t].second<<endl;
                QueryDebug(s,t);
            }else if(algoChoice==5){
                cout<<"Correct Test. InCorrect! Algorithm "<<algoQuery<<", "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<" ; Partition Tag: "<< PartiTags[s].first<<" "<<PartiTags[t].first<<endl;
                QueryPostMHLDebug(s,t);
            }else if(algoChoice==0){
                cout<<"InCorrect! "<<i<<": "<<s<<" "<<t<<": "<<d2<<" "<<d1<<endl;
            }
            exit(1);
        }
    }
    cout<<" Average Query Time: "<<1000*runT/runtimes<<" ms."<<endl;
}

//function for Query processing, debug version
int Graph::QueryDebug(int ID1, int ID2){
    int dis=INF;

    if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in core
        cout<<"Core-Core"<<endl;
//        dis=QueryCoreDebug(ID1, ID2);

    }else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
        cout<<"Core-Parti"<<endl;
        dis=QueryPartiCoreDebug(ID2, ID1);

    }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
        cout<<"Parti-Core"<<endl;
        dis=QueryPartiCoreDebug(ID1, ID2);

    }else if(!PartiTag[ID1].second && !PartiTag[ID2].second){//both in partition

        if(PartiTag[ID1].first != PartiTag[ID2].first){//Case 3: in different peripheries
            cout<<"Parti-Parti"<<endl;
            int d=INF;
            int b1,b2,d1,d2;//final results
            int pid1=PartiTag[ID1].first;
            int pid2=PartiTag[ID2].first;

            vector<int> B1=BoundVertex[pid1];
            vector<int> B2=BoundVertex[pid2];

            map<int,int> m1,m2;
            m1.clear();
            m2.clear();
            int bID1, bID2, tempdis;
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];

//            m1.insert(make_pair(bID1,Tree[rank[ID1]].disInf[i]));
                m1.insert(make_pair(bID1, QueryH2HPartition(ID1,bID1,pid1)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
//            m2.insert(make_pair(bID2,Tree[rank[ID2]].disInf[j]));
                m2.insert(make_pair(bID2, QueryH2HPartition(ID2,bID2,pid2)));
            }

            for(int k=0;k<B1.size();k++){
                bID1=B1[k];

                if(m1[bID1]>d)
                    continue;

                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];

                    if(m2[bID2]>d)
                        continue;

                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d){
                        d=tempdis;
                        b1=bID1; b2=bID2; d1=m1[bID1]; d2=m2[bID2];
//                        cout<<b1<<" "<<b2<<" "<<d<<endl;
                    }
                }
            }
            dis=d;
            int d_12=QueryCore(b1,b2), dDijk_s=Dijkstra(ID1,b1,Neighbor), dDijk_12=Dijkstra(b1,b2,Neighbor), dDijk_t=Dijkstra(b2,ID2,Neighbor);
            cout<<ID1<<" "<<b1<<"("<<NodeOrder[b1]<<","<<PartiTag[b1].first<<") "<<b2<<"("<<NodeOrder[b2]<<","<<PartiTag[b2].first<<") "<<ID2<<" : "<<d1<<" "<<d_12<<" "<<d2<<" ; "<<dDijk_s<<" "<<dDijk_12<<"("<<DijkstraCore(b1,b2)<<") "<<dDijk_t<<endl;

            QueryCoreDebug(b1,b2);
            QueryPartiPartiExtLCADebug(ID1,ID2);

//                if(d1!=dDijk_s){
//                    DijkstraPath(ID1,b1);
//                }
//                if(d_12!=dDijk_12){
//                    DijkstraPath(b1,b2);
//                }
//                if(d2!=dDijk_t){
//                    DijkstraPath(b2,ID2);
//                }

        }
        else{//Case 4: in the same periphery
            cout<<"Same-Parti"<<endl;
//                dis= QuerySameParti(ID1,ID2);
            int d=INF;
            int b1,b2,df1,df2;
            int pid1=PartiTag[ID1].first;
            int pid2=PartiTag[ID2].first;

            int temp_dis = QueryH2HPartition(ID1, ID2, pid1);/// d2 may be wrong sometimes
            if (temp_dis < d){
                d = temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
                b1=b2=-1;
                df1=df2=-1;
            }

            vector<int> B = BoundVertex[pid1];
            map<int, int> m1, m2;
            m1.clear();
            m2.clear();
            vector<int> B1, B2;
            B1.clear();
            B2.clear();
            int bID, d1, d2;
            for (int i = 0; i < B.size(); i++) {
                bID = B[i];

                d1 = QueryH2HPartition(ID1,bID,pid1);
                d2 = QueryH2HPartition(ID2,bID,pid1);

                if (d1 < d) {
                    B1.push_back(bID);
                    m1.insert(make_pair(bID, d1));
                }
                if (d2 < d) {
                    B2.push_back(bID);
                    m2.insert(make_pair(bID, d2));
                }
            }

            int bID1, bID2, tempdis;
            if (!B1.empty() && !B2.empty()) {
                for (int k = 0; k < B1.size(); k++) {
                    bID1 = B1[k];
                    if (m1[bID1] > d)
                        continue;
                    for (int z = 0; z < B2.size(); z++) {
                        bID2 = B2[z];
                        if (m2[bID2] > d)
                            continue;
                        tempdis = m1[bID1] + QueryCore(bID1, bID2) + m2[bID2];
                        if (tempdis < d){
                            d = tempdis;
                            b1=bID1;b2=bID2;
                            df1=m1[bID1];df2=m2[bID2];
                        }
                    }
                }
            }

            if(b1!=-1){
                cout<<"d4: "<<ID1<<" "<<b1<<" "<<b2<<" "<<ID2<<" : "<<df1<<" "<<QueryCore(b1,b2)<<" "<<df2<<" ; "<<Dijkstra(ID1,b1,Neighbor)<<" "<<Dijkstra(b1,b2,Neighbor)<<" "<<Dijkstra(b2,ID2,Neighbor)<<endl;
            }else{
                int dDijk2 = Dijkstra(ID1,ID2,Neighbor);
                cout<<"d2: "<<d<<"; "<<dDijk2<<endl;
                if(d!=dDijk2){
//                        DijkstraPath(ID1,ID2);
                }
            }

            dis = d;

        }

    }
    return dis;
}

//function for core index correctness check
void Graph::CorrectnessCheckCore(int runtimes){
    srand (time(NULL));
    int s, t, d1, d2, d3;
    vector<int> coreVertex;
    for(int i=0;i<node_num;++i){
        if(PartiTag[i].second){
            coreVertex.emplace_back(i);
        }
    }
    int corenum=coreVertex.size();
    cout<<"Core graph correctness check ("<<runtimes<<" rounds)..."<<endl;
    for(int i=0;i<runtimes;i++){
        s=coreVertex[rand()%corenum];
        t=coreVertex[rand()%corenum];
        if(PartiTag[s].second && PartiTag[t].second){//for core vertex
            d1=QueryCore(s,t);
            d2=DijkstraCore(s,t);

            if(d1!=d2){
                cout<<"InCorrect! "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<"): "<<d1<<" "<<d2<<endl;
//				DijkstraPath(s,t);
//				DijkstraCorePath(s,t);
                exit(1);
            }
        }else
            i--;
    }
}

//function for efficiency test
void Graph::EffiCheck(string filename,int runtimes){
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Query file: "<<filename<<endl;
    int num, ID1, ID2;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.push_back(make_pair(ID1, ID2));
    }
    if(runtimes > num){
        runtimes = num;
    }
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;

    double runT=0;
    int d1,d2;
    bool ifDebug=false;
//    ifDebug=true;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }
    clock_t start = clock();
    vector<int> results(runtimes,-1);
    for(int i=0;i<runtimes;i++){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
//        if(PartiTag[ID1].first!=PartiTag[ID2].first){
//            cout<<"Different Partition: "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<endl;
//        }
        if(algoChoice==1){
            tt.start();
            d1=QueryNP(ID1,ID2);
            tt.stop();
        }else if(algoChoice==2){
            tt.start();
            d1=Query(ID1,ID2);
            tt.stop();
        }else if(algoChoice==3){
            tt.start();
            d1=QueryPMHL(ID1,ID2);
            tt.stop();
        }else if(algoChoice==4){
            tt.start();
            d1=QueryPMHLOpt(ID1,ID2);
            tt.stop();
        }else if(algoChoice==5){
            tt.start();
            d1=QueryPostMHL(ID1,ID2);
            tt.stop();
        }else if(algoChoice==0){
            tt.start();
            d1=Astar(ID1,ID2,Neighbor);
//            d1=Dijkstra(ID1,ID2,Neighbor);
            tt.stop();
        }

        runT+=tt.GetRuntime();
        results[i]=d1;
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                if(algoChoice==5){
                    cout<<"Incorrect! "<<i<<": "<<ID1<<"("<<PartiTags[ID1].first<<") "<<ID2<<"("<<PartiTags[ID2].first<<"): "<<d1<<" "<<d2<<endl;
                    QueryPostMHLDebug(ID1,ID2);
                    exit(1);
                }else if(algoChoice>=1 && algoChoice<=4){
                    cout<<"Incorrect! "<<i<<": "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<"): "<<d1<<" "<<d2<<endl; exit(1);
                }

            }
        }

    }


    cout<<"Average Query Time: "<<(double)runT*1000/runtimes<<" ms. "<<1000*(double)(clock() - start) / (CLOCKS_PER_SEC*runtimes)<<" ms."<<endl;
}

//function for efficiency test
void Graph::EffiCheckStages(vector<pair<int,int>> & ODpair, int runtimes, int intervalT, unsigned long long & throughputNum, vector<double>& stageUpdateT, vector<double>& stageQueryT){
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;

    double runT=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    double updateT=0;//overall update time
    unsigned long long qt1=0,qt2=0,qt3=0,qt4=0,qt5=0,qtAll=0;
    int d1,d2;
    bool ifDebug=false;
//    ifDebug=true;

    if(ifDebug){
        cout<<"With correctness check."<<endl;
    }
    clock_t start = clock();
    vector<int> results(runtimes,-1);

    if(algoChoice==1){
        runT1=EffiMHLStage(ODpair,1000,Dijk);
        runT2=EffiMHLStage(ODpair,runtimes/10,CH);
        runT3=EffiMHLStage(ODpair,runtimes,H2H);
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT){
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: CH
            updateT+=stageDurations[CH];
            if(updateT<intervalT){
                dt2=stageDurations[CH]; qt2=dt2/runT2;
                //Stage 3: CH
                stageDurations[H2H]=intervalT-updateT;
                dt3=stageDurations[H2H]; qt3=dt3/runT3;
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }else{
            dt1=intervalT; qt1=dt1/runT1;
        }


        qtAll=qt1+qt2+qt3;
        cout<<"Stage 1 (BiDijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (CH): Duration: "<<dt2<<" ("<<stageDurations[CH]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (H2H): Duration: "<<dt3<<" ("<<stageDurations[H2H]<<") s; Throughput number: "<<qt3<<" ; query time: "<<1000 * runT3 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[CH]+=stageDurations[CH];
        stageUpdateT[H2H]+=stageDurations[H2H];
        stageQueryT[Dijk]+=runT1;
        stageQueryT[CH]+=runT2;
        stageQueryT[H2H]+=runT3;
    }
    else if(algoChoice==3){
        runT1=EffiPMHLStage(ODpair,1000,Dijk);
        runT2=EffiPMHLStage(ODpair,runtimes/10,PCH_No);
        runT3=EffiPMHLStage(ODpair,runtimes,PH2H_No);
        runT4=EffiPMHLStage(ODpair,runtimes,PH2H_Post);
        runT5=EffiPMHLStage(ODpair,runtimes,PH2H_Cross);
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT){
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT){
                dt2=stageDurations[PCH_No]; qt2=dt2/runT2;
                //Stage 3: No-boundary PMHL
                updateT+=stageDurations[PH2H_No];
                if(updateT<intervalT){
                    dt3=stageDurations[PH2H_No]; qt3=dt3/runT3;
                    //Stage 4: Post-boundary PMHL
                    updateT+=stageDurations[PH2H_Post];
                    if(updateT<intervalT){
                        dt4=stageDurations[PH2H_Post]; qt4=dt4/runT4;
                        //Stage 5: Cross-boundary PMHL
                        stageDurations[PH2H_Cross]=intervalT-updateT;
                        dt5=stageDurations[PH2H_Cross]; qt5=dt5/runT5;
                    }else{
                        dt4=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]-stageDurations[PH2H_No]; qt4=dt4/runT4;
                    }
                }else{
                    dt3=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]; qt3=dt3/runT3;
                }
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }
        else{
            dt1=intervalT; qt1=dt1/runT1;
        }

        qtAll=qt1+qt2+qt3+qt4+qt5;
        cout<<"Stage 1 (Dijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (PCH): Duration: "<<dt2<<" ("<<stageDurations[PCH_No]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-No): Duration: "<<dt3<<" ("<<stageDurations[PH2H_No]<<") s; Throughput number: "<<qt3<<" ; query time: "<<1000 * runT3 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Post): Duration: "<<dt4<<" ("<<stageDurations[PH2H_Post]<<") s; Throughput number: "<<qt4<<" ; query time: "<<1000 * runT4 << " ms."<<endl;
        cout<<"Stage 5 (PH2H-Extend): Duration: "<<dt5<<" ("<<stageDurations[PH2H_Cross]<<") s; Throughput number: "<<qt5<<" ; query time: "<<1000 * runT5 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
        stageUpdateT[PH2H_No]+=stageDurations[PH2H_No];
        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
        stageQueryT[Dijk]+=runT1;
        stageQueryT[PCH_No]+=runT2;
        stageQueryT[PH2H_No]+=runT3;
        stageQueryT[PH2H_Post]+=runT4;
        stageQueryT[PH2H_Cross]+=runT5;
    }
    else if(algoChoice==5){
        runT1 = EffiPostMHLStage(ODpair,1000,Dijk);
        runT2 = EffiPostMHLStage(ODpair, runtimes / 10, PCH_No);
        runT4 = EffiPostMHLStage(ODpair, runtimes, PH2H_Post);
        runT5 = EffiPostMHLStage(ODpair, runtimes, PH2H_Cross);
        //Stage 1: Dijkstra
        updateT=stageDurations[Dijk];
        if(updateT<intervalT) {
            dt1=stageDurations[Dijk]; qt1=dt1/runT1;
            //Stage 2: PCH
            updateT+=stageDurations[PCH_No];
            if(updateT<intervalT) {
                dt2=stageDurations[PCH_No]; qt2=dt2/runT2;
                //Stage 4: Post-boundary PMHL
                updateT+=stageDurations[PH2H_Post];
                if(updateT<intervalT) {
                    dt4=stageDurations[PH2H_Post]; qt4=dt4/runT4;
                    //Stage 5: Cross-boundary PMHL
                    stageDurations[PH2H_Cross]=intervalT-updateT;
                    dt5=stageDurations[PH2H_Cross]; qt5=dt5/runT5;
                }else{
                    dt4=intervalT-stageDurations[Dijk]-stageDurations[PCH_No]; qt4=dt4/runT4;
                }
            }else{
                dt2=intervalT-stageDurations[Dijk]; qt2=dt2/runT2;
            }
        }else{
            dt1=intervalT; qt1=dt1/runT1;
        }

        qtAll=qt1+qt2+qt4+qt5;
        cout<<"Stage 1 (Dijkstra): Duration: "<<dt1<<" ("<<stageDurations[Dijk]<<") s; Throughput number: "<<qt1<<" ; query time: "<<1000 * runT1 << " ms."<<endl;
        cout<<"Stage 2 (PCH): Duration: "<<dt2<<" ("<<stageDurations[PCH_No]<<") s; Throughput number: "<<qt2<<" ; query time: "<<1000 * runT2 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-Post): Duration: "<<dt4<<" ("<<stageDurations[PH2H_Post]<<") s; Throughput number: "<<qt4<<" ; query time: "<<1000 * runT4 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Extend): Duration: "<<dt5<<" ("<<stageDurations[PH2H_Cross]<<") s; Throughput number: "<<qt5<<" ; query time: "<<1000 * runT5 << " ms."<<endl;
        cout<<"Throughput number: "<<qtAll<<endl;
        throughputNum = throughputNum + qtAll;
        stageUpdateT[Dijk]+=stageDurations[Dijk];
        stageUpdateT[PCH_No]+=stageDurations[PCH_No];
        stageUpdateT[PH2H_Post]+=stageDurations[PH2H_Post];
        stageUpdateT[PH2H_Cross]+=stageDurations[PH2H_Cross];
        stageQueryT[Dijk]+=runT1;
        stageQueryT[PCH_No]+=runT2;
        stageQueryT[PH2H_Post]+=runT4;
        stageQueryT[PH2H_Cross]+=runT5;
    }else {
        cout<<"Wrong query type! "<<algoChoice<<endl; exit(1);
    }

}

double Graph::EffiMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType) {
    algoQuery=queryType;
    Timer tt;
    int ID1, ID2, d1, d2;
    double runT=0;
    vector<int> results(runtimes,0);
    bool ifDebug=false;
    ifDebug=true;

    for(int i=0;i<runtimes;++i){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
        tt.start();
        d1=QueryNP(ID1,ID2);
        results[i]=d1;
        tt.stop();
        runT+=tt.GetRuntime();
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                cout<<"Wrong result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
            }
        }
    }


    cout<<"Efficiency Test of Stage "<<queryType<<" : "<<runT<<" s."<<endl;
    return runT/runtimes;
}

double Graph::EffiPMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType) {
    algoQuery=queryType;
    Timer tt;
    int ID1, ID2, d1, d2;
    double runT=0;
    vector<int> results(runtimes,0);
    bool ifDebug=false;
//    ifDebug=true;

    for(int i=0;i<runtimes;++i){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
        tt.start();
        d1=QueryPMHL(ID1,ID2);
        results[i]=d1;
        tt.stop();
        runT+=tt.GetRuntime();
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                cout<<"Wrong result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
            }
        }
    }

    cout<<"Efficiency Test of Stage "<<queryType<<" : "<<runT<<" s."<<endl;
    return runT/runtimes;
}

double Graph::EffiPostMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType) {
    algoQuery=queryType;
    Timer tt;
    int ID1, ID2, d1, d2;
    double runT=0;
    vector<int> results(runtimes,0);
    bool ifDebug=false;
//    ifDebug=true;

    for(int i=0;i<runtimes;++i){
        ID1=ODpair[i].first, ID2=ODpair[i].second;
        tt.start();
        d1=QueryPostMHL(ID1,ID2);
        results[i]=d1;
        tt.stop();
        runT+=tt.GetRuntime();
        if(ifDebug){
            d2= Dijkstra(ID1,ID2,Neighbor);
            if(d1!=d2){
                cout<<"Wrong result! "<<ID1<<" "<<ID2<<" "<<d1<<" "<<d2<<endl; exit(1);
            }
        }
    }

    cout<<"Efficiency Test of Stage "<<queryType<<" : "<<runT<<" s."<<endl;
    return runT/runtimes;
}

void Graph::DFSTree(vector<int>& tNodes, int id){
    tNodes.push_back(Tree[id].uniqueVertex);
    for(int i=0;i<Tree[id].ch.size();++i){
        DFSTree(tNodes,Tree[id].ch[i]);
    }
}





/// Index Maintenance


//function of testing the throughput of path-finding system, batchInterval is the time interval between two adjacent update batch (in seconds)
void Graph::SPThroughputTest(int updateType, bool ifBatch, int batchNum, int batchSize, int batchInterval, int runtimes) {
    cout<<"Shortest path query throughput test..."<<endl;
    // read updates
    string file = sourcePath+dataset + ".update";
    bool ifDebug=false;
//    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand(0);
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);

    if(batchNum*batchSize>updateData.size()){
        batchNum=floor(updateData.size()/batchSize);
    }
    cout<<"Update batch: "<<batchNum<<" ; Batch size: "<<batchSize<<" ; Batch interval: "<< batchInterval<<endl;

    string queryF = sourcePath+dataset + ".query";
    if(samePartiPortion!=-1){
        queryF=sourcePath + "tmp/"+dataset+".PostMHL_"+to_string(bandWidth)+".sameParti_"+to_string(samePartiPortion)+".query";
    }
    ifstream IF(queryF);
    if(!IF){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF.close();
    unsigned long long throughputNum=0;
    if(algoChoice==5){
        ProBeginVertexSetParti.assign(partiNum,vector<int>());
        vertexIDChLParti.assign(partiNum,set<int>());
        ProBeginVertexSetPartiExtend.assign(partiNum,vector<int>());
    }

    Timer tt;
    Timer tRecord;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 0:{
            break;
        }
        case 1:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;

            if(ifBatch){//for batch update
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
//                Graph g2=*this;
                vector<double> stageUpdateT(5,0);
                vector<double> stageQueryT(5,0);

                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    stageDurations.assign(5,0);

                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[u*batchSize+i].first.first;
                        ID2 = updateData[u*batchSize+i].first.second;
                        oldW = updateData[u*batchSize+i].second;
                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
                        if(ifDebug){
                            //cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

//                    int s,t;
//                    s=112485, t=148515;//core-core
//                    s=41522, t=47616;//core-parti
//                    cout<<"Before update. "<<QueryPostMHLDebug(s,t)<<endl;

                    tRecord.start();
//                    boost::thread_group thread;

                    if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::DecBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT1)));
                        DecBatchThroughputNP(wBatch, u, runT1);
                    }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                        PMHLBatchUpdateDec(wBatch, u, runT1);
                    }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateDec, this, boost::ref(wBatch), u, boost::ref(runT1)));
                        PostMHLBatchUpdateDec(wBatch, u, runT1);
                    }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

                    EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);

                    if(ifDebug){
//                        g2.CorrectnessCheckH2H(100);
                        CorrectnessCheck(100);
                    }
                }
                StagePerformanceShow(batchNum,stageUpdateT,stageQueryT);
                cout<<"\nPartiNum: "<<partiNum<<". Overall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average Decrease batch update Time: "<<runT1/batchNum<<" s."<<endl;
            }
            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            if(ifBatch){//for batch update
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                vector<double> stageUpdateT(5,0);
                vector<double> stageQueryT(5,0);

//                for(int u=5;u<6;++u){
//                for(int u=117;u<118;++u){
                for(int u=0;u<batchNum;++u){
                    wBatch.clear();
                    stageDurations.assign(5,0);
                    for(int i=0;i<batchSize;++i) {
                        ID1 = updateData[u * batchSize + i].first.first;
                        ID2 = updateData[u * batchSize + i].first.second;
                        oldW = updateData[u * batchSize + i].second;
                        newW = oldW * 1.5;
                        if (ifDebug) {
                            //cout << "Batch " << u << ": " << ID1 << " " << ID2 << " " << oldW << " " << newW << endl;
                        }
                        wBatch.emplace_back(make_pair(ID1, ID2), make_pair(oldW, newW));
                    }
//                    int s,t,hub;
//                    s=67432, t=67975;//same-parti, 13076
////                    s=115541, t=116033;//same-parti//115999
//                    cout<<"Before update. "<<QueryPostMHLDebug(s,t)<<endl;
//                    hub=s;
//                    for(int i=0;i<Tree[rank[t]].disPost.size();++i){
//                        if(Tree[rank[t]].vAncestor[i]==hub){
//                            cout<<"Cnt: "<<Tree[rank[t]].disPost[i]<<" "<<Tree[rank[t]].cntPost[i]<<endl;
//                        }
//                    }
//                    for(auto it=Tree[rank[t]].vert.begin();it!=Tree[rank[t]].vert.end();++it){
//                        cout<<"sc: "<<t<<" "<<it->first<<" "<<it->second.first<<"("<<it->second.second<<")"<<endl;
//                    }
//                    int bv=11922;
//                    cout<<"disB: "<<s<<" "<<bv<<" "<<Tree[rank[s]].disInf[BoundVertexMap[3][bv]]<<endl;

                    tRecord.start();
                    boost::thread_group thread;
                    if(algoChoice==1){
//                        thread.add_thread(new boost::thread(&Graph::IncBatchThroughputNP, this, boost::ref(wBatch), u, boost::ref(runT2)));
                        IncBatchThroughputNP(wBatch, u, runT2);
                    }else if(algoChoice==3){
//                        thread.add_thread(new boost::thread(&Graph::PMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));
                        PMHLBatchUpdateInc(wBatch, u, runT2);
                    }else if(algoChoice==5){
//                        thread.add_thread(new boost::thread(&Graph::PostMHLBatchUpdateInc, this, boost::ref(wBatch), u, boost::ref(runT2)));//increase update is incorrect
                        PostMHLBatchUpdateInc(wBatch, u, runT2);
                    }

//                    thread.add_thread(new boost::thread(&Graph::EffiCheckThroughput, this, boost::ref(ODpair),  boost::ref(tRecord), batchInterval, boost::ref(throughputNum)));
//                    thread.join_all();

                    EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);

//                    cout<<"After update: "<<endl;
//                    for(int i=0;i<Tree[rank[t]].disPost.size();++i){
//                        if(Tree[rank[t]].vAncestor[i]==hub){
//                            cout<<"Cnt: "<<Tree[rank[t]].disPost[i]<<" "<<Tree[rank[t]].cntPost[i]<<endl;
//                        }
//                    }
//                    for(auto it=Tree[rank[t]].vert.begin();it!=Tree[rank[t]].vert.end();++it){
//                        cout<<"sc: "<<t<<" "<<it->first<<" "<<it->second.first<<"("<<it->second.second<<")"<<endl;
//                    }
//                    cout<<"disB: "<<s<<" "<<bv<<" "<<Tree[rank[s]].disInf[BoundVertexMap[3][bv]]<<endl;

                    if(ifDebug){
                        CorrectnessCheck(100);
                    }
//                    exit(0);
                }
                StagePerformanceShow(batchNum,stageUpdateT,stageQueryT);
                cout<<"\nPartiNum: "<<partiNum<<". Overall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average Increase batch update Time: "<<runT2/batchNum<<" s."<<endl;
            }


            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

void Graph::StagePerformanceShow(int batchNum, vector<double>& stageUpdateT, vector<double>& stageQueryT){
    if(algoChoice==1){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]/batchNum<<" ms."<<endl;
        cout<<"Stage 2 (CH). Average update time: "<<stageUpdateT[CH]/batchNum<<" s; average query time: "<<1000*stageQueryT[CH]/batchNum<<" ms."<<endl;
        cout<<"Stage 3 (H2H). Average duration: "<<stageUpdateT[H2H]/batchNum<<" s; average query time: "<<1000*stageQueryT[H2H]/batchNum<<" ms."<<endl;
    }else if(algoChoice==3){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]/batchNum<<" ms."<<endl;
        cout<<"Stage 2 (PCH). Average update time: "<<stageUpdateT[PCH_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PCH_No]/batchNum<<" ms."<<endl;
        cout<<"Stage 3 (No-boundary PMHL). Average update time: "<<stageUpdateT[PH2H_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_No]/batchNum<<" ms."<<endl;
        cout<<"Stage 4 (Post-boundary PMHL). Average update time: "<<stageUpdateT[PH2H_Post]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Post]/batchNum<<" ms."<<endl;
        cout<<"Stage 5 (Cross-boundary PMHL). Average duration: "<<stageUpdateT[PH2H_Cross]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Cross]/batchNum<<" ms."<<endl;
    }else if(algoChoice==5){
        cout<<"\nStage 1 (BiDijkstra). Average update time: "<<stageUpdateT[Dijk]/batchNum<<" s; average query time: "<<1000*stageQueryT[Dijk]/batchNum<<" ms."<<endl;
        cout<<"Stage 2 (PCH). Average update time: "<<stageUpdateT[PCH_No]/batchNum<<" s; average query time: "<<1000*stageQueryT[PCH_No]/batchNum<<" ms."<<endl;
        cout<<"Stage 3 (Post-boundary PostMHL). Average update time: "<<stageUpdateT[PH2H_Post]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Post]/batchNum<<" ms."<<endl;
        cout<<"Stage 4 (Cross-boundary PostMHL). Average duration: "<<stageUpdateT[PH2H_Cross]/batchNum<<" s; average query time: "<<1000*stageQueryT[PH2H_Cross]/batchNum<<" ms."<<endl;
    }
}

//function for throughput test of decrease update
void Graph::PMHLBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false, boost::ref(updatedSCs[pid]) ));
        }
        thread.join_all();
    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }

            if(newdis<olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                sm->notify();
            }
        }
    }
    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<<updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);//weightOverlay collect the changed edges on overlay graph
    }

    cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,false);
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]=tt2.GetRuntime();
    tt2.start();
//    cout<<"algoQuery: PCH-No"<<endl;

    DecreaseOverlayBatchLabel(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchLabel, this, boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], boost::ref(ProBeginVertexSetParti[pid]), boost::ref(vertexIDChLParti[pid]) ));
        }
        thread.join_all();
    }

    algoQuery=PH2H_No;
    tt2.stop();
    stageDurations[PCH_No]=tt2.GetRuntime();
    // repair the partition index
    if(algoUpdate>=PH2H_Post){
        tt2.start();
        Repair_PartiIndex(true, false, partiBatch);//post
        algoQuery=PH2H_Post;
        tt2.stop();
        stageDurations[PH2H_No]=tt2.GetRuntime();
        if(algoUpdate==PH2H_Cross){
            tt2.start();
//            RefreshExtensionLabels(partiBatch);//extend
            RefreshExtensionLabelsNoAllPair(partiBatch);//extend
            algoQuery=PH2H_Cross;
            tt2.stop();
            stageDurations[PH2H_Post]=tt2.GetRuntime();
            tt2.start();
        }
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_No]+stageDurations[PH2H_Post];
    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

//function for throughput test of decrease update
void Graph::PMHLBatchUpdateDecOpt(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
//            cout<<"Parti update: "<<pid1<<". "<<a<<" "<<b<<" "<<newW<<endl;
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), true , boost::ref(updatedSCs[pid])));
        }
        thread.join_all();
    }

    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }

            if(newdis<olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                sm->notify();
            }
        }
    }
    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<<updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);//weightOverlay collect the changed edges on overlay graph
    }
    cout<<"overlayBatch size: "<<overlayBatch.size()<<endl;
    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,false);
    algoQuery=PCH_No;
//    cout<<"algoQuery: PCH-No"<<endl;
    if(algoUpdate>PCH_No){
        DecreaseOverlayBatchLabel(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL
    }
    // repair the partition index
    if(algoUpdate>=PH2H_Post){

        Repair_PartiIndex(true, false, partiBatch);//post
//        Repair_PartiIndexForOpt(false, false, partiBatch);//post


        algoQuery=PH2H_Post;
        if(algoUpdate==PH2H_Cross){
//            RefreshExtensionLabels(partiBatch);//extend
            RefreshExtensionLabelsNoAllPair(partiBatch);//extend
            algoQuery=PH2H_Cross;
        }
    }

    tt.stop();
    runT1+=tt.GetRuntime();
    cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
}

//function for throughput test of decrease update
void Graph::PostMHLBatchUpdateDec(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTags[a].first, pid2=PartiTags[b].first;
        if(pid1==-1 && pid2==-1){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(pid1!=-1 && pid2==-1){
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid1) == partiBatch.end()){
                    partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid1].emplace_back(wBatch[k]);
            }else if(pid1==-1 && pid2!=-1){
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid2) == partiBatch.end()){
                    partiBatch.insert({pid2,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid2].emplace_back(wBatch[k]);
            }else
            {
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    if(partiBatch.find(pid1) == partiBatch.end()){
                        partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                    }
                    partiBatch[pid1].emplace_back(wBatch[k]);
                }else{
                    cout<<"Wrong for this edge update. "<<pid1<<" "<<pid2<<endl; exit(1);
                }

            }
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();
    for(int i=0;i<partiNum;++i){
        ProBeginVertexSetParti[i].clear(); vertexIDChLParti[i].clear();
        ProBeginVertexSetPartiExtend[i].clear();
    }
    for(int i=0;i<Tree.size();++i){
        Tree[i].DisRe.clear();
        Tree[i].DisRePost.clear();
    }
    // Step 1: Update shortcuts
    // Approach 1: partitioned version
    DecreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, true);//multi-thread
//    DecreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, false);//single-thread
    cout<<"overlayBatch size: "<<overlayBatch.size()<<endl;
    DecreaseOverlayBatchPostMHL(overlayBatch,Tree,rank,heightMax,false);

    // Approach 2: non-partitioned version
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,true);//direct bottom-up, with label construction
//    DecreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,false);//direct bottom-up, without label construction
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]=tt2.GetRuntime();
    cout<<"algoQuery: PCH-No"<<endl;

    //unparalleled version
//    if(algoUpdate>=PH2H_Post){
//        double tOverlay=0,tPost=0,tCross=0;
//        tt2.start();
//        // Step 2: Update overlay index
//        DecreaseOverlayBatchLabelPostMHL(Tree,rank,heightMax,ProBeginVertexSetOverlay,vertexIDChLOverlay);//top-down update for overlay
//        tt2.stop();
//        tOverlay=tt2.GetRuntime();
//        cout<<"Overlay index update: "<<tOverlay<<endl;
//        tt2.start();
//        // Step 3: Update partition index
//        Repair_PartiIndexPostMHLPost(true, false, partiBatch, tPost);//post
////        Repair_PartiIndexPostMHLPost(false, false, partiBatch);//post
//        tt2.stop();
//        stageDurations[PCH_No]=tt2.GetRuntime();
//        if(algoUpdate==PH2H_Cross){
//            tt2.start();
//            RefreshExtensionLabelsPostMHL(true, false, tCross);//extend, single-thread
////            RefreshExtensionLabelsPostMHL(false, false);//extend, single-thread
//            tt2.stop();
//            stageDurations[PH2H_Post]=tt2.GetRuntime();
//        }
//    }
    //parallel version
    double tOverlay=0,tPost=0,tCross=0;
    if(algoUpdate==PH2H_Post) {
        tt2.start();
        IncreaseOverlayBatchLabelPostMHL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
        tt2.stop();
        tOverlay = tt2.GetRuntime();
        cout << "Overlay label update time: " << tOverlay << " s." << endl;
        tt2.start();
        Repair_PartiIndexPostMHLPost(true, false, partiBatch, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, false, partiBatch);//post
        tt2.stop();
        stageDurations[PCH_No] = tt2.GetRuntime() + tOverlay;
    }else if(algoUpdate==PH2H_Cross){
        tt2.start();
        IncreaseOverlayBatchLabelPostMHL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
        tt2.stop();
        tOverlay = tt2.GetRuntime();
        cout << "Overlay label update time: " << tOverlay << " s." << endl;

        boost::thread_group thread;
        thread.add_thread(new boost::thread(&Graph::Repair_PartiIndexPostMHLPost, this, true, false, boost::ref(partiBatch), boost::ref(tPost) ));
        thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelsPostMHL, this, true, false, boost::ref(tCross)));
        thread.join_all();
        if(tPost<tCross){
            stageDurations[PCH_No] = tOverlay+tPost;
            stageDurations[PH2H_Post]=tCross-tPost;
        }else{
            cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
            stageDurations[PCH_No] = tOverlay+tCross;
            stageDurations[PH2H_Post]=0;
        }
    }


    tt.stop();
//    runT1+=tt.GetRuntime();
    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

//function for throughput test of increase update, optimized version
void Graph::PostMHLBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    tAncestor=0, tBoundary=0;
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTags[a].first, pid2=PartiTags[b].first;
        if(pid1==-1 && pid2==-1){
//            cout<<"Overlay update: "<<a<<" "<<b<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(pid1!=-1 && pid2==-1){
//                cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid1) == partiBatch.end()){
                    partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid1].emplace_back(wBatch[k]);
            }else if(pid1==-1 && pid2!=-1){
//                cout<<"Parti update: "<<pid2<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                if(partiBatch.find(pid2) == partiBatch.end()){
                    partiBatch.insert({pid2,vector<pair<pair<int, int>, pair<int, int>>>()});
                }
                partiBatch[pid2].emplace_back(wBatch[k]);
            }else
            {
                if(pid1==pid2){
//                    cout<<"Parti update: "<<pid1<<". "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<newW<<endl;
                    if(partiBatch.find(pid1) == partiBatch.end()){
                        partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
                    }
                    partiBatch[pid1].emplace_back(wBatch[k]);
                }else{
                    cout<<"Wrong for this edge update. "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<")"<<endl; exit(1);
                }

            }
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    ProBeginVertexSetOverlay.clear();
    for(int i=0;i<partiNum;++i){
        ProBeginVertexSetParti[i].clear();
        ProBeginVertexSetPartiExtend[i].clear();
    }
    for(int i=0;i<Tree.size();++i){
        Tree[i].DisRe.clear();
        Tree[i].DisRePost.clear();
    }
    // Step 1: Update shortcuts
    // Approach 1: partitioned version
    IncreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, true);//multi-thread, bottom-up shortcut update
//    IncreasePartiBatchUpdateCheckPostMHL(partiBatch, overlayBatch, false);//single-thread
    cout<<"overlayBatch size: "<<overlayBatch.size()<<endl;
    IncreaseOverlayBatchPostMHL(overlayBatch,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);

    // Approach 2: non-partitioned version
//    IncreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);//direct bottom-up, with label construction
//    IncreaseH2HBatch(wBatch,Neighbor,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);//direct bottom-up, without label construction
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]=tt2.GetRuntime();
//    cout<<"algoQuery: PCH-No"<<endl;


    // repair the partition index, unparalleled version
//    if(algoUpdate>=PH2H_Post){
//        tt2.start();
//        IncreaseOverlayBatchLabelPostMHL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
//        tt2.stop();
//        double toverlay=tt2.GetRuntime();
//        cout<<"Overlay label update time: "<<toverlay<<" s."<<endl;
////        for(int i=0;i<Tree[rank[67975]].disPost.size();++i){
////            if(Tree[rank[67975]].vAncestor[i]==67432){
////                cout<<"B Cnt: "<<Tree[rank[67975]].disPost[i]<<" "<<Tree[rank[67975]].cntPost[i]<<endl;
////            }
////        }
//        double tPost=0,tCross=0;
//        tt2.start();
//        Repair_PartiIndexPostMHLPost(true, true, partiBatch, tPost);//post
////        Repair_PartiIndexPostMHLPost(false, true, partiBatch);//post
//        tt2.stop();
//        cout<<"Post-boundary index update time: "<<tt2.GetRuntime()<<" s; Boundary array update: "<<tBoundary<<" s; Ancestor array update: "<<tAncestor<<" s."<<endl;
//        stageDurations[PCH_No]=tt2.GetRuntime()+toverlay;
//        if(algoUpdate==PH2H_Cross){
//            tt2.start();
//            RefreshExtensionLabelsPostMHL(true, true, tCross);//extend, single-thread
////            RefreshExtensionLabelsPostMHL(false, true);//extend, single-thread
//            algoQuery=PH2H_Cross;
//            tt2.stop();
//            stageDurations[PH2H_Post]=tt2.GetRuntime();
//        }
//    }
    double tOverlay=0,tPost=0,tCross=0;
    if(algoUpdate==PH2H_Post) {
        tt2.start();
        IncreaseOverlayBatchLabelPostMHL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
        tt2.stop();
        tOverlay = tt2.GetRuntime();
        cout << "Overlay label update time: " << tOverlay << " s." << endl;
        tt2.start();
        Repair_PartiIndexPostMHLPost(true, true, partiBatch, tPost);//post
//        Repair_PartiIndexPostMHLPost(false, true, partiBatch);//post
        tt2.stop();
        stageDurations[PCH_No] = tt2.GetRuntime() + tOverlay;
    }else if(algoUpdate==PH2H_Cross){
        tt2.start();
        IncreaseOverlayBatchLabelPostMHL(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
        tt2.stop();
        tOverlay = tt2.GetRuntime();
        cout << "Overlay label update time: " << tOverlay << " s." << endl;

        boost::thread_group thread;
        thread.add_thread(new boost::thread(&Graph::Repair_PartiIndexPostMHLPost, this, true, true, boost::ref(partiBatch), boost::ref(tPost) ));
        thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelsPostMHL, this, true, true, boost::ref(tCross)));
        thread.join_all();
        if(tPost<tCross){
            stageDurations[PCH_No] = tOverlay+tPost;
            stageDurations[PH2H_Post]=tCross-tPost;
        }else{
            cout<<"!!! Cross-boundary update is faster than post-boundary update! "<<tPost<<" "<<tCross<<endl;
            stageDurations[PCH_No] = tOverlay+tCross;
            stageDurations[PH2H_Post]=0;
        }
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post];
    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;

//    CorrectnessCheck(100);
}

//function for throughput test of increase update
void Graph::PMHLBatchUpdateInc(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    Timer tt2;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    tt2.start();
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), false, boost::ref(updatedSCs[pid]) ));
        }
        thread.join_all();

    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;

        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }
//        cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
            if(newdis>olddis){//if '=', not problem; if '<', problem
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<< updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);
    }

    cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;
    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    algoQuery=PCH_No;
    tt2.stop();
    stageDurations[Dijk]=tt2.GetRuntime();
    tt2.start();
    IncreaseOverlayBatchLabel(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
//        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchLabel, this, boost::ref(Trees[pid]), boost::ref(ranks[pid]), heightMaxs[pid], boost::ref(ProBeginVertexSetParti[pid]), boost::ref(VidtoTNidP) ));//vector<Node> &Tree, vector<int> &rank, int heightMax, int& checknum, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid
        }
        thread.join_all();
    }

    algoQuery=PH2H_No;
    tt2.stop();
    stageDurations[PCH_No]=tt2.GetRuntime();
    // repair the partition index
    if(algoUpdate>=PH2H_Post){
//        Trees=TreesNo;
        tt2.start();
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);
        algoQuery=PH2H_Post;
        tt2.stop();
        stageDurations[PH2H_No]=tt2.GetRuntime();
        if(algoUpdate==PH2H_Cross){
            tt2.start();
//            RefreshExtensionLabels(partiBatch);
            RefreshExtensionLabelsNoAllPair(partiBatch);
            algoQuery=PH2H_Cross;
            tt2.stop();
            stageDurations[PH2H_Post]=tt2.GetRuntime();
        }
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
    runT1=runT1+stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_No]+stageDurations[PH2H_Post];
    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[PCH_No]+stageDurations[PH2H_No]+stageDurations[PH2H_Post]<<" s; "<<tt.GetRuntime()<<" s."<<endl;

//    CorrectnessCheck(100);
}

//function for throughput test of increase update, optimized version
void Graph::PMHLBatchUpdateIncOpt(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double &runT1) {
    Timer tt;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
//            cout<<"Overlay update: "<<a<<"("<<PartiTag[a].first<<","<<PartiTag[a].second<<") "<<b<<"("<<PartiTag[b].first<<","<<PartiTag[b].second<<") "<<oldW<<" "<<NeighborsOverlay[a][b]<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
//            cout<<"Parti update: "<<pid1<<". "<<a<<"("<<PartiTag[a].first<<","<<PartiTag[a].second<<") "<<" "<<b<<"("<<PartiTag[b].first<<","<<PartiTag[b].second<<") "<<" "<<oldW<<" "<<newW<<endl;
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    vUpdated.assign(node_num, false);
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
//    cout<<"Before parti CH, size of overlayBatch: "<<overlayBatch.size()<<endl;
    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheckCH, this, pid, boost::ref(it->second), boost::ref(overlayBatch), true, boost::ref(updatedSCs[pid]) ));
//            IncreasePartiBatchUpdateCheckCH(pid, it->second, overlayBatch, true );
        }
        thread.join_all();

    }
    map<pair<int,int>,pair<int,int>> updateSCTrue;
    int updateSCSize=0;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        vector<int> Bid=BoundVertex[pid];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;

        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second; newdis=it->second;
            updateSCSize++;
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
                cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
            }
//        cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
            if(newdis>olddis){//if '=', not problem; if '<', problem
//                sm->wait();
                if(updateSCTrue.find(make_pair(bid1,bid2))==updateSCTrue.end()){
                    updateSCTrue.insert({make_pair(bid1,bid2), make_pair(olddis,newdis)});
                }else if(updateSCTrue[make_pair(bid1,bid2)].second>newdis){//if found and newdis is smaller
//                    cout<<"More than one supportive vertices. "<<bid1<<" "<<bid2<<" "<<updateSCTrue[make_pair(bid1,bid2)].second<<" "<<newdis<<endl;
                    updateSCTrue[make_pair(bid1,bid2)].second=newdis;
                }
//                overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
    cout<<"updateSCTrue size: "<<updateSCTrue.size()<<" "<<updateSCSize<<endl;
    for(auto it=updateSCTrue.begin();it!=updateSCTrue.end();++it){
        overlayBatch.emplace_back(it->first,it->second);
    }

    cout<<"OverlayBatch size: "<<overlayBatch.size()<<endl;

//    cout<<"After parti CH, size of overlayBatch: "<<overlayBatch.size()<<endl;
    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,false);
    algoQuery=PCH_No;
    if(algoUpdate>PCH_No) {
        IncreaseOverlayBatchLabel(Tree, rank, heightMax, ProBeginVertexSetOverlay, VidtoTNid);
    }

    // repair the partition index
    if(algoUpdate>=PH2H_Post){
        Repair_PartiIndex(true, true, partiBatch);//post
//        Repair_PartiIndexForOpt(false, true, partiBatch);//post

        algoQuery=PH2H_Post;
        if(algoUpdate==PH2H_Cross){
//            RefreshExtensionLabels(partiBatch);
            RefreshExtensionLabelsNoAllPair(partiBatch);
            algoQuery=PH2H_Cross;
        }
    }

    tt.stop();
    runT1+=tt.GetRuntime();
    cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;

//    CorrectnessCheck(100);
}

//function for throughput test of decrease update
void Graph::DecBatchThroughput(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    vUpdated.assign(node_num, false);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();
    }

    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,true);

    algoQuery=PH2H_No;
    // repair the partition index
    if(algoUpdate>=PH2H_Post){
        Repair_PartiIndex(true, false, partiBatch);//post
        algoQuery=PH2H_Post;
        if(algoUpdate==PH2H_Cross){
            RefreshExtensionLabels(partiBatch);//extend
            algoQuery=PH2H_Cross;
        }
    }

    tt.stop();
    runT1+=tt.GetRuntime();
    cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
}
//function for throughput test of increase update
void Graph::IncBatchThroughput(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    tt.start();
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int oldW = wBatch[k].second.first;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }

        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
//            cout<<"Overlay update. "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<oldW<<" "<<newW<<endl;
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
//            cout<<"Parti update. "<<a<<"("<<pid1<<") "<<b<<"("<<pid2<<") "<<oldW<<" "<<newW<<endl;
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }
    algoQuery=Dijk;
    vUpdated.assign(node_num, false);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();

    }

    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);
    algoQuery=PH2H_No;

    // repair the partition index
    if(algoUpdate>=PH2H_Post){
//        Trees=TreesNo;
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);
        algoQuery=PH2H_Post;
        if(algoUpdate==PH2H_Cross){
            RefreshExtensionLabels(partiBatch);
            algoQuery=PH2H_Cross;
        }
    }

    tt.stop();
    runT1+=tt.GetRuntime();
    cout<<"Batch "<<batch_i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
}


//function for throughput test of decrease update
void Graph::DecBatchThroughputNP(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    Timer tt2;
    tt.start();

    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompp> OC; //OC.clear();//vertexID in decreasing node order

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    int a,b,oldW,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second; oldW=wBatch[k].second.first;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
                Neighbor[b][i].second=newW;
                break;
            }
        }
    }

    algoQuery=Dijk;

    tt2.start();
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second; oldW=wBatch[k].second.first;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }


        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompp(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }
    }

    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw;
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }else{
                    Cw=Vert[j].second.first;
                }
            }

            if(Tree[rank[ProID]].dis[cidH]>Cw){
                Tree[rank[ProID]].dis[cidH]=Cw;
                Tree[rank[ProID]].FN[cidH]=true;
                ProIDdisCha=true;
                Tree[rank[ProID]].DisRe.insert(Cid);
            }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                Tree[rank[ProID]].FN[cidH]=true;
            }

            int hid,hidHeight,lid,lidHeight,wsum;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                if(Hnei.find(hid)!=Hnei.end()){
                    wsum=Cw+Hnei[hid];
                    if(wsum<Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.first=wsum;
                        Tree[rank[Cid]].vert[j].second.second=1;
                        SCre[Cid].insert(hid);
                        OC.insert(OrderCompp(Cid));
                    }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;lidHeight=Tree[rank[lid]].height-1;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid]].vert[k].second.first>wsum){
                            Tree[rank[lid]].vert[k].second.first=wsum;
                            Tree[rank[lid]].vert[k].second.second=1;
                            SCre[lid].insert(Cid);
                            OC.insert(OrderCompp(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            Tree[rank[lid]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChL.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQuery(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;
        }
    }

    algoQuery=CH;
    tt2.stop();
    stageDurations[Dijk]=tt2.GetRuntime();
    tt2.start();
    //cout<<"Finish bottom-up refresh"<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5H2H(rank[ProBeginVertexID], linee, vertexIDChL,checkedDis);
    }
    //return checkedDis.size();
    algoQuery=H2H;
    tt2.stop();
    stageDurations[CH]=tt2.GetRuntime();
    tt.stop();
//    runT1+=tt.GetRuntime();
    runT1=runT1+stageDurations[Dijk]+stageDurations[CH];
    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[CH]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

void Graph::EachNodeProBDis5H2H(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
    bool ProIDdisCha=false;

    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }
            }
        }
    }
    else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                        Tree[child].FN[i]=false;
                        ProIDdisCha=true;
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5H2H(Tree[child].ch[i], line, vertexIDChL,checkedDis);
    }
    line.pop_back();

}

//function for throughput test of decrease update
void Graph::IncBatchThroughputNP(vector<pair<pair<int, int>, pair<int, int>>> &wBatch, int batch_i, double& runT1) {
    Timer tt;
    Timer tt2;
    tt.start();

    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]), original weight of the affected shortcut, maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    //NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear(); the affected shortcut pair
    set<int> ss; ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompp> OC; OC.clear();//the lower-order vertex of the affected shortcut, vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW<newW){
            for(int i=0;i<Neighbor[a].size();i++){
                if(Neighbor[a][i].first==b){
                    Neighbor[a][i].second=newW;
                    break;
                }
            }
            for(int i=0;i<Neighbor[b].size();i++){
                if(Neighbor[b][i].first==a){
                    Neighbor[b][i].second=newW;
                    break;
                }
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }

        }
    }

    algoQuery=Dijk;
    tt2.start();
    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW<newW){

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }

            for(int i=0;i<Tree[rank[lid]].vert.size();i++){
                if(Tree[rank[lid]].vert[i].first==hid){
                    if(Tree[rank[lid]].vert[i].second.first==oldW){
                        Tree[rank[lid]].vert[i].second.second-=1;
                        if(Tree[rank[lid]].vert[i].second.second<1){//the shortcut needs update, should be increased
                            OCdis[make_pair(lid,hid)]=oldW;//original weight of the affected shortcut
                            SCre[lid].insert(hid);//the affected shortcut pair
                            OC.insert(OrderCompp(lid));//the lower-order vertex of the affected shortcut
                        }
                    }
                    break;
                }
            }
        }
    }

    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID; vector<int> line;
    /// Shortcut update
    while(!OC.empty()){
        ProID=(*OC.begin()).x;//from the lowest-order vertex
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        // get the ancestors of ProID, each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[pachid]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.emplace_back(Vert[j].first,Vert[j].second.first);
                }
            }
            //check the affected shortcuts
            int hid,lid;
            for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                hid=Tree[rank[Cid]].vert[j].first;
                if(Hnei.find(hid)!=Hnei.end()){
                    if(Cw+Hnei[hid]==Tree[rank[Cid]].vert[j].second.first){
                        Tree[rank[Cid]].vert[j].second.second-=1;
                        if(Tree[rank[Cid]].vert[j].second.second<1){
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompp(Cid));
                            OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid=Lnei[j].first;
                for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                    if(Tree[rank[lid]].vert[k].first==Cid){
                        if(Tree[rank[lid]].vert[k].second.first==Cw+Lnei[j].second){
                            Tree[rank[lid]].vert[k].second.second-=1;
                            if(Tree[rank[lid]].vert[k].second.second<1){
                                SCre[lid].insert(Cid);
                                OC.insert(OrderCompp(lid));
                                OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if(Tree[rank[ProID]].FN[cidH]){//if the distance label is from shortcut, then the label may be affected.
                influence=true;
                //higher than Cid
                for(int i=0;i<cidH;i++){
                    if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                        Tree[rank[ProID]].cnt[i]-=1;
                    }
                }

                //equal to Cid
                Tree[rank[ProID]].FN[cidH]=false;//? may still be the source
                Tree[rank[ProID]].cnt[cidH]-=1;

                //lower than Cid
                for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                    if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                        Tree[rank[ProID]].cnt[i]-=1;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            Cw=INF; int countwt=0;

            for(int i=0;i<Neighbor[ProID].size();i++){
                if(Neighbor[ProID][i].first==Cid){
                    Cw=Neighbor[ProID][i].second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(!Wnodes.empty()){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first;
                    for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                        if(Tree[rank[wid]].vert[j].first==ProID){
                            ssw=Tree[rank[wid]].vert[j].second.first;
                        }
                        if(Tree[rank[wid]].vert[j].first==Cid){
                            wtt=Tree[rank[wid]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                if(Tree[rank[ProID]].vert[i].first==Cid){
                    Tree[rank[ProID]].vert[i].second.first=Cw;
                    Tree[rank[ProID]].vert[i].second.second=countwt;
                    break;
                }
            }
        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSet.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSet.size();i++){
                r=rank[ProBeginVertexSet[i]];
                if(LCAQuery(rnew,r)!=rnew){//if they are in different branches
                    ProBeginVertexSetNew.push_back(ProBeginVertexSet[i]);
                }
            }
            ProBeginVertexSet=ProBeginVertexSetNew;//identify the roots
        }

    }
    algoQuery=CH;
    tt2.stop();
    stageDurations[Dijk]=tt2.GetRuntime();

    if(algoUpdate==H2H){
        tt2.start();
        int ProBeginVertexID;
//    cout<<"Root number: "<<ProBeginVertexSet.size()<<endl;
        for(int i=0;i<ProBeginVertexSet.size();i++){//for each root
            ProBeginVertexID=ProBeginVertexSet[i];
//        cout<<i<<" "<<ProBeginVertexID<<" "<<Tree[rank[ProBeginVertexID]].height<<endl;
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1H2H(rank[ProBeginVertexID], linee,checknum);
        }
        //return checknum;
        algoQuery=H2H;
        tt2.stop();
        stageDurations[CH]=tt2.GetRuntime();
    }

    tt.stop();
//    runT1+=tt.GetRuntime();
    runT1=runT1+stageDurations[Dijk]+stageDurations[CH];
    cout<<"Batch "<<batch_i<<". Update time: "<<stageDurations[Dijk]+stageDurations[CH]<<" s; "<<tt.GetRuntime()<<" s."<<endl;
}

void Graph::eachNodeProcessIncrease1H2H(int children, vector<int>& line, int& changelabel){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
    for(int i=0;i<Tree[children].dis.size();i++){
        if(Tree[children].cnt[i]==0){//if the distance label to i-th ancestor should be maintained
//        if(true){
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int PID;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){//for the tree node that contains childID as vert element
                PID=VidtoTNid[childID][k];
                if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){//if label is from shortcut
                    Tree[PID].cnt[i]-=1;
                }
            }

            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                PID=VidtoTNid[line[i]][k];
//                if(Tree[PID].height>Tree[children].height){///modified for correctness, PID may not be the descendant of children
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){///modified for correctness, PID may not be the descendant of children
                    if(PID>Tree.size()){
                        cout<<"PID error! "<<PID<<" "<<Tree.size()<<endl; exit(1);
                    }
                    if(childH>Tree[PID].dis.size()){
                        cout<<"childH error! "<<childH<<" "<<Tree[PID].dis.size()<<": "<<children<<"("<<Tree[children].height<<") "<<PID<<"("<<Tree[PID].height<<")"<<endl; exit(1);
                    }
                    if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){///
                        Tree[PID].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[b]].height-1;
                if(bH<i){
                    if(Dvb+Tree[rank[line[i]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[line[i]]].dis[bH];
                        count=1;
                    }else if(Dvb+Tree[rank[line[i]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{
                    if(Dvb+Tree[rank[b]].dis[i]<dis){
                        dis=Dvb+Tree[rank[b]].dis[i];
                        count=1;
                    }else if(Dvb+Tree[rank[b]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) {
                Tree[children].FN[i]=true;
            }
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1H2H(Tree[children].ch[i],line,changelabel);
    }
    line.pop_back();
}

//function for efficiency test
unsigned long long Graph::EffiCheckThroughput(vector<pair<int,int>>& ODpair, Timer& tRecord, int batchInterval, unsigned long long& throughputNum) {
    bool ifDebug=false;
//    ifDebug=true;
    int s, t;
    double runT = 0, runT0=0, runT1=0, runT2=0, runT3=0, runT4=0, runT5=0;
    int d1, d2;
    Timer tt;
    unsigned long long runtimes = 0, runtimes0=0, runtimes1=0, runtimes2=0, runtimes3=0, runtimes4=0, runtimes5=0;
    vector<int> results(ODpair.size(), -1);
    int i = 0;
    double tNow = 0;
    if(algoChoice==1){//CH+H2H
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == CH) {//CH
                tt.start();
                d2 = QueryCHWP(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == H2H) {//H2H
                tt.start();
                d2 = QueryH2H(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT2 += tt.GetRuntime();
                ++runtimes2;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery==Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
                }
            }

        }

        cout<<"CH+H2H. Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (CH): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (H2H): Duration: "<<runT2<<" s; Throughput number: "<<runtimes2<<" ; query time: "<<1000 * runT2 / runtimes2 << " ms."<<endl;
    }
    else if(algoChoice==2){//PH2H
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == PH2H_No) {//PH2H-No
                tt.start();
                d2 = Query(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == PH2H_Post) {//PH2H-Post
                tt.start();
                d2 = Query(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT2 += tt.GetRuntime();
                ++runtimes2;
            }
            else if (algoQuery == PH2H_Cross) {//PH2H-Extend
                tt.start();
                d2 = Query(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT3 += tt.GetRuntime();
                ++runtimes3;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery==Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" ("<<PartiTag[s].first<<", "<<PartiTag[s].second<<") "<<t<<" ("<<PartiTag[t].first<<", "<<PartiTag[t].second<<") "<<d2<<" "<<d1<<endl;
                    exit(1);
                }
            }
        }

        cout<<"PH2H. Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (No-boundary): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (Post-boundary): Duration: "<<runT2<<" s; Throughput number: "<<runtimes2<<" ; query time: "<<1000 * runT2 / runtimes2 << " ms."<<endl;
        cout<<"Stage 4 (Extension): Duration: "<<runT3<<" s; Throughput number: "<<runtimes3<<" ; query time: "<<1000 * runT3 / runtimes3 << " ms."<<endl;
    }
    else if(algoChoice==3){//PCH+PH2H
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == PCH_No) {//PCH-No
                tt.start();
                d2 = QueryPMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == PH2H_No) {//PH2H-No
                tt.start();
                d2 = QueryPMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT2 += tt.GetRuntime();
                ++runtimes2;
            }
            else if (algoQuery == PH2H_Post) {//PH2H-Post
                tt.start();
                d2 = QueryPMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT4 += tt.GetRuntime();
                ++runtimes4;
            }
            else if (algoQuery == PH2H_Cross) {//PH2H-Extend
                tt.start();
                d2 = QueryPMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT5 += tt.GetRuntime();
                ++runtimes5;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery== Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" ("<<PartiTag[s].first<<", "<<PartiTag[s].second<<") "<<t<<" ("<<PartiTag[t].first<<", "<<PartiTag[t].second<<") "<<d2<<" "<<d1<<endl;
                    exit(1);
                }
            }
        }

        cout<<"PCH+PH2H. Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (PCH-No): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-No): Duration: "<<runT2<<" s; Throughput number: "<<runtimes2<<" ; query time: "<<1000 * runT2 / runtimes2 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Post): Duration: "<<runT4<<" s; Throughput number: "<<runtimes4<<" ; query time: "<<1000 * runT4 / runtimes4 << " ms."<<endl;
        cout<<"Stage 5 (PH2H-Extend): Duration: "<<runT5<<" s; Throughput number: "<<runtimes5<<" ; query time: "<<1000 * runT5 / runtimes5 << " ms."<<endl;
    }
    else if(algoChoice==4){//PCH+PH2H with optimization
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == PCH_No) {//PCH-No
                tt.start();
                d2 = QueryPMHLOpt(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == PH2H_Post) {//PCH-Post
                tt.start();
                d2 = QueryPMHLOpt(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT3 += tt.GetRuntime();
                ++runtimes3;
            }
            else if (algoQuery == PH2H_Post) {//PH2H-Post
                tt.start();
                d2 = QueryPMHLOpt(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT4 += tt.GetRuntime();
                ++runtimes4;
            }
            else if (algoQuery == PH2H_Cross) {//PH2H-Extend
                tt.start();
                d2 = QueryPMHLOpt(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT5 += tt.GetRuntime();
                ++runtimes5;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery== Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" ("<<PartiTag[s].first<<", "<<PartiTag[s].second<<") "<<t<<" ("<<PartiTag[t].first<<", "<<PartiTag[t].second<<") "<<d2<<" "<<d1<<endl;
                    exit(1);
                }
            }
        }

        cout<<"PCH+PH2H (Opt). Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (PCH-No): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-Post): Duration: "<<runT4<<" s; Throughput number: "<<runtimes4<<" ; query time: "<<1000 * runT4 / runtimes4 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Extend): Duration: "<<runT5<<" s; Throughput number: "<<runtimes5<<" ; query time: "<<1000 * runT5 / runtimes5 << " ms."<<endl;
    }
    else if(algoChoice==5){//PostMHL
        while (true) {
            tRecord.stop();
            tNow = tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if (tNow > batchInterval) {
                break;
            }

            s = ODpair[i].first;
            t = ODpair[i].second;

            if (algoQuery == Dijk) {//Dijkstra
                tt.start();
                d2 = Dijkstra(s, t, Neighbor);
//                d2 = Astar(s, t, Neighbor);
                tt.stop();
                runT += tt.GetRuntime();
                runT0 +=tt.GetRuntime();
                ++runtimes0;
            } else if (algoQuery == PCH_No) {//PCH-No
                tt.start();
                d2 = QueryPostMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT1 += tt.GetRuntime();
                ++runtimes1;
            }
            else if (algoQuery == PH2H_Post) {//PH2H-Post
                tt.start();
                d2 = QueryPostMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT4 += tt.GetRuntime();
                ++runtimes4;
            }
            else if (algoQuery == PH2H_Cross) {//PH2H-Extend
                tt.start();
                d2 = QueryPostMHL(s, t);
                tt.stop();
                runT += tt.GetRuntime();
                runT5 += tt.GetRuntime();
                ++runtimes5;
            }

            results[i] = d2;
            ++runtimes;
            ++i;
            if (i == ODpair.size()) {
                i = 0;
            }
            if(ifDebug){
                if(algoQuery== Dijk){
                    continue;
                }
                d1=Dijkstra(s,t,Neighbor);
                if(d1!=d2){
                    cout<<"Effi Test. Incorrect! Algorithm "<<algoQuery<<": "<<s<<" ("<<PartiTags[s].first<<") "<<t<<" ("<<PartiTags[t].first<<") "<<d2<<" "<<d1<<endl;
                    QueryPostMHLDebug(s,t);
                    exit(1);
                }
            }
        }

        cout<<"PostMHL. Throughput number: "<<runtimes;
        if(runtimes>0) {
            cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
        }
        cout<<endl;
        cout<<"Stage 1 (Dijkstra): Duration: "<<runT0<<" s; Throughput number: "<<runtimes0<<" ; query time: "<<1000 * runT0 / runtimes0 << " ms."<<endl;
        cout<<"Stage 2 (PCH-No): Duration: "<<runT1<<" s; Throughput number: "<<runtimes1<<" ; query time: "<<1000 * runT1 / runtimes1 << " ms."<<endl;
        cout<<"Stage 3 (PH2H-Post): Duration: "<<runT4<<" s; Throughput number: "<<runtimes4<<" ; query time: "<<1000 * runT4 / runtimes4 << " ms."<<endl;
        cout<<"Stage 4 (PH2H-Extend): Duration: "<<runT5<<" s; Throughput number: "<<runtimes5<<" ; query time: "<<1000 * runT5 / runtimes5 << " ms."<<endl;
    }

    throughputNum+=runtimes;
    return runtimes;
}

void Graph::IndexMaintenance(int updateType, int updateSize, bool ifBatch, int batchSize) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = sourcePath+dataset + ".update";
    bool ifDebug=false;
    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand(0);
    int updateBatch=1;
    updateBatch=max(updateBatch,updateSize/batchSize);
    cout<<"Update batch: "<<updateBatch<<" ; Batch size: "<<batchSize<<endl;
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);

    Timer tt;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 0:{
            break;
        }
        case 1:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;
            Graph g2=*this;
            if(ifBatch){//for batch update
                if(updateBatch*batchSize>updateData.size()){
                    updateBatch=floor(updateData.size()/batchSize);
                }
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<updateBatch;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        update_i=u*batchSize+i;
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;
                        oldW = updateData[update_i].second;
                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    }

//                    g2.QueryPartiPartiExtDebug(3321,212184);

                    tt.start();
                    g2.DecreaseBatch(wBatch);
//                    DecreaseBatch(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease update Time: "<<runT1/(updateBatch*batchSize)<<" s.\n"<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<updateBatch;u++){
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*0.5;
                    if(newW < 1) {
                        cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        exit(1);
                    }
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                    }
                    tt.start();
                    g2.DecreaseSingle(ID1,ID2,oldW,newW);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        g2.CorrectnessCheck(100);
                    }

                }

                cout<<"Average Decrease update Time: "<<runT1/updateBatch<<" s.\n"<<endl;
            }

//            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            ifIncrease=true;
            if(ifBatch){//for batch update
                if(updateBatch*batchSize>updateData.size()){
                    updateBatch=floor(updateData.size()/batchSize);
                }
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=1;u<updateBatch;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        update_i=u*batchSize+i;
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;
                        oldW = updateData[update_i].second;
                        newW=oldW*1.5;
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
//                        ++update_i;
                    }

//                    QueryCoreDebug(111009,104090);
//                    QueryPartiPartiExtDebug(107830,104726);

                    tt.start();
                    IncreaseBatch(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }
                cout<<"Average Increase update Time: "<<runT2/(updateBatch*batchSize)<<" s.\n"<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<updateBatch;u++){
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*2;
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<"("<<PartiTag[ID1].first<<") "<<ID2<<"("<<PartiTag[ID2].first<<") "<<oldW<<" "<<newW<<endl;
                    }

                    tt.start();
                    IncreaseSingle(ID1,ID2,oldW,newW);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }

                cout<<"Average Increase update Time: "<<runT2/updateBatch<<" s.\n"<<endl;
            }

            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

void Graph::DecreaseBatch(vector<pair<pair<int, int>, pair<int, int>>> &wBatch) {
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }
        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;
        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();
    }

    DecreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,true);

    // repair the partition index
    if(algoUpdate>=PH2H_Post){
        Repair_PartiIndex(true, false, partiBatch);//post
        if(algoUpdate==PH2H_Cross){
            RefreshExtensionLabels(partiBatch);//extend
        }
    }
}
//partition update of PH2H
void Graph::DecreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay){
    //partition batch decrease update
    vector<pair<pair<int,int>,int>> updatedSC;

    DecreasePartiBatch(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, true);

    vector<int> Bid=BoundVertex[pid];
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(int i=0;i<Bid.size();i++){
        bid1=Bid[i];
        for(int j=i+1;j<Bid.size();j++){
            bid2=Bid[j];
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                continue;//exit(1);
            }

            newdis=QueryH2HPartition(bid1,bid2,pid);
            if(newdis<olddis){
//                    cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
                sm->wait();
                weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
                sm->notify();
            }
        }
    }
}
//partition update PCH
void Graph::DecreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC){
    //partition batch decrease update
//    vector<pair<pair<int,int>,int>> updatedSC;


    if(ifOpt){
        DecreasePartiBatchForOpt(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, false, false);
    }else{
        DecreasePartiBatch(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, false);
    }


//    vector<int> Bid=BoundVertex[pid];
//    //check the boundary edge within partition
//    int bid1,bid2,olddis,newdis;
//    for(auto it=updatedSC.begin();it!=updatedSC.end();++it){
//        bid1=it->first.first, bid2=it->first.second; newdis=it->second;
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }
//
//        if(newdis<olddis){
////            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
////                NeighborsOverlay[bid1][bid2]=newdis;
////                NeighborsOverlay[bid2][bid1]=newdis;
//            sm->wait();
//            weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
//        }
//    }

}
//shortcut update for PostMHL

void Graph::DecreasePartiBatchUpdateCheckPostMHL(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int,int>,pair<int,int>>>& overlayBatch, bool ifParallel){
    //partition batch decrease update
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());

    if(ifParallel){
        if(!partiBatch.empty()){
            if(partiBatch.size()>threadnum){
                cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
            }
            cout<<"Update Partition number: "<<partiBatch.size()<<endl;
            boost::thread_group thread;
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
                thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchPostMHLShortcut, this, pid, boost::ref(it->second), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCs[pid]) ));
            }
            thread.join_all();
        }
    }else{
        cout<<"Single thread computation."<<endl;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            cout<<"Shortcut update of partition "<<pid<<endl;
            DecreasePartiBatchPostMHLShortcut(pid, it->second, Tree, rank, heightMax, updatedSCs[pid]);
        }
    }


//    cout<<"updatedSC size: "<<updatedSC.size()<<endl;
    map<pair<int,int>, int> updatedSCTrue;
    int updateSCSize=0;
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second, newdis=it->second;
            updateSCSize++;
            if(updatedSCTrue.find(make_pair(bid1,bid2))==updatedSCTrue.end()){//if not found
                updatedSCTrue.insert({make_pair(bid1,bid2),newdis});
            }else {//if found
//            cout<<bid1<<" "<<bid2<<": "<<updatedSCTrue[make_pair(bid1,bid2)]<<" "<<newdis<<endl;
                if(updatedSCTrue[make_pair(bid1,bid2)]>newdis){//if found
                    updatedSCTrue[make_pair(bid1,bid2)]=newdis;
                }
            }
        }
    }

    cout<<"updatedSCTrue size: "<<updatedSCTrue.size()<<" "<<updateSCSize<< endl;

    int sVertex;
    for(auto it=updatedSCTrue.begin();it!=updatedSCTrue.end();++it){
        bid1=it->first.first, bid2=it->first.second, newdis=it->second;
//        for(int i=0;i<SCconNodesMT[bid1][bid2].size();i++){
//            sVertex=SCconNodesMT[bid1][bid2][i].first;
//            if(NeighborsOverlay[bid1][bid2] > Tree[rank[sVertex]].)
//            SCconNodesMT[ID1][ID2].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
//        }

        olddis=-1;
        for(auto it2=Tree[rank[bid1]].vert.begin();it2!=Tree[rank[bid1]].vert.end();++it2){
            if(it2->first==bid2){
                olddis=it2->second.first;
                break;
            }
        }
        if(olddis==-1){
            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
        }
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }

        if(newdis<olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//            sm->wait();
            overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
        }else if(newdis>olddis){
            cout<<"Something wrong happens. "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl; exit(1);
        }
    }

}

void Graph::IncreaseBatch(vector<pair<pair<int, int>, pair<int, int>>> &wBatch) {
    map<int, vector<pair<pair<int, int>, pair<int, int>>>> partiBatch; partiBatch.clear();
    vector<pair<pair<int, int>, pair<int, int>>> overlayBatch;
    for(int k=0;k<wBatch.size();k++) {
        int a = wBatch[k].first.first;
        int b = wBatch[k].first.second;
        int newW = wBatch[k].second.second;
        for(int i=0;i<Neighbor[a].size();i++){
            if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
                Neighbor[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbor[b].size();i++){
            if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
                Neighbor[b][i].second=newW;
                break;
            }
        }
        int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
        if(pid1 != pid2){
            overlayBatch.emplace_back(wBatch[k]);
        }else{
            if(partiBatch.find(pid1) == partiBatch.end()){
                partiBatch.insert({pid1,vector<pair<pair<int, int>, pair<int, int>>>()});
            }
            partiBatch[pid1].emplace_back(wBatch[k]);
        }
    }

    vUpdated.assign(node_num, false);

    if(!partiBatch.empty()){
        if(partiBatch.size()>threadnum){
            cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
        }
        cout<<"Update Partition number: "<<partiBatch.size()<<endl;

        boost::thread_group thread;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchUpdateCheck, this, pid, boost::ref(it->second), boost::ref(overlayBatch) ));
        }
        thread.join_all();

    }

    IncreaseOverlayBatch(overlayBatch,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);

    // repair the partition index
    if(algoUpdate>=PH2H_Post){
//        Trees=TreesNo;
        Repair_PartiIndex(true,true, partiBatch);
//        Repair_PartiIndex(false,true, partiBatch);
        if(algoUpdate==PH2H_Cross){
            RefreshExtensionLabels(partiBatch);
        }
    }
}
//for PH2H
void Graph::IncreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay){
    //partition batch Increase update
    vector<pair<pair<int,int>,int>> updatedSC;

    IncreasePartiBatch(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,true);

    vector<int> Bid=BoundVertex[pid];
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis;
    for(int i=0;i<Bid.size();i++){
        bid1=Bid[i];
        for(int j=i+1;j<Bid.size();j++){
            bid2=Bid[j];
            if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                olddis=NeighborsOverlay[bid1][bid2];
            }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                continue;//exit(1);
            }

            newdis=QueryH2HPartition(bid1,bid2,pid);
            if(newdis>olddis){//if '=', not problem; if '<', problem
                sm->wait();
                weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                sm->notify();
            } else if(newdis<olddis){
                cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                exit(1);
            }
        }
    }
}
//for PCH
void Graph::IncreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC){
    //partition batch Increase update
//    vector<pair<pair<int,int>,int>> updatedSC;

    if(ifOpt){
        IncreasePartiBatchForOpt(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,false);
    }else{
        IncreasePartiBatch(pid, wBatch,NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP, updatedSC,false);
    }

//    cout<<pid<<". Size of updatedSC: "<<updatedSC.size()<<endl;
//    vector<int> Bid=BoundVertex[pid];
//    //check the boundary edge within partition
//    int bid1,bid2,olddis,newdis;
//
//    for(auto it=updatedSC.begin();it!=updatedSC.end();++it){
//        bid1=it->first.first, bid2=it->first.second; newdis=it->second;
//        if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
//            olddis=NeighborsOverlay[bid1][bid2];
//        }else{//if not found
//            cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//        }
////        cout<<pid<<": "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//        if(newdis>olddis){//if '=', not problem; if '<', problem
//            sm->wait();
//            weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
//            sm->notify();
//        } else if(newdis<olddis){
//            cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
//            exit(1);
//        }
//    }

}
//bottom-up shortcut update for post-boundary partition index
void Graph::IncreasePartiBatchUpdateCheckPostMHL(map<int, vector<pair<pair<int, int>, pair<int, int>>>> &partiBatch, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatch, bool ifParallel) {
    //partition batch Increase update
    vector<vector<pair<pair<int,int>,int>>> updatedSCs;//the overlay shortcuts that may be affected
    updatedSCs.assign(partiNum,vector<pair<pair<int,int>,int>>());
    vector<int> PropagateOverlay;

    if(ifParallel){
        if(!partiBatch.empty()){
            if(partiBatch.size()>threadnum){
                cout<<"partiBatch is larger than thread number! "<<partiBatch.size()<<" "<<threadnum<< endl;
            }
            cout<<"Update Partition number: "<<partiBatch.size()<<endl;
            boost::thread_group thread;
            for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
                int pid=it->first;
                thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchPostMHLShortcut, this, pid, boost::ref(it->second), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(updatedSCs[pid]), boost::ref(PropagateOverlay) ));
            }
            thread.join_all();
        }
    }else{
        cout<<"Single thread computation."<<endl;
        for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
            int pid=it->first;
            cout<<"Shortcut update of partition "<<pid<<endl;
            IncreasePartiBatchPostMHLShortcut(pid, it->second, Tree, rank, heightMax, updatedSCs[pid], PropagateOverlay);
        }
    }

//    cout<<"updatedSC size: "<<updatedSC.size()<<endl;
    map<pair<int,int>, pair<int,int>> updatedSCTrue;
    int updateSCSize=0;
    //check the boundary edge within partition
    int bid1,bid2,olddis,newdis,cnt, wsum;
    for(auto it1=partiBatch.begin();it1!=partiBatch.end();++it1){
        int pid=it1->first;
        for(auto it=updatedSCs[pid].begin();it!=updatedSCs[pid].end();++it){
            bid1=it->first.first, bid2=it->first.second, wsum=it->second;
            updateSCSize++;
//            cout<<bid1<<" "<<bid2<<": "<<newdis<<endl;
            if(updatedSCTrue.find(make_pair(bid1,bid2))==updatedSCTrue.end()){//if not found
                updatedSCTrue.insert({make_pair(bid1,bid2),make_pair(wsum,1)});
            }else {//if found
//            cout<<bid1<<" "<<bid2<<": "<<updatedSCTrue[make_pair(bid1,bid2)]<<" "<<newdis<<endl;
                if(updatedSCTrue[make_pair(bid1,bid2)].first!=wsum){
                    cout<<"1. Inconsistent wsum! "<<bid1<<" "<<bid2<<" "<<updatedSCTrue[make_pair(bid1,bid2)].first<<" "<<wsum<<endl; exit(1);
                }
                updatedSCTrue[make_pair(bid1,bid2)].second+=1;
            }
        }
    }

    cout<<"updatedSCTrue size: "<<updatedSCTrue.size()<<" "<<updateSCSize<< endl;
//    cout<<"Initial overlayBatch size: "<<weightOverlay.size()<<endl;
    int sVertex;
    bool ifFound;
    for(auto it=updatedSCTrue.begin();it!=updatedSCTrue.end();++it){
        bid1=it->first.first, bid2=it->first.second, wsum=it->second.first, cnt=it->second.second;
        ifFound=false;
        for(int i=0;i<Tree[rank[bid1]].vert.size();i++){
            if(Tree[rank[bid1]].vert[i].first==bid2){
                ifFound=true;
                if(Tree[rank[bid1]].vert[i].second.second < 1){//need update
                    if(cnt>1){//which means this overlay shortcut is supported by multiple partitions
//                        cout<<"cnt is larger than one! "<<bid1<<" "<<bid2<<" "<<cnt<<endl;
//                        Tree[rank[bid1]].vert[i].second.second-=(cnt-1);//
                    }
                    newdis=INF; int countwt=0;
                    for(auto it2=Neighbor[bid1].begin();it2!=Neighbor[bid1].end();++it2){
                        if(it2->first==bid2){
                            newdis=it2->second;//the weight value in the original graph
                            countwt=1;
                            break;
                        }
                    }
                    int ssw,wtt,wid;
                    vector<pair<int,int>> Wnodes;
                    Wnodes.clear();
                    if(bid1<bid2)
                        Wnodes=SCconNodesMT[bid1][bid2]; //cout<<"wid num "<<Wnodes.size()<<endl;
                    else
                        Wnodes=SCconNodesMT[bid2][bid1];
                    if(Wnodes.size()>0){
                        for(int i=0;i<Wnodes.size();i++){
                            wid=Wnodes[i].first;
                            for(int j=0;j<Tree[rank[wid]].vert.size();j++){
                                if(Tree[rank[wid]].vert[j].first==bid1){
                                    ssw=Tree[rank[wid]].vert[j].second.first;
                                }
                                if(Tree[rank[wid]].vert[j].first==bid2){
                                    wtt=Tree[rank[wid]].vert[j].second.first;
                                }
                            }

                            if(ssw+wtt<newdis){
                                newdis=ssw+wtt;
                                countwt=1;
                            }else if(ssw+wtt==newdis){
                                countwt+=1;
                            }
                        }
                    }

//                    olddis=-1;
//                    for(auto it2=NeighborsOverlayV[bid1].begin();it2!=NeighborsOverlayV[bid1].end();++it2){
//                        if(it2->first==bid2){
//                            olddis=it2->second;
//                            break;
//                        }
//                    }
//                    if(olddis==-1){
//                        cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl; exit(1);
//                    }
//                    if(olddis!=wsum){
//                        cout<<"wsum is inconsitent with olddis! "<<bid1<<" "<<bid2<<" "<<wsum<<" "<<olddis<<endl; exit(1);//？
//                    }
                    olddis=wsum;
                    if(newdis>olddis){
//            cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//            sm->wait();
                        overlayBatch.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));//weightOverlay collect the changed edges on overlay graph
//            sm->notify();
                    }else if(newdis<olddis){
                        cout<<"Something wrong happens. "<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl; exit(1);
                    }
                }
                else{
//                    cout<<"No update for this boundary shortcut! "<<bid1<<" "<<bid2<<" "<<Tree[rank[bid1]].vert[i].second.first<<" "<<Tree[rank[bid1]].vert[i].second.second<<" "<<cnt<<endl;
                }
                break;
            }
        }
        if(!ifFound){
            cout<<"Not found "<<bid2<<" in "<<bid1<<" 's vert. "<<endl ;exit(1);
        }
    }
    //get the LCA of PropagateOverlay
//    if(!PropagateOverlay.empty()){
////        cout<<"PropagateOverlay size: "<<PropagateOverlay.size()<<endl;
//        int ProID=PropagateOverlay[0];
//        int rLCA;
//        int VID;
//        for(int i=1;i<PropagateOverlay.size();++i){
//            VID=PropagateOverlay[i];
//            rLCA = LCAQueryOverlay(rank[ProID], rank[VID]);
//            ProID = Tree[rLCA].uniqueVertex;
//        }
//        ProBeginVertexSetOverlay.push_back(ProID);
//    }


}

//Function for single-edge decrease update
void Graph::DecreaseSingle(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
//            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
//            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
            Neighbor[b][i].second=newW;
            break;
        }
    }
    for(int i=0;i<node_num;++i){
        vUpdated[i]=false;
    }
    map<int, vector<pair<pair<int,int>,pair<int,int>>>> partiBatch; partiBatch.clear();
    int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
    if(pid1!=pid2){//for cut edge
        DecreaseOverlay(a, b, newW,NeighborsOverlay,Tree,rank,heightMax);
    }else{//for intra edge
        vector<pair<pair<int,int>,pair<int,int>>> tempV;
        tempV.emplace_back(make_pair(a,b), make_pair(oldW,newW));
        partiBatch.insert({pid1, tempV});
        vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
        weightOverlay.clear();
        DecreaseParti(a,b,newW,NeighborsParti,Trees[pid1],ranks[pid1],heightMaxs[pid1]);

        //weightOverlay collect the changed edges on overlay graph
        vector<int> Bid=BoundVertex[pid1];
        //check the boundary edge within partition
        int bid1,bid2,olddis,newdis;
        for(int i=0;i<Bid.size();i++){
            bid1=Bid[i];
            for(int j=i+1;j<Bid.size();j++){
                bid2=Bid[j];
                if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                    olddis=NeighborsOverlay[bid1][bid2];
                }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                    continue;//exit(1);
                }

                newdis=QueryH2HPartition(bid1,bid2,pid1);
                if(newdis<olddis){
//                    cout<<bid1<<" "<<bid2<<" "<<olddis<<" "<<newdis<<endl;
                    weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                }
            }
        }

        DecreaseOverlayBatch(weightOverlay,NeighborsOverlay,Tree,rank,heightMax,true);
        //cout<<"Overlay update number "<<weightOverlay.size()<<endl;
        //update the overlay graph index, after partition index update
        /*for(int l=0;l<weightOverlay.size();l++){
            Decrease(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay);
        }*/
    }
    // repair the partition index
    if(algoUpdate==PH2H_Post){
        Repair_PartiIndex(true, false, partiBatch);
    }

}

//Function for single-edge increase update
void Graph::IncreaseSingle(int a, int b, int oldW, int newW){
    for(int i=0;i<Neighbor[a].size();i++){
        if(Neighbor[a][i].first==b){
            Neighbor[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbor[b].size();i++){
        if(Neighbor[b][i].first==a){
            Neighbor[b][i].second=newW;
            break;
        }
    }

    map<int, vector<pair<pair<int,int>,pair<int,int>>>> partiBatch; partiBatch.clear();

    for(int i=0;i<node_num;++i){
        vUpdated[i]=false;
    }

    int pid1=PartiTag[a].first, pid2=PartiTag[b].first;
    //cout<<"increase edge in partition "<<pid<<endl;
    if(pid1!=pid2) {//for cut edge
        cout<<"Inter edge update"<<endl;
        IncreaseOverlay(a, b, oldW, newW,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid);
    }
    else {
        vector<pair<pair<int,int>,pair<int,int>>> tempV;
        tempV.emplace_back(make_pair(a,b), make_pair(oldW,newW));
        partiBatch.insert({pid1, tempV});
//        cout<<"zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"<<endl;
        vector<pair<pair<int,int>,pair<int,int>>> weightOverlay;//collect the changed edges on overlay graph
        weightOverlay.clear();

        IncreaseParti(a,b,oldW,newW,NeighborsParti,Trees[pid1],ranks[pid1],heightMaxs[pid1],SCconNodesMTP,VidtoTNidP);

        //cout<<"/////////////////////////////////////////"<<endl;

        //cout<<"boundary edge checkkkkkkkkkkkkkkkkkkkkkkkkkkk"<<endl;
        //boundary edges check
        vector<int> Bid=BoundVertex[pid1];
        int bid1,bid2,olddis,newdis;
        for(int i=0;i<Bid.size();i++){
            bid1=Bid[i];
            for(int j=i+1;j<Bid.size();j++){
                bid2=Bid[j];
                if(NeighborsOverlay[bid1].find(bid2) != NeighborsOverlay[bid1].end()){//if found
                    olddis=NeighborsOverlay[bid1][bid2];//only works for no-boundary
//                    olddis= QueryCore(bid1,bid2);
                }else{//if not found
//                    cout<<"Not found edge e("<<bid1<<","<<bid2<<") in overlay graph!"<<endl;
                    continue;//exit(1);
                }
                newdis=QueryH2HPartition(bid1,bid2,pid1);
//                newdis=ShortcutDisCheck(bid1,bid2);
//                NeighborsOverlay[bid1][bid2]=newdis;
//                NeighborsOverlay[bid2][bid1]=newdis;
//                int overlaydis=QueryCore(bid1,bid2);
                if(newdis>olddis)//if '=', not problem; if '<', problem
                    weightOverlay.emplace_back(make_pair(bid1,bid2),make_pair(olddis,newdis));
                else if(newdis<olddis){
                    cout<<"Something wrong happens. "<<bid1<<"("<<PartiTag[bid1].first<<") "<<bid2<<"("<<PartiTag[bid2].first<<") : "<<newdis<<" "<<olddis<< endl;
                    exit(1);
                }


            }
        }

        IncreaseOverlayBatch(weightOverlay,NeighborsOverlay,Tree,rank,heightMax,SCconNodesMT,VidtoTNid,true);
//        CorrectnessCheckCore(100);
        //update the overlay graph index, after partition index update
        /*cout<<"Overlay update number "<<weightOverlay.size()<<endl;
        for(int l=0;l<weightOverlay.size();l++){
            Increase(weightOverlay[l].first.first,weightOverlay[l].first.second,weightOverlay[l].second.first,weightOverlay[l].second.second,NeighborsOverlay,TreeOverlay,rankOverlay,heightMaxOverlay,SCconNodesOverlayMT,VidtoTNidOverlay);
        }*/
        //cout<<"''''''''''''''''''''''''''''''''''''''''''"<<endl;
    }
    if(algoUpdate==PH2H_Post){
//        Trees=TreesNo;
        Repair_PartiIndex(true, true, partiBatch);
    }
}



