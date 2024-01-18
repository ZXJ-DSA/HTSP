/*
 * TreeIndex.cpp
 *
 *  Created on: 16 June 2023
 *      Author: Xinjie ZHOU
 */
#include "head.h"

vector<int> NodeOrder_;//nodeID order
vector<int> _DD_;//true degree, temporal degree ,_DD2_
vector<int> _DD2_;//count

//// Index Construction
//Function of constructing tree index for partitions
void Graph::Construct_PartiIndex(bool ifParallel, bool ifLabelC){
    //for H2H update
    SCconNodesMTP.assign(node_num, map<int, vector<pair<int,int>>>());
    VidtoTNidP.assign(node_num,vector<int>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
//    DD.assign(node_num,0); //DD2.assign(node_num,0);
    Trees.assign(partiNum,vector<Node>());
    toRMQs.assign(partiNum,vector<int>());
    RMQIndexs.assign(partiNum,vector<vector<int>>());
    ranks.assign(partiNum,vector<int>());
    heightMaxs.assign(partiNum,0);
    ProBeginVertexSetParti.assign(partiNum,vector<int>());
    vertexIDChLParti.assign(partiNum,set<int>());
//    BoundEdges.assign(node_num,map<int,pair<int,int>>());

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsParti.size();i++){
        for(auto it=NeighborsParti[i].begin();it!=NeighborsParti[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    for(int pid=0;pid<partiNum;++pid){
        ranks[pid].assign(PartiVertex[pid].size(),-1);
    }

    if(ifParallel){
        cout<<"Multiple thread computation for partition index construction!"<<endl;
        //multi-thread
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
            if(ifLabelC){
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiV, this, boost::ref(processID[j]), boost::ref(Trees), boost::ref(ranks), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs), boost::ref(toRMQs), boost::ref(RMQIndexs)));
                }
            }else{
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiVCH, this, boost::ref(processID[j]), boost::ref(Trees), boost::ref(ranks), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs), boost::ref(toRMQs), boost::ref(RMQIndexs)));
                }
            }

            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                if(ifLabelC){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_Parti, this, j, boost::ref(Trees), boost::ref(ranks), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs), boost::ref(toRMQs), boost::ref(RMQIndexs)));
                }else{
                    thread.add_thread(new boost::thread(&Graph::H2HCreateTree_Parti, this, j, boost::ref(Trees[j]), boost::ref(ranks[j]), boost::ref(SCconNodesMTP), boost::ref(VidtoTNidP), boost::ref(heightMaxs)));
                }

            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            if(ifLabelC){
                ConstructPH2H_Parti(pid, Trees, ranks, SCconNodesMTP, VidtoTNidP, heightMaxs, toRMQs, RMQIndexs);
            }else{
                H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
            }

        }
    }

    vector<int> treeSize;
    int aveHeight=0;
    for(int i=0;i<partiNum;++i){
        treeSize.emplace_back(Trees[i].size());
        aveHeight+=heightMaxs[i];
    }
    cout<<"Partition graph! Maximum tree node number: "<< *max_element(treeSize.begin(),treeSize.end()) <<" ; Maximum tree height: "<< *max_element(heightMaxs.begin(),heightMaxs.end())<<" ; Average tree height: "<< aveHeight/partiNum<< endl;
}

//function of vertex allocation
void Graph::ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID){
//        processID.assign(threadnum, vector<NodeId>());
    int pid=0;
    for(int i=0;i<vertices.size();++i){
        pid=i%processID.size();
        processID[pid].emplace_back(vertices[i]);
    }
}
void Graph::ConstructBoundaryShortcut(int pid){
    //boundary edges
    int ID1,ID2,weight;
    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
            weight=QueryH2HPartition(ID1,ID2,pid);
            NeighborsOverlay[ID1][ID2]=weight;
            NeighborsOverlay[ID2][ID1]=weight;
        }
    }
}
void Graph::ConstructBoundaryShortcutNoAllPair(int pid){
    //boundary edges
    int ID1,ID2,weight;
    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        if(!NeighborCon[ID1].empty()){
            for(auto it=NeighborCon[ID1].begin();it!=NeighborCon[ID1].end();++it){
                ID2=it->first; weight=it->second.first;
                if(!PartiTag[ID2].second){//if ID2 is not boundary vertex
                    cout<<"Wrong for this shortcut! "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<")"<<endl; exit(1);//?
                }
                if(NeighborsOverlay[ID1].find(ID2)==NeighborsOverlay[ID1].end()){//if not found
                    NeighborsOverlay[ID1][ID2]=weight;
                    NeighborsOverlay[ID2][ID1]=weight;
                }else if(NeighborsOverlay[ID1][ID2]>weight){
                    NeighborsOverlay[ID1][ID2]=weight;
                    NeighborsOverlay[ID2][ID1]=weight;
                }
            }
        }
    }
}
void Graph::ConstructBoundaryShortcutV(vector<int> & p, bool ifAllPair){
    if(ifAllPair){
        for(int i=0;i<p.size();++i){
            ConstructBoundaryShortcut(p[i]);
        }
    }else{
        for(int i=0;i<p.size();++i){
            ConstructBoundaryShortcutNoAllPair(p[i]);
        }
    }

}

void Graph::Construct_OverlayGraph(bool ifParallel){
    if(ifParallel){
        //multiple threads
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
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutV, this, boost::ref(processID[j]) ,true ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcut, this, j));
            }
            thread.join_all();
        }

    }
    else{
        //single thread
        for(int k=0;k<partiNum;k++){
            ConstructBoundaryShortcut(k);
        }
    }


}

void Graph::Construct_OverlayGraphNoAllPair(bool ifParallel){
    if(ifParallel){
        //multiple threads
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
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutV, this, boost::ref(processID[j]), false ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructBoundaryShortcutNoAllPair, this, j));
            }
            thread.join_all();
        }

    }
    else{
        //single thread
        for(int k=0;k<partiNum;k++){
            ConstructBoundaryShortcutNoAllPair(k);
        }
    }


}

//Function of constructing tree index for overlay grpah
void Graph::Construct_OverlayIndex(bool ifLabelC){
    //Create tree for partition
    H2HCreateTree_Overlay();
    if(ifLabelC){
        //Create labels for partition
        H2HCreateIndex_Overlay();
    }
}

void Graph::PreConstructAllPairs(bool ifParallel){
    if(ifParallel){
        cout<<"Multiple thread computation for partition index construction!"<<endl;
        //multi-thread
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
                thread.add_thread(new boost::thread(&Graph::PreConstructAllPairsPartiV, this, boost::ref(processID[j]) ));
            }

            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::PreConstructAllPairsParti, this, j));
            }
            thread.join_all();
        }

    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            PreConstructAllPairsParti(pid);
        }
    }
}

void Graph::PreConstructAllPairsPartiV(vector<int> & p){
    for(int i=0;i<p.size();++i){
        PreConstructAllPairsParti(p[i]);
    }
}

void Graph::PreConstructAllPairsParti(int pid){
    int ID1,ID2;
    bool flagEdge;
    int newENum=0;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            ID2=BoundVertex[pid][j];
            flagEdge=false;
            for(auto it=NeighborsParti[ID1].begin();it!=NeighborsParti[ID1].end();++it){
                if(ID2==it->first){
                    flagEdge=true;
                    break;
                }
            }
            if(!flagEdge){//if not exist edge
                NeighborsParti[ID1].emplace_back(ID2,INF/3);//insert edge with weight equals INF
                NeighborsParti[ID2].emplace_back(ID1,INF/3);//insert edge with weight equals INF
                newENum+=2;
            }
        }
    }
//    cout<<"New edge number of partition "<<pid<<" : "<<newENum<<endl;
}


void Graph::ConstructPartitionPostIndexOpt(bool ifParallel){

    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){

                thread.add_thread(new boost::thread(&Graph::RefreshBoundaryEdgesAndLabelingPartiV, this, boost::ref(processID[j])));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::RefreshBoundaryEdgesAndLabelingParti, this, j));
            }
            thread.join_all();
        }
    }
    else{
        // single thread
        for(int k=0;k<partiNum;k++){
            cout<<"Repairing partition "<<k<<endl;
            RefreshBoundaryEdgesAndLabelingParti(k);
        }
    }

}

void Graph::RefreshBoundaryEdgesAndLabelingPartiV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        RefreshBoundaryEdgesAndLabelingParti(p[i]);
    }
}

void Graph::RefreshBoundaryEdgesAndLabelingParti(int pid) {
    int ID;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;

    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        NeighborsPartiPost[ID].insert(NeighborsParti[ID].begin(),NeighborsParti[ID].end());
    }

    int ID1,ID2;
    int wlocal, woverlay;
    int k;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            ID2=BoundVertex[pid][j];
            if(NeighborsPartiPost[ID1].find(ID2)==NeighborsPartiPost[ID1].end()){//if not found
                cout<<pid<<": boundary edge between "<<ID1<<" and "<<ID2<<" does not exist! "<<endl; exit(1);
            }
            wlocal=NeighborsPartiPost[ID1][ID2];
//            wlocal=-1;
//            for(k=0;k<NeighborsParti[ID1].size();++k){
//                if(NeighborsParti[ID1][k].first==ID2){
//                    wlocal=NeighborsParti[ID1][k].second;
//                    break;
//                }
//            }
            woverlay= QueryCore(ID1,ID2);
//            int d2= Dijkstra(ID1,ID2,Neighbor);
//            if(d2!=weight){
//                cout<<"Incorrect! "<<ID1<<" "<<ID2<<" "<<weight<<" "<<d2<<endl; exit(1);
//            }
            if(woverlay<wlocal){
//                BoundEdges[ID1].insert({ID2, make_pair(woverlay,-1)});
//                BoundEdges[ID2].insert({ID1, make_pair(woverlay,-1)});
//                NeighborsParti[ID1][k].second=woverlay;
                wBatch.emplace_back(make_pair(ID1,ID2), make_pair(wlocal,woverlay));
            }else if(woverlay>wlocal){
                cout<<"Something wrong happens. woverlay is larger! "<<wlocal<<" "<<woverlay<<endl; exit(1);
            }
//            else{
//                BoundEdges[ID1].insert({ID2, make_pair(woverlay,pid)});
//                BoundEdges[ID2].insert({ID1, make_pair(woverlay,pid)});
//            }

        }
    }
//    cout<<"wBatch size: "<<wBatch.size()<<endl;

//    E.assign(node_num,map<int,pair<int,int>>());
//    for(int i=0;i<PartiVertex[pid].size();i++){
//        int id=PartiVertex[pid][i];
//        for(auto it=NeighborsParti[id].begin();it!=NeighborsParti[id].end();++it){
//            E[id].insert(make_pair(it->first,make_pair(it->second,1)));
//        }
//    }
//    H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
    vector<pair<pair<int,int>,int>> updatedSC;
    //update shortcuts
    DecreasePartiBatchForOpt(pid, wBatch, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, false, true);
    //Create LCA index
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    //Create labels for partition
    H2HCreateIndex_Parti(pid, Trees[pid], ranks[pid]);
}

void Graph::ConstructPartitionPost(bool ifParallel){
    NeighborsPartiPost.assign(node_num,unordered_map<int,int>());

    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPostPartiV, this, boost::ref(processID[j])));

            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructPostParti, this, j));

            }
            thread.join_all();
        }
    }
    else{
        // single thread
        for(int k=0;k<partiNum;k++){
            cout<<"Repairing partition "<<k<<endl;
            ConstructPostParti(k);

        }
    }
}

void Graph::ConstructPostParti(int pid){
    int ID;
    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        NeighborsPartiPost[ID].insert(NeighborsParti[ID].begin(),NeighborsParti[ID].end());
    }
    int ID1,ID2,weight=-1;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();++j){
            ID2=BoundVertex[pid][j];
//            assert(NeighborsOverlay[ID1].find(ID2)!=NeighborsOverlay[ID1].end());//may not be true for no all-pair strategy
//            weight=NeighborsOverlay[ID1][ID2];
//            if(NeighborsPartiPost[ID1].find(ID2)==NeighborsPartiPost[ID1].end()){//if not found
//                NeighborsPartiPost[ID1].insert({ID2,weight});
//                NeighborsPartiPost[ID2].insert({ID1,weight});
//            }
            weight= QueryCore(ID1,ID2);
            NeighborsPartiPost[ID1][ID2]=weight;
            NeighborsPartiPost[ID2][ID1]=weight;
        }
    }
}

void Graph::ConstructPostPartiV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        ConstructPostParti(p[i]);
    }
}
//old version
void Graph::ConstructPartitionPostIndex(bool ifParallel, bool ifLabelC){
    SCconNodesMTPost.assign(node_num, map<int, vector<pair<int,int>>>());
    VidtoTNidPost.assign(node_num,vector<int>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
//    DD.assign(node_num,0); //DD2.assign(node_num,0);
    TreesPost.assign(partiNum,vector<Node>());
    toRMQsPost.assign(partiNum,vector<int>());
    RMQIndexsPost.assign(partiNum,vector<vector<int>>());
    ranksPost.assign(partiNum,vector<int>());
    heightMaxsPost.assign(partiNum,0);

    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsPartiPost.size();i++){
        for(auto it=NeighborsPartiPost[i].begin();it!=NeighborsPartiPost[i].end();++it){
            E[i].insert(make_pair(it->first,make_pair(it->second,1)));
        }
    }
    for(int pid=0;pid<partiNum;++pid){
        ranksPost[pid].assign(PartiVertex[pid].size(),-1);
    }

    if(ifParallel){
        cout<<"Multiple thread computation for partition index construction!"<<endl;
        //multi-thread
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
                if(ifLabelC){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiV, this, boost::ref(processID[j]), boost::ref(TreesPost), boost::ref(ranksPost), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost), boost::ref(toRMQsPost), boost::ref(RMQIndexsPost) ));
                }else{
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_PartiVCH, this, boost::ref(processID[j]), boost::ref(TreesPost), boost::ref(ranksPost), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost), boost::ref(toRMQsPost), boost::ref(RMQIndexsPost) ));
                }

            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                if(ifLabelC){
                    thread.add_thread(new boost::thread(&Graph::ConstructPH2H_Parti, this, j, boost::ref(TreesPost), boost::ref(ranksPost), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost), boost::ref(toRMQsPost), boost::ref(RMQIndexsPost)));
                }else{
                    thread.add_thread(new boost::thread(&Graph::H2HCreateTree_Parti, this, j, boost::ref(TreesPost[j]), boost::ref(ranksPost[j]), boost::ref(SCconNodesMTPost), boost::ref(VidtoTNidPost), boost::ref(heightMaxsPost)));
                }

            }
            thread.join_all();
        }
    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
//            cout<<"Partition "<<pid<<endl;
            if(ifLabelC){
                ConstructPH2H_Parti(pid,TreesPost, ranksPost, SCconNodesMTPost, VidtoTNidPost, heightMaxsPost, toRMQsPost, RMQIndexsPost);
            }else{
                H2HCreateTree_Parti(pid,TreesPost[pid], ranksPost[pid], SCconNodesMTPost, VidtoTNidPost, heightMaxsPost);
            }
        }
    }
    vector<int> treeSize;
    int aveHeight=0;
    for(int i=0;i<partiNum;++i){
        treeSize.emplace_back(TreesPost[i].size());
        aveHeight+=heightMaxsPost[i];
    }
    cout<<"Post-boundary partition graph! Maximum tree node number: "<< *max_element(treeSize.begin(),treeSize.end()) <<" ; Maximum tree height: "<< *max_element(heightMaxsPost.begin(),heightMaxsPost.end())<<" ; Average tree height: "<< aveHeight/partiNum<< endl;
}

//Function of constructing tree index for overlay grpah
void Graph::PostMHLIndexConstructOverlay() {
    //Create LCA index
    makeRMQCore();
    //Create labels for overlay graph
    PostMHLCreateIndex_Overlay();

}

//Function of tree-label index construction for partition
void Graph::PostMHLCreateIndex_Overlay(){
    //initialize
    vector<int> list; //list.clear();
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);
    Tree[0].dis.push_back(0);
    Tree[0].vAncestor.push_back(Tree[0].uniqueVertex);
    list.push_back(Tree[0].uniqueVertex);

    for(int i=0;i<Tree[0].ch.size();i++){
        PostMHLmakeIndexDFSOverlay(Tree[0].ch[i],list);
    }

}

void Graph::PostMHLmakeIndexDFSOverlay(int p, vector<int>& ancestors){
    int vid=Tree[p].uniqueVertex;
    if(PartiTags[vid].first==-1){//if it is overlay vertex
        //initialize
        int NeiNum=Tree[p].vert.size();
        Tree[p].pos.assign(NeiNum+1,-1);
        Tree[p].dis.assign(ancestors.size(),INF);
        Tree[p].cnt.assign(ancestors.size(),0);
        Tree[p].FN.assign(ancestors.size(),true);
        Tree[p].vAncestor.assign(ancestors.begin(),ancestors.end());
        Tree[p].vAncestor.push_back(vid);
        Tree[p].dis.push_back(0);
        //pos
        //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
        for(int i=0;i<NeiNum;i++){
            for(int j=0;j<ancestors.size();j++){
                if(Tree[p].vert[i].first==ancestors[j]){
                    Tree[p].pos[i]=j;
                    Tree[p].dis[j]=Tree[p].vert[i].second.first;
                    Tree[p].cnt[j]=1;
                    break;
                }
            }
        }
        Tree[p].pos[NeiNum]=ancestors.size();

        //dis
        for(int i=0;i<NeiNum;i++){
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;
            int k=Tree[p].pos[i];//the kth ancestor is x
            if(k==-1){
                cout<<"k="<<k<<endl; exit(1);
            }

            for(int j=0;j<ancestors.size();j++){//check the distance to the j-th ancestor could be updated by neighbors, including the valley path and peak path
                int y=ancestors[j];//the jth ancestor is y

                int z;//the distance from x to y
                if(k!=j){
                    if(k<j)//x is the ancestor of y, peak path
                        z=Tree[rank[y]].dis[k];
                    else if(k>j)//y is the ancestor of x, valley path
                        z=Tree[rank[x]].dis[j];

                    if(Tree[p].dis[j]>z+disvb){
                        Tree[p].dis[j]=z+disvb;
                        Tree[p].FN[j]=false;
                        Tree[p].cnt[j]=1;
                    }else if(Tree[p].dis[j]==z+disvb){
                        Tree[p].cnt[j]+=1;
                    }
                }
            }
        }

        //nested loop
        ancestors.push_back(Tree[p].uniqueVertex);
        for(int i=0;i<Tree[p].ch.size();i++){
            PostMHLmakeIndexDFSOverlay(Tree[p].ch[i],ancestors);
        }
        ancestors.pop_back();
    }
}

void Graph::PostMHLIndexConstructPost(bool ifParallel, double & t) {
    Timer tt;
    tt.start();
    BoundShortcuts.assign(node_num,map<int,int>());
    if(ifParallel){
        cout<<"Multiple thread computation for partition index construction!"<<endl;
        //multi-thread
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
                thread.add_thread(new boost::thread(&Graph::PostMHLPartitionIndexConstructV, this, boost::ref(processID[j]), true ));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::PostMHLPartitionIndexConstruct, this, j, true));
            }
            thread.join_all();
        }
    }
    else{
        cout<<"Single thread computation!"<<endl;
        //single-thread
        for(int pid=0;pid<partiNum;++pid){
            cout<<"Partition "<<pid<<endl;
            PostMHLPartitionIndexConstruct(pid,true);
        }
    }
    tt.stop();
    t=tt.GetRuntime();
    if(algoQuery<PH2H_Post){
        algoQuery=PH2H_Post;
    }
}

void Graph::PostMHLPartitionIndexConstructV(vector<int>& p, bool ifPost) {
    int pid;
    for(int i=0;i<p.size();++i){
        pid=p[i];
        PostMHLPartitionIndexConstruct(pid, ifPost);
    }
}

void Graph::PostMHLPartitionIndexConstruct(int pid, bool ifPost){
    vector<int> ancestors; //ancestors, the virtual root node is omited
    vector<int> interfaces;//interfaces
    int rootID=partiRoots[pid];
    int rankRoot=rank[rootID];
    ancestors.clear();
    ancestors=Tree[Tree[rankRoot].pa].vAncestor;
//    cout<<"Ancestor size of partition "<<pid<<" is "<<ancestors.size()<<endl;
    int ID1,ID2;
    /// interface distance
//    for(int i=0;i<Tree[rankRoot].vert.size();++i){
//        ID2 = Tree[rankRoot].vert[i].first;
//        interfaces.push_back(ID2);
//    }
    interfaces=BoundVertex[pid];
    if(ifPost){
        //get the all-pair distances among the interface vertices, which is equal to the shortcuts among them
        map<int,unordered_map<int,int>> disInfs;//distance from interface vertex u to another interface vertex v
        for(int i=0;i<Tree[rankRoot].vert.size();++i){
            ID1 = Tree[rankRoot].vert[i].first;
            for(int j=i+1;j<Tree[rankRoot].vert.size();++j){
                ID2 = Tree[rankRoot].vert[j].first;
                int dis = INF;// dis might be INF

                dis = QueryOverlay(ID1,ID2);//use the information from core
//                int d2= Dijkstra(ID1,ID2,Neighbor);
//                if(dis!=d2){
//                    cout<<"Wrong for boundary shortcut "<<ID1<<" "<<ID2<<" "<<dis<<" "<<d2<<endl; exit(1);
//                }
                if(dis == INF || dis == 0){
                    cout<<"Wrong distance between boundary vertices! "<<ID1<<" "<<ID2<<" "<<dis<<endl; exit(1);
                }
                if(ID1 <= ID2){
                    disInfs[ID1].insert({ID2, dis});
                    vSm[ID1]->wait();
                    BoundShortcuts[ID1][ID2]=dis;
                    vSm[ID1]->notify();
                }else{
                    disInfs[ID2].insert({ID1, dis});
                    vSm[ID2]->wait();
                    BoundShortcuts[ID2][ID1]=dis;
                    vSm[ID2]->notify();
                }
            }
        }
        PostMHLmakeIndexDFSPartiPost(rankRoot,ancestors,interfaces,disInfs);//query-orient version
    }else{
//        PostMHLmakeIndexDFSParti(rankRoot,ancestors,interfaces);
    }
}

void Graph::PostMHLmakeIndexDFSPartiPost(int p, vector<int> &ancestor, vector<int> &interface, map<int, unordered_map<int,int>> &disInfs) {
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,-1);
    Tree[p].disPost.assign(ancestor.size(),INF);
    Tree[p].cntPost.assign(ancestor.size(),0);
    Tree[p].FNPost.assign(ancestor.size(),true);
    Tree[p].FNInf.assign(interface.size(),true);
//    for(int i=0;i<interface.size();++i){
//        Tree[p].FNInf.insert({interface[i],false});//true
//    }
    Tree[p].vAncestor.assign(ancestor.begin(), ancestor.end());
    Tree[p].vAncestor.push_back(Tree[p].uniqueVertex);//the last vertex is the tree node
    Tree[p].disPost.push_back(0);
    //for interface
    int InfNum = interface.size();
    Tree[p].disInf.assign(InfNum,INF);
    int ID1=Tree[p].uniqueVertex;
    int ID2;
    int pid=PartiTags[ID1].first;

//    bool flagDebug=false;

    /// interface
    for(int j=0;j<interface.size();j++){//for each interface vertex
        ID2=interface[j];
//        Tree[p].disInf.insert({ID2,INF});
        for(int i=0;i<NeiNum;i++){//for each neighbor
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;//
            if(x == ID2){//if it is current interface vertex
                if(Tree[p].disInf[j] > disvb){
                    Tree[p].disInf[j]=disvb;//
                    Tree[p].FNInf[j]=true;//obtain from vert
                }
                continue;
            }
            int z;
            if(PartiTags[x].first != -1 ){//if x is not boundary vertex, i.e., if it is ancestor
                z = Tree[rank[x]].disInf[j];
                if(Tree[p].disInf[j]>z+disvb){
                    Tree[p].disInf[j]=z+disvb;
                    Tree[p].FNInf[j]=false;
                }
            }else{//if it is other interface vertex
                if(x<ID2){
                    if(disInfs[x].find(ID2) != disInfs[x].end()){
                        z = disInfs[x][ID2];
                    }else{
                        cout<<"Not found boundary shortcut "<<x<<"("<<PartiTags[x].first<<","<<disInfs[x].size()<<") "<<ID2<<"("<<PartiTags[ID2].first<<","<<disInfs[ID2].size()<<")"<<endl; exit(1);
                    }

                }else{
                    if(disInfs[ID2].find(x) != disInfs[ID2].end()){
                        z = disInfs[ID2][x];
                    }else{
                        cout<<"Not found boundary shortcut "<<ID2<<"("<<PartiTags[ID2].first<<","<<disInfs[ID2].size()<<") "<<x<<"("<<PartiTags[x].first<<","<<disInfs[x].size()<<")"<<endl; exit(1);
                    }
                }
                assert(z>0);
                if(Tree[p].disInf[j]>z+disvb){
                    Tree[p].disInf[j]=z+disvb;
                    Tree[p].FNInf[j]=false;
                }
            }
        }
//        int d1= Tree[p].disInf[ID2];
//        int d2= Dijkstra(ID1,ID2,Neighbor);
//        if(d1!=d2){
//            cout<<"Interface Incorrect! "<<p<<" "<<ID1<<" "<<ID2<<": "<<d1<<" "<<d2<<endl;
//        }
    }


    /// ancestor
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<ancestor.size();j++){
            if(Tree[p].vert[i].first==ancestor[j]){
                Tree[p].pos[i]=j;
                if(PartiTags[ancestor[j]].first!=-1){//if partition vertex
                    Tree[p].disPost[j]=Tree[p].vert[i].second.first;//
                    Tree[p].cntPost[j]=1;
                }
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=ancestor.size();

    bool flag= false;
    for(int j=0;j<ancestor.size();j++){
        int y=ancestor[j];//the jth ancestor is y
        int z;//the distance from x to y
//        if(p == rank[210695] && y == 208312)
//            cout<<y<<endl;
        if(PartiTags[y].first==-1){//skip the overlay vertex
//            cout<<"This ancestor is boundary vertex! "<<ID1<<" "<<y<<endl; exit(1);
            continue;
        }

//        if(ID1==116089 && y==115541){//115999
//        if(ID1==67975 && y==67432){//23022
//            cout<<ID1<<" "<<y<<": "<<j<<endl;
//            flagDebug=true;
//        }else{
//            flagDebug=false;
//        }

        for(int i=0;i<NeiNum;i++){
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;//shortcut
            int k=Tree[p].pos[i];//the kth ancestor is x
            int inf_i;
            if(PartiTags[x].first==-1){//if the neighbor is interface vertex
                flag = false;
//                for(inf_i=0;inf_i<interface.size();++inf_i){
//                    if(interface[inf_i]==x){
//                        flag=true;
//                        break;
//                    }
//                }
                if(BoundVertexMap[pid].find(x)!=BoundVertexMap[pid].end()){//if found
                    flag=true;
                    inf_i=BoundVertexMap[pid][x];
                }

                if(!flag){
                    cout<<"Not found this boundary vertex! "<<x<<endl; exit(1);
                }

                z=Tree[rank[y]].disInf[inf_i];//from the ancestor to interface
                if(Tree[p].disPost[j]>z+disvb){
                    Tree[p].disPost[j]=z+disvb;
                    Tree[p].FNPost[j]=false;
                    Tree[p].cntPost[j]=1;
//                    if(flagDebug){
//                        cout<<Tree[p].disPost[j]<<"("<<disvb<<","<<z<<") "<<Tree[p].cntPost[j]<<" "<<x<<"("<<PartiTags[x].first<<","<<k<<")"<<endl;
//                    }
                }else if(Tree[p].disPost[j]==z+disvb){
                    Tree[p].cntPost[j]+=1;//record how many path has the shortest distance
//                    if(flagDebug){
//                        cout<<Tree[p].disPost[j]<<"("<<disvb<<","<<z<<") "<<Tree[p].cntPost[j]<<" "<<x<<"("<<PartiTags[x].first<<","<<k<<")"<<endl;
//                    }
                }

            }
            else{//if the neighbor is ancestor
                assert(k!=-1);
                if(k!=j){
                    if(k<j){//if this neighbor is the ancestor of y
                        z=Tree[rank[y]].disPost[k];//get the distance to the lower-order vertex
                    }
                    else if(k>j){//if this neighbor is the descendant of y
                        z=Tree[rank[x]].disPost[j];
                    }

                    if(Tree[p].disPost[j]>z+disvb){
                        Tree[p].disPost[j]=z+disvb;
                        Tree[p].FNPost[j]=false;
                        Tree[p].cntPost[j]=1;
//                        if(flagDebug){
//                            cout<<Tree[p].disPost[j]<<"("<<disvb<<","<<z<<") "<<Tree[p].cntPost[j]<<" "<<x<<"("<<PartiTags[x].first<<","<<k<<")"<<endl;
//                        }
                    }else if(Tree[p].disPost[j]==z+disvb){
                        Tree[p].cntPost[j]+=1;
//                        if(flagDebug){
//                            cout<<Tree[p].disPost[j]<<"("<<disvb<<","<<z<<") "<<Tree[p].cntPost[j]<<" "<<x<<"("<<PartiTags[x].first<<","<<k<<")"<<endl;
//                        }
                    }
                }
            }

        }
//        int d1= Tree[p].dis[j];
//        int d2= Dijkstra(ID1,y,Neighbor);
//        if(d1!=d2){
//            cout<<"Ancestor Incorrect! "<<p<<" "<<ID1<<" "<<y<<": "<<d1<<" "<<d2<<endl;
//        }
    }
    assert(Tree[p].disInf.size()==InfNum);
    //nested loop
    ancestor.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        PostMHLmakeIndexDFSPartiPost(Tree[p].ch[i],ancestor,interface,disInfs);
    }
    ancestor.pop_back();
}


void Graph::PostMHLIndexConstructExtend(bool ifParallel, double &t){
    Timer tt;
    tt.start();
    if(ifParallel){//use multi-thread
        //multiple thread
        if(partiNum>threadnum){
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
                thread.add_thread(new boost::thread(&Graph::PostMHLPartitionIndexConstructExtendV, this, boost::ref(processID[j]) ));
            }
            thread.join_all();
        }else{
            boost::thread_group threadf;
            for(int pid=0;pid<partiNum;++pid) {
                threadf.add_thread(new boost::thread(&Graph::PostMHLPartitionIndexConstructExtend, this, pid));
            }
            threadf.join_all();
        }
    }
    else{//single thread
        for(int i=0;i<partiNum;++i){
            PostMHLPartitionIndexConstructExtend(i);
        }
    }
    tt.stop();
    t=tt.GetRuntime();
    algoQuery=PH2H_Cross;
}

void Graph::PostMHLPartitionIndexConstructExtendV(vector<int>& p) {
    int pid;
    for(int i=0;i<p.size();++i){
        pid=p[i];
        PostMHLPartitionIndexConstructExtend(pid);
    }
}

void Graph::PostMHLPartitionIndexConstructExtend(int pid){
    vector<int> list; //list.clear();
    int ID=partiRoots[pid];
    int rootRank=rank[ID];
    int pa=Tree[rootRank].pa;
    list=Tree[pa].vAncestor;
//    list.push_back(Tree[rootRank].uniqueVertex);

    PostMHLmakeIndexDFSPartiExtend(rootRank,list);
}

void Graph::PostMHLmakeIndexDFSPartiExtend(int p, vector<int>& ancestors){

    //initialize
    int NeiNum=Tree[p].vert.size();
//    Tree[p].pos.assign(NeiNum+1,0);
    Tree[p].dis.assign(ancestors.size(),INF);
    Tree[p].cnt.assign(ancestors.size(),0);
    Tree[p].FN.assign(ancestors.size(),true);
//    Tree[p].vAncestor.assign(ancestors.begin(), ancestors.end());
//    Tree[p].vAncestor.push_back(Tree[p].uniqueVertex);//the last vertex is the tree node
    Tree[p].dis.push_back(0);
    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<ancestors.size();j++){
            if(Tree[p].vert[i].first==ancestors[j]){
//                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
//    Tree[p].pos[NeiNum]=ancestors.size();

    //dis
    for(int i=0;i<NeiNum;i++){
        int x=Tree[p].vert[i].first;
        int disvb=Tree[p].vert[i].second.first;
        int k=Tree[p].pos[i];//the kth ancestor is x

        for(int j=0;j<ancestors.size();j++){//check the distance to the j-th ancestor could be updated by neighbors, including the valley path and peak path
            int y=ancestors[j];//the jth ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//x is the ancestor of y, peak path
                    z=Tree[rank[y]].dis[k];
                else if(k>j)//y is the ancestor of x, valley path
                    z=Tree[rank[x]].dis[j];

                if(Tree[p].dis[j]>z+disvb){
                    Tree[p].dis[j]=z+disvb;
                    Tree[p].FN[j]=false;
                    Tree[p].cnt[j]=1;
                }else if(Tree[p].dis[j]==z+disvb){
                    Tree[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    ancestors.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        PostMHLmakeIndexDFSPartiExtend(Tree[p].ch[i],ancestors);
    }
    ancestors.pop_back();

}

//Function of repair the partition index
void Graph::Repair_PartiIndex(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
//    repairShortcuts.assign(node_num, unordered_map<vertex,pair<int,int>>());
    ifRepaired.assign(partiNum, false);
    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexV, this, boost::ref(processID[j]), ifIncrease, boost::ref(partiBatch), false));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            if(ifIncrease){//increase update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexIncrease, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }
            else{//decrease update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexDecrease, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }

        }
    }
    else{
        // single thread
        if(ifIncrease){
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexIncrease(k,partiBatch);
            }
        }
        else{
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexDecrease(k,partiBatch);
            }
        }

    }

//    int pNum=0;
//    for(auto it=ifRepaired.begin();it!=ifRepaired.end();++it){
//        if(*it){
//            ++pNum;
//        }
//    }
//    cout<<"Repaired partition number: "<<pNum<<endl;

}

//Function of repair the partition index
void Graph::Repair_PartiIndexPostMHLPost(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, double & runT){
//    repairShortcuts.assign(node_num, unordered_map<vertex,pair<int,int>>());
    Timer tt,tt2;
    tt2.start();
    tt.start();
    BoundShortcutsCheck(true, ifIncrease);// check which partition is affected and update the boundary shortcuts
//    BoundShortcutsCheck(false, ifIncrease);
    tt.stop();
//    int pid;
//    for(auto it=partiBatch.begin();it!=partiBatch.end();++it){
//        pid=it->first;
//        if(affectedParti.find(pid)==affectedParti.end()){
//            affectedParti.insert(pid);
//        }
//    }
    cout<<"All-pair boundary check time: "<<tt.GetRuntime()<<" s. affectedParti size: "<<affectedParti.size()<<endl;
    tt.start();
    ifRepaired.assign(partiNum, false);
    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        vector<int> vertices;
        for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
            vertices.emplace_back(*it);
        }
        if(threadnum<vertices.size()){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());

            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::PostMHLIndexUpdatePostPartV, this, boost::ref(processID[j]), ifIncrease));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                thread.add_thread(new boost::thread(&Graph::PostMHLIndexUpdatePostPart, this, *it, ifIncrease));
            }
            thread.join_all();
        }
    }
    else{
        for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
//            if(*it != 0){
//                continue;
//            }
//            cout<<"Repairing partition "<<*it<<endl;
//                AncestorEntryDecreaseUpdateParti(*it);
            PostMHLIndexUpdatePostPart(*it, ifIncrease);
        }

    }
    tt.stop();
    cout<<"Top-down label update time: "<<tt.GetRuntime()<<" s."<<endl;
    tt2.stop();
    runT=tt2.GetRuntime();
    if(algoQuery<PH2H_Post){
        algoQuery = PH2H_Post;
    }

}

void Graph::PostMHLIndexUpdatePostPartV(vector<int>& p, bool ifIncrease){
    int pid;
    for(int i=0;i<p.size();++i){
        pid = p[i];
        PostMHLIndexUpdatePostPart(pid, ifIncrease);
    }
}

void Graph::PostMHLIndexUpdatePostPart(int pid, bool ifIncrease){
    map<int,int> checkedDis;//map<tree node ID, distance index>
    vector<int> interface;

    int rootVertex=partiRoots[pid];
    int rootRank=rank[rootVertex];
    int rootHeight=Tree[rootRank].height-1;

    interface=BoundVertex[pid];

    vector<int> ancestors; //linee.clear();
    ancestors=Tree[Tree[rootRank].pa].vAncestor;
//    cout<<"Post parti update: "<<pid<<" ; root vertex: "<<rootVertex<< endl;
    Timer tt;
//    PostMHLIndexUpdatePostPartDFS(pid, rootRank, ancestors, interface, BoundShortcuts, rootHeight , ifIncrease);
    PostMHLIndexUpdatePostPartDFS(pid, rootRank, ancestors, interface, BoundShortcuts, vertexIDChLParti[pid],  rootHeight , ifIncrease, tt);
//    PostMHLIndexUpdatePostPartDFS(pid, rootRank, ancestors, interface, BoundShortcuts, vertexIDChLOverlay,  rootHeight , ifIncrease);

}
//correct version
void Graph::PostMHLIndexUpdatePostPartDFS(int pid, int p, vector<int>& ancestor, vector<int> & interface, vector<map<int,int>>& disInfs, int rootHeight, bool ifIncrease){
    int ID1,ID2;
    int NeiNum=Tree[p].vert.size();
    ID1=Tree[p].uniqueVertex;
    bool flagDebug=false;
//    flagDebug=true;
//    if(ID1==141527){
//        cout<<pid<<" "<<ID1<<endl;
////        flagDebug=true;
//    }
    int finalhub,dh1,dh2;
    /// interface
    for(int j=0;j<interface.size();j++){//for each interface vertex
        ID2=interface[j];
//        if(ID2==141623){
//            cout<<ID1<<" "<<ID2<<" "<<Tree[p].disInf[j]<< endl;
//            flagDebug=true;
//        }
        if(ifIncrease){
            Tree[p].disInf[j] = INF;
        }
//        Tree[p].disInf.insert({ID2,INF});
        for(int i=0;i<NeiNum;i++){//for each neighbor
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;//

//            if(flagDebug && ID1==141527){
//                cout<<"Neighbor "<<x<<"("<<PartiTags[x].first<<")"<<endl;
//            }


            if(x == ID2){//if it is current interface vertex
                if(Tree[p].disInf[j] > disvb){
                    Tree[p].disInf[j]=disvb;//
                    Tree[p].FNInf[j]=true;//obtain from vert
                    finalhub=x;
                    dh1=disvb, dh2=0;
                }
                continue;
            }
            int z=INF;
            if(PartiTags[x].first != -1 ){//if x is not boundary vertex, i.e., if it is ancestor
                z = Tree[rank[x]].disInf[j];
//                if(flagDebug && ID1==141527){
//                    int d2= Dijkstra(x,ID2,Neighbor);
//                    if(z!=d2){
//                        cout<<"Wrong for non-boundary "<<x<<" "<<z<<" "<<d2<<endl; exit(1);
//                    }
//                }
                if(Tree[p].disInf[j]>z+disvb){
                    Tree[p].disInf[j]=z+disvb;
                    Tree[p].FNInf[j]=false;
                    finalhub=x;
                    dh1=disvb, dh2=z;
                }
            }else{//if it is other interface vertex
                if(x<ID2){
                    if(disInfs[x].find(ID2) != disInfs[x].end()){
                        z = disInfs[x][ID2];
                    }else{
                        cout<<"Not found boundary shortcut "<<x<<"("<<PartiTags[x].first<<","<<disInfs[x].size()<<") "<<ID2<<"("<<PartiTags[ID2].first<<","<<disInfs[ID2].size()<<")"<<endl; exit(1);
                    }

                }else{
                    if(disInfs[ID2].find(x) != disInfs[ID2].end()){
                        z = disInfs[ID2][x];
                    }else{
                        cout<<"Not found boundary shortcut "<<ID2<<"("<<PartiTags[ID2].first<<","<<disInfs[ID2].size()<<") "<<x<<"("<<PartiTags[x].first<<","<<disInfs[x].size()<<")"<<endl; exit(1);
                    }
                }
                assert(z>0);
//                if(flagDebug && ID1==141527){
//                    int d2= Dijkstra(x,ID2,Neighbor);
//                    if(z!=d2){
//                        cout<<"Wrong for boundary "<<x<<" "<<z<<" "<<d2<< endl; exit(1);
//                    }
//                }
                if(Tree[p].disInf[j]>z+disvb){
                    Tree[p].disInf[j]=z+disvb;
                    Tree[p].FNInf[j]=false;
                    finalhub=x;
                    dh1=disvb, dh2=z;
                }
            }
        }
//        if(flagDebug){
//            int d1= Tree[p].disInf[j];
//            int d2= Dijkstra(ID1,ID2,Neighbor);
//            if(d1!=d2){
//                cout<<"Interface Incorrect! "<<pid<<". "<<ID1<<" "<<ID2<<": "<<d1<<" "<<d2<<endl;
//                cout<<ID1<<" "<<finalhub<<" "<<ID2<<": "<<dh1<<" "<<dh2<<" "<<d1<<" ; "<<Dijkstra(ID1,finalhub,Neighbor)<<" "<<Dijkstra(finalhub,ID2,Neighbor)<<" "<<d2<<endl;
//                exit(1);
//            }
//        }

    }

    /// ancestor
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    bool flag= false;

    if(!Tree[p].DisRe.empty()) {//Tree[child].DisRe.size()!=0
        for(int j=rootHeight;j<ancestor.size();j++){
            int y=ancestor[j];//the jth ancestor is y
            int z;//the distance from x to y
//        if(p == rank[210695] && y == 208312)
//            cout<<y<<endl;
            if(PartiTags[y].first==-1){//skip the overlay vertex
                cout<<"This ancestor is boundary vertex! "<<ID1<<" "<<y<<endl; exit(1);
//                continue;
            }

            for(int i=0;i<NeiNum;i++){
                int x=Tree[p].vert[i].first;
                int disvb=Tree[p].vert[i].second.first;//shortcut
                int k=Tree[p].pos[i];//the kth ancestor is x
                int inf_i;
                if(PartiTags[x].first==-1){//if the neighbor is interface vertex
                    flag = false;
//                for(inf_i=0;inf_i<interface.size();++inf_i){
//                    if(interface[inf_i]==x){
//                        flag=true;
//                        break;
//                    }
//                }
                    if(BoundVertexMap[pid].find(x)!=BoundVertexMap[pid].end()){//if found
                        flag=true;
                        inf_i=BoundVertexMap[pid][x];
                    }

                    if(!flag){
                        cout<<"Not found this boundary vertex! "<<x<<endl; exit(1);
                    }

                    z=Tree[rank[y]].disInf[inf_i];//from the ancestor to interface
                    if(Tree[p].disPost[j]>z+disvb){
                        Tree[p].disPost[j]=z+disvb;
                        Tree[p].FNPost[j]=false;
                        Tree[p].cntPost[j]=1;
                    }else if(Tree[p].disPost[j]==z+disvb){
                        Tree[p].cntPost[j]+=1;//record how many path has the shortest distance
                    }

                }
                else{//if the neighbor is ancestor
                    assert(k!=-1);
                    if(k!=j){
                        if(k<j){//if this neighbor is the ancestor of y
                            z=Tree[rank[y]].disPost[k];//get the distance to the lower-order vertex
                        }
                        else if(k>j){//if this neighbor is the descendant of y
                            z=Tree[rank[x]].disPost[j];
                        }

                        if(Tree[p].disPost[j]>z+disvb){
                            Tree[p].disPost[j]=z+disvb;
                            Tree[p].FNPost[j]=false;
                            Tree[p].cntPost[j]=1;
                        }else if(Tree[p].disPost[j]==z+disvb){
                            Tree[p].cntPost[j]+=1;
                        }
                    }
                }

            }
//        int d1= Tree[p].dis[j];
//        int d2= Dijkstra(ID1,y,Neighbor);
//        if(d1!=d2){
//            cout<<"Ancestor Incorrect! "<<p<<" "<<ID1<<" "<<y<<": "<<d1<<" "<<d2<<endl;
//        }
        }
    }
    else{
        for(int j=rootHeight;j<ancestor.size();j++){
            int y=ancestor[j];//the jth ancestor is y
            int z;//the distance from x to y
//        if(p == rank[210695] && y == 208312)
//            cout<<y<<endl;
            if(PartiTags[y].first==-1){//skip the overlay vertex
                cout<<"This ancestor is boundary vertex! "<<ID1<<" "<<y<<endl; exit(1);
//                continue;
            }

            for(int i=0;i<NeiNum;i++){
                int x=Tree[p].vert[i].first;
                int disvb=Tree[p].vert[i].second.first;//shortcut
                int k=Tree[p].pos[i];//the kth ancestor is x
                int inf_i;
                if(PartiTags[x].first==-1){//if the neighbor is interface vertex
                    flag = false;
//                for(inf_i=0;inf_i<interface.size();++inf_i){
//                    if(interface[inf_i]==x){
//                        flag=true;
//                        break;
//                    }
//                }
                    if(BoundVertexMap[pid].find(x)!=BoundVertexMap[pid].end()){//if found
                        flag=true;
                        inf_i=BoundVertexMap[pid][x];
                    }

                    if(!flag){
                        cout<<"Not found this boundary vertex! "<<x<<endl; exit(1);
                    }

                    z=Tree[rank[y]].disInf[inf_i];//from the ancestor to interface
                    if(Tree[p].disPost[j]>z+disvb){
                        Tree[p].disPost[j]=z+disvb;
                        Tree[p].FNPost[j]=false;
                        Tree[p].cntPost[j]=1;
                    }else if(Tree[p].disPost[j]==z+disvb){
                        Tree[p].cntPost[j]+=1;//record how many path has the shortest distance
                    }

                }
                else{//if the neighbor is ancestor
                    assert(k!=-1);
                    if(k!=j){
                        if(k<j){//if this neighbor is the ancestor of y
                            z=Tree[rank[y]].disPost[k];//get the distance to the lower-order vertex
                        }
                        else if(k>j){//if this neighbor is the descendant of y
                            z=Tree[rank[x]].disPost[j];
                        }

                        if(Tree[p].disPost[j]>z+disvb){
                            Tree[p].disPost[j]=z+disvb;
                            Tree[p].FNPost[j]=false;
                            Tree[p].cntPost[j]=1;
                        }else if(Tree[p].disPost[j]==z+disvb){
                            Tree[p].cntPost[j]+=1;
                        }
                    }
                }

            }
//        int d1= Tree[p].dis[j];
//        int d2= Dijkstra(ID1,y,Neighbor);
//        if(d1!=d2){
//            cout<<"Ancestor Incorrect! "<<p<<" "<<ID1<<" "<<y<<": "<<d1<<" "<<d2<<endl;
//        }
        }
    }

    //nested loop
    ancestor.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        PostMHLIndexUpdatePostPartDFS(pid, Tree[p].ch[i],ancestor,interface,disInfs, rootHeight, ifIncrease);
    }
    ancestor.pop_back();
}
//new version
void Graph::PostMHLIndexUpdatePostPartDFS(int pid, int p, vector<int>& ancestor, vector<int> & interface, vector<map<int,int>>& disInfs, set<int>& vertexIDChL, int rootHeight, bool ifIncrease, Timer& tt){
    int ID1,ID2;
    int NeiNum=Tree[p].vert.size();
    ID1=Tree[p].uniqueVertex;
    int childH=Tree[p].height-1;
    bool flagDebug=false;
//    flagDebug=true;
//    if(ID1==115541){
//        cout<<pid<<" "<<ID1<<" "<<Tree[p].DisRePost.size()<<" "<<ancestor.size()<<" "<<rootHeight<< " "<<Tree[p].vert.size()<< endl;
//        flagDebug=true;
//    }
    tt.start();
    int finalhub,dh1,dh2;
    /// interface
    for(int j=0;j<interface.size();j++){//for each interface vertex
        ID2=interface[j];
//        int boundaryH=Tree[rank[ID2]].height-1;
//        if(ID2==141623){
//            cout<<ID1<<" "<<ID2<<" "<<Tree[p].disInf[j]<< endl;
//            flagDebug=true;
//        }
//        if(ifIncrease){
//            Tree[p].disInf[j] = INF;
//        }
//        Tree[p].disInf.insert({ID2,INF});
        int dis=INF; bool FnInf=false;
        if(!ifIncrease){
            dis=Tree[p].disInf[j];
        }
        for(int i=0;i<NeiNum;i++){//for each neighbor
            int x=Tree[p].vert[i].first;
            int disvb=Tree[p].vert[i].second.first;//

//            if(flagDebug && ID2==116005){
//                cout<<"Neighbor "<<x<<"("<<PartiTags[x].first<<") "<<Tree[p].disInf[j]<<" "<<dis<<endl;
//            }


            if(x == ID2){//if it is current interface vertex
                if(dis > disvb){//Tree[p].disInf[j]
                    dis=disvb;//
                    FnInf=true;//obtain from vert. Tree[p].FNInf[j]
                    finalhub=x;
                    dh1=disvb, dh2=0;
                }
                continue;
            }
            int z=INF;
            if(PartiTags[x].first != -1 ){//if x is not boundary vertex, i.e., if it is ancestor
                z = Tree[rank[x]].disInf[j];
//                if(flagDebug && ID1==141527){
//                    int d2= Dijkstra(x,ID2,Neighbor);
//                    if(z!=d2){
//                        cout<<"Wrong for non-boundary "<<x<<" "<<z<<" "<<d2<<endl; exit(1);
//                    }
//                }
                if(dis>z+disvb){//Tree[p].disInf[j]
                    dis=z+disvb;
                    FnInf=false;//Tree[p].FNInf[j]
                    finalhub=x;
                    dh1=disvb, dh2=z;
                }
            }
            else{//if it is other interface vertex
                if(x<ID2){
                    if(disInfs[x].find(ID2) != disInfs[x].end()){
                        z = disInfs[x][ID2];
                    }else{
                        cout<<"Not found boundary shortcut "<<x<<"("<<PartiTags[x].first<<","<<disInfs[x].size()<<") "<<ID2<<"("<<PartiTags[ID2].first<<","<<disInfs[ID2].size()<<")"<<endl; exit(1);
                    }

                }else{
                    if(disInfs[ID2].find(x) != disInfs[ID2].end()){
                        z = disInfs[ID2][x];
                    }else{
                        cout<<"Not found boundary shortcut "<<ID2<<"("<<PartiTags[ID2].first<<","<<disInfs[ID2].size()<<") "<<x<<"("<<PartiTags[x].first<<","<<disInfs[x].size()<<")"<<endl; exit(1);
                    }
                }
//                assert(z>0);
                if(z<=0){
                    cout<<"Wrong! "<<pid<<" "<<ID2<<" "<<x<<" "<<z<<" "<<Dijkstra(x,ID2,Neighbor) <<endl; exit(1);
                }
//                if(flagDebug && ID1==141527){
//                    int d2= Dijkstra(x,ID2,Neighbor);
//                    if(z!=d2){
//                        cout<<"Wrong for boundary "<<x<<" "<<z<<" "<<d2<< endl; exit(1);
//                    }
//                }
                if(dis>z+disvb){//Tree[p].disInf[j]
                    dis=z+disvb;//Tree[p].disInf[j]
                    FnInf=false;//Tree[p].FNInf[j]
                    finalhub=x;
                    dh1=disvb, dh2=z;
                }
            }
        }

        if(ifIncrease){
//            if(flagDebug && ID2==116005){
//                cout<<ID1<<" "<<ID2<<" "<<Tree[p].disInf[j]<<" "<<dis<<endl;
//            }
            if(dis>Tree[p].disInf[j]){
                int PID;
                int disBF=Tree[p].disInf[j];
                //identify the affected post-boundary labels
//            for(int k=0;k<VidtoTNid[ID1].size();k++){//for the tree node that contains childID as vert element
//                PID=VidtoTNid[ID1][k];
//                if(Tree[PID].FNPost[childH] && Tree[PID].disPost[boundaryH]==Tree[p].disInf[j]+Tree[PID].disPost[childH]){//if label is from shortcut
//                    Tree[PID].cntPost[boundaryH]-=1;
//                    cout<<"enter 1."<<endl;
//                }
//            }

                //line[i]
                for(int k=0;k<VidtoTNid[ID2].size();k++){
                    PID=VidtoTNid[ID2][k];
//                    if(flagDebug && ID2==116005){
////                        cout<<VidtoTNid[ID2].size()<<",k: "<<Tree[PID].uniqueVertex<<"("<<PartiTags[Tree[PID].uniqueVertex].first<<")"<<endl;
//                        if(Tree[PID].uniqueVertex==116033){
//                            cout<<"Find! "<<Tree[PID].uniqueVertex<<": "<<Tree[PID].disPost[childH]<<" "<<disBF<<" "<<Tree[PID].disInf[j]<<" "<<Tree[PID].cntPost[childH]<<endl;
//                        }
//                    }
                    if(PartiTags[Tree[PID].uniqueVertex].first==-1){//if PID is overlay tree node
                        continue;
                    }
//                if(Tree[PID].height>Tree[children].height){///modified for correctness, PID may not be the descendant of children
                    if(Tree[PID].height>Tree[p].height && Tree[PID].vAncestor[childH] == ID1){///modified for correctness, PID may not be the descendant of children
                        if(Tree[PID].FNInf[j] && Tree[PID].disPost[childH]==disBF+Tree[PID].disInf[j]){// Tree[PID].FNPost[boundaryH] &&
                            Tree[PID].cntPost[childH]-=1;
//                        cout<<"enter 2. "<<Tree[PID].uniqueVertex<<" "<<ID1<<" "<<Tree[PID].cntPost[childH]<< endl;

//                            if(Tree[PID].uniqueVertex==67975 && ID1==67432){
//                                cout<<"2 Find!! "<<Tree[PID].uniqueVertex<<" "<<ID1<<" "<<Tree[PID].cntPost[childH]<<endl;
//                            }
                        }
                    }
                }

                Tree[p].disInf[j]=dis;
                Tree[p].FNInf[j]=FnInf;
            }
            else if(dis<Tree[p].disInf[j]){
                cout<<"Wrong dis for disInf. "<<ID1<<" "<<ID2<<" "<<Tree[p].disInf[j]<<" "<<dis<<endl; exit(1);
            }
        }
        else{
            if(dis<Tree[p].disInf[j]){
                Tree[p].disInf[j]=dis;
                Tree[p].FNInf[j]=FnInf;
            }else if(dis>Tree[p].disInf[j]){
                cout<<"Wrong dis for disInf. "<<ID1<<" "<<ID2<<" "<<Tree[p].disInf[j]<<" "<<dis<<endl; exit(1);
            }
        }

//        if(flagDebug){
//            int d1= Tree[p].disInf[j];
//            int d2= Dijkstra(ID1,ID2,Neighbor);
//            if(d1!=d2){
//                cout<<"Interface Incorrect! "<<pid<<". "<<ID1<<" "<<ID2<<": "<<d1<<" "<<d2<<endl;
//                cout<<ID1<<" "<<finalhub<<" "<<ID2<<": "<<dh1<<" "<<dh2<<" "<<d1<<" ; "<<Dijkstra(ID1,finalhub,Neighbor)<<" "<<Dijkstra(finalhub,ID2,Neighbor)<<" "<<d2<<endl;
//                exit(1);
//            }
//        }

    }
    tt.stop();
//    tBoundary+=tt.GetRuntime();
    /// ancestor
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    bool ProIDdisCha=false;
//    cout<<pid<<": "<<ID<<" "<<PartiTags[ID].first<<" "<<Tree[rank[ID]].height-1<<endl;

    int inf_i;
    tt.start();
    /// update ancestor entry
    if(ifIncrease){
//        int childID=Tree[p].uniqueVertex;
        for(int i=0;i<Tree[p].disPost.size()-1;i++){
            ID2=Tree[p].vAncestor[i];
            if(PartiTags[ID2].first==-1){
                continue;
            }
//            if(flagDebug && ID2==68428 ){
//                cout<<ID1<<" "<<ID2<<" "<<Tree[p].disPost[i]<<endl;
//                for(int j=0;j<Tree[p].vert.size();j++){
//                    if(Tree[p].vert[j].first==ID2){
//                        cout<<ID2<<" is in the vert of "<<ID1<<endl;
//                    }
//                }
//            }

            if(Tree[p].cntPost[i]<=0){//if the distance label to i-th ancestor should be maintained, cntPost may not be correctly maintained, some cntPost may be -1 due to twice update triggers (shortcut and boundary array)
//            if(true){//correct
//                changelabel+=1;
                //firstly, check which dis can be infected
                int disBF=Tree[p].disPost[i];
                int PID;
                //chidlID
                for(int k=0;k<VidtoTNid[ID1].size();k++){//for the tree node that contains childID as vert element
                    PID=VidtoTNid[ID1][k];
                    if(Tree[PID].FNPost[childH] && Tree[PID].disPost[i]==disBF+Tree[PID].disPost[childH]){//if label is from shortcut
                        Tree[PID].cntPost[i]-=1;
//                        cout<<"enter 3. "<<Tree[PID].uniqueVertex<<" "<<ID2<<" "<<Tree[PID].cntPost[i]<< endl;
//                        if(Tree[PID].uniqueVertex==116033 && ID2==115541){
//                            cout<<"3 Find!! "<<Tree[PID].uniqueVertex<<" "<<ID2<<" "<<Tree[PID].cntPost[i]<<endl;
//                        }
                    }
                }

                //line[i]
                for(int k=0;k<VidtoTNid[ancestor[i]].size();k++){
                    PID=VidtoTNid[ancestor[i]][k];
//                if(Tree[PID].height>Tree[children].height){///modified for correctness, PID may not be the descendant of children
                    if(Tree[PID].height>Tree[p].height && Tree[PID].vAncestor[childH] == ID1){///modified for correctness, PID may not be the descendant of children
//                        if(PID>Tree.size()){
//                            cout<<"PID error! "<<PID<<" "<<Tree.size()<<endl; exit(1);
//                        }
//                        if(childH>Tree[PID].disPost.size()){
//                            cout<<"childH error! "<<childH<<" "<<Tree[PID].disPost.size()<<": "<<ID1<<"("<<Tree[p].height<<") "<<PID<<"("<<Tree[PID].height<<")"<<endl; exit(1);
//                        }
                        if(Tree[PID].FNPost[i] && Tree[PID].disPost[childH]==disBF+Tree[PID].disPost[i]){///
                            Tree[PID].cntPost[childH]-=1;
//                            cout<<"enter 4. "<<Tree[PID].uniqueVertex<<" "<<ID1<<" "<<Tree[PID].cntPost[childH]<< endl;
//                            if(Tree[PID].uniqueVertex==116033 && ID1==115541){
//                                cout<<"4 Find!! "<<Tree[PID].uniqueVertex<<" "<<ID1<<" "<<Tree[PID].cntPost[childH]<<endl;
//                            }
                        }
                    }
                }

                //secondly, calculate the actual distance
                int dis=INF; int count=0;
                int Dvb; int b,bH; int DDvb=INF;
                for(int j=0;j<Tree[p].vert.size();j++){
                    Dvb=Tree[p].vert[j].second.first;
                    b=Tree[p].vert[j].first;
                    bH=Tree[rank[b]].height-1;
                    if(PartiTags[b].first==-1){//if boundary vertex
                        assert(BoundVertexMap[pid].find(b)!=BoundVertexMap[pid].end());
                        inf_i=BoundVertexMap[pid][b];
//                        if(flagDebug && ID2==204850){
//                            int d2= Dijkstra(ID2,b,Neighbor);
//                            if(d2!=Tree[rank[ancestor[i]]].disInf[inf_i]){
//                                cout<<"Incorrect! "<<ID1<<" "<<ID2<<" "<<b<<" "<<Tree[rank[ancestor[i]]].disInf[inf_i]<<" "<<d2<<endl; exit(1);
//                            }
//                        }
                        if(Dvb+Tree[rank[ancestor[i]]].disInf[inf_i]<dis){
//                            if(flagDebug && ID2==204850){
//                                cout<<ID1<<" "<<ID2<<": "<<b<<" "<<dis<<" "<<Dvb+Tree[rank[ancestor[i]]].disInf[inf_i]<<endl;
//                            }
                            dis=Dvb+Tree[rank[ancestor[i]]].disInf[inf_i];
                            count=1;
                        }else if(Dvb+Tree[rank[ancestor[i]]].disInf[inf_i]==dis){
                            count+=1;
                        }
                    }else{//if partition vertex
                        if(bH<i){
//                            if(flagDebug && ID2==204850){
//                                int d2= Dijkstra(ID2,b,Neighbor);
//                                if(d2!=Tree[rank[ancestor[i]]].disPost[bH]){
//                                    cout<<"Incorrect! "<<ID1<<" "<<ID2<<" "<<b<<" "<<Tree[rank[ancestor[i]]].disPost[bH]<<" "<<d2<<endl; exit(1);
//                                }
//                            }
                            if(Dvb+Tree[rank[ancestor[i]]].disPost[bH]<dis){
//                                if(flagDebug && ID2==204850){
//                                    cout<<ID1<<" "<<ID2<<": "<<b<<" "<<dis<<" "<<Dvb+Tree[rank[ancestor[i]]].disPost[bH]<<endl;
//                                }
                                dis=Dvb+Tree[rank[ancestor[i]]].disPost[bH];
                                count=1;
                            }else if(Dvb+Tree[rank[ancestor[i]]].disPost[bH]==dis){
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
//                            if(flagDebug && ID2==204850){
//                                int d2= Dijkstra(ID2,b,Neighbor);
//                                if(d2!=Tree[rank[b]].disPost[i]){
//                                    cout<<"Incorrect! "<<ID1<<" "<<ID2<<" "<<b<<" "<<Tree[rank[b]].disPost[i]<<" "<<d2<<endl; exit(1);
//                                }
//                            }
                            if(Dvb+Tree[rank[b]].disPost[i]<dis){
//                                if(flagDebug && ID2==204850){
//                                    cout<<ID1<<" "<<ID2<<": "<<b<<" "<<dis<<" "<<Dvb+Tree[rank[b]].disPost[i]<<endl;
//                                }
                                dis=Dvb+Tree[rank[b]].disPost[i];
                                count=1;
                            }else if(Dvb+Tree[rank[b]].disPost[i]==dis){
                                count+=1;
                            }
                        }
                    }

                }
                if(DDvb==dis) {
                    Tree[p].FNPost[i]=true;
                }
                if(Tree[p].disPost[i]<dis){
                    Tree[p].disPost[i]=dis;
                    Tree[p].cntPost[i]=count;
//                    ProIDdisCha=true;
                }
                else{
//                    cout<<"Ineffective update. "<<ID1<<" "<<ancestor[i]<<" "<<Tree[p].disPost[i]<<" "<<dis<<endl; exit(1);
//                    Tree[p].disPost[i]=dis;
//                    Tree[p].cntPost[i]=count;
                }

            }
            else if(Tree[p].cntPost[i]<0){
                cout<<"cntPost less than 0! "<<ID1<<"("<<Tree[p].height<<") "<<ID2<<" "<<Tree[p].cntPost[i]<<endl;
            }
        }
    }
    else{//if decrease
        if(!Tree[p].DisRePost.empty()){//Tree[child].DisRe.size()!=0
            /// update ancestor entries
            for(int k=0;k<Tree[p].vert.size();k++){
                int b=Tree[p].vert[k].first;
                if(PartiTags[b].first == -1){//if it is interface vertex, check whether the distance label may be affected by interface
                    int vbW=Tree[p].vert[k].second.first;
                    assert(BoundVertexMap[pid].find(b)!=BoundVertexMap[pid].end());
                    inf_i=BoundVertexMap[pid][b];
                    for(int i=rootHeight;i<ancestor.size();i++){
//                    if(Tree[rank[line[i]]].disInf.find(b) == Tree[rank[line[i]]].disInf.end()){
//                        cout<<"DisRePost. Not found this boundary vertex! "<<i<<": "<<Tree[child].uniqueVertex<<"("<< NodeOrder[Tree[child].uniqueVertex]<<","<<PartiTags[Tree[child].uniqueVertex].first<<") "<<line[i]<<"("<< NodeOrder[line[i]]<<","<<PartiTags[line[i]].first<<") "<<b<<"("<<NodeOrder[b]<<","<<PartiTags[b].first<<")"<<endl; exit(1);
//                    }

                        if(Tree[p].disPost[i]>vbW+Tree[rank[ancestor[i]]].disInf[inf_i]){
                            Tree[p].disPost[i]=vbW+Tree[rank[ancestor[i]]].disInf[inf_i];
                            Tree[p].FNPost[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                }
                else{//if not interface vertex
                    int bH=Tree[rank[b]].height-1,vbW=Tree[p].vert[k].second.first;
                    if(Tree[p].FNPost[bH]){//if distance label from child to b is directly sourced from their shortcut
                        if(Tree[p].DisRePost.find(b)!=Tree[p].DisRePost.end()){//if found b, all ancestor check
                            for(int i=rootHeight;i<bH;i++){
                                if(Tree[p].disPost[i]>vbW+Tree[rank[b]].disPost[i]){//update ancestor distance label
                                    Tree[p].disPost[i]=vbW+Tree[rank[b]].disPost[i];
                                    Tree[p].FNPost[i]=false;
                                    ProIDdisCha=true;
                                }
                            }
                            for(int i=bH+1;i<ancestor.size();i++){
                                if(Tree[p].disPost[i]>vbW+Tree[rank[ancestor[i]]].disPost[bH]){
                                    Tree[p].disPost[i]=vbW+Tree[rank[ancestor[i]]].disPost[bH];
                                    Tree[p].FNPost[i]=false;
                                    ProIDdisCha=true;
                                }
                            }

                        }else{//partial ancestor check if we cannot find b

                            if(vertexIDChL.find(b)!=vertexIDChL.end()){
                                for(int i=rootHeight;i<bH;i++){
                                    if(Tree[p].disPost[i]>vbW+Tree[rank[b]].disPost[i]){
                                        Tree[p].disPost[i]=vbW+Tree[rank[b]].disPost[i];
                                        Tree[p].FNPost[i]=false;
                                        ProIDdisCha=true;
                                    }
                                }
                            }
                            for(int i=bH+1;i<ancestor.size();i++){
                                if(Tree[p].disPost[i]>vbW+Tree[rank[ancestor[i]]].disPost[bH]){
                                    Tree[p].disPost[i]=vbW+Tree[rank[ancestor[i]]].disPost[bH];
                                    Tree[p].FNPost[i]=false;
                                    ProIDdisCha=true;
                                }
                            }

                        }
                    }
                }


            }
        }
        else{// if theisRePostre is no label for checking
            for(int k=0;k<Tree[p].vert.size();k++){
                int b=Tree[p].vert[k].first;
//            if(flagDebug ){//&& ancestor[i]==12771
//                cout<<"Neighbor. "<<ID1<<" "<<b<<"("<<PartiTags[b].first<<") "<<Tree[rank[b]].height-1<<endl;
//            }
                if(PartiTags[b].first == -1){//if it is interface vertex, check whether the distance label may be affected by interface
                    assert(BoundVertexMap[pid].find(b)!=BoundVertexMap[pid].end());
                    inf_i=BoundVertexMap[pid][b];
                    int vbW=Tree[p].vert[k].second.first;
                    for(int i=rootHeight;i<ancestor.size();i++){
//                    if(Tree[rank[ancestor[i]]].disInf.find(b) == Tree[rank[ancestor[i]]].disInf.end()){
//                        cout<<"Not found this boundary vertex! "<<i<<": "<<Tree[child].uniqueVertex<<"("<< NodeOrder[Tree[child].uniqueVertex]<<","<<PartiTags[Tree[child].uniqueVertex].first<<") "<<ancestor[i]<<"("<< NodeOrder[ancestor[i]]<<","<<PartiTags[ancestor[i]].first<<") "<<b<<"("<<NodeOrder[b]<<","<<PartiTags[b].first<<")"<<endl; exit(1);
//                    }
//                    if(flagDebug ){//&& ancestor[i]==12771
//                        cout<<"Boundary neighbor. "<<ID1<<" "<<b<<" "<<ancestor[i]<< " "<<Tree[p].disPost[i]<<endl;
//                    }
                        if(Tree[p].disPost[i]>vbW+Tree[rank[ancestor[i]]].disInf[inf_i]){
                            Tree[p].disPost[i]=vbW+Tree[rank[ancestor[i]]].disInf[inf_i];
                            Tree[p].FNPost[i]=false;
                            ProIDdisCha=true;
                        }
                    }


                }
                else {//if not interface vertex
                    int bH = Tree[rank[b]].height - 1, vbW = Tree[p].vert[k].second.first;
                    if (Tree[p].FNPost[bH]) {//Property 5, how to ensure the value of FNPost is correct?
                        if (vertexIDChL.find(b) != vertexIDChL.end()) {//if the distance label of b is changed
                            for (int i = rootHeight; i < bH; i++) {//check ancestor from 0 to bH
//                            if(flagDebug && ancestor[i]==12771){
//                                cout<<"Non-boundary neighbor. "<<ID1<<" "<<b<<" "<<ancestor[i]<< " "<<Tree[p].disPost[i]<<" "<<vbW + Tree[rank[b]].disPost[i]<<endl;
//                            }
                                if (Tree[p].disPost[i] > vbW + Tree[rank[b]].disPost[i]) {
                                    Tree[p].disPost[i] = vbW + Tree[rank[b]].disPost[i];
                                    Tree[p].FNPost[i] = false;
                                    ProIDdisCha = true;
                                }
                            }
                        }
                        for (int i = bH + 1; i < ancestor.size(); i++) {
//                        if(flagDebug && ancestor[i]==12771){
//                            cout<<"Non-boundary neighbor. "<<ID1<<" "<<b<<" "<<ancestor[i]<< " "<<Tree[p].disPost[i]<<endl;
//                        }
                            if (Tree[p].disPost[i] > vbW + Tree[rank[ancestor[i]]].disPost[bH]) {
                                Tree[p].disPost[i] = vbW + Tree[rank[ancestor[i]]].disPost[bH];
                                Tree[p].FNPost[i] = false;
                                ProIDdisCha = true;
                            }
                        }
                    }
                }
            }
        }
    }
    tt.stop();
//    tAncestor+=tt.GetRuntime();

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[p].uniqueVertex);
    }


    //nested loop
    ancestor.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        PostMHLIndexUpdatePostPartDFS(pid, Tree[p].ch[i],ancestor,interface,disInfs, vertexIDChL,  rootHeight, ifIncrease, tt);
    }
    ancestor.pop_back();
}

void Graph::Repair_PartiIndexForOpt(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
//    repairShortcuts.assign(node_num, unordered_map<vertex,pair<int,int>>());
    ifRepaired.assign(partiNum, false);
    if(ifParallel){
        // multi-thread
//        cout<<"Multi-thread computation!"<<endl;
        //multi-thread
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexV, this, boost::ref(processID[j]), ifIncrease, boost::ref(partiBatch), true));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            if(ifIncrease){//increase update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexIncreaseForOpt, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }
            else{//decrease update
                boost::thread_group thread;
                for(auto j=0;j<partiNum;++j){
                    thread.add_thread(new boost::thread(&Graph::RepairPartitionIndexDecreaseForOpt, this, j, boost::ref(partiBatch)));
                }
                thread.join_all();
            }

        }
    }
    else{
        // single thread
        if(ifIncrease){
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexIncreaseForOpt(k,partiBatch);
            }
        }
        else{
            for(int k=0;k<partiNum;k++){
//                cout<<"Repairing partition "<<k<<endl;
                RepairPartitionIndexDecreaseForOpt(k,partiBatch);
            }
        }

    }

//    int pNum=0;
//    for(auto it=ifRepaired.begin();it!=ifRepaired.end();++it){
//        if(*it){
//            ++pNum;
//        }
//    }
//    cout<<"Repaired partition number: "<<pNum<<endl;

}

//Function of repair the partition index
void Graph::RepairPartitionIndexV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, bool ifOpt) {
    if(ifOpt){
        if(ifIncrease){
            for(int i=0;i<p.size();++i){
                RepairPartitionIndexIncreaseForOpt(p[i], partiBatch);
            }
        }
        else{
            for(int i=0;i<p.size();++i){
                RepairPartitionIndexDecreaseForOpt(p[i], partiBatch);
            }
        }
    }else{
        if(ifIncrease){
            for(int i=0;i<p.size();++i){
                RepairPartitionIndexIncrease(p[i], partiBatch);
            }
        }
        else{
            for(int i=0;i<p.size();++i){
                RepairPartitionIndexDecrease(p[i], partiBatch);
            }
        }
    }


}
void Graph::RepairPartitionIndexDecrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();

    if(partiBatch.find(pid)!=partiBatch.end()){//if found
        weightsParti=partiBatch[pid];
    }

    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];

            if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                wlocal=NeighborsPartiPost[ID1][ID2];
            }else{
                cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryCore(ID1,ID2);
            bool found=false;//whether the boundary edge exist or not
            int wei;
            if(woverlay<wlocal){
                if(ID1<ID2){
                    weightsParti.emplace_back(make_pair(ID1,ID2),make_pair(wlocal,woverlay));
                }else{
                    weightsParti.emplace_back(make_pair(ID2,ID1),make_pair(wlocal,woverlay));
                }

//                DecreaseParti(ID1,ID2, woverlay, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);
            }else if(woverlay>wlocal){
                cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
            }

        }
    }
    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        bool ifPost=false;
        if(algoUpdate>=PH2H_Post){
            ifPost=true;
        }
        DecreasePartiBatch(pid,weightsParti, NeighborsPartiPost, TreesPost[pid], ranksPost[pid], heightMaxsPost[pid], ifPost);
        ifRepaired[pid]=true;
    }

}

void Graph::RepairPartitionIndexDecreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();

//    if(pid==12){
//        cout<<pid<<endl;
//    }

//    if(partiBatch.find(pid)!=partiBatch.end()){//if found
//        weightsParti=partiBatch[pid];
//    }

    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
//            wlocal=-1;
//            for(auto it=NeighborsParti[ID1].begin();it!=NeighborsParti[ID1].end();++it){
//                if(it->first==ID2){
//                    wlocal=it->second;
//                    break;
//                }
//            }
            if(NeighborsPartiPost[ID1].find(ID2)==NeighborsPartiPost[ID1].end()){
                cout<<"Not found! "<<ID1<<" "<<ID2<<endl; exit(1);
            }
            wlocal=NeighborsPartiPost[ID1][ID2];
//            wlocal= QueryH2HPartition(ID1,ID2,pid);
            if(wlocal==-1){//if found
                cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryCore(ID1,ID2);
            bool found=false;//whether the boundary edge exist or not
            int wei;
            if(woverlay<wlocal){
                if(ID1<ID2){
                    weightsParti.emplace_back(make_pair(ID1,ID2),make_pair(wlocal,woverlay));
                }else{
                    weightsParti.emplace_back(make_pair(ID2,ID1),make_pair(wlocal,woverlay));
                }

//                DecreaseParti(ID1,ID2, woverlay, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);
            }else if(woverlay>wlocal){
                cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
            }

        }
    }
    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti (with new) of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        bool ifPost=false;
        if(algoUpdate>=PH2H_Post){
            ifPost=true;
        }
        vector<pair<pair<int,int>,int>> updatedSC;
        DecreasePartiBatchForOpt(pid,weightsParti, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid], updatedSC, ifPost, false);
        ifRepaired[pid]=true;
    }else if(partiBatch.find(pid)!=partiBatch.end()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        DecreasePartiBatchLabel(TreesPost[pid], ranks[pid], heightMaxs[pid], ProBeginVertexSetParti[pid], vertexIDChLParti[pid]);
        ifRepaired[pid]=true;
    }

}

void Graph::RepairPartitionIndexIncrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();
    map<pair<int,int>,pair<int,int>> updatesSet;

    if(partiBatch.find(pid)!=partiBatch.end()){//if found
//        cout<<"Find for "<<pid<<" "<<partiBatch[pid].size()<<endl;
//        weightsParti=partiBatch[pid];
        for(auto it=partiBatch[pid].begin();it!=partiBatch[pid].end();++it){
            ID1=it->first.first, ID2=it->first.second;
            if(ID1<ID2){
                if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
                    updatesSet.insert(*it);
                }else{
                    cout<<"Original edge. Already exists! "<<ID1<<" "<<ID2<<" "<<updatesSet[make_pair(ID1,ID2)].first<<" "<<updatesSet[make_pair(ID1,ID2)].second<<endl; exit(1);
                }
//                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
            }else{
                if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
                    updatesSet.insert(*it);
                }else{
                    cout<<"Original edge. Already exists! "<<ID2<<" "<<ID1<<" "<<updatesSet[make_pair(ID2,ID1)].first<<" "<<updatesSet[make_pair(ID2,ID1)].second<<endl; exit(1);
                }
//                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
            }
        }
    }
    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
            if(NeighborsPartiPost[ID1].find(ID2)!=NeighborsPartiPost[ID1].end()){//if found
                wlocal=NeighborsPartiPost[ID1][ID2];
            }else{
                cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryCore(ID1,ID2);
            bool found=false;//whether the boundary edge exist or not
            int wei;
//            if(pid==49){
//                cout<<pid<<": "<<ID1<<" "<<ID2<<" "<<wlocal<<" "<<woverlay<<endl;
//            }

            if(woverlay>wlocal){
                // update partition index
                if(ID1<ID2){
                    if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID1,ID2),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"All-pair check. Already exists! "<<ID1<<" "<<ID2<<" "<<updatesSet[make_pair(ID1,ID2)].first<<" "<<updatesSet[make_pair(ID1,ID2)].second<<endl;
                        updatesSet[make_pair(ID1,ID2)].first=wlocal,updatesSet[make_pair(ID1,ID2)].second=woverlay;
//                        exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
                }else{
                    if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID2,ID1),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"All-pair check. Already exists! "<<ID2<<" "<<ID1<<" "<<updatesSet[make_pair(ID2,ID1)].first<<" "<<updatesSet[make_pair(ID2,ID1)].second<<endl;
                        updatesSet[make_pair(ID2,ID1)].first=wlocal,updatesSet[make_pair(ID1,ID2)].second=woverlay;
//                        exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
                }

//                DecreaseParti(ID1,ID2, woverlay, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);
            }else if(woverlay<wlocal){
                cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
            }

        }
    }
    for(auto it=updatesSet.begin();it!=updatesSet.end();++it){
        weightsParti.emplace_back(*it);
    }

    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        bool ifPost=false;
        if(algoUpdate>=PH2H_Post){
            ifPost=true;
        }
        IncreasePartiBatch(pid, weightsParti, NeighborsPartiPost, TreesPost[pid], ranksPost[pid], heightMaxsPost[pid],SCconNodesMTPost,VidtoTNidPost, ifPost);
        ifRepaired[pid]=true;
    }

}

void Graph::RepairPartitionIndexIncreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch) {
    int ID1,ID2,woverlay;
    int wlocal=INF;
    //boundary edges
    vector<pair<pair<int,int>,pair<int,int>>> weightsParti;//collect the changed edges on overlay graph
    weightsParti.clear();
    map<pair<int,int>,pair<int,int>> updatesSet;

//    if(partiBatch.find(pid)!=partiBatch.end()){//if found
////        cout<<"Find for "<<pid<<" "<<partiBatch[pid].size()<<endl;
////        weightsParti=partiBatch[pid];
//        for(auto it=partiBatch[pid].begin();it!=partiBatch[pid].end();++it){
//            ID1=it->first.first, ID2=it->first.second;
//            if(ID1<ID2){
//                if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
//                    updatesSet.insert(*it);
//                }else{
//                    cout<<"Already exists! "<<ID1<<" "<<ID2<<endl; exit(1);
//                }
////                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
//            }else{
//                if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
//                    updatesSet.insert(*it);
//                }else{
//                    cout<<"Already exists! "<<ID1<<" "<<ID2<<endl; exit(1);
//                }
////                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
//            }
//        }
//    }
    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];
            wlocal=-1;
            for(auto it=NeighborsParti[ID1].begin();it!=NeighborsParti[ID1].end();++it){
                if(it->first==ID2){
                    wlocal=it->second;
                    break;
                }
            }
            if(wlocal==-1){//if found
                cout<<"Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
            }
            woverlay=QueryCore(ID1,ID2);
            bool found=false;//whether the boundary edge exist or not
            int wei;
//            if(pid==49){
//                cout<<pid<<": "<<ID1<<" "<<ID2<<" "<<wlocal<<" "<<woverlay<<endl;
//            }
            if(woverlay>wlocal){
                // update partition index
                if(ID1<ID2){
                    if(updatesSet.find(make_pair(ID1,ID2))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID1,ID2),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"Already exists! "<<ID1<<" "<<ID2<<" "<<updatesSet[make_pair(ID1,ID2)].first<<" "<<updatesSet[make_pair(ID1,ID2)].second<<endl; exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID1,ID2), make_pair(wlocal, woverlay));//oldW,newW
                }else{
                    if(updatesSet.find(make_pair(ID2,ID1))==updatesSet.end()){//if not found
                        updatesSet.insert({make_pair(ID2,ID1),make_pair(wlocal, woverlay)});
                    }else{
                        cout<<"Already exists! "<<ID2<<" "<<ID1<<" "<<updatesSet[make_pair(ID2,ID1)].first<<" "<<updatesSet[make_pair(ID2,ID1)].second<<endl; exit(1);
                    }
//                    weightsParti.emplace_back(make_pair(ID2,ID1), make_pair(wlocal, woverlay));//oldW,newW
                }

//                DecreaseParti(ID1,ID2, woverlay, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid]);
            }else if(woverlay<wlocal){
                cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
            }

        }
    }
    for(auto it=updatesSet.begin();it!=updatesSet.end();++it){
        weightsParti.emplace_back(*it);
    }


    if(!weightsParti.empty()){
//        cout<<"Size of weightsParti (with new) of partition "<<pid<<" : "<<weightsParti.size()<<endl;
        bool ifPost=false;
        if(algoUpdate>=PH2H_Post){
            ifPost=true;
        }
        vector<pair<pair<int,int>,int>> updatedSC;
        IncreasePartiBatchForOpt(pid, weightsParti, NeighborsParti, Trees[pid], ranks[pid], heightMaxs[pid],SCconNodesMTP,VidtoTNidP,updatedSC, ifPost);
        ifRepaired[pid]=true;
    }else if(partiBatch.find(pid)!=partiBatch.end()){
//        cout<<"Size of weightsParti of partition "<<pid<<" : "<<weightsParti.size()<<endl;
//        if(pid==49){
//            cout<<pid<<endl;
//            for(auto it=ProBeginVertexSetParti[pid].begin();it!=ProBeginVertexSetParti[pid].end();++it){
//                cout<<*it<<"("<<NodeOrder[*it]<<")"<<endl;
//            }
//        }
        IncreasePartiBatchLabel(Trees[pid], ranks[pid], heightMaxs[pid], ProBeginVertexSetParti[pid], VidtoTNidP);
        ifRepaired[pid]=true;
    }

}

/// For extension labels
void Graph::ConstructExtensionLabelsNoAllPairTopDown(){
    TreeExt.clear();
    rankExt.assign(node_num,-1);
    vUpdated.assign(node_num, false);

    Node root;//root node
    int ID=Tree[0].uniqueVertex;
    root.height=1;
    root.uniqueVertex=ID;
//    root.dis.push_back(0);
    TreeExt.push_back(root);
    rankExt[ID] = 0;
    heightMaxExt=1;
//    cout<<PartiTag[ID].first<<" "<<ID<<" "<<NodeOrder[ID]<<endl;

    vector<int> list;

    makeTreeIndexExtDFS(ID,list);

    cout<<"TreeExt size: "<<TreeExt.size()<<" "<<node_num<<" ; tree height: "<<heightMaxExt<<" ; treewidth: "<<treewidth<< endl;
    for(int i=0;i<node_num;++i){
        if(rankExt[i]==-1){
            cout<<"Wrong! rankExt: "<<i<<" "<<rankExt[i]<<endl; exit(1);
        }
    }

    makeRMQ(toRMQExt, RMQIndexExt, TreeExt);//build LCA index
}

void Graph::ConstructExtensionLabelsNoAllPair(){
    TreeExt.clear();
    rankExt.assign(node_num,-1);
    vUpdated.assign(node_num, false);

    Node root;//root node
    int ID=Tree[0].uniqueVertex;
    root.height=1;
    root.uniqueVertex=ID;
//    root.dis.push_back(0);
    TreeExt.push_back(root);
    rankExt[ID] = 0;
    heightMaxExt=1;
//    cout<<PartiTag[ID].first<<" "<<ID<<" "<<NodeOrder[ID]<<endl;

    vector<int> list;

    makeTreeExtDFS(ID,list);

    cout<<"TreeExt size: "<<TreeExt.size()<<" "<<node_num<<" ; tree height: "<<heightMaxExt<<endl;
    for(int i=0;i<node_num;++i){
        if(rankExt[i]==-1){
            cout<<"Wrong! rankExt: "<<i<<" "<<rankExt[i]<<endl; exit(1);
        }
    }

    makeRMQ(toRMQExt, RMQIndexExt, TreeExt);//build LCA index

    bool ifParallel=true;

    if(ifParallel){//multi-thread
        cout<<"Multi-thread computation."<<endl;
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::ConstructExtensionLabelPartiV, this, boost::ref(processID[j]), false));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructExtensionLabelPartiNoAllPair, this, j));
            }
            thread.join_all();
        }
    }else{
        cout<<"Single-thread computation."<<endl;
        for(int pid=0;pid<node_num;++pid){
            ConstructExtensionLabelPartiNoAllPair(pid);
        }
    }


//    exit(0);
}

void Graph::ConstructExtensionLabels(){
    TreeExt.clear();
    rankExt.assign(node_num,-1);
    vUpdated.assign(node_num, false);

    Node root;//root node
    int ID=Tree[0].uniqueVertex;
    root.height=1;
    root.uniqueVertex=ID;
    TreeExt.push_back(root);
    rankExt[ID] = 0;
    heightMaxExt=1;
//    cout<<PartiTag[ID].first<<" "<<ID<<" "<<NodeOrder[ID]<<endl;

    vector<int> list;

    makeTreeExtDFS(ID,list);

    cout<<"TreeExt size: "<<TreeExt.size()<<" "<<node_num<<" ; tree height: "<<heightMaxExt<<endl;
    for(int i=0;i<node_num;++i){
        if(rankExt[i]==-1){
            cout<<"Wrong! rankExt: "<<i<<" "<<rankExt[i]<<endl; exit(1);
        }
    }

    makeRMQ(toRMQExt, RMQIndexExt, TreeExt);//build LCA index

    bool ifParallel=true;

    if(ifParallel){//multi-thread
        cout<<"Multi-thread computation."<<endl;
        if(threadnum<partiNum){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            vector<int> vertices;
            for(int pid=0;pid<partiNum;++pid){
                vertices.emplace_back(pid);
            }
            ThreadDistribute(vertices, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::ConstructExtensionLabelPartiV, this, boost::ref(processID[j]), true));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<partiNum;++j){
                thread.add_thread(new boost::thread(&Graph::ConstructExtensionLabelParti, this, j));
            }
            thread.join_all();
        }
    }else{
        cout<<"Single-thread computation."<<endl;
        for(int pid=0;pid<node_num;++pid){
            ConstructExtensionLabelParti(pid);
        }
    }


//    exit(0);
}

void Graph::makeTreeExtDFS(int ID, vector<int> &list){
    int PID=PartiTag[ID].first;
    unordered_set<int> children;//vertex id of children
    children.clear();
    int p=rankExt[ID];
//    int pa=list[list.size()-1];
    int tempRank;
    if(PartiTag[ID].second){//if boundary vertex
        for(int i=0;i<Tree[rank[ID]].ch.size();++i){
            tempRank=Tree[rank[ID]].ch[i];
            children.insert(Tree[tempRank].uniqueVertex);
        }
    }

//    cout<<"TreeExt size: "<<TreeExt.size()<<", "<<ID<<", Children size of Tree: "<<children.size()<<" , children size of Trees: "<<Trees[PID][ranks[PID][IDMap[ID]]].ch.size();
    bool ifCoexist=false;
    for(auto it=Trees[PID][ranks[PID][IDMap[ID]]].ch.begin();it!=Trees[PID][ranks[PID][IDMap[ID]]].ch.end();++it){
        tempRank=*it;
        if(PartiTag[Trees[PID][tempRank].uniqueVertex].second){//if boundary vertex
            ifCoexist=true;
            continue;
        }
        if(children.find(Trees[PID][tempRank].uniqueVertex) == children.end()){
            children.insert(Trees[PID][tempRank].uniqueVertex);
//            if(ifCoexist){
//                cout<<"!!! Vertex "<<ID<<"("<<PartiTag[ID].first<<","<< PartiTag[ID].second<<") has both boundary children and non-boundary children! "<<Trees[PID][tempRank].uniqueVertex<<endl;
//            }
        }
    }
//    cout<< ", overall children size: "<<children.size()<<" , height: "<<TreeExt[p].height<<endl;
//    if(TreeExt.size()%10000==0){
//        cout<<"TreeExt size: "<<TreeExt.size()<<", "<<ID<<", Children size of Tree: "<<children.size()<<" , children size of Trees: "<<Trees[PID][ranks[PID][IDMap[ID]]].ch.size()<< ", overall children size: "<<children.size()<<" , height: "<<TreeExt[p].height<<endl;
//    }
    Node node;
    if(!children.empty()){
        for(auto it=children.begin();it!=children.end();++it){
            node.uniqueVertex=*it;
            node.pa=ID;
            node.height=TreeExt[p].height+1;
            if(rankExt[*it]==-1){
                rankExt[*it]=TreeExt.size();
            }else{
//                cout<<"Wrong! "<<*it<<" "<<rankExt[*it]<<" "<<TreeExt.size()<<endl;
                continue;
                //exit(1);
            }

            TreeExt[p].ch.emplace_back(TreeExt.size());
            TreeExt.push_back(node);
            if(heightMaxExt<node.height){
                heightMaxExt=node.height;
            }
        }
    }


    //nested loop
    TreeExt[p].vAncestor=list;
    if(PartiTag[ID].second){//if boundary vertex
        TreeExt[p].pos=Tree[rank[ID]].pos;
    }

    list.push_back(TreeExt[p].uniqueVertex);
    for(int i=0;i<TreeExt[p].ch.size();i++){
        makeTreeExtDFS(TreeExt[TreeExt[p].ch[i]].uniqueVertex,list);
    }
    list.pop_back();
}

void Graph::makeTreeIndexExtDFS(int ID, vector<int> &list){
    int PID=PartiTag[ID].first;
    unordered_set<int> children;//vertex id of children
    children.clear();
    int p=rankExt[ID];
//    int pa=list[list.size()-1];
    int tempRank;
    if(PartiTag[ID].second){//if boundary vertex
        for(int i=0;i<Tree[rank[ID]].ch.size();++i){
            tempRank=Tree[rank[ID]].ch[i];
            children.insert(Tree[tempRank].uniqueVertex);
        }
    }

//    cout<<"TreeExt size: "<<TreeExt.size()<<", "<<ID<<", Children size of Tree: "<<children.size()<<" , children size of Trees: "<<Trees[PID][ranks[PID][IDMap[ID]]].ch.size();
    bool ifCoexist=false;
    for(auto it=Trees[PID][ranks[PID][IDMap[ID]]].ch.begin();it!=Trees[PID][ranks[PID][IDMap[ID]]].ch.end();++it){
        tempRank=*it;
        if(PartiTag[Trees[PID][tempRank].uniqueVertex].second){//if boundary vertex
            ifCoexist=true;
            continue;
        }
        if(children.find(Trees[PID][tempRank].uniqueVertex) == children.end()){
            children.insert(Trees[PID][tempRank].uniqueVertex);
//            if(ifCoexist){
//                cout<<"!!! Vertex "<<ID<<"("<<PartiTag[ID].first<<","<< PartiTag[ID].second<<") has both boundary children and non-boundary children! "<<Trees[PID][tempRank].uniqueVertex<<endl;
//            }
        }
    }
//    cout<< ", overall children size: "<<children.size()<<" , height: "<<TreeExt[p].height<<endl;
//    if(TreeExt.size()%10000==0){
//        cout<<"TreeExt size: "<<TreeExt.size()<<", "<<ID<<", Children size of Tree: "<<children.size()<<" , children size of Trees: "<<Trees[PID][ranks[PID][IDMap[ID]]].ch.size()<< ", overall children size: "<<children.size()<<" , height: "<<TreeExt[p].height<<endl;
//    }
    Node node;
    if(!children.empty()){
        for(auto it=children.begin();it!=children.end();++it){
            node.uniqueVertex=*it;
            node.pa=ID;
            node.height=TreeExt[p].height+1;
            if(rankExt[*it]==-1){
                rankExt[*it]=TreeExt.size();
            }else{
//                cout<<"Wrong! "<<*it<<" "<<rankExt[*it]<<" "<<TreeExt.size()<<endl;
                continue;
                //exit(1);
            }

            TreeExt[p].ch.emplace_back(TreeExt.size());
            TreeExt.push_back(node);
            if(heightMaxExt<node.height){
                heightMaxExt=node.height;
            }
        }
    }


    //nested loop
    TreeExt[p].vAncestor=list;
    int ancestor, neighbor, d1, d2;
    int dis=INF, vid;
    if(PartiTag[ID].second){//if boundary vertex
        TreeExt[p].pos=Tree[rank[ID]].pos;
        if(treewidth<TreeExt[p].pos.size()){
            treewidth=TreeExt[p].pos.size();
        }
        for(int i=0;i<list.size();++i){
            ancestor = list[i];
            dis=INF;
            for(int j=0;j<Tree[rank[ID]].vert.size();++j){
                neighbor = Tree[rank[ID]].vert[j].first;
                d1=Tree[rank[ID]].vert[j].second.first;
//                if(TreeExt[rankExt[neighbor]].dis.size()<list.size()){
//                    cout<<ID<<": "<<neighbor<<"("<<NodeOrder[neighbor]<<","<<PartiTag[neighbor].first<<","<<PartiTag[neighbor].second<<") "<<TreeExt[rankExt[neighbor]].dis.size()<<" "<<list.size()<<endl;
//                }
                if(TreeExt[rankExt[neighbor]].height<TreeExt[rankExt[ancestor]].height){//if neighbor is ancestor of ancestor
//                    d2= QueryCore(neighbor,ancestor);
                    d2=TreeExt[rankExt[ancestor]].dis[TreeExt[rankExt[neighbor]].height-1];;
                }else{
                    d2=TreeExt[rankExt[neighbor]].dis[i];
                }

                if(dis>d1+d2){
                    dis=d1+d2;
                }
            }
            TreeExt[p].dis.push_back(dis);
//            int disD= Dijkstra(ancestor,ID,Neighbor);
//            if(disD!=dis){
//                cout<<"Boundary. Incorrect! "<<ID<<" "<<ancestor<<" "<<dis<<" "<<disD<<endl;
//            }
        }
    }
    else{//if non-boundary vertex
        if(treewidth<Trees[PID][ranks[PID][IDMap[ID]]].vert.size()+1){
            treewidth=Trees[PID][ranks[PID][IDMap[ID]]].vert.size()+1;
        }
        for(int i=0;i<list.size();++i){
            ancestor = list[i];
            dis=INF;
            for(int j=0;j<Trees[PID][ranks[PID][IDMap[ID]]].vert.size();++j){
                neighbor = Trees[PID][ranks[PID][IDMap[ID]]].vert[j].first;
                d1=Trees[PID][ranks[PID][IDMap[ID]]].vert[j].second.first;
//                if(TreeExt[rankExt[neighbor]].dis.size()<list.size()){
//                    cout<<ID<<": "<<neighbor<<"("<<NodeOrder[neighbor]<<","<<PartiTag[neighbor].first<<","<<PartiTag[neighbor].second<<") "<<TreeExt[rankExt[neighbor]].dis.size()<<" "<<list.size()<<endl;
//                }

                if(TreeExt[rankExt[neighbor]].height<TreeExt[rankExt[ancestor]].height){//if neighbor is ancestor of ancestor
                    d2= TreeExt[rankExt[ancestor]].dis[TreeExt[rankExt[neighbor]].height-1];
                    if(dis>d1+d2){
                        dis=d1+d2;
                        vid=ancestor;
                    }
                }else{
                    d2=TreeExt[rankExt[neighbor]].dis[i];
                    if(dis>d1+d2){
                        dis=d1+d2;
                        vid=neighbor;
                    }
                }
//                d2=TreeExt[rankExt[neighbor]].dis[i];

//                if(ID==25991){
//                    cout<<ID<<": "<<neighbor<<" "<<ancestor<<" "<<d1<<" "<<d2<<" "<<dis<<endl;
//                }
            }
            TreeExt[p].dis.push_back(dis);
            TreeExt[p].cnt.push_back(vid);
//            int disD= Dijkstra(ancestor,ID,Neighbor);
//            if(disD!=dis){
//                cout<<"Non-boundary. Incorrect! "<<ID<<" "<<ancestor<<" "<<dis<<" "<<disD<<endl;
//            }
        }
    }
    TreeExt[p].dis.push_back(0);
    TreeExt[p].cnt.push_back(ID);

    list.push_back(TreeExt[p].uniqueVertex);
    for(int i=0;i<TreeExt[p].ch.size();i++){
        makeTreeIndexExtDFS(TreeExt[TreeExt[p].ch[i]].uniqueVertex,list);
    }
    list.pop_back();
}

void Graph::ConstructExtensionLabelPartiV(vector<int>& p, bool ifAllPair){
    if(ifAllPair){
        for(int i=0;i<p.size();++i){
            ConstructExtensionLabelParti(p[i]);
        }
    }else{
        for(int i=0;i<p.size();++i){
            ConstructExtensionLabelPartiNoAllPair(p[i]);
        }
    }

}
//top-down propagation
void Graph::ConstructExtensionLabelPartiNoAllPair(int pid) {
    int ID;
    vector<int> ancestorsV;
//    for(auto it=ancestors.begin();it!=ancestors.end();++it){
//        ancestorsV.emplace_back(it->ID);
////        cout<<"Ancestor "<<ancestorsV.size()<<": "<<it->ID<<" "<<NodeOrder[it->ID]<<" "<<PartiTag[it->ID].second<<endl;
//    }
    int bid,d1,d2,dis=INF;
    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        ancestorsV=TreeExt[rankExt[ID]].vAncestor;
        if(PartiTag[ID].second){
            continue;
        }
//        cout<<ID<<" "<<TreeExt[rankExt[ID]].vAncestor.size()<<": "<<endl;
//        for(auto it=TreeExt[rankExt[ID]].vAncestor.begin();it!=TreeExt[rankExt[ID]].vAncestor.end();++it){
//            cout<<*it<<" "<<NodeOrder[*it]<<" "<<PartiTag[*it].second<<endl;
//        }
//        cout<<endl;
//        TreeExt[rankExt[ID]].vAncestor=ancestorsV;
        for(auto it=ancestorsV.begin();it!=ancestorsV.end();++it){
            if(!PartiTag[*it].second){//if non-boundary vertex
//                if(ID==38175){
//                    cout<<*it<<" "<<PartiTag[*it].first<<","<<PartiTag[*it].second<<endl;
//                }
                continue;
            }
            dis = INF;
            for(auto it2=BoundVertex[pid].begin();it2!=BoundVertex[pid].end();++it2){
//                d1 = QueryH2HPartitionPost(ID,*it2,pid);//original
                d1 = QueryH2HPartition(ID,*it2,pid);
                d2 = QueryCore(*it2,*it);
                if(d1+d2<dis){
                    dis=d1+d2;
                    bid=*it2;
                }
            }
            TreeExt[rankExt[ID]].dis.push_back(dis);//the distance to ancestors
            TreeExt[rankExt[ID]].cnt.push_back(bid);//the corresponding boundary vertex for the label
        }
    }
}
//4-hop concatenation
/*void Graph::ConstructExtensionLabelPartiNoAllPair(int pid) {
    int ID;
    vector<int> ancestorsV;
//    for(auto it=ancestors.begin();it!=ancestors.end();++it){
//        ancestorsV.emplace_back(it->ID);
////        cout<<"Ancestor "<<ancestorsV.size()<<": "<<it->ID<<" "<<NodeOrder[it->ID]<<" "<<PartiTag[it->ID].second<<endl;
//    }
    int bid,d1,d2,dis=INF;
    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        ancestorsV=TreeExt[rankExt[ID]].vAncestor;
        if(PartiTag[ID].second){
            continue;
        }
//        cout<<ID<<" "<<TreeExt[rankExt[ID]].vAncestor.size()<<": "<<endl;
//        for(auto it=TreeExt[rankExt[ID]].vAncestor.begin();it!=TreeExt[rankExt[ID]].vAncestor.end();++it){
//            cout<<*it<<" "<<NodeOrder[*it]<<" "<<PartiTag[*it].second<<endl;
//        }
//        cout<<endl;
//        TreeExt[rankExt[ID]].vAncestor=ancestorsV;
        for(auto it=ancestorsV.begin();it!=ancestorsV.end();++it){
            if(!PartiTag[*it].second){//if non-boundary vertex
//                if(ID==38175){
//                    cout<<*it<<" "<<PartiTag[*it].first<<","<<PartiTag[*it].second<<endl;
//                }
                continue;
            }
            dis = INF;
            for(auto it2=BoundVertex[pid].begin();it2!=BoundVertex[pid].end();++it2){
//                d1 = QueryH2HPartitionPost(ID,*it2,pid);//original
                d1 = QueryH2HPartition(ID,*it2,pid);
                d2 = QueryCore(*it2,*it);
                if(d1+d2<dis){
                    dis=d1+d2;
                    bid=*it2;
                }
            }
            TreeExt[rankExt[ID]].dis.push_back(dis);//the distance to ancestors
            TreeExt[rankExt[ID]].cnt.push_back(bid);//the corresponding boundary vertex for the label
        }
    }
}*/
//old
void Graph::ConstructExtensionLabelParti(int pid){
    int ID;
    set<OrderCompMax> ancestors;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID=BoundVertex[pid][i];
//        cout<<pid<<" "<<ID<<" "<<NodeOrder[ID]<<endl;
        for(auto it=Tree[rank[ID]].vAncestor.begin();it!=Tree[rank[ID]].vAncestor.end();++it){
            if(ancestors.find(OrderCompMax(*it))==ancestors.end()){
                ancestors.insert(OrderCompMax(*it));
            }
        }
        if(ancestors.find(OrderCompMax(ID))==ancestors.end()){
            ancestors.insert(OrderCompMax(ID));
        }
    }
//    cout<<"Partition "<<pid<<", Size of ancestors: "<<ancestors.size()<<endl;
    vector<int> ancestorsV;
    for(auto it=ancestors.begin();it!=ancestors.end();++it){
        ancestorsV.emplace_back(it->ID);
//        cout<<"Ancestor "<<ancestorsV.size()<<": "<<it->ID<<" "<<NodeOrder[it->ID]<<" "<<PartiTag[it->ID].second<<endl;
    }
    int bid,d1,d2,dis=INF;
    for(int i=0;i<PartiVertex[pid].size();++i){
        ID=PartiVertex[pid][i];
        if(PartiTag[ID].second){
            continue;
        }
//        cout<<ID<<" "<<TreeExt[rankExt[ID]].vAncestor.size()<<": "<<endl;
//        for(auto it=TreeExt[rankExt[ID]].vAncestor.begin();it!=TreeExt[rankExt[ID]].vAncestor.end();++it){
//            cout<<*it<<" "<<NodeOrder[*it]<<" "<<PartiTag[*it].second<<endl;
//        }
//        cout<<endl;
        TreeExt[rankExt[ID]].vAncestor=ancestorsV;
        for(auto it=ancestorsV.begin();it!=ancestorsV.end();++it){
            dis = INF;
            for(auto it2=BoundVertex[pid].begin();it2!=BoundVertex[pid].end();++it2){
//                d1 = QueryH2HPartitionPost(ID,*it2,pid);
                d1 = QueryH2HPartition(ID,*it2,pid);
                d2 = QueryCore(*it2,*it);
                if(d1+d2<dis){
                    dis=d1+d2;
                    bid=*it2;
                }
            }
            TreeExt[rankExt[ID]].dis.push_back(dis);//the distance to ancestors
            TreeExt[rankExt[ID]].cnt.push_back(bid);//the corresponding boundary vertex for the label
        }
    }
}
//Top down mechanism
void Graph::RefreshExtensionLabelsNoAllPair(map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
//    updated.assign(node_num, false);
    cout<<"Top-down extension label update!"<<endl;
    int bid;
    vector<int> vParti;
    for(int pid=0;pid<partiNum;++pid){
        if(partiBatch.find(pid)!=partiBatch.end()){//if found
            vParti.emplace_back(pid);
            continue;
        }
        for(int i=0;i<BoundVertex[pid].size();++i){
            bid = BoundVertex[pid][i];
            if(vUpdated[bid]){
                vParti.emplace_back(pid);
                break;
            }
        }
    }
//    vParti.emplace_back(69);
    cout<<"vParti size: "<<vParti.size()<<endl;

    int pid, ProID;
    int rnew,r;
    vector<int> rootsTemp, rootsFinal;
    rootsFinal.push_back(PartiVertex[vParti[0]][0]);
    for(int i=1;i<vParti.size();++i){
        pid=vParti[i];
        ProID=PartiVertex[pid][0];
        rnew=rankExt[ProID];
        rootsTemp.clear();
        bool ifbelong=false;
        for(int j=0;j<rootsFinal.size();++j){
            r=rankExt[rootsFinal[j]];
//            if(!ifbelong){
                int LCA=LCAQuery(rnew,r,toRMQExt,RMQIndexExt,TreeExt);
                if(LCA==rnew){
                    rootsTemp.push_back(ProID);
                    ifbelong=true;
                }else if(LCA==r){
                    rootsTemp.push_back(rootsFinal[j]);
                    ifbelong= true;
                }else{
//                    rootsTemp.push_back(ProID);
                    rootsTemp.push_back(rootsFinal[j]);
                }
//            }else{//if r belongs to existing branch
//                rootsTemp.push_back(rootsFinal[j]);
//            }

        }
        if(!ifbelong){
            rootsTemp.push_back(ProID);
        }
        rootsFinal=rootsTemp;
    }
    cout<<"rootsFinal size: "<<rootsFinal.size()<<endl;


    bool ifParallel=true;
//    ifParallel=false;

    if(ifParallel){//multi-thread
        cout<<"Multi-thread computation."<<endl;
        if(threadnum<vParti.size()){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            ThreadDistribute(rootsFinal, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelTopDownV, this, boost::ref(processID[j])));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<rootsFinal.size();++j){
                thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelTopDown, this, rootsFinal[j]));
            }
            thread.join_all();
        }
    }else{//single thread
        cout<<"Single-thread computation."<<endl;
//        for(int i=0;i<vParti.size();++i) {
//            RefreshExtensionLabelPartiTopDown(vParti[i]);
//        }

        for(int i=0;i<rootsFinal.size();++i) {
            RefreshExtensionLabelTopDown(rootsFinal[i]);
        }
    }

}
//old version, no LCA computation
/*void Graph::RefreshExtensionLabelsNoAllPair(map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
//    updated.assign(node_num, false);
    cout<<"Top-down extension label update!"<<endl;
    int bid;
    vector<int> vParti;
    for(int pid=0;pid<partiNum;++pid){
        if(partiBatch.find(pid)!=partiBatch.end()){
            vParti.emplace_back(pid);
            continue;
        }
        for(int i=0;i<BoundVertex[pid].size();++i){
            bid = BoundVertex[pid][i];
            if(vUpdated[bid]){
                vParti.emplace_back(pid);
                break;
            }
        }
    }
//    vParti.emplace_back(69);
    cout<<"vParti size: "<<vParti.size()<<endl;


    bool ifParallel=true;
    //ifParallel=false;

    if(ifParallel){//multi-thread
        cout<<"Multi-thread computation."<<endl;
        if(threadnum<vParti.size()){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            ThreadDistribute(vParti, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelPartiV, this, boost::ref(processID[j]), true));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<vParti.size();++j){
//                thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelParti, this, vParti[j]));
                thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelPartiTopDown, this, vParti[j]));
            }
            thread.join_all();
        }
    }else{//single thread
        cout<<"Single-thread computation."<<endl;
        for(int i=0;i<vParti.size();++i) {
            RefreshExtensionLabelPartiTopDown(vParti[i]);
        }
    }

}*/

void Graph::RefreshExtensionLabelsPostMHL(bool ifParallel, bool ifIncrease, double & runT){
//    updated.assign(node_num, false);
    Timer tt;
    tt.start();
    int bid;
    vector<int> vParti;
    for(int pid=0;pid<partiNum;++pid){
        if(!ProBeginVertexSetPartiExtend[pid].empty()){
            vParti.push_back(pid);
        }
    }
//    vParti.emplace_back(69);
    cout<<"vParti size: "<<vParti.size()<<endl;

//    for(int i=0;i<vParti.size();++i){
//        cout<<vParti[i]<<" ";
//    }
//    cout<<endl;

    if(ifParallel){//multi-thread
        cout<<"Multi-thread computation."<<endl;
        if(threadnum<vParti.size()){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            ThreadDistribute(vParti, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            if(ifIncrease){
                boost::thread_group thread;
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchLabelPostMHLExtendV, this, boost::ref(processID[j]), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(ProBeginVertexSetPartiExtend), boost::ref(VidtoTNid)));
                }
                thread.join_all();
            }else{
                boost::thread_group thread;
                for(auto j=0;j<processID.size();++j){
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchLabelPostMHLExtendV, this, boost::ref(processID[j]), boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(ProBeginVertexSetPartiExtend), boost::ref(vertexIDChLOverlay)));
                }
                thread.join_all();
            }

        }
        else{//thread number is sufficient
            if(ifIncrease){
                boost::thread_group thread;
                for(auto j=0;j<vParti.size();++j){
                    thread.add_thread(new boost::thread(&Graph::IncreasePartiBatchLabelPostMHLExtend, this, vParti[j], boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(ProBeginVertexSetPartiExtend[vParti[j]]), boost::ref(VidtoTNid)));
                }
                thread.join_all();
            }else{
                boost::thread_group thread;
                for(auto j=0;j<vParti.size();++j){
                    thread.add_thread(new boost::thread(&Graph::DecreasePartiBatchLabelPostMHLExtend, this, vParti[j], boost::ref(Tree), boost::ref(rank), heightMax, boost::ref(ProBeginVertexSetPartiExtend[vParti[j]]), boost::ref(vertexIDChLOverlay)));
                }
                thread.join_all();
            }

        }
    }else{//single thread
        cout<<"Single-thread computation."<<endl;
        if(ifIncrease){
            for(int i=0;i<vParti.size();++i){
                int pid=vParti[i];
//            cout<<"Extend update of partition "<<pid<<endl;
                IncreasePartiBatchLabelPostMHLExtend(pid, Tree, rank, heightMax, ProBeginVertexSetPartiExtend[pid], VidtoTNid);
            }
        }else{
            for(int i=0;i<vParti.size();++i){
                int pid=vParti[i];
//            cout<<"Extend update of partition "<<pid<<endl;
                DecreasePartiBatchLabelPostMHLExtend(pid, Tree, rank, heightMax, ProBeginVertexSetPartiExtend[pid], vertexIDChLOverlay);
            }
        }

    }
    tt.stop();
    runT=tt.GetRuntime();
    algoQuery=PH2H_Cross;
}

void Graph::RefreshExtensionLabels(map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch){
//    updated.assign(node_num, false);
    int bid;
    vector<int> vParti;
    for(int pid=0;pid<partiNum;++pid){
        if(partiBatch.find(pid)!=partiBatch.end()){
            vParti.emplace_back(pid);
            continue;
        }
        for(int i=0;i<BoundVertex[pid].size();++i){
            bid = BoundVertex[pid][i];
            if(vUpdated[bid]){
                vParti.emplace_back(pid);
                break;
            }
        }
    }
//    vParti.emplace_back(69);
    cout<<"vParti size: "<<vParti.size()<<endl;

//    for(int i=0;i<vParti.size();++i){
//        cout<<vParti[i]<<" ";
//    }
//    cout<<endl;
    bool ifParallel=true;

    if(ifParallel){//multi-thread
        cout<<"Multi-thread computation."<<endl;
        if(threadnum<vParti.size()){
            vector<vector<int>> processID;
            processID.assign(threadnum, vector<int>());
            ThreadDistribute(vParti, processID);
//            cout<<"Batch number: "<<processID[0].size()<<endl;
            boost::thread_group thread;
            for(auto j=0;j<processID.size();++j){
                thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelPartiV, this, boost::ref(processID[j]), false));
            }
            thread.join_all();
        }
        else{//thread number is sufficient
            boost::thread_group thread;
            for(auto j=0;j<vParti.size();++j){
                thread.add_thread(new boost::thread(&Graph::RefreshExtensionLabelParti, this, vParti[j]));
            }
            thread.join_all();
        }
    }else{//single thread
        cout<<"Single-thread computation."<<endl;
        for(int i=0;i<vParti.size();++i){
            RefreshExtensionLabelParti(vParti[i]);
        }
    }

}

pair<int,int> Graph::ExtensionLabelCompute(int pid, int ID, int ancestor){
    int dis = INF;
    int d1,d2,bid=-1;
    for(auto it2=BoundVertex[pid].begin();it2!=BoundVertex[pid].end();++it2){
//        d1 = QueryH2HPartitionPost(ID,*it2,pid);
        d1 = QueryH2HPartition(ID,*it2,pid);
        d2 = QueryCore(*it2,ancestor);
        if(d1+d2<dis){
            dis=d1+d2;
            bid=*it2;
        }
    }
//    TreeExt[rankExt[ID]].dis[j]=dis;//the distance to ancestors
//    TreeExt[rankExt[ID]].cnt[j]=bid;//the corresponding boundary vertex for the label
    return make_pair(dis, bid);
}

/// PostMHL

void Graph::TreeDecompositionPartitioningNaive() {
    // MDE contraction
    TDPContract();
    // Create tree and partitions
    TDPCreateTreeAndPartiNaive();
}

void Graph::TreeDecompositionPartitioning(int pNum, double bRatioUpper, double bRatioLower) {
    // MDE contraction
    TDPContract();
    // Create tree and partitions
    TDPCreateTreeAndParti(pNum, bRatioUpper, bRatioLower);
}

//Function of MDE-based contracting for H2H
void Graph::TDPContract(){
    HighestOrder = node_num;
    //for H2H update
    SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());

    //initialize E
    E.assign(node_num,map<int,pair<int,int>>());
    for(int i=0;i<Neighbor.size();i++){
        for(int j=0;j<Neighbor[i].size();j++)
            E[i].insert(make_pair(Neighbor[i][j].first,make_pair(Neighbor[i][j].second,1)));
    }
//    EOrder.assign(node_num,map<OrderCompMin,pair<int,int>>());
//    for(int i=0;i<Neighbor.size();i++){
//        for(int j=0;j<Neighbor[i].size();j++)
//            EOrder[i].insert(make_pair(OrderCompMin(Neighbor[i][j].first),make_pair(Neighbor[i][j].second,1)));
//    }


    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
    DD.assign(node_num,0); //DD2.assign(node_num,0);

    set<DegComp> Deg;//min first
//    vector<bool> active(node_num,false);//flag that indicate whether a vertex is active for contraction
    int degree;
    for(int i=0;i<node_num;i++){
        degree=Neighbor[i].size();

        if(degree > 0){//get degree
            _DD_[i]=degree;
            DD[i]=degree;
            Deg.insert(DegComp(i));
//            active[i] = true;
        }else{
            cout<<"Wrong!! Degree of "<<i<<" is "<<degree<<endl;
            exit(1);
        }
    }

    vNodeOrder.clear();
    //vector<bool> exist;
    existOverlay.assign(node_num,true);//if in the core, all vertices is originally in core
    vector<bool> change;
    change.assign(node_num,false);//whether the neighbor (degree) has changed

    vector<pair<int,pair<int,int>>> vect;
    NeighborCon.assign(node_num,vect);//temporal graph to store Neighbors in the core, for graph contraction

    bool CutLabel=false;
    int count=0;
    int ID1,ID2;

    //Get the order of all vertices by MDE
    while(!Deg.empty()){
        count+=1;
        int x=(*Deg.begin()).x;//minimum degree first

        while(change[x]){//update the degree if it is changed
            Deg.erase(DegComp(x));
            _DD_[x]=DD[x];
            Deg.insert(DegComp(x));
            change[x]=false;
            x=(*Deg.begin()).x;
        }

        vNodeOrder.push_back(x);//least important vertex first
        Deg.erase(Deg.begin());

        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[x].begin();it!=E[x].end();it++){
            if(existOverlay[(*it).first]){
                Neigh.push_back(*it);
            }else{
                cout<<"Strange! "<<endl; exit(1);
            }
        }
        NeighborCon[x].assign(Neigh.begin(),Neigh.end());

        /// if still need to contract
        if(!CutLabel){
            if(Neigh.size()<bandWidth){//if the neighbor is smaller than tree width threshold, the vertex will be a part of tree
                existOverlay[x]=false;
                //delete the star
                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteECore(x,y);//delete x from y's adjacency list
                    change[y]=true;
                }
                //add all-pair neighbors
                for(int i=0;i<Neigh.size();i++){
                    ID1=Neigh[i].first;
                    for(int j=i+1;j<Neigh.size();j++){
                        ID2=Neigh[j].first;
                        insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                        /// For TD update
                        if(ID1<ID2){
                            SCconNodesMT[ID1][ID2].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                        }
                        else if(ID1>ID2){
                            SCconNodesMT[ID2][ID1].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);
                        }

                        change[ID1]=true;
                        change[ID2]=true;
                    }
                }
            }
            else{//else vertices in the core, Neigh.size()>=Width
                HighestOrder=vNodeOrder.size()-1;
                cout<<"Highest Order for vertex in periphery "<<HighestOrder<<", x = "<<x<<endl;
                CutLabel=true;
                //delete the star
                for(int i=0;i<Neigh.size();i++){
                    int y=Neigh[i].first;
                    deleteECore(x,y);//delete x from y's adjacency list
                    change[y]=true;
                }
                //add all-pair neighbors
                for(int i=0;i<Neigh.size();i++){
                    ID1=Neigh[i].first;
                    for(int j=i+1;j<Neigh.size();j++){
                        ID2=Neigh[j].first;
                        insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                        /// For TD update
                        if(ID1<ID2){
                            SCconNodesMT[ID1][ID2].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                        }
                        else if(ID1>ID2){
                            SCconNodesMT[ID2][ID1].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);
                        }

                        change[ID1]=true;
                        change[ID2]=true;
                    }
                }
            }
        }else{
            //delete the star
            for(int i=0;i<Neigh.size();i++){
                int y=Neigh[i].first;
                deleteECore(x,y);//delete x from y's adjacency list
                change[y]=true;
            }
            //add all-pair neighbors
            for(int i=0;i<Neigh.size();i++){
                ID1=Neigh[i].first;
                for(int j=i+1;j<Neigh.size();j++){
                    ID2=Neigh[j].first;
                    insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                    /// For TD update
                    if(ID1<ID2){
                        SCconNodesMT[ID1][ID2].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                    }
                    else if(ID1>ID2){
                        SCconNodesMT[ID2][ID1].emplace_back(x,Neigh[i].second.first+Neigh[j].second.first);
                    }

                    change[ID1]=true;
                    change[ID2]=true;
                }
            }

        }

    }

//    cout<<existCore[46862]<<" "<<existCore[47489]<<endl;

    NodeOrder.assign(node_num,-1);
    for(int k=0;k<vNodeOrder.size();k++){
        NodeOrder[vNodeOrder[k]]=k;
    }

    NodeOrder_ = NodeOrder;

//    string orderF = graphfile + ".orderH2H" + to_string(bandWidth);
//    ifstream IF(orderF);
//    if(!IF){//if the label file does not exist, construct it
//        WriteOrder(orderF);
//    }
//    IF.close();
//    WriteOrder(graphfile+".order1");
//    WriteOrder(graphfile+".order2");
//    CompareOrder(graphfile+".order1", graphfile+".order2");
//    exit(0);

}
//Function of creating tree
void Graph::TDPCreateTreeAndPartiNaive(){
    cout<<"Creating tree..."<<endl;
    //// Get tree
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(node_num,vecemp);

    //rank.assign(HighestOrder+2,0);
    rank.clear();
    rank.assign(node_num,0);//the vector index of tree nodes, map from vertex to tree node
    int len=HighestOrder-1;
    len=node_num-1;
    heightMax=0;
    set<OrderCompMax> orderOverlayV; orderOverlayV.clear();
    Node root;//virtual root node
    int x = vNodeOrder[len];
    if(NeighborCon[x].empty()){
        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=x;
        orderOverlayV.insert(OrderCompMax(x));
        len--;
    }else{
        root.uniqueVertex=-1;
    }
    root.height=1;
    Tree.push_back(root);


    int nn;
    for(;len>=0;len--){//check the vertices with order lower than HighestOrder
        x=vNodeOrder[len];
        if(existOverlay[x]){
            orderOverlayV.insert(OrderCompMax(x));
//            OverlayVertex.push_back(x);
//            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[x];//
        nod.uniqueVertex=x;
        int pa=matchCore(x,NeighborCon[x],rank);

        //cout<<"pa "<<pa<<endl;

        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        /// for update
        nod.hdepth=Tree[pa].height+1;
        for(int i=0;i<NeighborCon[x].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[x][i].first;
            VidtoTNid[nn].push_back(Tree.size());//record the child tree node rank who has direct super edge to nn
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1){
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
            }
        }
        if(nod.height>heightMax){
            heightMax=nod.height;
        }

        rank[x]=Tree.size();//the position of tree, higher-order vertex has lower rank
        Tree.push_back(nod);
        if(existOverlay[Tree[pa].uniqueVertex] && !existOverlay[x]){//if the father is in overlay while x is not in overlay
            partiRoots.push_back(x);
        }
    }
    partiNum=partiRoots.size();
    cout<<"Tree node number: "<<Tree.size()<<endl;



    //// Get partitions
//    vector<int> cr;//the vertex rank of root's children
//    cr.clear();
//    int temp_ch=0;
//    int temp_vert=0;
//    int sum_singular=0;
    BoundVertex.assign(partiNum,vector<int>());
    BoundVertexMap.assign(partiNum, unordered_map<int,int>());
    PartiTags.assign(node_num, make_pair(-1,set<int>()));
    PartiVertex.assign(partiNum, vector<int>());
    vector<set<OrderCompMax>> orderPartiV(partiNum,set<OrderCompMax>());

//    vector<int> tRoots;
    for(int k=0;k<partiRoots.size();k++){//Tree[0] is the core
        int id=partiRoots[k];
        //get boundary vertices
        for(int j=0;j<Tree[rank[id]].vert.size();++j){
            BoundVertex[k].push_back(Tree[rank[id]].vert[j].first);
            BoundVertexMap[k].insert({Tree[rank[id]].vert[j].first,j});
            if(orderOverlayV.find(OrderCompMax(Tree[rank[id]].vert[j].first))==orderOverlayV.end()){
                cout<<"Not found this boundary vertex! "<<Tree[rank[id]].vert[j].first<<" "<< NodeOrder[Tree[rank[id]].vert[j].first] <<endl; exit(1);
            }
            if( PartiTags[Tree[rank[id]].vert[j].first].first!=-1){
                cout<<"!!! Strange! "<<Tree[rank[id]].vert[j].first<<" "<<PartiTags[Tree[rank[id]].vert[j].first].first<<endl; exit(1);
//                PartiTags[it->first].first=-1;
            }

            PartiTags[Tree[rank[id]].vert[j].first].second.insert(k);
        }
        //get partition vertex and graph
        makePartitionDFS(orderPartiV,k,id);
    }
    int noBoundayNum=0;
    for(auto it=orderOverlayV.begin();it!=orderOverlayV.end();++it){
//        cout<<PartiTag[it->ID].first<<": "<<it->ID<<" "<<NodeOrder[it->ID]<<" "<<NeighborsOverlay[it->ID].size()<<endl;
        OverlayVertex.emplace_back(it->ID);
        if(!PartiTags[it->ID].second.empty()){
//            for(auto it2=PartiTags[it->ID].second.begin();it2!=PartiTags[it->ID].second.end();++it2){
//                BoundVertex[*it2].emplace_back(it->ID);
//            }
        }else{
            noBoundayNum++;
        }
    }
//    IDMap.assign(node_num,-1);
//    for(int pid=0;pid<partiNum;++pid){
//        int id=0;
//        if(BoundVertex[pid].empty()){
//            cout<<"Boundary vertex number of partition "<<pid<<" is empty!"<<endl; exit(1);
//        }
//        for(auto it=orderPartiV[pid].begin();it!=orderPartiV[pid].end();++it){
////            cout<<pid<<": "<<it->ID<<" "<<NodeOrder[it->ID]<<" "<<PartiTag[it->ID].second<<endl;
//            PartiVertex[pid].emplace_back(it->ID);
////            IDMap[it->ID] = id;
//            ++id;
//        }
//    }

    int pmin=INF,pmax=0,pave=0;
    for(int i=0;i<partiNum;++i){
        if(pmin>PartiVertex[i].size()){
            pmin=PartiVertex[i].size();
        }
        if(pmax<PartiVertex[i].size()){
            pmax=PartiVertex[i].size();
        }
        pave+=PartiVertex[i].size();
    }

    cout<<"Partition number: "<<partiNum<<" . Average vertex number: "<<pave/partiNum<<" ; Maximal vertex number: "<<pmax<<" ; Minimal vertex number: "<<pmin<< endl;
    cout<<"Overlay vertex number: "<<OverlayVertex.size()<<" . Non-boundary overlay vertex number: "<<noBoundayNum<<endl;

//    ifstream inFile(graphfile + ".queryST",ios::in);
//    if(!inFile){//if cannot open
//        cout << "Cannot open file " << graphfile + ".queryST"<<endl;
//        SameTreeQueryGen(tRoots, 1000);//generate same-tree queries
//    }
//    inFile.close();

//    exit(0);


    /// for AdjaCore
    NeighborsOverlay.assign(node_num,unordered_map<int,int>());
    NeighborsOverlayV.assign(node_num,vector<pair<vertex,int>>());
    cout<<"Generating overlay graph..."<<endl;
    int coreVNum=0;
    unsigned long long coreENum=0;
    int neiID,weight;
    for(int NodeID=0;NodeID<node_num;NodeID++){
        if(PartiTags[NodeID].first==-1){//core vertex
            ++coreVNum;
            for(int i=0;i<NeighborCon[NodeID].size();i++) {//for each neighbor
                neiID = NeighborCon[NodeID][i].first;
                weight = NeighborCon[NodeID][i].second.first;
                if(PartiTags[neiID].first != -1){
                    cout<<"Wrong! The contracted neighbor "<<neiID<<" of "<<NodeID<<" is not core vertex!!!"<<endl;
                    exit(1);
                }
                if(NeighborsOverlay[NodeID].find(neiID) == NeighborsOverlay[NodeID].end()){//if not found
                    coreENum+=2;
                    NeighborsOverlay[NodeID][neiID]=weight;//the root vertex is regarded as core vertex
                    NeighborsOverlay[neiID][NodeID]=weight;//the root vertex is regarded as core vertex
                    NeighborsOverlayV[NodeID].emplace_back(neiID, weight);
                    NeighborsOverlayV[neiID].emplace_back(NodeID, weight);
                }else{//if found
                    cout<<"Redundant edge! "<<NodeID<<" "<<neiID<<" "<<NeighborsOverlay[NodeID][neiID]<<" "<<weight<<endl; exit(1);
                }

            }
        }
    }
    cout<<"Core vertex number: "<<coreVNum<<"; core edge number: "<<coreENum<< endl;

    SuppPartiID.assign(node_num, map<int,map<int,int>>());//record the partition and its supportive value for a given interface edge
    SuppPartiIDReal.assign(node_num, map<int,pair<int,set<int>>>());//record the partitions that really support a given interface edge
    int ID1,ID2;
    for(int PID=0;PID<partiNum;PID++){
        int childrank=rank[PID];
        for(int i=0;i<BoundVertex[PID].size();++i){//for each root vertex
            ID1 = BoundVertex[PID][i];
            for(int j=i+1;j<BoundVertex[PID].size();++j){
                ID2 = BoundVertex[PID][j];
                if(NeighborsOverlay[ID1].find(ID2)!=NeighborsOverlay[ID1].end()) {//if found
                    if(ID1<ID2){
                        if(SuppPartiID[ID1].find(ID2)==SuppPartiID[ID1].end()){//if we cannot find ID2
                            SuppPartiID[ID1][ID2]=map<int,int>();
                        }
                        SuppPartiID[ID1][ID2].insert({PID,INF});
                        if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                            cout<<"No supportive vertex for shortcut sc("<<ID1<<","<<ID2<<") !"<<endl; exit(1);
                        }
                        for(auto it=SCconNodesMT[ID1][ID2].begin();it!=SCconNodesMT[ID1][ID2].end();++it){//for each supported vertex
                            if(NeighborsOverlay[ID1].find(ID2)==NeighborsOverlay[ID1].end()){
                                cout<<"Boundary edge e("<<ID1<<","<<ID2<<") does not exist!"<<endl; exit(1);
                            }//must exist
                            if(NeighborsOverlay[ID2].find(ID1)==NeighborsOverlay[ID2].end()){
                                cout<<"Boundary edge e("<<ID2<<","<<ID1<<") does not exist!"<<endl; exit(1);
                            }
                            if(PartiTags[it->first].first == PID){//if the contracted vertex belongs to PID
                                if(SuppPartiID[ID1][ID2][PID] > it->second){
                                    SuppPartiID[ID1][ID2][PID] = it->second;
                                }

                                if(SuppPartiIDReal[ID1].find(ID2) == SuppPartiIDReal[ID1].end()){//if not found
                                    SuppPartiIDReal[ID1][ID2]=make_pair(NeighborsOverlay[ID1][ID2],set<int>());
                                }
                                if(NeighborsOverlay[ID1][ID2] == it->second){
                                    SuppPartiIDReal[ID1][ID2].second.insert(PartiTags[it->first].first);
                                }
                            }

                        }

                    }
                    else{
                        if(SuppPartiID[ID2].find(ID1)==SuppPartiID[ID2].end()){//if we cannot find ID2
                            SuppPartiID[ID2][ID1]=map<int,int>();
                        }
                        SuppPartiID[ID2][ID1].insert({PID,INF});
                        if(SCconNodesMT[ID2].find(ID1)==SCconNodesMT[ID2].end()){//if not found
                            cout<<"No supportive vertex for shortcut sc("<<ID2<<","<<ID1<<") !"<<endl; exit(1);
                        }
                        for(auto it=SCconNodesMT[ID2][ID1].begin();it!=SCconNodesMT[ID2][ID1].end();++it){//for each supported vertex
                            if(NeighborsOverlay[ID1].find(ID2)==NeighborsOverlay[ID1].end()){
                                cout<<"Boundary edge e("<<ID1<<","<<ID2<<") does not exist!"<<endl; exit(1);
                            }//must exist
                            if(NeighborsOverlay[ID2].find(ID1)==NeighborsOverlay[ID2].end()){
                                cout<<"Boundary edge e("<<ID2<<","<<ID1<<") does not exist!"<<endl; exit(1);
                            }
                            if(PartiTags[it->first].first == PID){//if the contracted vertex belongs to PID
                                if(SuppPartiID[ID2][ID1][PID] > it->second){
                                    SuppPartiID[ID2][ID1][PID] = it->second;
                                }

                                if(SuppPartiIDReal[ID2].find(ID1) == SuppPartiIDReal[ID2].end()){
                                    SuppPartiIDReal[ID2][ID1]=make_pair(NeighborsOverlay[ID1][ID2],set<int>());
                                }
                                if(NeighborsOverlay[ID1][ID2] == it->second){
                                    SuppPartiIDReal[ID2][ID1].second.insert(PartiTags[it->first].first);
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    //clear useless variables
    existOverlay.clear();
    NeighborCon.clear();
}

void Graph::GetChildrenNumDFS(int p, vector<int>& ancestors){
    int ID=Tree[p].uniqueVertex;
    Tree[p].vAncestor=ancestors;
    if(!Tree[p].ch.empty()){
        ancestors.push_back(ID);
        for(int i=0;i<Tree[p].ch.size();++i){
            GetChildrenNumDFS(Tree[p].ch[i],ancestors);
        }
        ancestors.pop_back();
        int childID;
        for(auto it=Tree[p].ch.begin();it!=Tree[p].ch.end();++it){
            childID = Tree[*it].uniqueVertex;
            childNums[ID] += childNums[childID];
        }
    }
//    else{
//        childNums[ID]=1;
//    }


}

//Function of creating tree
void Graph::TDPCreateTreeAndParti(int pNum, double bRatioUpper, double bRatioLower){
    cout<<"Creating tree..."<<endl;
    //// Get tree
    vector<int> vecemp; //vecemp.clear();
    VidtoTNid.assign(node_num,vecemp);

    //rank.assign(HighestOrder+2,0);
    rank.clear();
    rank.assign(node_num,0);//the vector index of tree nodes, map from vertex to tree node
    int len=HighestOrder-1;
    len=node_num-1;
    heightMax=0;
    set<OrderCompMax> orderOverlayV; orderOverlayV.clear();
    Node root;//virtual root node
    int x = vNodeOrder[len];
    if(NeighborCon[x].empty()){
        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=x;
//        orderOverlayV.insert(OrderCompMax(x));
        len--;
    }else{
        root.uniqueVertex=-1;
    }
    root.height=1;
    Tree.push_back(root);


    int nn;
    for(;len>=0;len--){//check the vertices with order lower than HighestOrder
        x=vNodeOrder[len];
//        if(existOverlay[x]){
//            orderOverlayV.insert(OrderCompMax(x));
////            OverlayVertex.push_back(x);
////            cout<<"Wrong: should be out of core"<<endl; exit(1);
//        }
        Node nod;
        nod.vert=NeighborCon[x];//
        nod.uniqueVertex=x;
        if(treewidth<NeighborCon[x].size()+1){
            treewidth=NeighborCon[x].size()+1;
        }
        int pa=matchCore(x,NeighborCon[x],rank);

        //cout<<"pa "<<pa<<endl;

        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        /// for update
        nod.hdepth=Tree[pa].height+1;
        for(int i=0;i<NeighborCon[x].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[x][i].first;

            VidtoTNid[nn].push_back(Tree.size());//record the child tree node rank who has direct super edge to nn
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1){
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
            }
        }
        if(nod.height>heightMax){
            heightMax=nod.height;
        }

        rank[x]=Tree.size();//the position of tree, higher-order vertex has lower rank
        Tree.push_back(nod);
//        if(existOverlay[Tree[pa].uniqueVertex] && !existOverlay[x]){//if the father is in overlay while x is not in overlay
//            partiRoots.push_back(x);
//        }
    }
//    partiNum=partiRoots.size();
    cout<<"Tree node number: "<<Tree.size()<<" ; tree height: "<<heightMax<<" ; treewidth: "<<treewidth <<endl;

    //// Get partitions
    //get the children number of all tree nodes
    childNums.assign(node_num,1);
    vector<int> ancestors; ancestors.clear();
    GetChildrenNumDFS(0,ancestors);

//    cout<<Tree[0].uniqueVertex<<" "<<childNums[Tree[0].uniqueVertex]<<endl;
//    for(int i=0;i<Tree[0].ch.size();++i){
//        cout<<i<<": "<<Tree[Tree[0].ch[i]].uniqueVertex<<" "<<childNums[Tree[Tree[0].ch[i]].uniqueVertex]<<endl;
//    }

    int upperB, lowerB;
    upperB = bRatioUpper*node_num/pNum;
    lowerB = bRatioLower*node_num/pNum;
    cout<<"Imbalance ratio: "<<bRatioLower<<" "<<bRatioUpper<<". Lower bound of partition vertex number: "<<lowerB<<" ; upper bound: "<<upperB<<endl;

    MinimizeBoundaryNum(upperB, lowerB);

    BoundVertex.assign(partiNum,vector<int>());
    BoundVertexMap.assign(partiNum, unordered_map<int,int>());
    PartiTags.assign(node_num, make_pair(-1,set<int>()));
    PartiVertex.assign(partiNum, vector<int>());
    vector<set<OrderCompMax>> orderPartiV(partiNum,set<OrderCompMax>());

//    vector<int> tRoots;
    for(int k=0;k<partiRoots.size();k++){//Tree[0] is the core
        int id=partiRoots[k];
        //get boundary vertices
        for(int j=0;j<Tree[rank[id]].vert.size();++j){
            BoundVertex[k].push_back(Tree[rank[id]].vert[j].first);
            BoundVertexMap[k].insert({Tree[rank[id]].vert[j].first,j});
//            if(orderOverlayV.find(OrderCompMax(Tree[rank[id]].vert[j].first))==orderOverlayV.end()){
//                cout<<"Not found this boundary vertex! "<<Tree[rank[id]].vert[j].first<<" "<< NodeOrder[Tree[rank[id]].vert[j].first] <<endl; exit(1);
//            }
            if( PartiTags[Tree[rank[id]].vert[j].first].first!=-1){
                cout<<"!!! Strange! "<<Tree[rank[id]].vert[j].first<<" "<<PartiTags[Tree[rank[id]].vert[j].first].first<<endl; exit(1);
//                PartiTags[it->first].first=-1;
            }

            PartiTags[Tree[rank[id]].vert[j].first].second.insert(k);
        }
        //get partition vertex and graph
        makePartitionDFS(orderPartiV,k,id);
    }
    int noBoundayNum=0;
    orderOverlayV.clear();
    for(int i=0;i<node_num;++i){
        if(PartiTags[i].first==-1){
            orderOverlayV.insert(OrderCompMax(i));
        }
    }

    for(auto it=orderOverlayV.begin();it!=orderOverlayV.end();++it){
//        cout<<PartiTag[it->ID].first<<": "<<it->ID<<" "<<NodeOrder[it->ID]<<" "<<NeighborsOverlay[it->ID].size()<<endl;
        OverlayVertex.emplace_back(it->ID);
        if(!PartiTags[it->ID].second.empty()){
//            for(auto it2=PartiTags[it->ID].second.begin();it2!=PartiTags[it->ID].second.end();++it2){
//                BoundVertex[*it2].emplace_back(it->ID);
//            }
        }else{
            noBoundayNum++;
        }
    }


    int pmin=INF,pmax=0,pave=0;
    for(int i=0;i<partiNum;++i){
        if(pmin>PartiVertex[i].size()){
            pmin=PartiVertex[i].size();
        }
        if(pmax<PartiVertex[i].size()){
            pmax=PartiVertex[i].size();
        }
        pave+=PartiVertex[i].size();
    }

    cout<<"Partition number: "<<partiNum<<" . Average vertex number: "<<pave/partiNum<<" ; Maximal vertex number: "<<pmax<<" ; Minimal vertex number: "<<pmin<< endl;
    cout<<"Overlay vertex number: "<<OverlayVertex.size()<<" . Non-boundary overlay vertex number: "<<noBoundayNum<<endl;

//    ifstream inFile(graphfile + ".queryST",ios::in);
//    if(!inFile){//if cannot open
//        cout << "Cannot open file " << graphfile + ".queryST"<<endl;
//        SameTreeQueryGen(tRoots, 1000);//generate same-tree queries
//    }
//    inFile.close();

//    exit(0);


    /// for AdjaCore
    NeighborsOverlay.assign(node_num,unordered_map<int,int>());
    NeighborsOverlayV.assign(node_num,vector<pair<vertex,int>>());
//    cout<<"Generating overlay graph..."<<endl;
    int coreVNum=0;
    unsigned long long coreENum=0;
    int neiID,weight;
    for(int NodeID=0;NodeID<node_num;NodeID++){
        if(PartiTags[NodeID].first==-1){//core vertex
            ++coreVNum;
            for(int i=0;i<NeighborCon[NodeID].size();i++) {//for each neighbor
//                coreENum+=2;
                neiID = NeighborCon[NodeID][i].first;
                weight = NeighborCon[NodeID][i].second.first;
                if(PartiTags[neiID].first != -1){
                    cout<<"Wrong! The contracted neighbor "<<neiID<<" of "<<NodeID<<" is not core vertex!!!"<<endl;
                    exit(1);
                }
                if(NeighborsOverlay[NodeID].find(neiID) == NeighborsOverlay[NodeID].end()){//if not found
                    coreENum+=2;
                    NeighborsOverlay[NodeID][neiID]=weight;//the root vertex is regarded as core vertex
                    NeighborsOverlay[neiID][NodeID]=weight;//the root vertex is regarded as core vertex
                    NeighborsOverlayV[NodeID].emplace_back(neiID, weight);
                    NeighborsOverlayV[neiID].emplace_back(NodeID, weight);
                }else{//if found
                    cout<<"Redundant edge! "<<NodeID<<" "<<neiID<<" "<<NeighborsOverlay[NodeID][neiID]<<" "<<weight<<endl; exit(1);
                }

            }
        }
    }
    cout<<"Core vertex number: "<<coreVNum<<"; core edge number: "<<coreENum<< endl;

    SuppPartiID.assign(node_num, map<int,map<int,int>>());//record the partition and its supportive value for a given interface edge
    SuppPartiIDReal.assign(node_num, map<int,pair<int,set<int>>>());//record the partitions that really support a given interface edge
    int ID1,ID2;
    for(int PID=0;PID<partiNum;PID++){
        int childrank=rank[PID];
        for(int i=0;i<BoundVertex[PID].size();++i){//for each root vertex
            ID1 = BoundVertex[PID][i];
            for(int j=i+1;j<BoundVertex[PID].size();++j){
                ID2 = BoundVertex[PID][j];
                if(NeighborsOverlay[ID1].find(ID2)!=NeighborsOverlay[ID1].end()) {//if found
                    if(ID1<ID2){
                        if(SuppPartiID[ID1].find(ID2)==SuppPartiID[ID1].end()){//if we cannot find ID2
                            SuppPartiID[ID1][ID2]=map<int,int>();
                        }
                        SuppPartiID[ID1][ID2].insert({PID,INF});
                        if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                            cout<<"No supportive vertex for shortcut sc("<<ID1<<","<<ID2<<") !"<<endl; exit(1);
                        }
                        for(auto it=SCconNodesMT[ID1][ID2].begin();it!=SCconNodesMT[ID1][ID2].end();++it){//for each supported vertex
                            if(NeighborsOverlay[ID1].find(ID2)==NeighborsOverlay[ID1].end()){
                                cout<<"Boundary edge e("<<ID1<<","<<ID2<<") does not exist!"<<endl; exit(1);
                            }//must exist
                            if(NeighborsOverlay[ID2].find(ID1)==NeighborsOverlay[ID2].end()){
                                cout<<"Boundary edge e("<<ID2<<","<<ID1<<") does not exist!"<<endl; exit(1);
                            }
                            if(PartiTags[it->first].first == PID){//if the contracted vertex belongs to PID
                                if(SuppPartiID[ID1][ID2][PID] > it->second){
                                    SuppPartiID[ID1][ID2][PID] = it->second;
                                }

                                if(SuppPartiIDReal[ID1].find(ID2) == SuppPartiIDReal[ID1].end()){//if not found
                                    SuppPartiIDReal[ID1][ID2]=make_pair(NeighborsOverlay[ID1][ID2],set<int>());
                                }
                                if(NeighborsOverlay[ID1][ID2] == it->second){
                                    SuppPartiIDReal[ID1][ID2].second.insert(PartiTags[it->first].first);
                                }
                            }

                        }

                    }
                    else{
                        if(SuppPartiID[ID2].find(ID1)==SuppPartiID[ID2].end()){//if we cannot find ID2
                            SuppPartiID[ID2][ID1]=map<int,int>();
                        }
                        SuppPartiID[ID2][ID1].insert({PID,INF});
                        if(SCconNodesMT[ID2].find(ID1)==SCconNodesMT[ID2].end()){//if not found
                            cout<<"No supportive vertex for shortcut sc("<<ID2<<","<<ID1<<") !"<<endl; exit(1);
                        }
                        for(auto it=SCconNodesMT[ID2][ID1].begin();it!=SCconNodesMT[ID2][ID1].end();++it){//for each supported vertex
                            if(NeighborsOverlay[ID1].find(ID2)==NeighborsOverlay[ID1].end()){
                                cout<<"Boundary edge e("<<ID1<<","<<ID2<<") does not exist!"<<endl; exit(1);
                            }//must exist
                            if(NeighborsOverlay[ID2].find(ID1)==NeighborsOverlay[ID2].end()){
                                cout<<"Boundary edge e("<<ID2<<","<<ID1<<") does not exist!"<<endl; exit(1);
                            }
                            if(PartiTags[it->first].first == PID){//if the contracted vertex belongs to PID
                                if(SuppPartiID[ID2][ID1][PID] > it->second){
                                    SuppPartiID[ID2][ID1][PID] = it->second;
                                }

                                if(SuppPartiIDReal[ID2].find(ID1) == SuppPartiIDReal[ID2].end()){
                                    SuppPartiIDReal[ID2][ID1]=make_pair(NeighborsOverlay[ID1][ID2],set<int>());
                                }
                                if(NeighborsOverlay[ID1][ID2] == it->second){
                                    SuppPartiIDReal[ID2][ID1].second.insert(PartiTags[it->first].first);
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    //clear useless variables
    existOverlay.clear();
    NeighborCon.clear();
    NeighborsOverlay.clear();
    NeighborsOverlayV.clear();
}

void Graph::MinimizeBoundaryNum(int upperB, int lowerB){
    vector<pair<int,int>> candidates;
    map<int,vector<int>> candidatesMap;//from tree height to vertex IDs
    GetPartiRootCandidates(0,upperB,lowerB,candidates,candidatesMap);
    if(candidates.empty()){
        cout<<"Bandwidth is too small! Please enlarge it!"<<endl; exit(0);
    }
    cout<<"Partition root candidate size: "<<candidates.size()<<" ; tree height number: "<<candidatesMap.size() <<endl;
    vector<pair<int,vector<int>>> branches;//<root,children>
    int height,heightMin,ID1,ID2;
    heightMin=heightMax;
    for(auto it=candidatesMap.begin();it!=candidatesMap.end();++it){
        height=it->first;
        if(heightMin>height){
            heightMin=height;
        }
        for(auto it2=it->second.begin();it2!=it->second.end();++it2){
            ID1=*it2;
            if(height == heightMin){
                vector<int> tempV;
                tempV.push_back(ID1);
                branches.emplace_back(ID1,tempV);
            }else{
                bool ifAdd=true;
                for(int i=0;i<branches.size();++i){
                    ID2=branches[i].first;
                    int tempHeight=Tree[rank[ID2]].height-1;
                    if(Tree[rank[ID1]].vAncestor[tempHeight] == ID2){
                        ifAdd=false;
                        branches[i].second.push_back(ID1);
                        break;
                    }
//                    else{
//                        cout<<ID1<<"("<<Tree[rank[ID1]].height<<") "<<Tree[rank[ID1]].vAncestor[tempHeight]<<" "<< ID2<<"("<<tempHeight+1<<","<<Tree[rank[ID2]].height<<")"<<endl;
//                    }
                }
                if(ifAdd){
                    vector<int> tempV;
                    tempV.push_back(ID1);
                    branches.emplace_back(ID1,tempV);
                }
            }
        }
    }
    cout<<"Branch number: "<<branches.size()<<endl;
    partiRoots.clear();
    for(int i=0;i<branches.size();++i){
        // Goal 1: minimize the overlay vertex number
        partiRoots.push_back(branches[i].first);
        // Goal 2: minimize the boundary vertex number
//        pair<int,int> finalCut;
//        finalCut.first=-1 , finalCut.second=INF;
//        for(auto it=branches[i].second.begin();it!=branches[i].second.end();++it){
//            if(finalCut.second>Tree[rank[*it]].vert.size()){
//                finalCut.first=*it;
//                finalCut.second=Tree[rank[*it]].vert.size();
//            }
//        }
//        if(finalCut.first!=-1){
//            partiRoots.push_back(finalCut.first);
//        }else{
//            cout<<"Wrong! "<<endl; exit(1);
//        }

    }
    partiNum=partiRoots.size();
}

void Graph::GetPartiRootCandidates(int p, int upperB, int lowerB, vector<pair<int,int>>& candidates, map<int,vector<int>>& candidatesMap){
    int ID=Tree[p].uniqueVertex;
    if(childNums[ID]>=lowerB && childNums[ID]<=upperB && Tree[p].vert.size()<=bandWidth){
        candidates.emplace_back(ID,Tree[p].height);
        if(candidatesMap.find(Tree[p].height) == candidatesMap.end()){
            candidatesMap.insert({Tree[p].height,vector<int>()});
            candidatesMap[Tree[p].height].push_back(ID);
        }else{//if found
            candidatesMap[Tree[p].height].push_back(ID);
        }
    }
    if(!Tree[p].ch.empty()){
        for(int i=0;i<Tree[p].ch.size();++i){
            GetPartiRootCandidates(Tree[p].ch[i], upperB, lowerB, candidates, candidatesMap);
        }
    }

}

void Graph::makePartitionDFS(vector<set<OrderCompMax>>& orderPartiV, int pid, int vid){
    orderPartiV[pid].insert(OrderCompMax(vid));
    PartiVertex[pid].push_back(vid);
    if(PartiTags[vid].first==-1){
        PartiTags[vid].first=pid;
    }else{
        cout<<"Wrong! "<<pid<<": "<<vid<<" already belongs to partition "<<PartiTags[vid].first<<endl; exit(1);
    }

    for(auto it=Tree[rank[vid]].ch.begin();it!=Tree[rank[vid]].ch.end();++it){
        makePartitionDFS(orderPartiV , pid,Tree[*it].uniqueVertex);
    }
}


void Graph::RefreshExtensionLabelPartiV(vector<int>& p, bool ifTopDown){
    if(ifTopDown){
        for(int i=0;i<p.size();++i){
            RefreshExtensionLabelPartiTopDown(p[i]);
        }
    }
    else{
        for(int i=0;i<p.size();++i){
            RefreshExtensionLabelParti(p[i]);
        }
    }
}

void Graph::RefreshExtensionLabelTopDownV(vector<int>& p){
    for(int i=0;i<p.size();++i){
        RefreshExtensionLabelTopDown(p[i]);
    }
}

void Graph::RefreshExtensionLabelTopDown(int ID){
    int PID=PartiTag[ID].first;
    unordered_set<int> children;//vertex id of children
    children.clear();
    int p=rankExt[ID];
//    int pa=list[list.size()-1];
    int tempRank;
    //nested loop
    int ancestor, neighbor, d1, d2;
    int dis=INF, vid;
    if(PartiTag[ID].second){//if boundary vertex
        TreeExt[p].pos=Tree[rank[ID]].pos;
        for(int i=0;i<TreeExt[p].vAncestor.size();++i){
            ancestor = TreeExt[p].vAncestor[i];
            dis=INF;
//            for(int j=0;j<Tree[rank[ID]].vert.size();++j){
//                neighbor = Tree[rank[ID]].vert[j].first;
//                d1=Tree[rank[ID]].vert[j].second.first;
////                if(TreeExt[rankExt[neighbor]].dis.size()<list.size()){
////                    cout<<ID<<": "<<neighbor<<"("<<NodeOrder[neighbor]<<","<<PartiTag[neighbor].first<<","<<PartiTag[neighbor].second<<") "<<TreeExt[rankExt[neighbor]].dis.size()<<" "<<list.size()<<endl;
////                }
//                if(TreeExt[rankExt[neighbor]].height<TreeExt[rankExt[ancestor]].height){//if neighbor is ancestor of ancestor
////                    d2= QueryCore(neighbor,ancestor);
//                    d2=TreeExt[rankExt[ancestor]].dis[TreeExt[rankExt[neighbor]].height-1];;
//                }else{
//                    d2=TreeExt[rankExt[neighbor]].dis[i];
//                }
//
//                if(dis>d1+d2){
//                    dis=d1+d2;
//                }
//            }
            dis= QueryCore(ID,ancestor);
            TreeExt[p].dis[i]=dis;
//            int disD= Dijkstra(ancestor,ID,Neighbor);
//            if(disD!=dis){
//                cout<<"Boundary. Incorrect! "<<ID<<" "<<ancestor<<" "<<dis<<" "<<disD<<endl;
//            }
        }
    }
    else{//if non-boundary vertex
        for(int i=0;i<TreeExt[p].vAncestor.size();++i){
            ancestor = TreeExt[p].vAncestor[i];
            dis=INF;
            for(int j=0;j<Trees[PID][ranks[PID][IDMap[ID]]].vert.size();++j){
                neighbor = Trees[PID][ranks[PID][IDMap[ID]]].vert[j].first;
                d1=Trees[PID][ranks[PID][IDMap[ID]]].vert[j].second.first;
//                if(TreeExt[rankExt[neighbor]].dis.size()<list.size()){
//                    cout<<ID<<": "<<neighbor<<"("<<NodeOrder[neighbor]<<","<<PartiTag[neighbor].first<<","<<PartiTag[neighbor].second<<") "<<TreeExt[rankExt[neighbor]].dis.size()<<" "<<list.size()<<endl;
//                }

                if(TreeExt[rankExt[neighbor]].height<TreeExt[rankExt[ancestor]].height){//if neighbor is ancestor of ancestor
                    d2= TreeExt[rankExt[ancestor]].dis[TreeExt[rankExt[neighbor]].height-1];
                    if(dis>d1+d2){
                        dis=d1+d2;
                        vid=ancestor;
                    }
                }else{
                    d2=TreeExt[rankExt[neighbor]].dis[i];
                    if(dis>d1+d2){
                        dis=d1+d2;
                        vid=neighbor;
                    }
                }
//                d2=TreeExt[rankExt[neighbor]].dis[i];

//                if(ID==25991){
//                    cout<<ID<<": "<<neighbor<<" "<<ancestor<<" "<<d1<<" "<<d2<<" "<<dis<<endl;
//                }
            }
            TreeExt[p].dis[i]=dis;
            TreeExt[p].cnt[i]=vid;
//            int disD= Dijkstra(ancestor,ID,Neighbor);
//            if(disD!=dis){
//                cout<<"Non-boundary. Incorrect! "<<ID<<" "<<ancestor<<" "<<dis<<" "<<disD<<endl;
//            }
        }
    }

    for(int i=0;i<TreeExt[p].ch.size();i++){
        RefreshExtensionLabelTopDown(TreeExt[TreeExt[p].ch[i]].uniqueVertex);
    }
}

//top-down version
void Graph::RefreshExtensionLabelPartiTopDown(int pid){
    int ID,ancestor;

    ID=PartiVertex[pid][0];
    RefreshExtensionLabelTopDown(ID);
}
//boundary version
void Graph::RefreshExtensionLabelParti(int pid){
    int ID,ancestor;
    pair<int,int> result;
    bool flagBoundary=false;
    for(int i=0;i<BoundVertex[pid].size();++i){
        ID=BoundVertex[pid][i];
        if(vUpdated[ID]){
            flagBoundary = true;
        }
    }

    if(ifIncrease){//increase update
        for(int i=0;i<PartiVertex[pid].size();++i){
            ID = PartiVertex[pid][i];
            if(PartiTag[ID].second){
                continue;
            }
            for(int j=0;j<TreeExt[rankExt[ID]].dis.size()-1;++j){
                ancestor = TreeExt[rankExt[ID]].vAncestor[j];
                if(!PartiTag[ancestor].second){//if not boundary vertex, new
                    continue;
                }
                if(PartiTag[ancestor].first!=pid){//ancestor of other partitions
                    if(vUpdated[TreeExt[rankExt[ID]].cnt[j]] || vUpdated[ID]){//if the passing boundary's label has changed or this vertex's label has changed
                        //update
                        result=ExtensionLabelCompute(pid, ID, ancestor);
                        if(TreeExt[rankExt[ID]].dis[j]!=result.first){
//                        cout<<"Different parti: "<<ID<<" "<<ancestor<<": "<<TreeExt[rankExt[ID]].dis[j]<<"("<<TreeExt[rankExt[ID]].cnt[j]<<") "<<result.first<<"("<<result.second<<")"<<endl;
                            TreeExt[rankExt[ID]].dis[j]=result.first;//the distance to ancestors
                            TreeExt[rankExt[ID]].cnt[j]=result.second;//the corresponding boundary vertex for the label
                        }

                    }
                }
                else{//ancestor of this partition
                    if(vUpdated[ID]) {//if the passing boundary's label has changed or this vertex's label has changed
                        //update
                        result=ExtensionLabelCompute(pid, ID, ancestor);
                        if(TreeExt[rankExt[ID]].dis[j]!=result.first){
//                        cout<<"Same parti: "<<ID<<" "<<ancestor<<": "<<TreeExt[rankExt[ID]].dis[j]<<"("<<TreeExt[rankExt[ID]].cnt[j]<<") "<<result.first<<"("<<result.second<<")"<<endl;
                            TreeExt[rankExt[ID]].dis[j]=result.first;//the distance to ancestors
                            TreeExt[rankExt[ID]].cnt[j]=result.second;//the corresponding boundary vertex for the label
                        }

                    }
                }
            }
        }
    }
    else{//decrease update
        for(int i=0;i<PartiVertex[pid].size();++i){
            ID = PartiVertex[pid][i];
            if(PartiTag[ID].second){
                continue;
            }
            for(int j=0;j<TreeExt[rankExt[ID]].dis.size()-1;++j){
                ancestor = TreeExt[rankExt[ID]].vAncestor[j];
                if(!PartiTag[ancestor].second){//if not boundary vertex
                    continue;
                }
                if(PartiTag[ancestor].first!=pid){//ancestor of other partitions
                    if(flagBoundary || vUpdated[ID]){//if the passing boundary's label has changed or this vertex's label has changed
                        //update
                        result=ExtensionLabelCompute(pid, ID, ancestor);
                        if(TreeExt[rankExt[ID]].dis[j]!=result.first){
//                            cout<<"Different parti: "<<ID<<" "<<ancestor<<": "<<TreeExt[rankExt[ID]].dis[j]<<"("<<TreeExt[rankExt[ID]].cnt[j]<<") "<<result.first<<"("<<result.second<<")"<<endl;
                            TreeExt[rankExt[ID]].dis[j]=result.first;//the distance to ancestors
                            TreeExt[rankExt[ID]].cnt[j]=result.second;//the corresponding boundary vertex for the label
                        }
                    }
                }
                else{//ancestor of this partition
                    if(vUpdated[ID]) {//if the passing boundary's label has changed or this vertex's label has changed
                        //update
                        result=ExtensionLabelCompute(pid, ID, ancestor);
                        if(TreeExt[rankExt[ID]].dis[j]!=result.first){
//                            cout<<"Same parti: "<<ID<<" "<<ancestor<<": "<<TreeExt[rankExt[ID]].dis[j]<<"("<<TreeExt[rankExt[ID]].cnt[j]<<") "<<result.first<<"("<<result.second<<")"<<endl;
                            TreeExt[rankExt[ID]].dis[j]=result.first;//the distance to ancestors
                            TreeExt[rankExt[ID]].cnt[j]=result.second;//the corresponding boundary vertex for the label
                        }

                    }
                }
            }
        }

    }

}

void Graph::ConstructPH2H_PartiV(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    int PID;
    for(int i=0;i<P.size();++i){
        PID=P[i];
        ConstructPH2H_Parti(PID, Trees, ranks, SCconNodesMTP, VidtoTNidP, heightMaxs, toRMQs, RMQIndexs);
    }
}

void Graph::ConstructPH2H_PartiVCH(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    int pid;
    for(int i=0;i<P.size();++i){
        pid=P[i];
        H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
//        ConstructPH2H_Parti(PID, Trees, ranks, SCconNodesMTP, VidtoTNidP, heightMaxs, toRMQs, RMQIndexs, ifLabelC);
    }
}

void Graph::ConstructPH2H_Parti(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    //Create tree for partition
    H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
    //Create LCA index
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    //Create labels for partition
    H2HCreateIndex_Parti(pid, Trees[pid], ranks[pid]);
}
//void Graph::ConstructPH2H_PartiNoLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
//    //Create tree for partition
//    H2HCreateTree_Parti(pid, Trees[pid], ranks[pid], SCconNodesMTP, VidtoTNidP, heightMaxs);
//}
void Graph::ConstructPH2H_PartiLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs){
    //Create LCA index
    makeRMQCoreP(pid, toRMQs, RMQIndexs, Trees);
    //Create labels for partition
    H2HCreateIndex_Parti(pid, Trees[pid], ranks[pid]);
}

//Function of Creating tree for partition
void Graph::H2HCreateTree_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs) {
//    cout<<"Partition "<<pid<<": boundary number is "<<BoundVertex[pid].size()<<endl;
    /// Contraction
    int degree;
    int ID, ID1, ID2;
    unordered_map<vertex,bool> existCore; existCore.clear();
    for(auto it=PartiVertex[pid].begin();it!=PartiVertex[pid].end();++it){
        existCore.insert({*it,true});
    }

//    unordered_map<vertex,int> rankP; rankP.clear();
    for(int id=PartiVertex[pid].size()-1;id>=0;--id){
        ID = PartiVertex[pid][id];
//        rankP.insert({ID,-1});
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[ID].begin();it!=E[ID].end();it++){
            if(existCore.find(it->first)==existCore.end()){// not found
                cout<<"Wrong neighbor! "<<ID<<"("<<PartiTag[ID].first<<") "<< it->first<<"("<<PartiTag[it->first].first<<")"<<endl; exit(1);
            }
            else{// if found
                if(existCore[it->first]){
                    Neigh.emplace_back(*it);
                }else{
                    cout<<"Not in core!"<<it->first<<endl; exit(1);
                }
            }
        }
        NeighborCon[ID].assign(Neigh.begin(),Neigh.end());
//        cout<<ID<<" "<<NeighborCon[ID].size()<<" "<<PartiTag[ID].second<<endl;

        existCore[ID]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteECore(ID,y);//delete ID from y's adjacency list
        }
        //add all-pair neighbors
        for(int i=0;i<Neigh.size();i++){
            ID1=Neigh[i].first;
            for(int j=i+1;j<Neigh.size();j++){
                ID2=Neigh[j].first;
                insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                /// For TD update
                if(ID1<ID2){
                    SCconNodesMTP[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                }
                else{
                    SCconNodesMTP[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
                }

            }
        }


    }

    /// Create Tree
    ID=PartiVertex[pid][0];
    Node root;//virtual root node
    if(NeighborCon[ID].empty()){
//        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=ID;
    }else{
        cout<<"Wrong!"<<endl; exit(1);
    }
    root.height=1;
    TreeP.push_back(root);
//    cout<<"0 "<<ID<<" "<<IDMap[ID]<<endl;
    rankP[IDMap[ID]] = 0;
//    rankP[ID] = 0;

    for(int id=1;id<PartiVertex[pid].size();++id){
        ID = PartiVertex[pid][id];
//        cout<<id<<" "<<ID<<" "<<IDMap[ID]<<endl;
        int nn;
        if(existCore[ID]){
            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[ID];//
        nod.uniqueVertex=ID;
        int pa=matchCoreParti(ID,NeighborCon[ID], rankP);

//        cout<<"pa "<<pa<<" "<<TreeP[pa].height<<endl;

        TreeP[pa].ch.push_back(TreeP.size());
        nod.pa=pa;
        nod.height=TreeP[pa].height+1;
        /// for update
        nod.hdepth=TreeP[pa].height+1;
        for(int i=0;i<NeighborCon[ID].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[ID][i].first;
            if(PartiTag[nn].first != pid){
                cout<<"Wrong nn! "<<PartiTag[nn].first <<" "<< pid<<endl; exit(1);
            }
            VidtoTNidP[nn].emplace_back(TreeP.size());
//            if(rankP.find(nn)==rankP.end()){
//                cout<<"Not found "<<nn<<" in rank!"<<endl; exit(1);
//            }
            if(TreeP[rankP[IDMap[nn]]].hdepth<TreeP[pa].height+1){
                TreeP[rankP[IDMap[nn]]].hdepth=TreeP[pa].height+1;
            }

        }
        if(nod.height>heightMaxs[pid]){
            heightMaxs[pid]=nod.height;
        }

//        if(rankP.find(ID)==rankP.end()){
//            cout<<"Not found "<<ID<<" in rank!"<<endl; exit(1);
//        }
        rankP[IDMap[ID]]=TreeP.size();//the position of tree, higher-order vertex has lower rank
        TreeP.push_back(nod);

    }

//    cout<<pid<<"'s tree node number: "<<TreeP.size()<<" ; tree height: "<< heightMaxs[pid]<<endl;
}
//Function of tree-label index construction for partition
void Graph::H2HCreateIndex_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP){
//    cout<<"Computing Tree Label for partition "<<pid<<endl;
    //initialize
    vector<int> list; //list.clear();
    list.push_back(TreeP[0].uniqueVertex);
    TreeP[0].pos.clear();
    TreeP[0].pos.push_back(0);

    for(int i=0;i<TreeP[0].ch.size();i++){
        makeTreeIndexDFSP(TreeP[0].ch[i],list,TreeP, rankP);
    }

}

//Function of Creating tree for overlay graph
void Graph::H2HCreateTree_Overlay() {
    //for H2H update
    SCconNodesMT.assign(node_num, map<int, vector<pair<int,int>>>());
    NeighborCon.assign(node_num, vector<pair<int,pair<int,int>>>());
//    _DD_.assign(node_num,0); //_DD2_.assign(node_num,0);
    DD.assign(node_num,0); //DD2.assign(node_num,0);
    VidtoTNid.assign(node_num,vector<int>());
    //initialize E
    map<int,pair<int,int>> m;
    E.assign(node_num,m);
    for(int i=0;i<NeighborsOverlay.size();i++){
        if(!NeighborsOverlay[i].empty()){
            for(auto it=NeighborsOverlay[i].begin();it!=NeighborsOverlay[i].end();++it){
                E[i].insert(make_pair(it->first,make_pair(it->second,1)));
            }
        }
    }
    /// Contraction
    int degree;
    int ID, ID1, ID2;
    vector<bool> existCore(node_num,true);
    rank.assign(node_num,-1);

    bool flagAdd = true;
    for(int id=OverlayVertex.size()-1;id>=0;--id){
        ID = OverlayVertex[id];
//        cout<<ID<<" "<<NodeOrder[ID]<<endl;
        vector<pair<int,pair<int,int>>> Neigh; //Neigh.clear();
        for(auto it=E[ID].begin();it!=E[ID].end();it++){
            if(existCore[it->first]){
                Neigh.emplace_back(*it);
            }else{
                cout<<"Not in core!"<<it->first<<endl; exit(1);
            }
        }
        NeighborCon[ID].assign(Neigh.begin(),Neigh.end());

        existCore[ID]=false;
        //delete the star
        for(int i=0;i<Neigh.size();i++){
            int y=Neigh[i].first;
            deleteECore(ID,y);//delete ID from y's adjacency list
        }

        if(Neigh.size()<=100){
            //single thread
            for(int i=0;i<Neigh.size();i++){
                ID1=Neigh[i].first;
                for(int j=i+1;j<Neigh.size();j++){
                    ID2=Neigh[j].first;
                    insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                    /// For TD update
                    if(ID1<ID2){
                        SCconNodesMT[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
                    }else{
                        SCconNodesMT[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
                    }
                }
            }
        }else{
//            cout<<"Multiple thread for contraction. "<<ID<<" "<<Neigh.size()<<endl;
            //multiple thread
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
                    thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, ID));
                }
                thread.join_all();
            }else{
                boost::thread_group thread;
                for(int i=0;i<Neigh.size();i++){
                    pair<int,int> p;
                    p.first=i; p.second=(i+1);
                    thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, ID));
                }
                thread.join_all();
            }
        }

        //add all-pair neighbors
//        for(int i=0;i<Neigh.size();i++){
//            ID1=Neigh[i].first;
//            for(int j=i+1;j<Neigh.size();j++){
//                ID2=Neigh[j].first;
//                insertECore(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
//                /// For TD update
//                SCconNodesMT[ID1][ID2].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
//                SCconNodesMT[ID2][ID1].emplace_back(ID,Neigh[i].second.first+Neigh[j].second.first);
//            }
//        }


    }
//    cout<<"Flag 1"<<endl;
    /// Create Tree
    ID=OverlayVertex[0];
    Node root;//virtual root node
    if(NeighborCon[ID].empty()){
//        cout<<"There exist non-virtual root!"<<endl;
        root.uniqueVertex=ID;
    }else{
        cout<<"Wrong!"<<endl; exit(1);
    }
    root.height=1;
    Tree.push_back(root);
//    rankP[IDMap[ID]] = 0;
    rank[ID] = 0;
    heightMax=0;

    for(int id=1;id<OverlayVertex.size();++id){
        ID = OverlayVertex[id];
//        cout<<ID<<" "<<NodeOrder[ID]<<endl;
        int nn;
        if(existCore[ID]){
            cout<<"Wrong: should be out of core"<<endl; exit(1);
        }
        Node nod;
        nod.vert=NeighborCon[ID];//
        nod.uniqueVertex=ID;
        int pa=matchCore(ID,NeighborCon[ID], rank);

        //cout<<"pa "<<pa<<endl;

        Tree[pa].ch.push_back(Tree.size());
        nod.pa=pa;
        nod.height=Tree[pa].height+1;
        /// for update
        nod.hdepth=Tree[pa].height+1;
        for(int i=0;i<NeighborCon[ID].size();i++){//for the neighbors which have higher order
            nn=NeighborCon[ID][i].first;
            VidtoTNid[nn].emplace_back(Tree.size());
            if(Tree[rank[nn]].hdepth<Tree[pa].height+1){
                Tree[rank[nn]].hdepth=Tree[pa].height+1;
            }

        }
        if(nod.height>heightMax){
            heightMax=nod.height;
        }

        rank[ID]=Tree.size();//the position of tree, higher-order vertex has lower rank
        Tree.push_back(nod);

    }

    /// LCA index
    makeRMQCore();//build LCA index

    cout<<"Overlay graph! Tree node number: "<<Tree.size()<<" ; Tree height: "<< heightMax<<endl;
}

void Graph::NeighborComorderH2H(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        w1=Neighvec[k].second.first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
            w2=Neighvec[h].second.first;
            if(ID1==ID2){
                continue;
            }
            insertEMTorder(ID1, ID2, w1+w2);
            if(ID1<ID2){
                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                    SCconNodesMT[ID1].insert({ID2,vector<pair<int,int>>()});
                }
                SCconNodesMT[ID1][ID2].emplace_back(x,w1+w2);
            }

        }
    }
//    sm->notify();
}

void Graph::NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x){
//    sm->wait();
    int ID1, w1;
    int ID2, w2;
    for(int k=p.first;k<p.second;k++){
        ID1=Neighvec[k].first;
        w1=Neighvec[k].second.first;
        for(int h=0;h<Neighvec.size();h++){
            ID2=Neighvec[h].first;
            w2=Neighvec[h].second.first;
            if(ID1==ID2){
                continue;
            }
            insertECoreMT(ID1,ID2,w1+w2);
            /// For TD update
            if(ID1<ID2){
                SCconNodesMT[ID1][ID2].emplace_back(x,w1+w2);//supportive vertex, no direction, may contain the supportive vertices for shortcuts between interface vertices
            }
        }
    }
//    sm->notify();
}

//Function of tree-label index construction for partition
void Graph::H2HCreateIndex_Overlay(){
    //initialize
    vector<int> list; //list.clear();
    list.push_back(Tree[0].uniqueVertex);
    Tree[0].pos.clear();
    Tree[0].pos.push_back(0);

    for(int i=0;i<Tree[0].ch.size();i++){
        makeIndexDFS(Tree[0].ch[i],list);
    }

}



//function of computing the H2H label of peripheries: original version
void Graph::makeTreeIndexDFSP(int p, vector<int>& list,  vector<Node>& TreeP, vector<int>& rankP){
//initialize
//    cout<<"Map "<<p<<" "<<IDMap[p]<<endl;
//    p=IDMap[p];
    int NeiNum=TreeP[p].vert.size();
    TreeP[p].pos.assign(NeiNum+1,0);
    TreeP[p].dis.assign(list.size(),INF);
    TreeP[p].cnt.assign(list.size(),0);
    TreeP[p].FN.assign(list.size(),true);
    TreeP[p].vAncestor=list;

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(TreeP[p].vert[i].first==list[j]){
                TreeP[p].pos[i]=j;//record the position of neighbors
                TreeP[p].dis[j]=TreeP[p].vert[i].second.first;
                TreeP[p].cnt[j]=1;
                break;
            }
        }
    }
    TreeP[p].pos[NeiNum]=list.size();


    //dis
    for(int i=0;i<NeiNum;i++){
        int x=TreeP[p].vert[i].first;
        int disvb=TreeP[p].vert[i].second.first;
        int k=TreeP[p].pos[i];//the kth ancestor is x

        for(int j=0;j<list.size();j++){//check the distance to the j-th ancestor could be updated by neighbors, including the valley path and peak path
            int y=list[j];//the jth ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//x is the ancestor of y, peak path
                    z=TreeP[rankP[IDMap[y]]].dis[k];
                else if(k>j)//y is the ancestor of x, valley path
                    z=TreeP[rankP[IDMap[x]]].dis[j];

                if(TreeP[p].dis[j]>z+disvb){
                    TreeP[p].dis[j]=z+disvb;
                    TreeP[p].FN[j]=false;
                    TreeP[p].cnt[j]=1;
                }else if(TreeP[p].dis[j]==z+disvb){
                    TreeP[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    list.push_back(TreeP[p].uniqueVertex);
    for(int i=0;i<TreeP[p].ch.size();i++){
        makeTreeIndexDFSP(TreeP[p].ch[i],list, TreeP, rankP);
    }
    list.pop_back();
}

void Graph::makeIndexDFS(int p, vector<int>& list){
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,0);
    Tree[p].dis.assign(list.size(),INF);
    Tree[p].cnt.assign(list.size(),0);
    Tree[p].FN.assign(list.size(),true);
    Tree[p].vAncestor=list;

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(Tree[p].vert[i].first==list[j]){
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();

    //dis
    for(int i=0;i<NeiNum;i++){
        int x=Tree[p].vert[i].first;
        int disvb=Tree[p].vert[i].second.first;
        int k=Tree[p].pos[i];//the kth ancestor is x

        for(int j=0;j<list.size();j++){//check the distance to the j-th ancestor could be updated by neighbors, including the valley path and peak path
            int y=list[j];//the jth ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//x is the ancestor of y, peak path
                    z=Tree[rank[y]].dis[k];
                else if(k>j)//y is the ancestor of x, valley path
                    z=Tree[rank[x]].dis[j];

                if(Tree[p].dis[j]>z+disvb){
                    Tree[p].dis[j]=z+disvb;
                    Tree[p].FN[j]=false;
                    Tree[p].cnt[j]=1;
                }else if(Tree[p].dis[j]==z+disvb){
                    Tree[p].cnt[j]+=1;
                }
            }
        }
    }

    //nested loop
    list.push_back(Tree[p].uniqueVertex);
    for(int i=0;i<Tree[p].ch.size();i++){
        makeIndexDFS(Tree[p].ch[i],list);
    }
    list.pop_back();
}



//function of erasing edge (u,v), i.e., erase u from v's adjacency list.
void Graph::deleteECore(int u,int v){
//	if(Emap[u].find(v)!=Emap[u].end()){
//		Emap[u].erase(Emap[u].find(v));
//		DD[u]--;
//	}

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
    }
}
//function of inserting edge (u,v)
void Graph::insertECore(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){//if not found
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//		DD2[u]++;
    }
    else{//if found
        if(E[u][v].first>w)
            E[u][v]= make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

    if(E[v].find(u)==E[v].end()){
        E[v].insert(make_pair(u,make_pair(w,1)));
        DD[v]++;
//		DD2[v]++;
    }
    else{
        if(E[v][u].first>w)
            E[v][u]=make_pair(w,1);
        else if(E[v][u].first==w)
            E[v][u].second+=1;
    }
}
void Graph::insertECoreMT(int u,int v,int w){
    if(E[u].find(v)==E[u].end()){//if not found
        E[u].insert(make_pair(v,make_pair(w,1)));
        DD[u]++;
//		DD2[u]++;
    }
    else{//if found
        if(E[u][v].first>w)
            E[u][v]= make_pair(w,1);
        else if(E[u][v].first==w)
            E[u][v].second+=1;
    }

}

//compute the father tree node
int Graph::matchCore(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(NodeOrder[vert[i].first]<NodeOrder[nearest])//get the least node order
            nearest=vert[i].first;
    }
    //cout<<nearest<<" "<<rankCore[nearest]<<endl;
    return rank[nearest];
}


//compute the father tree node
int Graph::matchCoreParti(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(NodeOrder[vert[i].first]<NodeOrder[nearest])//get the least node order
            nearest=vert[i].first;
    }
    //cout<<nearest<<" "<<rankCore[nearest]<<endl;
//    if(rank.find(nearest)==rank.end()){//if not found
//        cout<<"Not found "<<nearest<<" in rank!"<<endl; exit(1);
//    }
    return rank[IDMap[nearest]];
}

//construct RMQ index
void Graph::makeRMQCoreP(int pid, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs, vector<vector<Node>>& Trees){
    //EulerSeq.clear();
    toRMQs[pid].assign(PartiVertex[pid].size(),0);
    vector<int> EulerSeqP;
    //RMQIndex.clear();
    makeRMQDFSCoreP(pid, 0, 1, EulerSeqP, toRMQs, Trees);
    RMQIndexs[pid].push_back(EulerSeqP);

    int m = EulerSeqP.size();
//    cout<<"m: "<<m<<endl;
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndexs[pid][k - 1][j], y = RMQIndexs[pid][k - 1][j + i / 2];
//            cout<<"x and y: "<<x<<" "<<y<<endl;
            if (Trees[pid][x].height < Trees[pid][y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndexs[pid].push_back(tmp);
    }
}

void Graph::makeRMQDFSCoreP(int pid, int p, int height, vector<int>& EulerSeqP, vector<vector<int>>& toRMQs, vector<vector<Node>>& Trees){
    toRMQs[pid][p] = EulerSeqP.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeqP.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Trees[pid][p].ch.size(); i++){
        makeRMQDFSCoreP(pid,Trees[pid][p].ch[i], height + 1, EulerSeqP, toRMQs, Trees);
        EulerSeqP.push_back(p);
    }
}

void Graph::makeRMQ(vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree){
    vector<int> EulerSeq;
    EulerSeq.clear();
    toRMQ.assign(node_num,0);
    //RMQIndex.clear();
    makeRMQDFS(0, 1, EulerSeq, toRMQ, Tree);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}

void Graph::makeRMQDFS(int p, int height, vector<int>& EulerSeq, vector<int>& toRMQ, vector<Node>& Tree){
    toRMQ[p] = EulerSeq.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeq.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFS(Tree[p].ch[i], height + 1, EulerSeq, toRMQ, Tree);
        EulerSeq.push_back(p);
    }
}

//construct RMQ index
void Graph::makeRMQCore(){
    vector<int> EulerSeq;
    EulerSeq.clear();
    toRMQ.assign(node_num,0);
    //RMQIndex.clear();
    makeRMQDFSCore(0, 1, EulerSeq);
    RMQIndex.push_back(EulerSeq);

    int m = EulerSeq.size();
    for (int i = 2, k = 1; i < m; i = i * 2, k++){
        vector<int> tmp;
        //tmp.clear();
        tmp.assign(m,0);
        for (int j = 0; j < m - i; j++){
            int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
            if (Tree[x].height < Tree[y].height)
                tmp[j] = x;
            else tmp[j] = y;
        }
        RMQIndex.push_back(tmp);
    }
}

void Graph::makeRMQDFSCore(int p, int height, vector<int>& EulerSeq){
    toRMQ[p] = EulerSeq.size();//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    EulerSeq.push_back(p);//EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFSCore(Tree[p].ch[i], height + 1, EulerSeq);
        EulerSeq.push_back(p);
    }
}


/// Query Processing
//function for Query processing, new
int Graph::QueryNP(int ID1, int ID2){
    int dis=INF;

    if(algoQuery==Dijk){
//        dis = Dijkstra(ID1,ID2,Neighbor);
        dis = BiDijkstra(ID1,ID2,Neighbor);
//        dis = Astar(ID1,ID2,Neighbor);
    }else if(algoQuery==CH){//CH
        dis = QueryCHWP(ID1,ID2);
    }else if(algoQuery==H2H){//H2H
        dis = QueryH2H(ID1,ID2);
    }

    return dis;
}

int	Graph::QueryCHWP(int ID1, int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

//            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){//
            for(auto out=Tree[rank[topNodeIDForward]].vert.begin();out!=Tree[rank[topNodeIDForward]].vert.end();out++){//
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
            for(auto in=Tree[rank[topNodeIDBackward]].vert.begin();in!=Tree[rank[topNodeIDBackward]].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

int Graph::LCAQuery(int _p, int _q){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}

int Graph::QueryH2H(int ID1,int ID2){
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQuery(r1,r2);
//    cout<<r1<<" "<<r2<<" "<<LCA<<" "<<Tree.size()<<endl;
    if(LCA==r1)
        return Tree[r2].dis[Tree[r1].pos.back()];
    else if(LCA==r2)
        return Tree[r1].dis[Tree[r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Tree[LCA].pos.size();i++){
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
                tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
        }
        return tmp;
    }
}

//function for Query processing of PostMHL
int Graph::QueryPostMHL(int ID1, int ID2){
    int dis=INF;
    if(algoQuery==Dijk){
//        dis = Dijkstra(ID1,ID2,Neighbor);
        dis = BiDijkstra(ID1,ID2,Neighbor);
//        dis = Astar(ID1,ID2,Neighbor);
        return dis;
    }
    int algoQuery=this->algoQuery;

    if(algoQuery>PH2H_Cross){
        cout<<"Wrong query type! "<<algoQuery<<endl; exit(1);
    }
    if(algoQuery==PH2H_Cross){
        dis = QueryPostMHLPartiPartiExt(ID1, ID2);
    }else{
        if(PartiTags[ID1].first==PartiTags[ID2].first && PartiTags[ID1].first!=-1){//if same partition
            if(algoQuery==PH2H_Post){
//            cout<<"Same-parti"<<endl;
                dis= QueryPostMHLSamePartiPost(ID1,ID2);
//            dis = QueryPostMHLPartiPartiExt(ID1, ID2);
            }else if(algoQuery==PCH_No){
                if(PartiTags[ID1].first==-1 && PartiTags[ID2].first==-1){//Case 1: both in overlay graph
//                cout<<"Same-parti: Core-Core"<<endl;
                    dis=QueryOverlayCH(ID1, ID2);
                }
                else if(PartiTags[ID1].first==-1 && PartiTags[ID2].first!=-1){//Case 2: ID2 in partition, ID1 in core
//                cout<<"Same-parti: Core-Parti"<<endl;
                    dis=QueryOverlayCH(ID2, ID1);
                }else if(PartiTags[ID1].first!=-1 && PartiTags[ID2].first==-1){//Case 2: ID1 in partition, ID2 in core
//                cout<<"Same-parti: Parti-Core"<<endl;
                    dis = QueryOverlayCH(ID1, ID2);
                }else{//Case 3: Same partition
//                cout<<"Same-parti: same-parti"<<endl;
                    dis= QueryOverlayCH(ID1,ID2);
                }
            }
        }
        else{//if different partition
            if(PartiTags[ID1].first==-1 && PartiTags[ID2].first==-1){//Case 1: both in overlay graph
//            cout<<"Different Parti: Core-Core"<<endl;
                if(algoQuery==PCH_No){
                    dis=QueryOverlayCH(ID1, ID2);
                }else if(algoQuery==PH2H_Post || algoQuery==PH2H_Cross){
                    dis= QueryOverlay(ID1,ID2);
                }
            }
            else if(PartiTags[ID1].first==-1 && PartiTags[ID2].first!=-1){//Case 2: ID2 in partition, ID1 in core
//            cout<<"Different Parti: Core-Parti"<<endl;
                if(algoQuery==PCH_No){
                    dis=QueryOverlayCH(ID2, ID1);
                }else if(algoQuery==PH2H_Post){
                    dis=QueryPostMHLPartiOverlay(ID2, ID1);
                }else if(algoQuery==PH2H_Cross){
                    dis= QueryPostMHLPartiPartiExt(ID2,ID1);
                }

            }else if(PartiTags[ID1].first!=-1 && PartiTags[ID2].first==-1){//Case 2: ID1 in partition, ID2 in core
//            cout<<"Different Parti: Parti-Core"<<endl;
                if(algoQuery==PCH_No){
                    dis=QueryOverlayCH(ID1, ID2);
                }else if(algoQuery==PH2H_Post){
                    dis=QueryPostMHLPartiOverlay(ID1, ID2);
                }else if(algoQuery==PH2H_Cross){
                    dis= QueryPostMHLPartiPartiExt(ID1, ID2);
                }
            }
            else{//Case 4: Different Partitions
//            cout<<"Different Parti: Parti-Parti"<<endl;
                if(algoQuery==PCH_No) {//PCH-No, PCH-Post
                    dis = QueryOverlayCH(ID1, ID2);
                }else if(algoQuery==PH2H_Post){//PH2H-No, PH2H-Post
                    dis = QueryPostMHLPartiParti(ID1, ID2);
                }else if(algoQuery==PH2H_Cross){//PH2H-Extend
                    dis = QueryPostMHLPartiPartiExt(ID1, ID2);
                }
            }

        }
    }


    return dis;
}

int Graph::QueryPostMHLDebug(int ID1, int ID2){
    int dis=INF;

    int algoQuery=this->algoQuery;

    if(algoQuery>PH2H_Cross){
        cout<<"Wrong query type! "<<algoQuery<<endl; exit(1);
    }

    if(PartiTags[ID1].first==PartiTags[ID2].first && PartiTags[ID1].first!=-1){//if same partition
        if(algoQuery>=PH2H_Post){
            cout<<"Same-parti"<<endl;
            dis= QueryPostMHLSamePartiPostDebug(ID1,ID2);
        }else if(algoQuery==PCH_No){
            if(PartiTags[ID1].first==-1 && PartiTags[ID2].first==-1){//Case 1: both in overlay graph
                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryOverlayCH(ID1, ID2);
            }
            else if(PartiTags[ID1].first==-1 && PartiTags[ID2].first!=-1){//Case 2: ID2 in partition, ID1 in core
                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryOverlayCH(ID2, ID1);
            }else if(PartiTags[ID1].first!=-1 && PartiTags[ID2].first==-1){//Case 2: ID1 in partition, ID2 in core
                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryOverlayCH(ID1, ID2);
            }else{//Case 3: Same partition
                cout<<"Same-parti: same-parti"<<endl;
                dis= QueryOverlayCH(ID1,ID2);
            }
        }
    }
    else{//if different partition
        if(PartiTags[ID1].first==-1 && PartiTags[ID2].first==-1){//Case 1: both in overlay graph
            cout<<"Different Parti: Core-Core"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryOverlayCH(ID1, ID2);
            }else if(algoQuery==PH2H_Post || algoQuery==PH2H_Cross){
                dis= QueryOverlayDebug(ID1,ID2);
            }
        }
        else if(PartiTags[ID1].first==-1 && PartiTags[ID2].first!=-1){//Case 2: ID2 in partition, ID1 in core
            cout<<"Different Parti: Core-Parti"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryOverlayCH(ID2, ID1);
            }else if(algoQuery==PH2H_Post){
                dis=QueryPostMHLPartiOverlayDebug(ID2, ID1);
            }else if(algoQuery==PH2H_Cross){
                dis= QueryPostMHLPartiPartiExt(ID2,ID1);
            }

        }else if(PartiTags[ID1].first!=-1 && PartiTags[ID2].first==-1){//Case 2: ID1 in partition, ID2 in core
            cout<<"Different Parti: Parti-Core"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryOverlayCH(ID1, ID2);
            }else if(algoQuery==PH2H_Post){
                dis=QueryPostMHLPartiOverlayDebug(ID1, ID2);
            }else if(algoQuery==PH2H_Cross){
                dis= QueryPostMHLPartiPartiExt(ID1, ID2);
            }
        }
        else{//Case 4: Different Partitions
            cout<<"Different Parti: Parti-Parti"<<endl;
            if(algoQuery==PCH_No) {//PCH-No, PCH-Post
                dis = QueryOverlayCH(ID1, ID2);
            }else if(algoQuery==PH2H_Post){//PH2H-No, PH2H-Post
                dis = QueryPostMHLPartiPartiDebug(ID1, ID2);
            }else if(algoQuery==PH2H_Cross){//PH2H-Extend
                dis = QueryPostMHLPartiPartiExtDebug(ID1, ID2);
            }
        }

    }

    return dis;
}

//function for Query processing, new
int Graph::QueryPMHLOpt(int ID1, int ID2){
    int dis=INF;
    if(algoQuery==Dijk){
        dis = Dijkstra(ID1,ID2,Neighbor);
//        dis = Astar(ID1,ID2,Neighbor);
        return dis;
    }
    int algoQuery=this->algoQuery;

    if(algoQuery>PH2H_Cross){
        cout<<"Wrong query type! "<<algoQuery<<endl; exit(1);
    }

    if(PartiTag[ID1].first==PartiTag[ID2].first){//if same partition
        if(algoQuery>=PH2H_Post){
//            cout<<"Same-parti"<<endl;
            dis= QuerySamePartiPost(ID1,ID2);
        }else if(algoQuery==PCH_No){
            if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryCoreCH(ID1, ID2);
            }
            else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryPartiCoreCH(ID2, ID1);
            }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryPartiCoreCH(ID1, ID2);
            }else{//Case 3: Same partition
//                cout<<"Same-parti: same-parti"<<endl;
                dis= QueryPartiPartiCH(ID1,ID2);
            }
        }
    }
    else{//if different partition
        if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//            cout<<"Different Parti: Core-Core"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryCoreCH(ID1, ID2);
            }else if(algoQuery==PH2H_Post || algoQuery==PH2H_Cross){
                dis= QueryCore(ID1,ID2);
            }
        }
        else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//            cout<<"Different Parti: Core-Parti"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryPartiCoreCH(ID2, ID1);
            }else if(algoQuery==PH2H_Post){
                dis=QueryPartiCore(ID2, ID1);
            }else if(algoQuery==PH2H_Cross){
                dis= QueryPartiCoreExtLCA(ID2,ID1);
            }

        }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//            cout<<"Different Parti: Parti-Core"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryPartiCoreCH(ID1, ID2);
            }else if(algoQuery==PH2H_Post){
                dis=QueryPartiCore(ID1, ID2);
            }else if(algoQuery==PH2H_Cross){
                dis= QueryPartiCoreExtLCA(ID1, ID2);
            }
        }
        else{//Case 4: Different Partitions
//            cout<<"Different Parti: Parti-Parti"<<endl;
            if(algoQuery==PCH_No) {//PCH-No, PCH-Post
                dis = QueryPartiPartiCH(ID1, ID2);
            }else if(algoQuery==PH2H_Post){//PH2H-No, PH2H-Post
                dis = QueryPartiParti(ID1, ID2);
            }else if(algoQuery==PH2H_Cross){//PH2H-Extend
                dis = QueryPartiPartiExtLCA(ID1, ID2);
            }
        }

    }

    return dis;
}

//function for Query processing, new
int Graph::QueryPMHL(int ID1, int ID2){
    int dis=INF;
    if(algoQuery==Dijk){
//        dis = Dijkstra(ID1,ID2,Neighbor);
        dis = BiDijkstra(ID1,ID2,Neighbor);
//        dis = Astar(ID1,ID2,Neighbor);
        return dis;
    }
    int algoQuery=this->algoQuery;

    if(algoQuery>PH2H_Cross){
        cout<<"Wrong query type! "<<algoQuery<<endl; exit(1);
    }

    if(PartiTag[ID1].first==PartiTag[ID2].first){//if same partition
        if(algoQuery>=PH2H_Post){
//            cout<<"Same-parti"<<endl;
            dis= QuerySamePartiPost(ID1,ID2);
        }else if(algoQuery==PH2H_No){//no-boundary
            if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryCore(ID1, ID2);
            }
            else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryPartiCore(ID2, ID1);
            }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryPartiCore(ID1, ID2);
            }else{//Case 3: Same partition
//                cout<<"Same-parti: same-parti"<<endl;
                dis= QuerySameParti(ID1,ID2);
            }
        }else if(algoQuery==PCH_No){
            if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryCoreCH(ID1, ID2);
            }
            else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryPartiCoreCH(ID2, ID1);
            }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryPartiCoreCH(ID1, ID2);
            }else{//Case 3: Same partition
//                cout<<"Same-parti: same-parti"<<endl;
                if(algoQuery==PCH_No){
                    dis= QueryPartiPartiCH(ID1,ID2);
                }else{
                    dis= QueryCHPartition(ID1,ID2,PartiTag[ID1].first);
                }
            }
        }
    }
    else{//if different partition
        if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//            cout<<"Different Parti: Core-Core"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryCoreCH(ID1, ID2);
            }else if(algoQuery==PH2H_No || algoQuery==PH2H_Post || algoQuery==PH2H_Cross){
                dis= QueryCore(ID1,ID2);
            }
        }
        else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//            cout<<"Different Parti: Core-Parti"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryPartiCoreCH(ID2, ID1);
            }else if(algoQuery==PH2H_No || algoQuery==PH2H_Post){
                dis=QueryPartiCore(ID2, ID1);
            }else if(algoQuery==PH2H_Cross){
                dis= QueryPartiCoreExtLCA(ID2,ID1);
            }

        }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//            cout<<"Different Parti: Parti-Core"<<endl;
            if(algoQuery==PCH_No){
                dis=QueryPartiCoreCH(ID1, ID2);
            }else if(algoQuery==PH2H_No || algoQuery==PH2H_Post){
                dis=QueryPartiCore(ID1, ID2);
            }else if(algoQuery==PH2H_Cross){
                dis= QueryPartiCoreExtLCA(ID1, ID2);
            }
        }
        else{//Case 4: Different Partitions
//            cout<<"Different Parti: Parti-Parti"<<endl;
            if(algoQuery==PCH_No) {//PCH-No, PCH-Post
                dis = QueryPartiPartiCH(ID1, ID2);
            }else if(algoQuery==PH2H_No || algoQuery==PH2H_Post){//PH2H-No, PH2H-Post
                dis = QueryPartiParti(ID1, ID2);
            }else if(algoQuery==PH2H_Cross){//PH2H-Extend
                dis = QueryPartiPartiExtLCA(ID1, ID2);
            }
        }

    }

    return dis;
}

//function for Query processing, new
int Graph::Query(int ID1, int ID2){
    int dis=INF;
    if(algoQuery==Dijk){
        dis = Dijkstra(ID1,ID2,Neighbor);
//        dis = Astar(ID1,ID2,Neighbor);
        return dis;
    }

    if(PartiTag[ID1].first==PartiTag[ID2].first){//if same partition
        if(algoQuery>=PH2H_Post){
//            cout<<"Same-parti"<<endl;
            dis= QuerySamePartiPost(ID1,ID2);
        }else if(algoQuery==PH2H_No){//no-boundary
            if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//                cout<<"Same-parti: Core-Core"<<endl;
                dis=QueryCore(ID1, ID2);
            }
            else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//                cout<<"Same-parti: Core-Parti"<<endl;
                dis=QueryPartiCore(ID2, ID1);
            }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//                cout<<"Same-parti: Parti-Core"<<endl;
                dis = QueryPartiCore(ID1, ID2);
            }else{//Case 3: Same partition
//                cout<<"Same-parti: Non-boundary"<<endl;
                dis= QuerySameParti(ID1,ID2);
            }
        }
    }
    else{//if different partition
        if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
//            cout<<"Different Parti: Core-Core"<<endl;
            dis=QueryCore(ID1, ID2);
        }
        else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
//            cout<<"Different Parti: Core-Parti"<<endl;
            if(algoQuery<=PH2H_Post){
                dis=QueryPartiCore(ID2, ID1);
            }else if(algoQuery==PH2H_Cross){
                dis= QueryPartiCoreExt(ID2,ID1);
            }

        }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
//            cout<<"Different Parti: Parti-Core"<<endl;
            if(algoQuery<=PH2H_Post) {
                dis = QueryPartiCore(ID1, ID2);
            }else if(algoQuery==PH2H_Cross){
                dis= QueryPartiCoreExt(ID1,ID2);
            }
        }
        else{//Case 4: Different Partitions
//            cout<<"Different Parti: Parti-Parti"<<endl;
            if(algoQuery<=PH2H_Post) {
                dis = QueryPartiParti(ID1, ID2);
            }else if(algoQuery==PH2H_Cross){
                dis = QueryPartiPartiExt(ID1, ID2);
            }
        }

    }

    return dis;
}
//function for Query processing, old
/*int Graph::Query(int ID1, int ID2){
    int dis=INF;

    if(PartiTag[ID1].second && PartiTag[ID2].second){//Case 1: both in overlay graph
        cout<<"Core-Core"<<endl;
        dis=QueryCore(ID1, ID2);
    }else if(PartiTag[ID1].second && !PartiTag[ID2].second){//Case 2: ID2 in partition, ID1 in core
        cout<<"Core-Parti"<<endl;
        if(algoQuery<=1){
            dis=QueryPartiCore(ID2, ID1);
        }else if(algoQuery==ExtendedLabel){
            dis= QueryPartiCoreExt(ID2,ID1);
        }

    }else if(!PartiTag[ID1].second && PartiTag[ID2].second){//Case 2: ID1 in partition, ID2 in core
        cout<<"Parti-Core"<<endl;
        if(algoQuery<=1) {
            dis = QueryPartiCore(ID1, ID2);
        }else if(algoQuery==ExtendedLabel){
            dis= QueryPartiCoreExt(ID1,ID2);
        }
    }else if(!PartiTag[ID1].second && !PartiTag[ID2].second){//both in partition

        if(PartiTag[ID1].first != PartiTag[ID2].first){//Case 3: in different peripheries
            cout<<"Parti-Parti!"<<endl;
            if(algoQuery<=1) {
                dis = QueryPartiParti(ID1, ID2);
            }else if(algoQuery==ExtendedLabel){
                dis = QueryPartiPartiExt(ID1, ID2);
            }

        }else{//Case 4: in the same periphery
//            cout<<"Same partition!"<<endl;
            if(algoQuery==NO_Boundary){//no-boundary
                dis= QuerySameParti(ID1,ID2);
            }else if(algoQuery==Post_Boundary){
                dis= QuerySamePartiPost(ID1,ID2);
            }


        }
    }
    return dis;
}*/
//Case 1: query on overlay graph
int Graph::QueryCore(int ID1, int ID2){
    if(!PartiTag[ID1].second || !PartiTag[ID2].second){
        cout<<"Not overlay vertex! "<<ID1<<"("<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].second<<")"<<endl; exit(1);
    }
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);

    if(LCA==r1)
        return Tree[r2].dis[Tree[r1].pos.back()];
    else if(LCA==r2)
        return Tree[r1].dis[Tree[r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Tree[LCA].pos.size();i++){
//            if(Tree[LCA].pos[i]>=Tree[r1].dis.size() || Tree[LCA].pos[i]>=Tree[r2].dis.size()){
//                cout<<ID1<<"("<<r1<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<r2<<","<<PartiTag[ID2].second<<") "<<LCA<<" "<<Tree.size()<<": "<<Tree[LCA].pos[i]<<" "<<Tree[r1].dis.size()<<" "<<Tree[r2].dis.size()<<endl;
//            }
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
                tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
        }
        return tmp;
    }
}
int Graph::QueryCoreDebug(int ID1, int ID2){
    if(!PartiTag[ID1].second || !PartiTag[ID2].second){
        cout<<"Not overlay vertex! "<<ID1<<"("<<PartiTag[ID1].second<<") "<<ID2<<"("<<PartiTag[ID2].second<<")"<<endl; exit(1);
    }
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);
    int d1,d2,ancestor1,ancestor2;

    cout<<"Overlay Query."<<endl;

    if(LCA==r1){
        cout<<"r1. Hub of QueryCore: "<<Tree[LCA].uniqueVertex<<endl;
        return Tree[r2].dis[Tree[r1].pos.back()];
    }
    else if(LCA==r2){
        cout<<"r2. Hub of QueryCore: "<<Tree[LCA].uniqueVertex<<endl;
        return Tree[r1].dis[Tree[r2].pos.back()];
    }
    else{
        int tmp=INF;
        cout<<"LCA="<<LCA<<" "<<Tree[LCA].uniqueVertex<<endl;
        for(int i=0;i<Tree[LCA].pos.size();i++){
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]]){
                d1=Tree[r1].dis[Tree[LCA].pos[i]], d2=Tree[r2].dis[Tree[LCA].pos[i]];
                ancestor1=Tree[r1].vAncestor[Tree[LCA].pos[i]], ancestor2=Tree[r2].vAncestor[Tree[LCA].pos[i]];
                tmp=d1+d2;
                cout<<ancestor1<<"("<<NodeOrder[ancestor1]<<") "<<ancestor2<<"("<<NodeOrder[ancestor2]<<") "<<tmp<<endl;

            }
        }
        cout<<"Hub of QueryCore: "<<ancestor1<<"("<<NodeOrder[ancestor1]<<") "<<ancestor2<<"("<<NodeOrder[ancestor2]<<")"<< endl;
        return tmp;
    }
}

//Case 2: one core, one tree
int Graph::QueryPartiCore(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTag[ID1].first;
    int bid;
    int dis1,dis2;
    if(algoQuery==PH2H_No){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartition(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            if(d>dis1+dis2)
                d=dis1+dis2;
        }
    }
    else if(algoQuery==PH2H_Post){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartitionPost(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            if(d>dis1+dis2)
                d=dis1+dis2;
        }
    }

    return d;
}

int Graph::QueryPartiCoreDebug(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTag[ID1].first;
    int bid;
    int dis1,dis2;
    int finalbid,finaldis1,finaldis2;
    if(algoQuery==PH2H_No){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartition(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            if(d>dis1+dis2){
                d=dis1+dis2;
                cout<<bid<<": "<<dis1<<" "<<dis2<<" "<<d<<endl;
                finalbid=bid, finaldis1=dis1, finaldis2=dis2;
            }

        }
    }
    else if(algoQuery==PH2H_Post){
        for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
            bid=*it;
            dis1= QueryH2HPartitionPost(ID1,bid,pid);
            dis2= QueryCore(bid,ID2);
            if(d>dis1+dis2){
                d=dis1+dis2;
                cout<<bid<<": "<<dis1<<" "<<dis2<<" "<<d<<endl;
                finalbid=bid, finaldis1=dis1, finaldis2=dis2;
            }
        }
    }

    int dDijk_s=Dijkstra(ID1,finalbid,Neighbor), dDijk_t=Dijkstra(finalbid,ID2,Neighbor);
    cout<<ID1<<" "<<finalbid<<"("<<NodeOrder[finalbid]<<","<<PartiTag[finalbid].first<<","<<PartiTag[finalbid].second<<") "<<ID2<<" : "<<finaldis1<<" "<<finaldis2<<" "<<d<<" ; "<<dDijk_s<<" "<<dDijk_t<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;

    return d;
}
//Case 2: one core, one tree, extension label, version 1, PLL, correct
int Graph::QueryPartiCoreExt(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTag[ID1].first;
    int ancestor;
    int dis1,dis2;
    unordered_map<int,int> disMap; disMap.clear();
    for(int i=0;i<TreeExt[rankExt[ID1]].dis.size();++i){
        ancestor=TreeExt[rankExt[ID1]].vAncestor[i];
        disMap.insert({ancestor,TreeExt[rankExt[ID1]].dis[i]});
    }
    if(disMap.find(ID2)!=disMap.end()){//if found
        d=disMap[ID2];
    }
    else{//if not found
        for(int i=0;i<Tree[rank[ID2]].dis.size();++i){
            ancestor=Tree[rank[ID2]].vAncestor[i];
            if(disMap.find(ancestor)!=disMap.end()){//if found
                dis1=Tree[rank[ID2]].dis[i];
                dis2=disMap[ancestor];
                if(d>dis1+dis2){
                    d=dis1+dis2;
                }
            }
        }
    }


    return d;
}
//Case 2: one core, one tree, extension label, version 2, LCA
int Graph::QueryPartiCoreExtLCA(int ID1, int ID2){//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTag[ID1].first;

    int r1=rankExt[ID1], r2=rankExt[ID2];
    int LCA=LCAQuery(r1,r2, toRMQExt, RMQIndexExt, TreeExt);
//        cout<<"LCA: "<<TreeExt[LCA].uniqueVertex<<"("<<PartiTag[TreeExt[LCA].uniqueVertex].first<<","<<PartiTag[TreeExt[LCA].uniqueVertex].second<<") "<<TreeExt[LCA].pos.size()<<" "<<TreeExt[LCA].height<<endl;
//        cout<<"ID1: "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<TreeExt[r1].dis.size()<<" "<<TreeExt[r1].height<<" ; ID2: "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<") "<<TreeExt[r2].dis.size()<<" "<<TreeExt[r2].height<< endl;
    if(LCA==r1){
        d = TreeExt[r2].dis[TreeExt[r1].pos.back()];
    }
    else if(LCA==r2){
        d = TreeExt[r1].dis[TreeExt[r2].pos.back()];
    }
    else{
        for(int i=0;i<TreeExt[LCA].pos.size();i++){
            if(d>TreeExt[r1].dis[TreeExt[LCA].pos[i]]+TreeExt[r2].dis[TreeExt[LCA].pos[i]]){
//                    cout<<i<<"("<<TreeExt[LCA].pos[i]<<"): "<<TreeExt[r1].dis[TreeExt[LCA].pos[i]]<<"("<<TreeExt[r1].vAncestor[TreeExt[LCA].pos[i]]<<") "<<TreeExt[r2].dis[TreeExt[LCA].pos[i]]<<"("<<TreeExt[r2].vAncestor[TreeExt[LCA].pos[i]]<<") "<<d<<endl;
                d=TreeExt[r1].dis[TreeExt[LCA].pos[i]]+TreeExt[r2].dis[TreeExt[LCA].pos[i]];
            }

        }
    }

    return d;
}
//Case 3: Different trees
int Graph::QueryPartiParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;

        if(algoQuery==PH2H_No){
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];
                m1.insert(make_pair(bID1, QueryH2HPartition(ID1,bID1,pid1)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
                m2.insert(make_pair(bID2,QueryH2HPartition(ID2,bID2,pid2)));
            }
        }else if(algoQuery==PH2H_Post){
            for(int i=0;i<B1.size();i++){
                bID1=B1[i];
                m1.insert(make_pair(bID1, QueryH2HPartitionPost(ID1,bID1,pid1)));
            }
            for(int j=0;j<B2.size();j++){
                bID2=B2[j];
                m2.insert(make_pair(bID2,QueryH2HPartitionPost(ID2,bID2,pid2)));
            }
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
                    d1=m1[bID1]; d2=m2[bID2];
                    b1=bID1; b2=bID2;
                }

            }
        }

//        cout<<"b1, b2, d1, d2: "<<b1<<" "<<b2<<" "<<d1<<" "<<d2<<endl;
    }

    return d;
}

//Case 3: Different trees by extension label, version 1
int Graph::QueryPartiPartiExt(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;

//        map<int,int> m1,m2;
//        m1.clear();
//        m2.clear();
        int bID1, bID2, tempdis;
        int ancestor1,ancestor2,d1,d2;
        int lid, hid;
        int hub;

        if(TreeExt[rankExt[ID1]].dis.size()<TreeExt[rankExt[ID2]].dis.size()){
            hid=ID1, lid=ID2;
        }else{
            hid=ID2, lid=ID1;
        }

        for(int i=0;i<TreeExt[rankExt[hid]].dis.size();++i){
            ancestor1=TreeExt[rankExt[hid]].vAncestor[i];
            d1=TreeExt[rankExt[hid]].dis[i];
            for(int j=0;j<TreeExt[rankExt[lid]].dis.size();++j){
                ancestor2=TreeExt[rankExt[lid]].vAncestor[j];
                d2=TreeExt[rankExt[lid]].dis[j];
                if(ancestor1==ancestor2){
                    if(d>d1+d2){
                        d=d1+d2;
                        hub=ancestor1;
//                        cout<<d<<" "<<d1<<" "<<d2<<": "<<TreeExt[rankExt[hid]].cnt[i]<<" "<<TreeExt[rankExt[lid]].cnt[j]<<" "<< hub<<endl;
                    }
                    break;
                }
            }
        }

//        cout<<"b1, b2, d1, d2: "<<b1<<" "<<b2<<" "<<d1<<" "<<d2<<endl;
    }

    return d;
}

//Case 3: Different trees by extension label, version 2, LCA
int Graph::QueryPartiPartiExtLCA(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;

//        map<int,int> m1,m2;
//        m1.clear();
//        m2.clear();
        int bID1, bID2, tempdis;
        int ancestor1,ancestor2,d1,d2;
        int lid, hid;
        int hub;

        int r1=rankExt[ID1], r2=rankExt[ID2];
        int LCA=LCAQuery(r1,r2, toRMQExt, RMQIndexExt, TreeExt);
//        cout<<"LCA: "<<TreeExt[LCA].uniqueVertex<<"("<<PartiTag[TreeExt[LCA].uniqueVertex].first<<","<<PartiTag[TreeExt[LCA].uniqueVertex].second<<") "<<TreeExt[LCA].pos.size()<<" "<<TreeExt[LCA].height<<endl;
//        cout<<"ID1: "<<ID1<<"("<<PartiTag[ID1].first<<","<<PartiTag[ID1].second<<") "<<TreeExt[r1].dis.size()<<" "<<TreeExt[r1].height<<" ; ID2: "<<ID2<<"("<<PartiTag[ID2].first<<","<<PartiTag[ID2].second<<") "<<TreeExt[r2].dis.size()<<" "<<TreeExt[r2].height<< endl;
        if(LCA==r1){
            d = TreeExt[r2].dis[TreeExt[r1].pos.back()];
        }
        else if(LCA==r2){
            d = TreeExt[r1].dis[TreeExt[r2].pos.back()];
        }
        else{
            for(int i=0;i<TreeExt[LCA].pos.size();i++){
                if(d>TreeExt[r1].dis[TreeExt[LCA].pos[i]]+TreeExt[r2].dis[TreeExt[LCA].pos[i]]){
//                    cout<<i<<"("<<TreeExt[LCA].pos[i]<<"): "<<TreeExt[r1].dis[TreeExt[LCA].pos[i]]<<"("<<TreeExt[r1].vAncestor[TreeExt[LCA].pos[i]]<<") "<<TreeExt[r2].dis[TreeExt[LCA].pos[i]]<<"("<<TreeExt[r2].vAncestor[TreeExt[LCA].pos[i]]<<") "<<d<<endl;
                    d=TreeExt[r1].dis[TreeExt[LCA].pos[i]]+TreeExt[r2].dis[TreeExt[LCA].pos[i]];
                }

            }
        }

    }

    return d;
}

int Graph::QueryPartiPartiExtLCADebug(int ID1, int ID2){//both are within the partition
    int d=INF;
    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);
    }else{//if in different partitions
        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;
//        map<int,int> m1,m2;
//        m1.clear();
//        m2.clear();
        int bID1, bID2, tempdis;
        int ancestor1,ancestor2,d1,d2;
        int lid, hid;
        int hub;
        int r1=rankExt[ID1], r2=rankExt[ID2];
        int LCA=LCAQuery(r1,r2, toRMQExt, RMQIndexExt, TreeExt);

        if(LCA==r1){
            cout<<"LCA=r1. "<<TreeExt[r2].vAncestor[TreeExt[r1].pos.back()]<<endl;
            d = TreeExt[r2].dis[TreeExt[r1].pos.back()];
        }
        else if(LCA==r2){
            cout<<"LCA=r2. "<<TreeExt[r1].vAncestor[TreeExt[r2].pos.back()]<<endl;
            d = TreeExt[r1].dis[TreeExt[r2].pos.back()];
        }
        else{
            cout<<"LCA="<<LCA<<" "<<TreeExt[LCA].uniqueVertex<<endl;
            for(int i=0;i<TreeExt[LCA].pos.size();i++){
                if(d>TreeExt[r1].dis[TreeExt[LCA].pos[i]]+TreeExt[r2].dis[TreeExt[LCA].pos[i]]){
                    d=TreeExt[r1].dis[TreeExt[LCA].pos[i]]+TreeExt[r2].dis[TreeExt[LCA].pos[i]];
                    ancestor1=TreeExt[r1].vAncestor[TreeExt[LCA].pos[i]], ancestor2=TreeExt[r2].vAncestor[TreeExt[LCA].pos[i]];
                    cout<<ancestor1<<"("<<NodeOrder[ancestor1]<<") "<<ancestor2<<"("<<NodeOrder[ancestor2]<<") "<<d<<endl;
                }
            }
        }
        cout<<"Hub of QueryPartiParti: "<<ancestor1<<"("<<NodeOrder[ancestor1]<<") "<<ancestor2<<"("<<NodeOrder[ancestor2]<<")"<< endl;
    }

    return d;
}

//Case 4: Same tree, for original version
int Graph::QuerySameParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        int temp_dis = QueryH2HPartition(ID1,ID2,pid1);/// d2 may be wrong sometimes
        if(temp_dis<d)//QueryH2HPartition(ID1,ID2,pid1)
            d=temp_dis;//QueryH2HPartition(ID1,ID2,pid1);
        vector<int> B=BoundVertex[pid1];
        map<int,int> m1,m2;
        m1.clear();
        m2.clear();
        vector<int> B1,B2;
        B1.clear();
        B2.clear();
        int bID,d1,d2;
        for(int i=0;i<B.size();i++){
            bID=B[i];
            d1=QueryH2HPartition(ID1,bID,pid1);
            d2=QueryH2HPartition(ID2,bID,pid1);

            if(d1<d){
                B1.push_back(bID);
                m1.insert(make_pair(bID,d1));
            }
            if(d2<d){
                B2.push_back(bID);
                m2.insert(make_pair(bID,d2));
            }
        }

        int bID1, bID2, tempdis;
        if(!B1.empty() && !B2.empty()){
            for(int k=0;k<B1.size();k++){
                bID1=B1[k];
                if(m1[bID1]>d)
                    continue;
                for(int z=0;z<B2.size();z++){
                    bID2=B2[z];
                    if(m2[bID2]>d)
                        continue;
                    tempdis=m1[bID1]+QueryCore(bID1,bID2)+m2[bID2];
                    if(tempdis<d)
                        d=tempdis;
                }
            }
        }

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

//Case 4: Same tree, for query-orient version
int Graph::QuerySamePartiPost(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        d = QueryH2HPartitionPost(ID1,ID2,pid1);

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

int Graph::QuerySamePartiPostOpt(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTag[ID1].first;
    int pid2=PartiTag[ID2].first;
    if(pid1==pid2){//if in the same partition
//        cout<<"Same-Parti"<<endl;
        d = QueryH2HPartition(ID1,ID2,pid1);

    }else{//if in different partitions
        cout<<"Wrong for same partition query!"<<endl;
        exit(1);
    }

    return d;
}

int Graph::QueryCoreCH(int ID1, int ID2){
    if(ID1==ID2) return 0;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);
    int rForward, rBackward;

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

//            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
            rForward=rank[topNodeIDForward];
            for(auto out=Tree[rForward].vert.begin();out!=Tree[rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
            rBackward=rank[topNodeIDBackward];
            for(auto in=Tree[rBackward].vert.begin();in!=Tree[rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

//Query within one partition, no-boundary of PCH
int Graph::QueryCHPartition(int ID1, int ID2, int PID){
    if(ID1==ID2) return 0;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);
    int rForward, rBackward;

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

//            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
            rForward=ranks[PID][IDMap[topNodeIDForward]];
            for(auto out=Trees[PID][rForward].vert.begin();out!=Trees[PID][rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
            rBackward=ranks[PID][IDMap[topNodeIDBackward]];
            for(auto in=Trees[PID][rBackward].vert.begin();in!=Trees[PID][rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

//Query from partition to core, no-boundary of PCH
int Graph::QueryPartiCoreCH(int ID1, int ID2){//ID1: partition vertex
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);

    int rForward, rBackward;
    int PID;
    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

            //if(VtoParID[topNodeIDForward]==PID1){
//            for(auto out=NeighborCons[PID1][topNodeIDForward].begin();out!=NeighborCons[PID1][topNodeIDForward].end();out++){
            PID=PartiTag[topNodeIDForward].first;
            rForward=ranks[PID][IDMap[topNodeIDForward]];
            for(auto out=Trees[PID][rForward].vert.begin();out!=Trees[PID][rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
            //}else{
//            for(auto out=NeighborConOverlay[topNodeIDForward].begin();out!=NeighborConOverlay[topNodeIDForward].end();out++){
            rForward=rank[topNodeIDForward];
            if(rForward!=-1){
                for(auto out=Tree[rForward].vert.begin();out!=Tree[rForward].vert.end();out++){
                    neighborNodeID = (*out).first;
                    neighborLength = (*out).second.first;

                    int df = vDistanceForward[topNodeIDForward] + neighborLength;
                    if(!vVisitedF[neighborNodeID]){
                        if(vDistanceForward[neighborNodeID] > df){
                            //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                            vDistanceForward[neighborNodeID] = df;
                            fHeapForward.update(neighborNodeID, df);
                        }
                    }
                }
            }

            //}

        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborConOverlay[topNodeIDBackward].begin();in!=NeighborConOverlay[topNodeIDBackward].end();in++){
            rBackward=rank[topNodeIDBackward];
            for(auto in=Tree[rBackward].vert.begin();in!=Tree[rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

//Query from partition to core, no-boundary of PCH
int Graph::QueryPartiPartiCH(int ID1, int ID2) {//ID1 and ID2 partition vertices
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);

    int rForward, rBackward;
    int PID;
    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

            //if(VtoParID[topNodeIDForward]==PID1){
//            for(auto out=NeighborCons[PID1][topNodeIDForward].begin();out!=NeighborCons[PID1][topNodeIDForward].end();out++){
            PID=PartiTag[topNodeIDForward].first;
            rForward=ranks[PID][IDMap[topNodeIDForward]];
            for(auto out=Trees[PID][rForward].vert.begin();out!=Trees[PID][rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
            //}else{
//            for(auto out=NeighborConOverlay[topNodeIDForward].begin();out!=NeighborConOverlay[topNodeIDForward].end();out++){
            rForward=rank[topNodeIDForward];
            if(rForward!=-1){
                for(auto out=Tree[rForward].vert.begin();out!=Tree[rForward].vert.end();out++){
                    neighborNodeID = (*out).first;
                    neighborLength = (*out).second.first;

                    int df = vDistanceForward[topNodeIDForward] + neighborLength;
                    if(!vVisitedF[neighborNodeID]){
                        if(vDistanceForward[neighborNodeID] > df){
                            //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                            vDistanceForward[neighborNodeID] = df;
                            fHeapForward.update(neighborNodeID, df);
                        }
                    }
                }
            }

            //}

        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborConOverlay[topNodeIDBackward].begin();in!=NeighborConOverlay[topNodeIDBackward].end();in++){
            PID=PartiTag[topNodeIDBackward].first;
            rBackward=ranks[PID][IDMap[topNodeIDBackward]];
            for(auto in=Trees[PID][rBackward].vert.begin();in!=Trees[PID][rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }

            rBackward=rank[topNodeIDBackward];
            if(rBackward!=-1) {
                for (auto in = Tree[rBackward].vert.begin(); in != Tree[rBackward].vert.end(); in++) {
                    neighborNodeID = (*in).first;
                    neighborLength = (*in).second.first;

                    int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                    if (!vVisitedB[neighborNodeID]) {
                        if (vDistanceBackward[neighborNodeID] > db) {
                            vDistanceBackward[neighborNodeID] = db;
                            fHeapBackward.update(neighborNodeID, db);
                        }
                    }
                }
            }
        }
    }
    return d;
}

int Graph::QueryOverlayCH(int ID1, int ID2) {
    if(ID1==ID2) return 0;
    int d=INF;
    benchmark::heap<2,int,int> fHeapForward(node_num);
    benchmark::heap<2, int, int> fHeapBackward(node_num);

    //closed or not
    vector<bool> vVisitedF(node_num, false);
    vector<bool> vVisitedB(node_num, false);
    //the existing shortest distance
    vector<int>	vDistanceForward(node_num, INF);
    vector<int>	vDistanceBackward(node_num, INF);
    //stop search or not
    bool bF = false;
    bool bB = false;
    vDistanceForward[ID1] = 0;
    vDistanceBackward[ID2] = 0;
    fHeapForward.update(ID1,0);
    fHeapBackward.update(ID2,0);
    int rForward, rBackward;

    int topNodeIDForward, topNodeIDBackward,topDisForward,topDisBackward, neighborNodeID, neighborLength;

    while(!fHeapForward.empty() || !fHeapBackward.empty() )
    {
        if(bF && bB)
            break;
        if(bF && fHeapBackward.empty())
            break;
        if(bB && fHeapForward.empty())
            break;
        //Forward Search
        if(!fHeapForward.empty() && !bF)
        {
            fHeapForward.extract_min(topNodeIDForward, topDisForward);
            //cout<<topNodeIDForward<<" "<<topDisForward<<" "<<(NodeOrder[topNodeIDForward]>NodeOrder[137099])<<endl;

            if(vDistanceForward[topNodeIDForward] > d)
                bF = true;

            vVisitedF[topNodeIDForward] = true;

            if(vVisitedB[topNodeIDForward]){
                int distTmp=topDisForward+vDistanceBackward[topNodeIDForward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"forwardtopvalue "<<topDisForward<<" "<<vDistanceBackward[topNodeIDForward]<<" "<<d<<"meet "<<topNodeIDForward<<endl;
                }
            }

//            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
            rForward=rank[topNodeIDForward];
            for(auto out=Tree[rForward].vert.begin();out!=Tree[rForward].vert.end();out++){
                neighborNodeID = (*out).first;
                neighborLength = (*out).second.first;

                int df = vDistanceForward[topNodeIDForward] + neighborLength;
                if(!vVisitedF[neighborNodeID]){
                    if(vDistanceForward[neighborNodeID] > df){
                        //if(neighborNodeID==37783) cout<<"neighdis "<<vDistanceForward[neighborNodeID]<<" adddis "<<df<<endl;
                        vDistanceForward[neighborNodeID] = df;
                        fHeapForward.update(neighborNodeID, df);
                    }
                }
            }
        }

        //Backward Search
        if(!fHeapBackward.empty() && !bB)
        {
            fHeapBackward.extract_min(topNodeIDBackward, topDisBackward);

            if(vDistanceBackward[topNodeIDBackward] > d)
                bB = true;

            vVisitedB[topNodeIDBackward] = true;

            if(vVisitedF[topNodeIDBackward]){
                int distTmp=topDisBackward+vDistanceForward[topNodeIDBackward];
                if(distTmp<d){
                    d=distTmp;
                    //cout<<"backtopvalue "<<topDisBackward<<" "<<vDistanceForward[topNodeIDBackward]<<" "<<d<<"meet "<<topNodeIDBackward<<endl;
                }
            }

//            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
            rBackward=rank[topNodeIDBackward];
            for(auto in=Tree[rBackward].vert.begin();in!=Tree[rBackward].vert.end();in++){
                neighborNodeID = (*in).first;
                neighborLength = (*in).second.first;

                int db = vDistanceBackward[topNodeIDBackward] + neighborLength;
                if(!vVisitedB[neighborNodeID]){
                    if(vDistanceBackward[neighborNodeID]>db){
                        vDistanceBackward[neighborNodeID] = db;
                        fHeapBackward.update(neighborNodeID, db);
                    }
                }
            }
        }
    }
    return d;
}

//Case 1: query on overlay graph
int Graph::QueryOverlay(int ID1, int ID2){
    if(PartiTags[ID1].first != -1 || PartiTags[ID2].first!=-1){
        cout<<"Not overlay vertex! "<<ID1<<"("<<PartiTags[ID1].first<<") "<<ID2<<"("<<PartiTags[ID2].first<<")"<<endl; exit(1);
    }
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);

    if(LCA==r1)
        return Tree[r2].dis[Tree[r1].pos.back()];
    else if(LCA==r2)
        return Tree[r1].dis[Tree[r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Tree[LCA].pos.size();i++){
//            if(Tree[LCA].pos[i]>=Tree[r1].dis.size() || Tree[LCA].pos[i]>=Tree[r2].dis.size()){
//                cout<<ID1<<"("<<r1<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<r2<<","<<PartiTag[ID2].second<<") "<<LCA<<" "<<Tree.size()<<": "<<Tree[LCA].pos[i]<<" "<<Tree[r1].dis.size()<<" "<<Tree[r2].dis.size()<<endl;
//            }
            if(tmp>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]])
                tmp=Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]];
        }
        return tmp;
    }
}

int Graph::QueryOverlayDebug(int ID1, int ID2){
    if(PartiTags[ID1].first != -1 || PartiTags[ID2].first!=-1){
        cout<<"Not overlay vertex! "<<ID1<<"("<<PartiTags[ID1].first<<") "<<ID2<<"("<<PartiTags[ID2].first<<")"<<endl; exit(1);
    }
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int d=INF;
    int dis1,dis2,hub;
    int LCA=LCAQueryOverlay(r1,r2);
    cout<<"LCA: "<<Tree[LCA].uniqueVertex<<" , PartiTag: "<<PartiTags[Tree[LCA].uniqueVertex].first<<endl;

    if(LCA==r1){
        d=Tree[r2].dis[Tree[r1].pos.back()];
        dis1=0, dis2=d;
        hub=ID1;
    }
    else if(LCA==r2){
        d=Tree[r1].dis[Tree[r2].pos.back()];
        dis1=d, dis2=0;
        hub=ID2;
    }
    else{
        for(int i=0;i<Tree[LCA].pos.size();i++){
//            if(Tree[LCA].pos[i]>=Tree[r1].dis.size() || Tree[LCA].pos[i]>=Tree[r2].dis.size()){
//                cout<<ID1<<"("<<r1<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<r2<<","<<PartiTag[ID2].second<<") "<<LCA<<" "<<Tree.size()<<": "<<Tree[LCA].pos[i]<<" "<<Tree[r1].dis.size()<<" "<<Tree[r2].dis.size()<<endl;
//            }
            if(d>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]]){
                dis1=Tree[r1].dis[Tree[LCA].pos[i]], dis2=Tree[r2].dis[Tree[LCA].pos[i]];
                if(i<Tree[LCA].vert.size()){
                    hub=Tree[LCA].vert[i].first;
                }else{
                    hub=Tree[LCA].uniqueVertex;
                }
                cout<<i<<": "<<hub<<" "<<dis1<<" "<<dis2<<" "<<d<<endl;
                d=dis1+dis2;
            }
        }
    }

    int dDijk_s=Dijkstra(ID1,hub,Neighbor), dDijk_t=Dijkstra(hub,ID2,Neighbor), dDijk_st=Dijkstra(ID1,ID2,Neighbor);
    cout<<ID1<<" "<<hub<<"("<<NodeOrder[hub]<<","<<PartiTags[hub].first<<") "<<ID2<<": "<<dis1<<" "<<dis2<<" "<<d<<" ; "<<dDijk_s<<" "<<dDijk_t<<" "<<dDijk_st<<endl;
    return d;
}

//Query within one partition
int Graph::QueryPostMHLSamePartiPost(int ID1, int ID2) {
    if(PartiTags[ID1].first!=PartiTags[ID2].first || PartiTags[ID1].first==-1){
        cout<<"Not the same partition query! "<<ID1<<"("<<PartiTags[ID1].first<<") "<<ID2<<"("<<PartiTags[ID2].first<<")"<<endl; exit(1);
    }
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);
    int d1,d2,hub;
    int d=INF;
    int inf_i;
    int pid=PartiTags[ID1].first;
//    cout<<"LCA: "<<Tree[LCA].uniqueVertex<<" , PartiTag: "<<PartiTags[Tree[LCA].uniqueVertex].first<<endl;
//    cout<<"Boundary size: "<<BoundVertexMap[pid].size()<<endl;
    if(LCA==r1){
        d=Tree[r2].disPost[Tree[r1].pos.back()];
    }
    else if(LCA==r2){
        d=Tree[r1].disPost[Tree[r2].pos.back()];
    }
    else{

        for(int i=0;i<Tree[LCA].pos.size();i++){
            if(i<Tree[LCA].vert.size() && PartiTags[Tree[LCA].vert[i].first].first==-1){
                hub=Tree[LCA].vert[i].first;
                assert(BoundVertexMap[pid].find(hub)!=BoundVertexMap[pid].end());
                inf_i=BoundVertexMap[pid][hub];
//                inf_i=-1;
//                for(int j=0;j<BoundVertex[pid].size();++j){
//                    if(BoundVertex[pid][j] == hub){
//                        inf_i=j;
//                        break;
//                    }
//                }
//                if(inf_i==-1){
//                    cout<<"Wrong! "<<inf_i<<endl; exit(1);
//                }

                if(d>Tree[r1].disInf[inf_i]+Tree[r2].disInf[inf_i]){
                    d1=Tree[r1].disInf[inf_i], d2=Tree[r2].disInf[inf_i];
//                    cout<<i<<": "<<d1<<"("<<Tree[r1].vAncestorPost[Tree[LCA].posPost[i]]<<") "<<d2<<"("<<Tree[r2].vAncestorPost[Tree[LCA].posPost[i]]<<") "<<d<<endl;
                    d=d1+d2;
                }
            }
            else{
//            if(Tree[LCA].pos[i]>=Tree[r1].dis.size() || Tree[LCA].pos[i]>=Tree[r2].dis.size()){
//                cout<<ID1<<"("<<r1<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<r2<<","<<PartiTag[ID2].second<<") "<<LCA<<" "<<Tree.size()<<": "<<Tree[LCA].pos[i]<<" "<<Tree[r1].dis.size()<<" "<<Tree[r2].dis.size()<<endl; exit(1);
//            }
                if(d>Tree[r1].disPost[Tree[LCA].pos[i]]+Tree[r2].disPost[Tree[LCA].pos[i]]){
                    d1=Tree[r1].disPost[Tree[LCA].pos[i]], d2=Tree[r2].disPost[Tree[LCA].pos[i]];
                    hub=Tree[r1].vAncestor[Tree[LCA].pos[i]];
//                    cout<<i<<": "<<d1<<"("<<Tree[r1].vAncestorPost[Tree[LCA].posPost[i]]<<") "<<d2<<"("<<Tree[r2].vAncestorPost[Tree[LCA].posPost[i]]<<") "<<d<<endl;
                    d=d1+d2;
                }
            }

        }
        return d;
    }
}

int Graph::QueryPostMHLSamePartiPostDebug(int ID1, int ID2) {
    if(PartiTags[ID1].first!=PartiTags[ID2].first || PartiTags[ID1].first==-1){
        cout<<"Not the same partition query! "<<ID1<<"("<<PartiTags[ID1].first<<") "<<ID2<<"("<<PartiTags[ID2].first<<")"<<endl; exit(1);
    }
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);
    int d=INF;
    bool flag;
    int pid=PartiTags[ID1].first;
    cout<<"LCA: "<<Tree[LCA].uniqueVertex<<" , PartiTag: "<<PartiTags[Tree[LCA].uniqueVertex].first<<endl;

    int d1,d2,hub,finalhub;
    int inf_i;
    if(LCA==r1){
        d=Tree[r2].disPost[Tree[r1].pos.back()];
    }
    else if(LCA==r2){
        d=Tree[r1].disPost[Tree[r2].pos.back()];
    }
    else{
        for(int i=0;i<Tree[LCA].pos.size()-1;i++){
            if(PartiTags[Tree[LCA].vert[i].first].first==-1){
                hub=Tree[LCA].vert[i].first;
                flag=false;

                if(BoundVertexMap[pid].find(hub)!=BoundVertexMap[pid].end()){//if found
                    inf_i=BoundVertexMap[pid][hub];
                    flag=true;
                }

                if(!flag){
                    cout<<"Wrong for this boundary vertex! "<<hub<<endl; exit(1);
                }

                if(d>Tree[r1].disInf[inf_i]+Tree[r2].disInf[inf_i]){
                    finalhub=hub;
                    d1=Tree[r1].disInf[inf_i], d2=Tree[r2].disInf[inf_i];
                    cout<<"Interface. "<<i<<": "<<d1<<"("<<hub<<") "<<d2<<"("<<hub<<") "<<d<<" "<<d1+d2<<endl;
                    d=d1+d2;
                }
            }else{
                if(d>Tree[r1].disPost[Tree[LCA].pos[i]]+Tree[r2].disPost[Tree[LCA].pos[i]]){
                    d1=Tree[r1].disPost[Tree[LCA].pos[i]], d2=Tree[r2].disPost[Tree[LCA].pos[i]];
                    finalhub=Tree[r1].vAncestor[Tree[LCA].pos[i]];
                    cout<<"Ancestor. "<<i<<": "<<d1<<"("<<Tree[r1].vAncestor[Tree[LCA].pos[i]]<<") "<<d2<<"("<<Tree[r2].vAncestor[Tree[LCA].pos[i]]<<") "<<d<<" "<<d1+d2<<endl;
                    d=d1+d2;
                }
            }
        }

        int i=Tree[LCA].pos.size()-1;
        if(d>Tree[r1].disPost[Tree[LCA].pos[i]]+Tree[r2].disPost[Tree[LCA].pos[i]]){
            d1=Tree[r1].disPost[Tree[LCA].pos[i]], d2=Tree[r2].disPost[Tree[LCA].pos[i]];
            finalhub=Tree[r1].vAncestor[Tree[LCA].pos[i]];
            cout<<"Ancestor. "<<i<<": "<<d1<<"("<<Tree[r1].vAncestor[Tree[LCA].pos[i]]<<") "<<d2<<"("<<Tree[r2].vAncestor[Tree[LCA].pos[i]]<<") "<<d<<" "<<d1+d2<<endl;
            d=d1+d2;
        }
        int dDijk_s=Dijkstra(ID1,finalhub,Neighbor), dDijk_t=Dijkstra(finalhub,ID2,Neighbor), dDijk_st=Dijkstra(ID1,ID2,Neighbor);
        cout<<ID1<<" "<<finalhub<<"("<<NodeOrder[finalhub]<<","<<PartiTags[finalhub].first<<") "<<ID2<<": "<<d1<<" "<<d2<<" "<<d<<" ; "<<dDijk_s<<" "<<dDijk_t<<" "<<dDijk_st<<endl;
    }

    return d;
}

//Case 3: Different trees
int Graph::QueryPostMHLPartiParti(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTags[ID1].first;
    int pid2=PartiTags[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

//        map<int,int> m1,m2;
//        m1.clear();
//        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2;
        int inf1,inf2;
        for(int k=0;k<B1.size();k++){
            bID1=B1[k];
            assert(BoundVertexMap[pid1].find(bID1)!=BoundVertexMap[pid1].end());
//            if(BoundVertexMap[pid1].find(bID1)!=BoundVertexMap[pid1].end()){//if found
//                inf_i=BoundVertexMap[pid1][bID1];
//            }else{
//                cout<<"Wrong for this boundary vertex! "<<ID1<<" "<<bID1<<endl; exit(1);
//            }
            inf1=BoundVertexMap[pid1][bID1];

            if(Tree[rank[ID1]].disInf[inf1]>d)
                continue;

            for(int z=0;z<B2.size();z++){
                bID2=B2[z];
                assert(BoundVertexMap[pid2].find(bID2)!=BoundVertexMap[pid2].end());
//                if(BoundVertexMap[pid2].find(bID2)!=BoundVertexMap[pid2].end()){//if found
//                    inf2=BoundVertexMap[pid2][bID2];
//                }else{
//                    cout<<"Wrong for this boundary vertex! "<<ID2<<" "<<bID2<<endl; exit(1);
//                }
                inf2=BoundVertexMap[pid2][bID2];
                if(Tree[rank[ID2]].disInf[inf2]>d)
                    continue;

                tempdis=Tree[rank[ID1]].disInf[inf1]+QueryOverlay(bID1,bID2)+Tree[rank[ID2]].disInf[inf2];
                if(tempdis<d){
                    d=tempdis;
                    d1=Tree[rank[ID1]].disInf[inf1]; d2=Tree[rank[ID2]].disInf[inf2];
                    b1=bID1; b2=bID2;
                }

            }
        }

//        cout<<"b1, b2, d1, d2: "<<b1<<" "<<b2<<" "<<d1<<" "<<d2<<endl;
    }

    return d;
}

int Graph::QueryPostMHLPartiPartiDebug(int ID1, int ID2){//both are within the partition
    int d=INF;

    int pid1=PartiTags[ID1].first;
    int pid2=PartiTags[ID2].first;
    if(pid1==pid2){//if in the same partition
        cout<<"Wrong for partition-partition query!"<<endl;
        exit(1);

    }else{//if in different partitions
//        cout<<"Parti-Parti: "<<pid1<<" "<<pid2<<endl;
        vector<int> B1=BoundVertex[pid1];
        vector<int> B2=BoundVertex[pid2];

//        map<int,int> m1,m2;
//        m1.clear();
//        m2.clear();
        int bID1, bID2, tempdis;
        int b1,b2,d1,d2,doverlay,doverlayfinal;


        int inf1,inf2;
        for(int k=0;k<B1.size();k++){
            bID1=B1[k];
//            assert(BoundVertexMap[pid1].find(bID1)!=BoundVertexMap[pid1].end());
            if(BoundVertexMap[pid1].find(bID1)!=BoundVertexMap[pid1].end()){//if found
                inf1=BoundVertexMap[pid1][bID1];
            }else{
                cout<<"Wrong for this boundary vertex! "<<ID1<<" "<<bID1<<endl; exit(1);
            }
            inf1=BoundVertexMap[pid1][bID1];

            if(Tree[rank[ID1]].disInf[inf1]>d)
                continue;

            for(int z=0;z<B2.size();z++){
                bID2=B2[z];
//                assert(BoundVertexMap[pid2].find(bID2)!=BoundVertexMap[pid2].end());
                if(BoundVertexMap[pid2].find(bID2)!=BoundVertexMap[pid2].end()){//if found
                    inf2=BoundVertexMap[pid2][bID2];
                }else{
                    cout<<"Wrong for this boundary vertex! "<<ID2<<" "<<bID2<<endl; exit(1);
                }
                inf2=BoundVertexMap[pid2][bID2];
                if(Tree[rank[ID2]].disInf[inf2]>d)
                    continue;

                doverlay=QueryOverlay(bID1,bID2);
                tempdis=Tree[rank[ID1]].disInf[inf1]+doverlay+Tree[rank[ID2]].disInf[inf2];
                if(tempdis<d){
                    d=tempdis;
                    d1=Tree[rank[ID1]].disInf[inf1]; d2=Tree[rank[ID2]].disInf[inf2];
                    b1=bID1; b2=bID2;
                    doverlayfinal=doverlay;
                }

            }
        }

//        cout<<"b1, b2, d1, d2: "<<b1<<" "<<b2<<" "<<d1<<" "<<d2<<endl;

        int dDijk_s=Dijkstra(ID1,b1,Neighbor), dDijk_bb=Dijkstra(b1,b2,Neighbor), dDijk_t=Dijkstra(b2,ID2,Neighbor);
        cout<<ID1<<" "<<b1<<"("<<NodeOrder[b1]<<","<<PartiTags[b1].first<<") "<<b2<<"("<<NodeOrder[b2]<<","<<PartiTags[b2].first<<") "<<ID2<<" : "<<d1<<" "<<doverlayfinal<<" "<<d2<<" "<<d<<" ; "<<dDijk_s<<" "<<dDijk_bb<<" "<<dDijk_t<<" "<<Dijkstra(ID1,ID2,Neighbor)<<endl;
    }

    return d;
}

//Case 3: Different trees
int Graph::QueryPostMHLPartiPartiExt(int ID1, int ID2) {//both are within the partition
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);
    int d1,d2,hub;
    int d=INF;

    if(LCA==r1){
        d= Tree[r2].dis[Tree[r1].pos.back()];
    }
    else if(LCA==r2){
        d= Tree[r1].dis[Tree[r2].pos.back()];
    }
    else{
        for(int i=0;i<Tree[LCA].pos.size();i++){
//            if(Tree[LCA].pos[i]>=Tree[r1].dis.size() || Tree[LCA].pos[i]>=Tree[r2].dis.size()){
//                cout<<ID1<<"("<<r1<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<r2<<","<<PartiTag[ID2].second<<") "<<LCA<<" "<<Tree.size()<<": "<<Tree[LCA].pos[i]<<" "<<Tree[r1].dis.size()<<" "<<Tree[r2].dis.size()<<endl;
//            }
            d1=Tree[r1].dis[Tree[LCA].pos[i]], d2=Tree[r2].dis[Tree[LCA].pos[i]];
            if(d>d1+d2){
                hub=Tree[LCA].vert[i].first;
                d=d1+d2;
            }
        }

    }
    return d;
}

int Graph::QueryPostMHLPartiPartiExtDebug(int ID1, int ID2) {//both are within the partition
    if(ID1==ID2) return 0;
    if(NodeOrder[ID1]==-1 || NodeOrder[ID2]==-1) return INF;
    int r1=rank[ID1], r2=rank[ID2];
    int LCA=LCAQueryOverlay(r1,r2);
    int d1,d2,hub;
    int d=INF;

    cout<<"LCA: "<<Tree[LCA].uniqueVertex<<" , PartiTag: "<<PartiTags[Tree[LCA].uniqueVertex].first<<" , pos size: "<<Tree[LCA].pos.size()<<endl;

    if(LCA==r1){
        d= Tree[r2].dis[Tree[r1].pos.back()];
    }
    else if(LCA==r2){
        d= Tree[r1].dis[Tree[r2].pos.back()];
    }
    else{
        for(int i=0;i<Tree[LCA].pos.size();i++){
//            if(Tree[LCA].pos[i]>=Tree[r1].dis.size() || Tree[LCA].pos[i]>=Tree[r2].dis.size()){
//                cout<<ID1<<"("<<r1<<","<<PartiTag[ID1].second<<") "<<ID2<<"("<<r2<<","<<PartiTag[ID2].second<<") "<<LCA<<" "<<Tree.size()<<": "<<Tree[LCA].pos[i]<<" "<<Tree[r1].dis.size()<<" "<<Tree[r2].dis.size()<<endl;
//            }
            if(d>Tree[r1].dis[Tree[LCA].pos[i]]+Tree[r2].dis[Tree[LCA].pos[i]]){
                d1=Tree[r1].dis[Tree[LCA].pos[i]], d2=Tree[r2].dis[Tree[LCA].pos[i]];
                cout<<i<<": "<<d1<<"("<<Tree[r1].vAncestor[Tree[LCA].pos[i]]<<") "<<d2<<"("<<Tree[r2].vAncestor[Tree[LCA].pos[i]]<<") "<<d<<endl;
                hub=Tree[LCA].vert[i].first;
                d=d1+d2;
            }
        }
    }

    int dDijk_s=Dijkstra(ID1,hub,Neighbor), dDijk_t=Dijkstra(hub,ID2,Neighbor), dDijk_st=Dijkstra(ID1,ID2,Neighbor);
    cout<<ID1<<" "<<hub<<"("<<NodeOrder[hub]<<","<<PartiTags[hub].first<<") "<<ID2<<": "<<d1<<" "<<d2<<" "<<d<<" ; "<<dDijk_s<<" "<<dDijk_t<<" "<<dDijk_st<<endl;

    return d;
}

//Case 2: one core, one tree
int Graph::QueryPostMHLPartiOverlay(int ID1, int ID2) {//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTags[ID1].first;
    int bid;
    int dis1,dis2;
    if(algoQuery==PH2H_Post){
        int inf_i;

        if(BoundVertexMap[pid].find(ID2)!=BoundVertexMap[pid].end()){//if found
            inf_i=BoundVertexMap[pid][ID2];
            d=Tree[rank[ID1]].disInf[inf_i];
        }
        else{
            for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
                bid=*it;
//                if(BoundVertexMap[pid].find(bid)==BoundVertexMap[pid].end()){
//                    cout<<"Not found this boundary vertex! "<<ID1<<" "<<bid<<endl; exit(1);
//                }
                assert(BoundVertexMap[pid].find(bid)!=BoundVertexMap[pid].end());
                inf_i=BoundVertexMap[pid][bid];
                dis1= Tree[rank[ID1]].disInf[inf_i];
                dis2= QueryOverlay(bid,ID2);
                if(d>dis1+dis2){
                    d=dis1+dis2;
                }

            }
        }

    }

    return d;
}

int Graph::QueryPostMHLPartiOverlayDebug(int ID1, int ID2) {//ID1 partition, ID2 core
    int d=INF;

    int pid=PartiTags[ID1].first;
    int bid;
    int dis1,dis2;
    int disf1,disf2,hubf;

    if(algoQuery==PH2H_Post){
        int inf_i;

        if(BoundVertexMap[pid].find(ID2)!=BoundVertexMap[pid].end()){//if found
            cout<<"ID2 is boundary vertex! "<<endl;
            inf_i=BoundVertexMap[pid][ID2];
            d=Tree[rank[ID1]].disInf[inf_i];
            int dDijk_st=Dijkstra(ID1,ID2,Neighbor);
            cout<<ID1<<"("<<NodeOrder[ID1]<<","<<PartiTags[ID1].first<<") "<<ID2<<"("<<NodeOrder[ID2]<<","<<PartiTags[ID2].first<<"): "<<d<<" ; "<<dDijk_st<<endl;
        }
        else{
            cout<<"Not boundary vertex! "<<endl;
            for(auto it=BoundVertex[pid].begin();it!=BoundVertex[pid].end();++it){
                bid=*it;
//                if(BoundVertexMap[pid].find(bid)==BoundVertexMap[pid].end()){
//                    cout<<"Not found this boundary vertex! "<<ID1<<" "<<bid<<endl; exit(1);
//                }
                assert(BoundVertexMap[pid].find(bid)!=BoundVertexMap[pid].end());
                inf_i=BoundVertexMap[pid][bid];
                dis1= Tree[rank[ID1]].disInf[inf_i];
                dis2= QueryOverlay(bid,ID2);
                if(d>dis1+dis2){
                    cout<<ID1<<" "<<bid<<" "<<ID2<<": "<<d<<" "<<dis1+dis2<<endl;
                    disf1=dis1,disf2=dis2;
                    hubf=bid;
                    d=dis1+dis2;
                }
            }
            int dDijk_s=Dijkstra(ID1,hubf,Neighbor), dDijk_t=Dijkstra(hubf,ID2,Neighbor), dDijk_st=Dijkstra(ID1,ID2,Neighbor);
            cout<<ID1<<" "<<hubf<<"("<<NodeOrder[hubf]<<","<<PartiTags[hubf].first<<") "<<ID2<<": "<<disf1<<" "<<disf2<<" "<<d<<" ; "<<dDijk_s<<" "<<dDijk_t<<" "<<dDijk_st<<endl;
        }

    }

    return d;
}

//Query within one partition, no-boundary
int Graph::QueryH2HPartition(int ID1, int ID2, int PID){
    if(ID1==ID2) return 0;
    if(PartiTag[ID1].first!=PID || PartiTag[ID2].first!=PID){
        cout<<"Wrong! ID1 and ID2 are not in the same partition! "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<" "<<PID<<endl; exit(1);
    }
    int r1=ranks[PID][IDMap[ID1]], r2=ranks[PID][IDMap[ID2]];
//    cout<<"ID1: "<<ID1<<" "<<IDMap[ID1]<<" "<<r1<<" ; ID2: "<<ID2<<" "<<IDMap[ID2]<<" "<<r2<<endl;
    int LCA=LCAQueryPartition(r1,r2,PID);
//    cout<<"LCA: "<<LCA<<endl;
    if(LCA==r1)
        return Trees[PID][r2].dis[Trees[PID][r1].pos.back()];
    else if(LCA==r2)
        return Trees[PID][r1].dis[Trees[PID][r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<Trees[PID][LCA].pos.size();i++){
            if(tmp>Trees[PID][r1].dis[Trees[PID][LCA].pos[i]]+Trees[PID][r2].dis[Trees[PID][LCA].pos[i]])
                tmp=Trees[PID][r1].dis[Trees[PID][LCA].pos[i]]+Trees[PID][r2].dis[Trees[PID][LCA].pos[i]];
        }
        return tmp;
    }
}

//Query within one partition, no-boundary
int Graph::QueryH2HPartitionPost(int ID1, int ID2, int PID){
    if(ID1==ID2) return 0;
    if(PartiTag[ID1].first!=PID || PartiTag[ID2].first!=PID){
        cout<<"Wrong! ID1 and ID2 are not in the same partition! "<<PartiTag[ID1].first<<" "<<PartiTag[ID2].first<<" "<<PID<<endl; exit(1);
    }
    int r1=ranksPost[PID][IDMap[ID1]], r2=ranksPost[PID][IDMap[ID2]];
//    cout<<"ID1: "<<ID1<<" "<<IDMap[ID1]<<" "<<r1<<" ; ID2: "<<ID2<<" "<<IDMap[ID2]<<" "<<r2<<endl;
    int LCA=LCAQueryPartitionPost(r1,r2,PID);
//    cout<<"LCA: "<<LCA<<endl;
    if(LCA==r1)
        return TreesPost[PID][r2].dis[TreesPost[PID][r1].pos.back()];
    else if(LCA==r2)
        return TreesPost[PID][r1].dis[TreesPost[PID][r2].pos.back()];
    else{
        int tmp=INF;
        for(int i=0;i<TreesPost[PID][LCA].pos.size();i++){
            if(tmp>TreesPost[PID][r1].dis[TreesPost[PID][LCA].pos[i]]+TreesPost[PID][r2].dis[TreesPost[PID][LCA].pos[i]])
                tmp=TreesPost[PID][r1].dis[TreesPost[PID][LCA].pos[i]]+TreesPost[PID][r2].dis[TreesPost[PID][LCA].pos[i]];
        }
        return tmp;
    }
}

int Graph::LCAQueryPartition(int _p, int _q, int PID){
    int p = toRMQs[PID][_p], q = toRMQs[PID][_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Trees[PID][RMQIndexs[PID][k][p]].height < Trees[PID][RMQIndexs[PID][k][q]].height)
        return RMQIndexs[PID][k][p];
    else return RMQIndexs[PID][k][q];
}

int Graph::LCAQueryPartitionPost(int _p, int _q, int PID){
    int p = toRMQsPost[PID][_p], q = toRMQsPost[PID][_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (TreesPost[PID][RMQIndexsPost[PID][k][p]].height < TreesPost[PID][RMQIndexsPost[PID][k][q]].height)
        return RMQIndexsPost[PID][k][p];
    else return RMQIndexsPost[PID][k][q];
}

int Graph::LCAQueryOverlay(int _p, int _q){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}

int Graph::LCAQuery(int _p, int _q, vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree){
    int p = toRMQ[_p], q = toRMQ[_q];
    if (p > q){
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    int i = 1, k = 0;
    while (i * 2 < len){
        i *= 2;
        k++;
    }
    q = q - i + 1;
    if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
        return RMQIndex[k][p];
    else return RMQIndex[k][q];
}


/// Index maintenance
//H2H index update
void Graph::DecreaseOverlay(int a,int b, int newW, vector<unordered_map<vertex,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
//    map<int,int> checkedDis;//map<tree node ID, distance index>

    if(Neighbors[a].find(b)!=Neighbors[a].end()){
        Neighbors[a][b]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }
    if(Neighbors[b].find(a)!=Neighbors[b].end()){
        Neighbors[b][a]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }

    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> ss;//ss.clear();
    SCre.assign(ProH+1,ss);

    //map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
    //DisRe.clear();

    int MinH;

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){
            if(Tree[rank[lid]].vert[i].second.first>newW){
                Tree[rank[lid]].vert[i].second.first=newW;
                Tree[rank[lid]].vert[i].second.second=1;
                tri=true;
                SCre[ProH].insert(hid);
                MinH=IniH;
            }else if(Tree[rank[lid]].vert[i].second.first==newW){
                Tree[rank[lid]].vert[i].second.second+=1;
            }
            break;
        }
    }

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);

    //int ProBeginH;
    int ProBeginID;
    if(tri){
        //cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
        while(ProH>=MinH){

            ProIDRecord[ProH]=ProID;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
            bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw;//=OCdis[make_pair(ProID,Cid)];
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

                if(Tree[rank[ProID]].dis[cidH]>=Cw){
                    Tree[rank[ProID]].dis[cidH]=Cw;
                    Tree[rank[ProID]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProID]].DisRe.insert(Cid);
                    //DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
                }

                int hid,hidHeight,lid,lidHeight,wsum;
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.first=wsum;
                            Tree[rank[Cid]].vert[j].second.second=1;
                            SCre[Tree[rank[Cid]].height].insert(hid);
                            if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;

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
                                SCre[Tree[rank[lid]].height].insert(Cid);
                                if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;

                            }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                                Tree[rank[lid]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

            if(ProIDdisCha){//if the distance labeling is detected changed
                vertexIDChL.insert(ProID);
                //ProBeginH=ProH;
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }

        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        //top-down process
        EachNodeProBDis5(rank[ProBeginID], linee, vertexIDChL,  Tree, rank);
    }
    //return checkedDis.size();
}

void Graph::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank){
    bool ProIDdisCha=false;

    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
//                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }

                }
            }
        }
    }else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
//                    checkedDis.insert(make_pair(child,i));
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
//        sm->wait();
        vertexIDChL.insert(Tree[child].uniqueVertex);
//        vUpdated[Tree[child].uniqueVertex] = true;
//        sm->notify();
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,Tree, rank);
    }
    line.pop_back();

}

void Graph::EachNodeProBDis5PostMHLOverlay(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank){
    int ID=Tree[child].uniqueVertex;
    if(PartiTags[ID].first == -1){

        bool flagDebug=false;
//        if(ID==148515){
//            cout<<ID<<" "<<Tree[child].DisRe.size()<<" "<<line.size()<<" "<<Tree[child].vert.size()<<endl;
//            flagDebug=true;
//        }
        bool ProIDdisCha=false;
        if(Tree[child].DisRe.size()!=0){
            for(int k=0;k<Tree[child].vert.size();k++){
                int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
                if(Tree[child].FN[bH]){
                    if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                        for(int i=0;i<bH;i++){
//                        checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                        for(int i=bH+1;i<line.size();i++){
//                        checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                                Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }

                    }else{//partial ancestor check

                        if(vertexIDChL.find(b)!=vertexIDChL.end()){
                            for(int i=0;i<bH;i++){
//                            checkedDis.insert(make_pair(child,i));
                                if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                    Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                    Tree[child].FN[i]=false;
                                    ProIDdisCha=true;
                                }
                            }
                        }
                        for(int i=bH+1;i<line.size();i++){
//                        checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
                                Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }

                    }
                }
            }
        }else{
            for(int k=0;k<Tree[child].vert.size();k++){
                int b=Tree[child].vert[k].first, bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
                if(Tree[child].FN[bH]){
                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
//                        checkedDis.insert(make_pair(child,i));
                            if(flagDebug){
                                if(Tree[child].vAncestor[i] == 112485){
                                    cout<<"1: "<<ID<<" "<<Tree[child].vAncestor[i]<<" "<<Tree[rank[b]].vAncestor[i]<<": "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[b]].dis[i]<<endl;
                                    int d2= Dijkstra(b,Tree[rank[b]].vAncestor[i],Neighbor);
                                    if(d2!=Tree[rank[b]].dis[i]){
                                        cout<<"1. Incorrect! "<<b<<" "<<Tree[rank[b]].vAncestor[i]<<" "<<Tree[rank[b]].dis[i]<<" "<<d2<<endl; exit(1);
                                    }
                                }
                            }
                            if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                    checkedDis.insert(make_pair(child,i));
                        if(flagDebug){
                            if(Tree[child].vAncestor[i] == 112485){
                                cout<<"2: "<<ID<<" "<<Tree[child].vAncestor[i]<<" "<<line[i]<<": "<<Tree[child].dis[i]<<" "<<vbW+Tree[rank[line[i]]].dis[bH]<<endl;
                                int d2= Dijkstra(b,line[i],Neighbor);
                                if(d2!=Tree[rank[line[i]]].dis[bH]){
                                    cout<<"2. Incorrect! "<<line[i]<<" "<<b<<" "<<Tree[rank[line[i]]].dis[bH]<<" "<<d2<<endl; exit(1);
                                }
                            }
                        }
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
            vUpdated[Tree[child].uniqueVertex] = true;
        }

        line.push_back(Tree[child].uniqueVertex);
        for(int i=0;i<Tree[child].ch.size();i++){
            EachNodeProBDis5PostMHLOverlay(Tree[child].ch[i], line, vertexIDChL,Tree, rank);
        }
        line.pop_back();
    }
    else{//if partition vertex
        int pid=PartiTags[ID].first;
        vector<int> ProBeginVertexSetNew;
        ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetPartiExtend[pid].size()+1);
        ProBeginVertexSetNew.push_back(ID);
        int rnew=rank[ID],r;
        for(int i=0;i<ProBeginVertexSetPartiExtend[pid].size();i++){
            r=rank[ProBeginVertexSetPartiExtend[pid][i]];
            if(LCAQueryOverlay(rnew,r)!=rnew){
                ProBeginVertexSetNew.push_back(ProBeginVertexSetPartiExtend[pid][i]);
            }
        }
        ProBeginVertexSetPartiExtend[pid]=ProBeginVertexSetNew;
    }

}

void Graph::DecreaseParti(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax){
    map<int,int> checkedDis;//map<tree node ID, distance index>

    for(int i=0;i<Neighbors[a].size();i++){
        if(Neighbors[a][i].first==b){
            Neighbors[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbors[b].size();i++){
        if(Neighbors[b][i].first==a){
            Neighbors[b][i].second=newW;
            break;
        }
    }

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }
    int lidM=IDMap[lid];
    int IniH=Tree[rank[lidM]].height;//the height where weight change begins
    int ProH=Tree[rank[lidM]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> ss;//ss.clear();
    SCre.assign(ProH+1,ss);

    //map<int,set<int>> DisRe;//rankid; record the distance change caused by the shortcut in each height
    //DisRe.clear();

    int MinH;

    bool tri=false;
    for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
        if(Tree[rank[lidM]].vert[i].first==hid){
            if(Tree[rank[lidM]].vert[i].second.first>newW){
                Tree[rank[lidM]].vert[i].second.first=newW;
                Tree[rank[lidM]].vert[i].second.second=1;
                tri=true;
                SCre[ProH].insert(hid);
                MinH=IniH;
            }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                Tree[rank[lidM]].vert[i].second.second+=1;
            }
            break;
        }
    }

    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed

    vector<int> ProIDRecord; ProIDRecord.assign(ProH+1,0);

    //int ProBeginH;
    int ProBeginID;
    if(tri){
//        cout<<"Bottom-up ;;;;;;;;;;;;;;;;;; "<<endl;
        while(ProH>=MinH){

            ProIDRecord[ProH]=ProID;
            int ProIDM=IDMap[ProID];
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
            bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw;//=OCdis[make_pair(ProID,Cid)];
                int CidM=IDMap[Cid];
                int cidH=Tree[rank[CidM]].height-1;

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

                if(Tree[rank[ProIDM]].dis[cidH]>=Cw){
                    Tree[rank[ProIDM]].dis[cidH]=Cw;
                    Tree[rank[ProIDM]].FN[cidH]=true;
                    ProIDdisCha=true;
                    Tree[rank[ProIDM]].DisRe.insert(Cid);
                    //DisRe[rank[ProID]].insert(Cid); //cout<<"dischange Cid "<<Cid<<endl;
                }

                int hid,hidHeight,lid,lidHeight,wsum;
                for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                    hid=Tree[rank[CidM]].vert[j].first;hidHeight=Tree[rank[IDMap[hid]]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.first=wsum;
                            Tree[rank[CidM]].vert[j].second.second=1;
                            SCre[Tree[rank[CidM]].height].insert(hid);
                            if(Tree[rank[CidM]].height<MinH) MinH=Tree[rank[CidM]].height;

                        }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.second+=1;
                        }

                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;
                    int lidM=IDMap[lid];
                    lidHeight=Tree[rank[lidM]].height-1;
                    for(int k=0;k<Tree[rank[lidM]].vert.size();k++){
                        if(Tree[rank[lidM]].vert[k].first==Cid){
                            wsum=Cw+Lnei[j].second;
                            if(Tree[rank[lidM]].vert[k].second.first>wsum){
                                Tree[rank[lidM]].vert[k].second.first=wsum;
                                Tree[rank[lidM]].vert[k].second.second=1;
                                SCre[Tree[rank[lidM]].height].insert(Cid);
                                if(Tree[rank[lidM]].height<MinH) MinH=Tree[rank[lidM]].height;

                            }else if(Tree[rank[lidM]].vert[k].second.first==wsum){
                                Tree[rank[lidM]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

            if(ProIDdisCha){//if the distance labeling is detected changed
                vertexIDChL.insert(ProID);
                //ProBeginH=ProH;
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProIDM]].pa].uniqueVertex;
        }

        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        //top-down process
        EachNodeProBDis5Parti(rank[IDMap[ProBeginID]], linee, vertexIDChL, Tree, rank);
    }
    else{
//        cout<<"Not trigger update! "<<lid<<" "<<hid<<endl;
    }
    //return checkedDis.size();
}

void Graph::EachNodeProBDis5Parti(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank){
    bool ProIDdisCha=false;
    int childID=Tree[child].uniqueVertex;
    int bHeight=BoundVertex[PartiTag[childID].first].size();

    if(Tree[child].DisRe.size()!=0){
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[IDMap[b]]].height-1,vbW=Tree[child].vert[k].second.first;
            int bM=IDMap[b];
            if(Tree[child].FN[bH]){
                if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//all ancestor check
                    for(int i=0;i<bH;i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
//                            if(i<bHeight){
//                                vUpdated[childID] = true;
//                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
//                            if(i<bHeight){
//                                vUpdated[childID] = true;
//                            }
                        }
                    }

                }else{//partial ancestor check

                    if(vertexIDChL.find(b)!=vertexIDChL.end()){
                        for(int i=0;i<bH;i++){
//                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                                Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                                Tree[child].FN[i]=false;
                                ProIDdisCha=true;
//                                if(i<bHeight){
//                                    vUpdated[childID] = true;
//                                }
                            }
                        }
                    }
                    for(int i=bH+1;i<line.size();i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                            Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
//                            if(i<bHeight){
//                                vUpdated[childID] = true;
//                            }
                        }
                    }

                }
            }
        }
    }else{
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first, bH=Tree[rank[IDMap[b]]].height-1,vbW=Tree[child].vert[k].second.first;
            int bM=IDMap[b];
            if(Tree[child].FN[bH]){
                if(vertexIDChL.find(b)!=vertexIDChL.end()){
                    for(int i=0;i<bH;i++){
//                        checkedDis.insert(make_pair(child,i));
                        if(Tree[child].dis[i]>vbW+Tree[rank[bM]].dis[i]){
                            Tree[child].dis[i]=vbW+Tree[rank[bM]].dis[i];
                            Tree[child].FN[i]=false;
                            ProIDdisCha=true;
//                            if(i<bHeight){
//                                vUpdated[childID] = true;
//                            }
                        }
                    }
                }
                for(int i=bH+1;i<line.size();i++){
//                    checkedDis.insert(make_pair(child,i));
                    if(Tree[child].dis[i]>vbW+Tree[rank[IDMap[line[i]]]].dis[bH]){
                        Tree[child].dis[i]=vbW+Tree[rank[IDMap[line[i]]]].dis[bH];
                        Tree[child].FN[i]=false;
                        ProIDdisCha=true;
//                        if(i<bHeight){
//                            vUpdated[childID] = true;
//                        }
                    }
                }
            }
        }
    }

    if(ProIDdisCha){
        vertexIDChL.insert(Tree[child].uniqueVertex);
        vUpdated[Tree[child].uniqueVertex]=true;
    }

    line.push_back(Tree[child].uniqueVertex);
    for(int i=0;i<Tree[child].ch.size();i++){
        EachNodeProBDis5Parti(Tree[child].ch[i], line, vertexIDChL, Tree, rank);
    }
    line.pop_back();

}

//batch update for original graph
void Graph::DecreaseH2HBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    //    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
//    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();
//    ProBeginVertexSetParti.assign(partiNum,vector<int>()); vertexIDChLParti.assign(partiNum,set<int>());

    int a,b,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

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

        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }


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

            if(algoUpdate!=PCH_No){
                if(PartiTags[ProID].first == -1){//if overlay vertex
                    if(Tree[rank[ProID]].dis[cidH]>Cw){
                        Tree[rank[ProID]].dis[cidH]=Cw;
                        Tree[rank[ProID]].FN[cidH]=true;
                        ProIDdisCha=true;
                        vUpdated[ProID]=true;
                        Tree[rank[ProID]].DisRe.insert(Cid);
                    }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                        Tree[rank[ProID]].FN[cidH]=true;
                    }
                }else{//if partition vertex
                    if(algoUpdate>=PH2H_Post){
                        if(PartiTags[Cid].first != -1){//if partition vertex
                            if(Tree[rank[ProID]].vAncestor[cidH] != Cid){
                                cout<<"Wrong! "<<Tree[rank[ProID]].vAncestor[cidH]<<" "<<Cid<<endl; exit(1);
                            }
                            if(Tree[rank[ProID]].disPost[cidH]>Cw){
                                Tree[rank[ProID]].disPost[cidH]=Cw;
                                Tree[rank[ProID]].FNPost[cidH]=true;
                                ProIDdisCha=true;
                                vUpdated[ProID]=true;
                                Tree[rank[ProID]].DisRePost.insert(Cid);
                            }else if(Tree[rank[ProID]].disPost[cidH]==Cw){
                                Tree[rank[ProID]].FNPost[cidH]=true;
                            }
                        }

                    }
                    if(algoUpdate==PH2H_Cross){
//                        cout<<"Extend"<<endl;
                        if(Tree[rank[ProID]].dis[cidH]>Cw){
                            Tree[rank[ProID]].dis[cidH]=Cw;
                            Tree[rank[ProID]].FN[cidH]=true;
                            ProIDdisCha=true;
                            vUpdated[ProID]=true;
                            Tree[rank[ProID]].DisRe.insert(Cid);
                        }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                            Tree[rank[ProID]].FN[cidH]=true;
                        }
                    }

                }

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
                        OC.insert(OrderCompMin(Cid));
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
                            OC.insert(OrderCompMin(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            Tree[rank[lid]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            if(PartiTags[ProID].first==-1){
                vertexIDChLOverlay.insert(ProID);
                ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
                ProBeginVertexSetNew.push_back(ProID);
                int rnew=rank[ProID],r;
                for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                    r=rank[ProBeginVertexSetOverlay[i]];
                    if(LCAQueryOverlay(rnew,r)!=rnew){
                        ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                    }
                }
                ProBeginVertexSetOverlay=ProBeginVertexSetNew;
            }else{
                int pid=PartiTags[ProID].first;
                affectedParti.insert(pid);
                vertexIDChLParti[pid].insert(ProID);
                ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
                ProBeginVertexSetNew.push_back(ProID);
                int rnew=rank[ProID],r;
                for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                    r=rank[ProBeginVertexSetParti[pid][i]];
                    if(LCAQueryOverlay(rnew,r)!=rnew){
                        ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                    }
                }
                ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
            }

        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;

    if(ifLabelU){
        int ProBeginVertexID;
        cout<<"ProBeginVertexSet size: "<<ProBeginVertexSetOverlay.size()<<endl;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            cout<<"ProBeginVertexID: "<<ProBeginVertexID<<"("<<PartiTags[ProBeginVertexID].first<<","<<Tree[rank[ProBeginVertexID]].height<<")"<<endl;
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChLOverlay,Tree,rank);
        }
    }

    //return checkedDis.size();
}

//batch update for overlay graph
void Graph::DecreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    //    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();

    int a,b,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

        if(Neighbor[a].find(b)!=Neighbor[a].end()){
            Neighbor[a][b]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }
        if(Neighbor[b].find(a)!=Neighbor[b].end()){
            Neighbor[b][a]=newW;
        }else{
            cout<<"Wrong for Neighbors!"<<endl; exit(1);
        }

        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }


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
                vUpdated[ProID]=true;
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
                        OC.insert(OrderCompMin(Cid));
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
                            OC.insert(OrderCompMin(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            Tree[rank[lid]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLOverlay.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                r=rank[ProBeginVertexSetOverlay[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChLOverlay,Tree,rank);
        }
    }

    //return checkedDis.size();
}

void Graph::DecreaseOverlayBatchPostMHL(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

//    for(int i=0;i<Tree.size();i++){
//        if(PartiTags[Tree[i].uniqueVertex].first==-1){
//            Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
//        }
//    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    //    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
//    ProBeginVertexSetOverlay.clear(); vertexIDChLOverlay.clear();

    int a,b,newW,lid,hid;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }

//        if(Neighbor[a].find(b)!=Neighbor[a].end()){
//            Neighbor[a][b]=newW;
//        }else{
//            cout<<"Wrong for Neighbors!"<<endl; exit(1);
//        }
//        if(Neighbor[b].find(a)!=Neighbor[b].end()){
//            Neighbor[b][a]=newW;
//        }else{
//            cout<<"Wrong for Neighbors!"<<endl; exit(1);
//        }
//
//        for(int i=0;i<NeighborsOverlayV[a].size();i++){
//            if(NeighborsOverlayV[a][i].first==b){
////            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
//                NeighborsOverlayV[a][i].second=newW;
//                break;
//            }
//        }
//        for(int i=0;i<NeighborsOverlayV[b].size();i++){
//            if(NeighborsOverlayV[b][i].first==a){
////            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
//                NeighborsOverlayV[b][i].second=newW;
//                break;
//            }
//        }

        for(int i=0;i<Tree[rank[lid]].vert.size();i++){
            if(Tree[rank[lid]].vert[i].first==hid){
                if(Tree[rank[lid]].vert[i].second.first>newW){
                    Tree[rank[lid]].vert[i].second.first=newW;
                    Tree[rank[lid]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
                    Tree[rank[lid]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }


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

            if(algoUpdate>PCH_No){
                if(Tree[rank[ProID]].dis[cidH]>Cw){
                    Tree[rank[ProID]].dis[cidH]=Cw;
                    Tree[rank[ProID]].FN[cidH]=true;
                    ProIDdisCha=true;
                    vUpdated[ProID]=true;
                    Tree[rank[ProID]].DisRe.insert(Cid);
                }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                    Tree[rank[ProID]].FN[cidH]=true;
                }
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
                        OC.insert(OrderCompMin(Cid));
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
                            OC.insert(OrderCompMin(lid));
                        }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                            Tree[rank[lid]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLOverlay.insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                r=rank[ProBeginVertexSetOverlay[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5PostMHLOverlay(rank[ProBeginVertexID], linee, vertexIDChLOverlay,Tree,rank);
        }
    }

    //return checkedDis.size();
}

void Graph::DecreaseOverlayBatchLabelPostMHL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL) {
    int ProBeginVertexID;
//    cout<<"ProBeginVertexSet size: "<<ProBeginVertexSet.size()<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
//        cout<<"ProBeginVertexID: "<<ProBeginVertexID<<"("<<PartiTags[ProBeginVertexID].first<<","<<Tree[rank[ProBeginVertexID]].height<<")"<<endl;
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5PostMHLOverlay(rank[ProBeginVertexID], linee, vertexIDChL,Tree,rank);
    }
}

void Graph::DecreasePartiBatchLabelPostMHLExtendV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>>& ProBeginVertexSetV, set<int> &vertexIDChL){
    int pid;
    for(int i=0;i<p.size();++i){
        pid=p[i];
        DecreasePartiBatchLabelPostMHLExtend(pid, Tree, rank, heightMax, ProBeginVertexSetV[pid], vertexIDChL);
    }
}

void Graph::DecreasePartiBatchLabelPostMHLExtend(int pid, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, set<int> &vertexIDChL) {
    int ProBeginVertexID=partiRoots[pid];
    int rootRank=rank[ProBeginVertexID];
//    vector<int> linee; //linee.clear();
//    linee=Tree[Tree[rootRank].pa].vAncestor;
//    EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,Tree,rank);
//    cout<<"ProBeginVertexSetPartiExtend of Parti "<<pid<<": "<<ProBeginVertexSetPartiExtend[pid].size()<<endl;

    set<int> vertexIDChLParti=vertexIDChL;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
//        cout<<"ProBeginVertexID: "<<ProBeginVertexID<<"("<<PartiTags[ProBeginVertexID].first<<","<<Tree[rank[ProBeginVertexID]].height<<")"<<endl;
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
//        EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,Tree,rank);
        EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChLParti,Tree,rank);
    }
}

void Graph::DecreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, set<int> &vertexIDChL) {
    int ProBeginVertexID;
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
        EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,Tree,rank);
    }
}


//batch update for partition graph of PH2H
void Graph::DecreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order



    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

//        if(pid==23){
//            cout<<pid<<": "<<a<<" "<<b<<" "<<newW<<endl;
//        }

        for(int i=0;i<Neighbors[a].size();i++){
            if(Neighbors[a][i].first==b){
                Neighbors[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbors[b].size();i++){
            if(Neighbors[b][i].first==a){
                Neighbors[b][i].second=newW;
                break;
            }
        }

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

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

            if(PartiTag[ProID].second){//if boundary vertex
                updatedSC.emplace_back(make_pair(ProID,Cid),Cw);
            }

            if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                Tree[rank[ProIDM]].dis[cidH]=Cw;
                Tree[rank[ProIDM]].FN[cidH]=true;
                ProIDdisCha=true;
                vUpdated[ProID]=true;
                Tree[rank[ProIDM]].DisRe.insert(Cid);
            }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                Tree[rank[ProIDM]].FN[cidH]=true;
            }

            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            Tree[rank[lid2M]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    if(ifLabelU){
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

void Graph::DecreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU, bool ifConstruct){
    map<int,int> checkedDis;

//    for(int i=0;i<Tree.size();i++){
//        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
//    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,oldW,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;oldW=wBatch[k].second.first;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

        for(int i=0;i<Neighbors[a].size();i++){
            if(Neighbors[a][i].first==b){
                if(Neighbors[a][i].second!=oldW){
                    cout<<"Old edge weight inconsistent! "<<a<<" "<<b<<" "<<oldW<<" "<<Neighbors[a][i].second<<endl; exit(1);
                }
                Neighbors[a][i].second=newW;
                break;
            }
        }
        for(int i=0;i<Neighbors[b].size();i++){
            if(Neighbors[b][i].first==a){
                if(Neighbors[b][i].second!=oldW){
                    cout<<"Old edge weight inconsistent! "<<b<<" "<<a<<" "<<oldW<<" "<<Neighbors[b][i].second<<endl; exit(1);
                }
                Neighbors[b][i].second=newW;
                break;
            }
        }

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
//                    if(pid==12){
//                        cout<<"Equal. "<<pid<<": "<<lid<<" "<<hid<<" "<<oldW<<" "<<newW<<" ("<<Tree[rank[lidM]].vert[i].second.second<<")"<<endl;
//                    }
                    if(!ifConstruct){//should not plus one if in index construction
//                        cout<<"Equal. "<<pid<<": "<<lid<<" "<<hid<<" "<<oldW<<" "<<newW<<" ("<<Tree[rank[lidM]].vert[i].second.second<<")"<<endl;
                        Tree[rank[lidM]].vert[i].second.second+=1;
                    }
                }else{
                    cout<<"Unexpected result! "<<Tree[rank[lidM]].vert[i].second.first<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }

    }

//    if(pid==12){
//        cout<<pid<<endl;
//    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
    if(!ifLabelU){//for PH2H_Post update
        ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
        for(int i=0;i<Tree.size();i++){
            Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
        }
    }

    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

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

            if(PartiTag[ProID].second){//if boundary vertex
                updatedSC.emplace_back(make_pair(ProID,Cid),Cw);
            }

//            if(algoUpdate!=PCH_No && !ifConstruct){
//                if(Tree[rank[ProIDM]].dis[cidH]>Cw){
//                    Tree[rank[ProIDM]].dis[cidH]=Cw;
//                    Tree[rank[ProIDM]].FN[cidH]=true;
//                    ProIDdisCha=true;
//                    vUpdated[ProID]=true;
//                    Tree[rank[ProIDM]].DisRe.insert(Cid);
//                }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
//                    Tree[rank[ProIDM]].FN[cidH]=true;
//                }
//            }


            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            Tree[rank[lid2M]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;

    if(ifLabelU){
//        cout<<pid<<": ProBeginVertexSet size: "<<ProBeginVertexSetParti[pid].size()<<endl;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

void Graph::DecreasePartiBatchPostMHLShortcut(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC){
    map<int,int> checkedDis;

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,oldW,newW,lid,hid;
    int hidHeight,lidHeight,wsum;
    int inf_i;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;oldW=wBatch[k].second.first;
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
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lid]].vert[i].second.first==newW){
//                    if(pid==12){
//                        cout<<"Equal. "<<pid<<": "<<lid<<" "<<hid<<" "<<oldW<<" "<<newW<<" ("<<Tree[rank[lidM]].vert[i].second.second<<")"<<endl;
//                    }
                    Tree[rank[lid]].vert[i].second.second+=1;
                }else{
                    cout<<"Unexpected result! "<<Tree[rank[lid]].vert[i].second.first<<" "<<newW<<endl; exit(1);
                }
                break;
            }
        }

    }

//    if(pid==12){
//        cout<<pid<<endl;
//    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distanc labeling has changed
//    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
//    ProBeginVertexSetPartiExtend[pid].clear(); vertexIDChLPartiExtend[pid].clear();
//    for(int i=0;i<Tree.size();i++){
//        if(PartiTags[Tree[i].uniqueVertex].first==pid){
//            Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
//            Tree[i].DisRePost.clear();//record the star weight change (causing the distance change)
//        }
//    }


    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID;
    int MinHInf=INF;
    set<int> ProBeginIDInf;
//    set<int> affectedParti;
    bool ifAffected=false;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        bool ProIDdisChaExtend=false;//to see if the distance labeling of proID change or not

        if(PartiTags[ProID].first==-1){
            cout<<"Wrong! Overlay vertex! "<<ProID<<endl; exit(1);
        }

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

            //if Cid is not boundary vertex
            if(algoUpdate!=PCH_No){
                if(algoUpdate>=PH2H_Post){
                    if(PartiTags[Cid].first != -1){
                        if(Tree[rank[ProID]].vAncestor[cidH] != Cid){
                            cout<<"Wrong! "<<Tree[rank[ProID]].vAncestor[cidH]<<" "<<Cid<<endl; exit(1);
                        }
                        if(Tree[rank[ProID]].disPost[cidH]>Cw){
                            Tree[rank[ProID]].disPost[cidH]=Cw;
                            Tree[rank[ProID]].FNPost[cidH]=true;
                            ProIDdisCha=true;
                            vUpdated[ProID]=true;
                            Tree[rank[ProID]].DisRePost.insert(Cid);
                        }else if(Tree[rank[ProID]].disPost[cidH]==Cw){
                            Tree[rank[ProID]].FNPost[cidH]=true;
                        }
                    }
                }
                if(algoUpdate==PH2H_Cross){
//                    cout<<"Extend"<<endl;
                    if(Tree[rank[ProID]].dis[cidH]>Cw){
                        Tree[rank[ProID]].dis[cidH]=Cw;
                        Tree[rank[ProID]].FN[cidH]=true;
                        ProIDdisChaExtend=true;
                        vUpdated[ProID]=true;
                        Tree[rank[ProID]].DisRe.insert(Cid);
                    }else if(Tree[rank[ProID]].dis[cidH]==Cw){
                        Tree[rank[ProID]].FN[cidH]=true;
                    }
                }
            }

            ///deal with the boundary vertex, if the higher endpoint Cid is boundary vertex
            if(PartiTags[Cid].first == -1){
//                if(Tree[rank[ProID]].height<MinHInf){
//                    MinHInf = Tree[rank[ProID]].height;
//                    ProBeginIDInf.clear();
//                    ProBeginIDInf.insert(ProID);
//                }else if(Tree[rank[ProID]].height==MinHInf){
//                    ProBeginIDInf.insert(ProID);
//                }
                int wsum,lidHeight;
                /// For Hnei, update the super edges between Cid and other interface vertex, i.e., the shortcuts in AdjaCore
                for(int i=0;i<Tree[rank[Cid]].vert.size();++i){//NeighborsOverlayV
                    hid=Tree[rank[Cid]].vert[i].first;
//                for(int i=0;i<NeighborsOverlayV[Cid].size();++i){
//                    hid=NeighborsOverlayV[Cid][i].first;
                    if(Hnei.find(hid)!=Hnei.end()){//if hid is also the higher-order neighbor of ProID (a)
                        wsum=Cw+Hnei[hid];
//                        if(NeighborsOverlay[Cid].find(hid) == NeighborsOverlay[Cid].end()){
//                            cout<<"Cannot find hid in Cid's adjacency list!!"<<endl;
//                            exit(1);
//                        }else if(NeighborsOverlay[hid].find(Cid) == NeighborsOverlay[hid].end()){
//                            cout<<"Cannot find Cid in hid's adjacency list!!"<<endl;
//                            exit(1);
//                        }

//                        if(wsum<NeighborsOverlay[Cid][hid]){//update neighbor distance of Cid
                        if(wsum<Tree[rank[Cid]].vert[i].second.first){
//                            Tree[rank[Cid]].vert[i].second.first=wsum;
//                            Tree[rank[Cid]].vert[i].second.second=1;
//                            SCre[Cid].insert(hid);
//                            OC.insert(OrderCompMin(Cid));
//                            vSm[Cid]->wait(); vSm[hid]->wait();
//                            NeighborsOverlay[Cid][hid]=wsum;
//                            NeighborsOverlay[hid][Cid]=wsum;
//                            vSm[Cid]->notify(); vSm[hid]->notify();
//                            sm->wait();
                            updatedSC.emplace_back(make_pair(Cid,hid),wsum);
//                            sm->notify();
                            ///identify the affected partitions
                            if(Cid<hid){
                                if(SuppPartiID[Cid][hid].size()>1){
                                    for(auto it=SuppPartiID[Cid][hid].begin();it!=SuppPartiID[Cid][hid].end();++it){
                                        if(it->first != pid){
                                            sm->wait();
                                            affectedParti.insert(it->first);
                                            sm->notify();
                                        }
                                    }
                                }
                            }else{
                                if(SuppPartiID[hid][Cid].size()>1){
                                    for(auto it=SuppPartiID[hid][Cid].begin();it!=SuppPartiID[hid][Cid].end();++it){
                                        if(it->first != pid){
                                            sm->wait();
                                            affectedParti.insert(it->first);
                                            sm->notify();
                                        }
                                    }
                                }
                            }

                        }
//                        else if(wsum==Tree[rank[Cid]].vert[i].second.first){
//                            Tree[rank[Cid]].vert[i].second.second+=1;
//                        }
                    }
                }
                /// For Lnei, update the super edges of Lnei vertex
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;
                    if(PartiTags[lid].first == -1){//if lid is interface vertex
                        for(int i=0;i<Tree[rank[lid]].vert.size();++i) {//NeighborsOverlayV
                            int vertid = Tree[rank[lid]].vert[i].first;
//                        for(int i=0;i<NeighborsOverlayV[lid].size();++i) {
//                            int vertid = NeighborsOverlayV[lid][i].first;
                            if (vertid == Cid) {//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                                wsum=Cw+Lnei[j].second;
//                                if(NeighborsOverlay[Cid].find(lid) == NeighborsOverlay[Cid].end()){
//                                    cout<<"Cannot find lid in Cid's adjacency list!!"<<endl;
//                                    exit(1);
//                                }else if(NeighborsOverlay[lid].find(Cid) == NeighborsOverlay[lid].end()){
//                                    cout<<"Cannot find Cid in lid's adjacency list!!"<<endl;
//                                    exit(1);
//                                }
//                                if(NeighborsOverlay[lid][Cid]>wsum){//update neighbor distance of Lnei
                                if(Tree[rank[lid]].vert[i].second.first>wsum){//update neighbor distance of Lnei
//                                    Tree[rank[lid]].vert[i].second.first=wsum;
//                                    Tree[rank[lid]].vert[i].second.second=1;
//                                    SCre[lid].insert(Cid);//if lid has lower order than ProID, the inserted Cid will not be processed
//                                    OC.insert(OrderCompMin(lid));
//                                    vSm[lid]->wait(); vSm[Cid]->wait();
//                                    NeighborsOverlay[lid][Cid]=wsum;
//                                    NeighborsOverlay[Cid][lid]=wsum;
//                                    vSm[lid]->notify(); vSm[Cid]->notify();
//                                    sm->wait();
                                    updatedSC.emplace_back(make_pair(lid,Cid),wsum);
//                                    sm->notify();
                                    ///identify the affected partitions
                                    if(lid<Cid){
                                        if(SuppPartiID[lid][Cid].size()>1){
                                            for(auto it=SuppPartiID[lid][Cid].begin();it!=SuppPartiID[lid][Cid].end();++it){
                                                if(it->first != pid){
                                                    sm->wait();
                                                    affectedParti.insert(it->first);
                                                    sm->notify();
                                                }
                                            }
                                        }
                                    }else {
                                        if(SuppPartiID[Cid][lid].size()>1){
                                            for(auto it=SuppPartiID[Cid][lid].begin();it!=SuppPartiID[Cid][lid].end();++it){
                                                if(it->first != pid){
                                                    sm->wait();
                                                    affectedParti.insert(it->first);
                                                    sm->notify();
                                                }
                                            }
                                        }
                                    }

                                }
//                                else if(Tree[rank[lid]].vert[i].second.first==wsum){
//                                    Tree[rank[lid]].vert[i].second.second+=1;
//                                }
                                break;
                            }
                        }
                    }else{//if lid is periphery vertex
                        lidHeight=Tree[rank[lid]].height-1;
                        for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                            int vertid=Tree[rank[lid]].vert[k].first;
                            if(vertid==Cid){//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                                wsum=Cw+Lnei[j].second;
                                if(Tree[rank[lid]].vert[k].second.first>wsum){//update neighbor distance of Lnei
                                    Tree[rank[lid]].vert[k].second.first=wsum;
                                    Tree[rank[lid]].vert[k].second.second=1;
                                    SCre[lid].insert(Cid);//if lid has lower order than ProID, the inserted Cid will not be processed
                                    OC.insert(OrderCompMin(lid));
//                                    if(Tree[rank[lid]].height<MinH)
//                                        MinH=Tree[rank[lid]].height;
                                    if(algoUpdate!=PCH_No) {
                                        if (Tree[rank[lid]].height < MinHInf) {
                                            MinHInf = Tree[rank[lid]].height;
                                            ProBeginIDInf.clear();
                                            ProBeginIDInf.insert(lid);
                                        } else if (Tree[rank[lid]].height == MinHInf) {
                                            ProBeginIDInf.insert(lid);
                                        }
                                        /// interface
                                        if (PartiTags[vertid].first == -1) {//if the neighbor is interface vertex, never enter this branch
//                                            assert(Tree[rank[lid]].disInf.find(vertid) != Tree[rank[lid]].disInf.end());
//                                            if (Tree[rank[lid]].disInf[vertid] > wsum) {
//                                                Tree[rank[lid]].DisReInf.insert(vertid);//record the vertex id that the interface label should be updated
//                                                Tree[rank[lid]].disInf[vertid] = wsum;
//                                            }
                                        }
                                    }
                                }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                                    Tree[rank[lid]].vert[k].second.second+=1;
                                }

                                break;
                            }
                        }
                    }
                }
            }
            else{
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid=Tree[rank[Cid]].vert[j].first;hidHeight=Tree[rank[hid]].height-1;
                    if(Hnei.find(hid)!=Hnei.end()){
                        wsum=Cw+Hnei[hid];
                        if(wsum<Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.first=wsum;
                            Tree[rank[Cid]].vert[j].second.second=1;
                            SCre[Cid].insert(hid);
                            OC.insert(OrderCompMin(Cid));

                            /// interface
                            if(algoUpdate!=PCH_No) {
                                if (PartiTags[hid].first == -1) {//if the neighbor is interface vertex
//                                assert(Tree[rank[Cid]].disInf.find(hid) != Tree[rank[Cid]].disInf.end());
//                                if (Tree[rank[Cid]].disInf[hid] > wsum) {
////                                    Tree[rank[Cid]].DisReInf.insert(hid);//record the vertex id that the interface label should be updated
//                                    if (Tree[rank[Cid]].height < MinHInf) {
//                                        MinHInf = Tree[rank[Cid]].height;
//                                        ProBeginIDInf.clear();
//                                        ProBeginIDInf.insert(Cid);
//                                    } else if (Tree[rank[Cid]].height == MinHInf) {
//                                        ProBeginIDInf.insert(Cid);
//                                    }
////                                    Tree[rank[Cid]].disInf[hid] = wsum;
//                                }
                                }
                            }

                        }else if(wsum==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second+=1;
                        }
                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid=Lnei[j].first;
                    lidHeight=Tree[rank[lid]].height-1;
                    for(int k=0;k<Tree[rank[lid]].vert.size();k++){
                        if(Tree[rank[lid]].vert[k].first==Cid){
                            wsum=Cw+Lnei[j].second;
                            if(Tree[rank[lid]].vert[k].second.first>wsum){
                                Tree[rank[lid]].vert[k].second.first=wsum;
                                Tree[rank[lid]].vert[k].second.second=1;
                                SCre[lid].insert(Cid);
                                OC.insert(OrderCompMin(lid));
                            }else if(Tree[rank[lid]].vert[k].second.first==wsum){
                                Tree[rank[lid]].vert[k].second.second+=1;
                            }

                            break;
                        }
                    }
                }
            }

        }



        if(ProIDdisCha){//if the distance labeling is dectected changed
            ifAffected=true;
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[ProBeginVertexSetParti[pid][i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
        if(ProIDdisChaExtend){//if the distance labeling is dectected changed
            sm->wait();
            vertexIDChLOverlay.insert(ProID);
            sm->notify();
        }
    }

    if(ifAffected){
        sm->wait();
        affectedParti.insert(pid);
        sm->notify();
    }

    //cout<<"Finish bottom-up refresh"<<endl;

    //return checkedDis.size();
}
//function of identifying the affected partitions and update the shortcut shortcuts
void Graph::BoundShortcutsCheck(bool ifParallel, bool ifIncrease){
//    map<pair<int,int>,int> weightsParti;
    vector<vector<pair<pair<int,int>,int>>> affectedShortcuts;
    affectedShortcuts.assign(partiNum,vector<pair<pair<int,int>,int>>());
//    cout<<"Before. affectedParti size: "<<affectedParti.size()<<endl;
    if(ifParallel){//use multi-thread
        boost::thread_group thread;
        vector<int> vParti;
        for(int pid=0;pid<partiNum;++pid){
//            if(affectedParti.find(pid)==affectedParti.end()){//if not found
                vParti.push_back(pid);
//            }
        }
//                    cout<<"vParti size: "<<vParti.size()<<endl;
        if(vParti.size()>threadnum){
            int step=vParti.size()/threadnum;
            boost::thread_group threadf;
            for(int i=0;i<threadnum;i++){
                pair<int,int> p;
                p.first=i*step;
                if(i==threadnum-1)
                    p.second=vParti.size();
                else
                    p.second=(i+1)*step;
                threadf.add_thread(new boost::thread(&Graph::BoundShortcutsCheckPartiV, this, p, boost::ref(vParti), boost::ref(affectedShortcuts), ifIncrease));
            }
            threadf.join_all();
        }else{
            boost::thread_group threadf;
            for(int i=0;i<vParti.size();i++){
                threadf.add_thread(new boost::thread(&Graph::BoundShortcutsCheckParti, this, vParti[i], boost::ref(affectedShortcuts[vParti[i]]), ifIncrease));
            }
            threadf.join_all();
        }
    }
    else{//use single-thread
        cout<<"Single thread for boundary shortcut check."<<endl;
        for(int pid=0;pid<partiNum;++pid){
//            if(affectedParti.find(pid)==affectedParti.end()){//if not found
//                cout<<"Boundary check of partition "<<pid<<endl;
                BoundShortcutsCheckParti(pid, affectedShortcuts[pid], ifIncrease);
//            }
        }
    }
//    cout<<"After. affectedParti size: "<<affectedParti.size()<<endl;
//    cout<<"Size of changed boundary shortcuts: "<<weightsParti.size()<<endl;
//    affectedParti.clear();
    for(int pid=0;pid<partiNum;++pid){
        if(!affectedShortcuts[pid].empty()){
            affectedParti.insert(pid);
            for(auto it=affectedShortcuts[pid].begin();it!=affectedShortcuts[pid].end();++it){
                BoundShortcuts[it->first.first][it->first.second]=it->second;
            }
        }
    }

}

void Graph::BoundShortcutsCheckPartiV(pair<int,int> pRange, vector<int>& pids, vector<vector<pair<pair<int,int>,int>>>& affectedShortcuts, bool ifIncrease){
    for(int i=pRange.first;i<pRange.second;++i){
        int pid = pids[i];
        BoundShortcutsCheckParti( pid, affectedShortcuts[pid], ifIncrease);
    }
}

void Graph::BoundShortcutsCheckParti(int pid, vector<pair<pair<int,int>,int>>& affectedShortcut, bool ifIncrease){
    int ID1, ID2, wlocal, woverlay;
    for(int i=0;i<BoundVertex[pid].size();i++){
        ID1=BoundVertex[pid][i];
        for(int j=i+1;j<BoundVertex[pid].size();j++){
            ID2=BoundVertex[pid][j];

            if(ID1<=ID2){
                if(BoundShortcuts[ID1].find(ID2)!=BoundShortcuts[ID1].end()){//if found
                    wlocal=BoundShortcuts[ID1][ID2];
                }else{
                    cout<<pid<<": Not found edge e("<<ID1<<","<<ID2<<") in overlay graph!"<<endl; exit(1);
                }
            }else{
                if(BoundShortcuts[ID2].find(ID1)!=BoundShortcuts[ID2].end()){//if found
                    wlocal=BoundShortcuts[ID2][ID1];
                }else{
                    cout<<pid<<": Not found edge e("<<ID2<<","<<ID1<<") in overlay graph!"<<endl; exit(1);
                }
            }
            woverlay=QueryOverlay(ID1,ID2);
//            int d2= Dijkstra(ID1,ID2,Neighbor);
//            if(d2!=woverlay){
//                cout<<"Wrong for boundary shortcut! "<<pid<<" "<<ID1<<" "<<ID2<<" "<<woverlay<<" "<<d2<< endl; exit(1);
//            }
            bool found=false;//whether the boundary edge exist or not
            int wei;
            if(ifIncrease){
                if(woverlay>wlocal){
//                cout<<"bound check. "<<pid<<": "<<ID1<<" "<<ID2<<" "<<wlocal<<" "<<woverlay<<endl;
//                    sm->wait();
//                    affectedParti.insert(pid);
                    if(ID1<ID2){
                        affectedShortcut.push_back({make_pair(ID1,ID2),woverlay});
                    }else{
                        affectedShortcut.push_back({make_pair(ID2,ID1),woverlay});
                    }
//                    sm->notify();
                }else if(woverlay<wlocal){
                    cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
                }
            }else{
                if(woverlay<wlocal){
//                cout<<"bound check. "<<pid<<": "<<ID1<<" "<<ID2<<" "<<wlocal<<" "<<woverlay<<endl;
//                    sm->wait();
//                    affectedParti.insert(pid);
                    if(ID1<ID2){
                        affectedShortcut.push_back({make_pair(ID1,ID2),woverlay});
                    }else{
                        affectedShortcut.push_back({make_pair(ID2,ID1),woverlay});
                    }
//                    sm->notify();
                }else if(woverlay>wlocal){
                    cout<<"Something wrong: shortest path in the overlay graph rather than in the subgraph. "<<ID1<<" "<<ID2<<" "<< wlocal<<" "<<woverlay<< endl; exit(1);
                }
            }


        }
    }

//    if(!affectedShortcut.empty()){
//        sm->wait();
//        affectedParti.insert(pid);
//        sm->notify();
//    }
}

void Graph::InterfaceEntryDecreaseUpdate(bool ifParallel){
    ///interface entry update
    if(affectedParti.size()>=1){
//            cout<<"Flag 3: affectedParti number: "<<affectedParti.size()<<endl;
        if(ifParallel){//use multi-thread
            boost::thread_group thread;
            vector<int> vParti;
            for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                vParti.push_back(*it);
//                PartiUpdateExt[*it] = true;//if the interface entry needs to be updated, we may also need to update the extended labels
            }
//                    cout<<"vParti size: "<<vParti.size()<<endl;
            if(vParti.size()>threadnum){
                int step=vParti.size()/threadnum;
                boost::thread_group threadf;
                for(int i=0;i<threadnum;i++){
                    pair<int,int> p;
                    p.first=i*step;
                    if(i==threadnum-1)
                        p.second=vParti.size();
                    else
                        p.second=(i+1)*step;
                    threadf.add_thread(new boost::thread(&Graph::InterfacePropagateParallel, this, p, boost::ref(vParti), false));
                }
                threadf.join_all();
            }else{
                boost::thread_group threadf;
                for(int i=0;i<vParti.size();i++){
                    threadf.add_thread(new boost::thread(&Graph::InterfacePropagateParallel, this, make_pair(i,i+1), boost::ref(vParti), false));
                }
                threadf.join_all();
            }

        }
        else{//use single-thread
            if(!affectedParti.empty()){//if we need to update the shortcut of other partition
                for(auto it=affectedParti.begin();it!=affectedParti.end();++it){
                    int Pid = *it;
                    int rootID = partiRoots[Pid];
                    vector<int> interfaceP;
//                    PartiUpdateExt[*it] = true;//if the interface entry needs to be updated, we may also need to update the extended labels
//                    for(auto it2=Tree[rank[rootID]].vert.begin();it2!=Tree[rank[rootID]].vert.end();++it2){
//                        interfaceP.emplace_back(it2->first);
//                    }
                    interfaceP=BoundVertex[Pid];
                    InterfacePropagate( rank[rootID], interfaceP,Tree,false);
                }
            }
        }
    }
}

void Graph::AncestorEntryDecreaseUpdateParti(int pid) {
//    int ProBeginVertexID;
    map<int,int> checkedDis;//map<tree node ID, distance index>
    vector<int> interface;
//    int ProBeginID;
//    ProBeginID=partiRoots[pid];//ProBeginIDV[pid];
    int rootVertex=partiRoots[pid];
    int rootRank=rank[rootVertex];
    int rootHeight=Tree[rootRank].height-1;
//    rootHeight=0;
//    for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
//        ProBeginVertexID=ProBeginVertexSetParti[pid][i];
//        vector<int> linee; //linee.clear();
//        linee.reserve(heightMax);
//        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
//        while(Tree[rank[pachidd]].height>1){
//            linee.insert(linee.begin(),pachidd);
//            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
//        }
//        linee.insert(linee.begin(),pachidd);
//        AncestorEntryDecreaseUpdate(pid, rank[ProBeginVertexID], linee, interface, vertexIDChLParti[pid], checkedDis, Tree, rank, rootHeight);
//    }

    vertexIDChLParti[pid].clear();
    interface=BoundVertex[pid];

    vector<int> linee; //linee.clear();
//    linee.reserve(heightMax);
//    int pachidd=ProBeginID;
//    while(Tree[rank[pachidd]].height>1){
//        linee.insert(linee.begin(),pachidd);
//        pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
//    }
//    linee.insert(linee.begin(),pachidd);
    linee=Tree[Tree[rootRank].pa].vAncestor;
    linee.push_back(rootVertex);
//    cout<<"Ancestor size of "<<pid<<" : "<<linee.size()<<" ; Root vertex: "<< rootVertex<<" ; Root height: "<<rootHeight<<endl;

    for(int i=0;i<Tree[rootRank].ch.size();++i){
        AncestorEntryDecreaseUpdate(pid, Tree[rootRank].ch[i], linee, interface, vertexIDChLParti[pid], checkedDis, Tree, rank, rootHeight);
    }


}

//function for edge decrease update
void Graph::AncestorEntryDecreaseUpdate(int pid, int child,vector<int>& line, vector<int>& interfaces, set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank, int rootHeight){
    bool ProIDdisCha=false;
    if(rootHeight >= line.size()){
        cout<<"Wrong root height! "<<rootHeight<<" "<<line.size()<<endl; exit(1);
    }
    int ID=Tree[child].uniqueVertex;
//    cout<<pid<<": "<<ID<<" "<<PartiTags[ID].first<<" "<<Tree[rank[ID]].height-1<<endl;
    if(PartiTags[ID].first!=pid){
        cout<<"This vertex does not belong to this partition! "<<pid<<": "<<ID<<"("<<NodeOrder[ID]<<","<<PartiTags[ID].first<<")"<<endl; exit(1);
    }
    int inf_i;
    /// update ancestor entry
    if(!Tree[child].DisRe.empty()){//Tree[child].DisRe.size()!=0
        /// update ancestor entries
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first;
            if(PartiTags[b].first == -1){//if it is interface vertex, check whether the distance label may be affected by interface
                int vbW=Tree[child].vert[k].second.first;
                for(int i=rootHeight;i<line.size();i++){
//                    if(Tree[rank[line[i]]].disInf.find(b) == Tree[rank[line[i]]].disInf.end()){
//                        cout<<"DisRePost. Not found this boundary vertex! "<<i<<": "<<Tree[child].uniqueVertex<<"("<< NodeOrder[Tree[child].uniqueVertex]<<","<<PartiTags[Tree[child].uniqueVertex].first<<") "<<line[i]<<"("<< NodeOrder[line[i]]<<","<<PartiTags[line[i]].first<<") "<<b<<"("<<NodeOrder[b]<<","<<PartiTags[b].first<<")"<<endl; exit(1);
//                    }
                    assert(BoundVertexMap[pid].find(b)!=BoundVertexMap[pid].end());
                    inf_i=BoundVertexMap[pid][b];
                    if(Tree[child].disPost[i]>vbW+Tree[rank[line[i]]].disInf[inf_i]){
                        Tree[child].disPost[i]=vbW+Tree[rank[line[i]]].disInf[inf_i];
                        Tree[child].FNPost[i]=false;
                        ProIDdisCha=true;
                    }
                }
            }
            else{//if not interface vertex
                int bH=Tree[rank[b]].height-1,vbW=Tree[child].vert[k].second.first;
                if(Tree[child].FNPost[bH]){//if distance label from child to b is directly sourced from their shortcut
                    if(Tree[child].DisRe.find(b)!=Tree[child].DisRe.end()){//if found b, all ancestor check
                        for(int i=rootHeight;i<bH;i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].disPost[i]>vbW+Tree[rank[b]].disPost[i]){//update ancestor distance label
                                Tree[child].disPost[i]=vbW+Tree[rank[b]].disPost[i];
                                Tree[child].FNPost[i]=false;
                                ProIDdisCha=true;
                            }
                        }
                        for(int i=bH+1;i<line.size();i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].disPost[i]>vbW+Tree[rank[line[i]]].disPost[bH]){
                                Tree[child].disPost[i]=vbW+Tree[rank[line[i]]].disPost[bH];
                                Tree[child].FNPost[i]=false;
                                ProIDdisCha=true;
                            }
                        }

                    }else{//partial ancestor check if we cannot find b

                        if(vertexIDChL.find(b)!=vertexIDChL.end()){
                            for(int i=rootHeight;i<bH;i++){
                                checkedDis.insert(make_pair(child,i));
                                if(Tree[child].disPost[i]>vbW+Tree[rank[b]].disPost[i]){
                                    Tree[child].disPost[i]=vbW+Tree[rank[b]].disPost[i];
                                    Tree[child].FNPost[i]=false;
                                    ProIDdisCha=true;
                                }
                            }
                        }
                        for(int i=bH+1;i<line.size();i++){
                            checkedDis.insert(make_pair(child,i));
                            if(Tree[child].disPost[i]>vbW+Tree[rank[line[i]]].disPost[bH]){
                                Tree[child].disPost[i]=vbW+Tree[rank[line[i]]].disPost[bH];
                                Tree[child].FNPost[i]=false;
                                ProIDdisCha=true;
                            }
                        }

                    }
                }
            }


        }
    }
    else{// if there is no label for checking
        for(int k=0;k<Tree[child].vert.size();k++){
            int b=Tree[child].vert[k].first;
            if(PartiTags[b].first == -1){//if it is interface vertex, check whether the distance label may be affected by interface

                int vbW=Tree[child].vert[k].second.first;
                for(int i=rootHeight;i<line.size();i++){
//                    if(Tree[rank[line[i]]].disInf.find(b) == Tree[rank[line[i]]].disInf.end()){
//                        cout<<"Not found this boundary vertex! "<<i<<": "<<Tree[child].uniqueVertex<<"("<< NodeOrder[Tree[child].uniqueVertex]<<","<<PartiTags[Tree[child].uniqueVertex].first<<") "<<line[i]<<"("<< NodeOrder[line[i]]<<","<<PartiTags[line[i]].first<<") "<<b<<"("<<NodeOrder[b]<<","<<PartiTags[b].first<<")"<<endl; exit(1);
//                    }
                    assert(BoundVertexMap[pid].find(b)!=BoundVertexMap[pid].end());
                    inf_i=BoundVertexMap[pid][b];
                    if(Tree[child].disPost[i]>vbW+Tree[rank[line[i]]].disInf[inf_i]){
                        Tree[child].disPost[i]=vbW+Tree[rank[line[i]]].disInf[inf_i];
                        Tree[child].FNPost[i]=false;
                        ProIDdisCha=true;
                    }
                }


            }else {//if not interface vertex
                int bH = Tree[rank[b]].height - 1, vbW = Tree[child].vert[k].second.first;
                if (Tree[child].FNPost[bH]) {//Property 5
                    if (vertexIDChL.find(b) != vertexIDChL.end()) {//if the distance label of b is changed
                        for (int i = rootHeight; i < bH; i++) {//check ancestor from 0 to bH
                            checkedDis.insert(make_pair(child, i));
                            if (Tree[child].disPost[i] > vbW + Tree[rank[b]].disPost[i]) {
                                Tree[child].disPost[i] = vbW + Tree[rank[b]].disPost[i];
                                Tree[child].FNPost[i] = false;
                                ProIDdisCha = true;
                            }
                        }
                    }
                    for (int i = bH + 1; i < line.size(); i++) {
                        checkedDis.insert(make_pair(child, i));
                        if (Tree[child].disPost[i] > vbW + Tree[rank[line[i]]].disPost[bH]) {
                            Tree[child].disPost[i] = vbW + Tree[rank[line[i]]].disPost[bH];
                            Tree[child].FNPost[i] = false;
                            ProIDdisCha = true;
                        }
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
        if(PartiTags[Tree[Tree[child].ch[i]].uniqueVertex].first != pid){
            cout<<"Wrong for this child! "<<pid<<": "<<ID<<" "<<Tree[Tree[child].ch[i]].uniqueVertex<<" "<<PartiTags[Tree[Tree[child].ch[i]].uniqueVertex].first<<endl;
        }
        AncestorEntryDecreaseUpdate(pid, Tree[child].ch[i], line, interfaces,vertexIDChL,checkedDis,Tree, rank, rootHeight);
    }
    line.pop_back();

}

//function for the interface entry propagation
void Graph::InterfacePropagate(int child, vector<int>& interfaces, vector<Node> &Tree, bool ifIncrease){
//    if(child == rank[210695]){
//        cout<<210695<<" "<<Tree[Tree[child].treeroot].uniqueVertex<<endl;
//    }
    // Solution 1: check all interface entries
    int inf_i;
    int pid=PartiTags[Tree[child].uniqueVertex].first;
    for(int j=0;j<interfaces.size();j++){//for each interface vertex
        int ID2=interfaces[j];
//        if(child == rank[197979] && ID2 == 144762)
//            cout<<ID2<<endl;
        assert(BoundVertexMap[pid].find(ID2)!=BoundVertexMap[pid].end());
        inf_i=BoundVertexMap[pid][ID2];
        if(ifIncrease){
            Tree[child].disInf[inf_i] = INF;
        }
        for(int k=0;k<Tree[child].vert.size();k++) {
            int b = Tree[child].vert[k].first;
            int vbW = Tree[child].vert[k].second.first;
            if(b == interfaces[j]){
                if(Tree[child].disInf[inf_i] > vbW){
                    Tree[child].disInf[inf_i]=vbW;//
                }
                continue;
            }
            int z = INF;
            if(PartiTags[b].first != -1 ){//if x is not boundary vertex, i.e., if it is ancestor
                z = Tree[rank[b]].disInf[inf_i];
                if(Tree[child].disInf[inf_i]>z+vbW){
                    Tree[child].disInf[inf_i]=z+vbW;
                }
            }else{//if it is interface vertex
//                    if(NeighborsOverlay[b].find(ID2) != NeighborsOverlay[b].end()){//if found
//                        z = NeighborsOverlay[b][ID2];
//                    }
                    if(b<=ID2){
                        if(BoundShortcuts[b].find(ID2) != BoundShortcuts[b].end()){//if found
                            z = BoundShortcuts[b][ID2];
                        }
                        else{
                            cout<<"Not found overlay shortcut "<<b<<" "<<ID2<<endl; exit(1);
                        }
                    }else{
                        if(BoundShortcuts[ID2].find(b) != BoundShortcuts[ID2].end()){//if found
                            z = BoundShortcuts[ID2][b];
                        }
                        else{
                            cout<<"Not found overlay shortcut "<<ID2<<" "<<b<<endl; exit(1);
                        }
                    }

                    assert(z>0);
                    if(Tree[child].disInf[inf_i]>z+vbW){// never enter this branch
                        Tree[child].disInf[inf_i]=z+vbW;
                    }


            }
        }
    }


    for(int i=0;i<Tree[child].ch.size();i++){
        InterfacePropagate(Tree[child].ch[i], interfaces,Tree,ifIncrease);
    }
}
//function for the interface entry propagation: parallel version
void Graph::InterfacePropagateParallel(pair<int,int> pRange, vector<int>& pids, bool ifIncrease){
    for(int i=pRange.first;i<pRange.second;++i){
        int Pid = pids[i];
        int rootID = partiRoots[Pid];
        vector<int> interfaceP;
//        for(auto it=Tree[rank[rootID]].vert.begin();it!=Tree[rank[rootID]].vert.end();++it){
//            interfaceP.emplace_back(it->first);
//        }
        interfaceP=BoundVertex[Pid];
        InterfacePropagate( rank[rootID], interfaceP,Tree,ifIncrease);
    }
}

void Graph::DecreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, set<int>& vertexIDChL){
    int ProBeginVertexID;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChL,Tree,rank);
    }
}

//batch update for partition graph
void Graph::DecreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU){
    map<int,int> checkedDis;

    for(int i=0;i<Tree.size();i++){
        Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
    }

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; //OC.clear();//vertexID in decreasing node order

    int a,b,newW,lid,hid,lidM,hidM;
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first; b=wBatch[k].first.second ;newW=wBatch[k].second.second;
        if(NodeOrder[a]<NodeOrder[b]){
            lid=a;hid=b;
        }else{
            lid=b;hid=a;
        }
        lidM = IDMap[lid]; hidM = IDMap[hid];

        if(Neighbors[a].find(b)!=Neighbors[a].end()){
            Neighbors[a][b]=newW;
        }else{
            cout<<"Not found edge! "<<endl; exit(1);
        }

        if(Neighbors[b].find(a)!=Neighbors[b].end()){
            Neighbors[b][a]=newW;
        }else{
            cout<<"Not found edge! "<<endl; exit(1);
        }

        for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
            if(Tree[rank[lidM]].vert[i].first==hid){
                if(Tree[rank[lidM]].vert[i].second.first>newW){
                    Tree[rank[lidM]].vert[i].second.first=newW;
                    Tree[rank[lidM]].vert[i].second.second=1;
                    SCre[lid].insert(hid);
                    OC.insert(OrderCompMin(lid));
                }else if(Tree[rank[lidM]].vert[i].second.first==newW){
                    Tree[rank[lidM]].vert[i].second.second+=1;
                }
                break;
            }
        }

    }

//    vector<int> ProBeginVertexSet; //ProBeginVertexSet.clear();
//    set<int> vertexIDChL; //vertexIDChL.clear();//record the vertex whose distance labeling has changed
    ProBeginVertexSetParti[pid].clear(); vertexIDChLParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    int ProBeginVertexID;
    int ProID, ProIDM;
    //processing the stars
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        bool ProIDdisCha=false;//to see if the distance labeling of proID change or not
        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw; int CidM=IDMap[Cid];
            int cidH=Tree[rank[CidM]].height-1;

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

            if(Tree[rank[ProIDM]].dis[cidH]>Cw){
                Tree[rank[ProIDM]].dis[cidH]=Cw;
                Tree[rank[ProIDM]].FN[cidH]=true;
                ProIDdisCha=true;
                Tree[rank[ProIDM]].DisRe.insert(Cid);
            }else if(Tree[rank[ProIDM]].dis[cidH]==Cw){
                Tree[rank[ProIDM]].FN[cidH]=true;
            }

            int hid2,hidHeight2,lid2,lidHeight2,wsum,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;hidHeight2=Tree[rank[IDMap[hid2]]].height-1;
                if(Hnei.find(hid2)!=Hnei.end()){
                    wsum=Cw+Hnei[hid2];
                    if(wsum<Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.first=wsum;
                        Tree[rank[CidM]].vert[j].second.second=1;
                        SCre[Cid].insert(hid2);
                        OC.insert(OrderCompMin(Cid));
                    }else if(wsum==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second+=1;
                    }

                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                lidHeight2=Tree[rank[lid2M]].height-1;
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        wsum=Cw+Lnei[j].second;
                        if(Tree[rank[lid2M]].vert[k].second.first>wsum){
                            Tree[rank[lid2M]].vert[k].second.first=wsum;
                            Tree[rank[lid2M]].vert[k].second.second=1;
                            SCre[lid2].insert(Cid);
                            OC.insert(OrderCompMin(lid2));
                        }else if(Tree[rank[lid2M]].vert[k].second.first==wsum){
                            Tree[rank[lid2M]].vert[k].second.second+=1;
                        }

                        break;
                    }
                }
            }
        }

        if(ProIDdisCha){//if the distance labeling is dectected changed
            vertexIDChLParti[pid].insert(ProID);
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
    }

    //cout<<"Finish bottom-up refresh"<<endl;
    if(ifLabelU){
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);
            EachNodeProBDis5Parti(rank[IDMap[ProBeginVertexID]], linee, vertexIDChLParti[pid],Tree,rank);
        }
    }

    //return checkedDis.size();
}

void Graph::IncreaseOverlay(int a,int b, int oldW, int newW, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int ChangeNum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    if(Neighbors[a].find(b)!=Neighbors[a].end()){
        Neighbors[a][b]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }
    if(Neighbors[b].find(a)!=Neighbors[b].end()){
        Neighbors[b][a]=newW;
    }else{
        cout<<"Wrong for Neighbors!"<<endl; exit(1);
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }

    int IniH=Tree[rank[lid]].height;//the height where weight change begins
    int ProH=Tree[rank[lid]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> vec; //vec.clear();
    SCre.assign(IniH+1,vec);
    int MinH;

    vector<int> line; //line.clear();
    line.reserve(heightMax);
    int pachid=ProID;
    while(Tree[rank[pachid]].height>1){
        line.insert(line.begin(),pachid);
        pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
    }
    line.insert(line.begin(),pachid);

    bool tri=false;
    for(int i=0;i<Tree[rank[lid]].vert.size();i++){
        if(Tree[rank[lid]].vert[i].first==hid){
            if(Tree[rank[lid]].vert[i].second.first==oldW){
                Tree[rank[lid]].vert[i].second.second-=1;
                if(Tree[rank[lid]].vert[i].second.second<1){
                    OCdis[make_pair(lid,hid)]=oldW;
                    SCre[ProH].insert(hid);
                    MinH=IniH;
                    tri=true;//cout<<"Trigger the Shortcut Change or not? "<<tri<<endl;
                }
            }
            break;
        }
    }

    bool influence; int ProBeginID;
    if(tri){
        //shortcut update
        while(ProH>=MinH){
            influence=false;
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[Cid]].height-1;

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
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
                                SCre[Tree[rank[Cid]].height].insert(hid);
                                if(Tree[rank[Cid]].height<MinH) MinH=Tree[rank[Cid]].height;
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
                                    SCre[Tree[rank[lid]].height].insert(Cid);
                                    if(Tree[rank[lid]].height<MinH) MinH=Tree[rank[lid]].height;
                                    OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                                }
                            }
                            break;
                        }
                    }
                }

                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProID]].FN[cidH]){
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProID]].FN[cidH]=false;
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

                for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                    if(it2->first==Cid){
                        Cw=it2->second;//the weight value in the original graph
                        countwt=1;
                        break;
                    }
                }

                int ssw,wtt,wid;
                vector<pair<int,int>> Wnodes; //Wnodes.clear();
                /*if(ProID<Cid)
                    Wnodes=SCconNodes[make_pair(ProID,Cid)]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodes[make_pair(Cid,ProID)];*/

                if(ProID<Cid)
                    Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodesMT[Cid][ProID];
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
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProID]].pa].uniqueVertex;
        }
    }

    vector<int> line1; //line1.clear();
    line1.reserve(heightMax);
    pachid=Tree[Tree[rank[ProBeginID]].pa].uniqueVertex;
    while(Tree[rank[pachid]].height>1){
        line1.insert(line1.begin(),pachid);
        pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
    }
    line1.insert(line1.begin(),pachid);

    eachNodeProcessIncrease1(rank[ProBeginID],line1,ChangeNum,Tree,rank,VidtoTNid);

    //return ChangeNum;
}

void Graph::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
    for(int i=0;i<Tree[children].dis.size();i++){
        if(Tree[children].cnt[i]==0){
            vUpdated[childID] = true;
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int PID;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){
                PID=VidtoTNid[childID][k];
                if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){
                    Tree[PID].cnt[i]-=1;
                }
            }

            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                PID=VidtoTNid[line[i]][k];
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){///
                    if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){
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
            if(DDvb==dis) Tree[children].FN[i]=true;
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
    }
    line.pop_back();
}

void Graph::eachNodeProcessIncrease1PostMHL(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
    for(int i=0;i<Tree[children].dis.size()-1;i++){//may not -1
        if(Tree[children].cnt[i]==0){
//            vUpdated[childID] = true;
            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int PID;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){
                PID=VidtoTNid[childID][k];
                if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){
                    Tree[PID].cnt[i]-=1;
                }
            }

            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                PID=VidtoTNid[line[i]][k];
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){///
                    if(Tree[PID].FN[i] && Tree[PID].dis[childH]==disBF+Tree[PID].dis[i]){
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
            if(DDvb==dis) Tree[children].FN[i]=true;
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1PostMHL(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
    }
    line.pop_back();
}

void Graph::eachNodeProcessIncrease1PostMHLOverlay(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
    int childID=Tree[children].uniqueVertex;
    if(PartiTags[childID].first==-1){
        int childH=Tree[children].height-1;
        for(int i=0;i<Tree[children].dis.size()-1;i++){
            if(Tree[children].cnt[i]==0){
                vUpdated[childID] = true;
                changelabel+=1;
                //firstly, check which dis can be infected
                int disBF=Tree[children].dis[i];
                int PID;
                //chidlID
                for(int k=0;k<VidtoTNid[childID].size();k++){
                    PID=VidtoTNid[childID][k];
                    if(PartiTags[Tree[PID].uniqueVertex].first==-1 || algoUpdate==PH2H_Cross){
                        if(Tree[PID].FN[childH] && Tree[PID].dis[i]==disBF+Tree[PID].dis[childH]){
                            Tree[PID].cnt[i]-=1;
                        }
                    }

                }

                //line[i]
                for(int k=0;k<VidtoTNid[line[i]].size();k++){
                    PID=VidtoTNid[line[i]][k];
                    if(PartiTags[Tree[PID].uniqueVertex].first==-1 || algoUpdate==PH2H_Cross) {
                        if (Tree[PID].height > Tree[children].height &&
                            Tree[PID].vAncestor[Tree[children].height - 1] == childID) {///
                            if (Tree[PID].FN[i] && Tree[PID].dis[childH] == disBF + Tree[PID].dis[i]) {
                                Tree[PID].cnt[childH] -= 1;
                            }
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
            eachNodeProcessIncrease1PostMHLOverlay(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
        }
        line.pop_back();
    }
    else{//if partition vertex
        int pid=PartiTags[childID].first;
        vector<int> ProBeginVertexSetNew;
        ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetPartiExtend[pid].size()+1);
        ProBeginVertexSetNew.push_back(childID);
        int rnew=rank[childID],r;
        for(int i=0;i<ProBeginVertexSetPartiExtend[pid].size();i++){
            r=rank[ProBeginVertexSetPartiExtend[pid][i]];
            if(LCAQueryOverlay(rnew,r)!=rnew){
                ProBeginVertexSetNew.push_back(ProBeginVertexSetPartiExtend[pid][i]);
            }
        }
        ProBeginVertexSetPartiExtend[pid]=ProBeginVertexSetNew;
    }
}

void Graph::IncreaseParti(int a,int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid){
    int ChangeNum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    for(int i=0;i<Neighbors[a].size();i++){
        if(Neighbors[a][i].first==b){
            Neighbors[a][i].second=newW;
            break;
        }
    }
    for(int i=0;i<Neighbors[b].size();i++){
        if(Neighbors[b][i].first==a){
            Neighbors[b][i].second=newW;
            break;
        }
    }

    int lid,hid;
    if(NodeOrder[a]<NodeOrder[b]){
        lid=a;hid=b;
    }else{
        lid=b;hid=a;
    }
    int lidM=IDMap[lid];
    int IniH=Tree[rank[lidM]].height;//the height where weight change begins
    int ProH=Tree[rank[lidM]].height; int ProID=lid;
    vector<set<int>> SCre;//record the shortcut change in each height
    set<int> vec; //vec.clear();
    SCre.assign(IniH+1,vec);
    int MinH;

    vector<int> line; //line.clear();
    line.reserve(heightMax);
    int pachid=ProID;
    while(Tree[rank[IDMap[pachid]]].height>1){
        line.insert(line.begin(),pachid);
        pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
    }
    line.insert(line.begin(),pachid);

    bool tri=false;
    for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
        if(Tree[rank[lidM]].vert[i].first==hid){
            if(Tree[rank[lidM]].vert[i].second.first==oldW){
                Tree[rank[lidM]].vert[i].second.second-=1;
                if(Tree[rank[lidM]].vert[i].second.second<1){
                    OCdis[make_pair(lid,hid)]=oldW;
                    SCre[ProH].insert(hid);
                    MinH=IniH;
                    tri=true;//cout<<"Trigger the Shortcut Change or not? "<<tri<<endl;
                }
            }
            break;
        }
    }

    bool influence; int ProBeginID;
    if(tri){
        //shortcut update
        while(ProH>=MinH){
            influence=false;
            int ProIDM=IDMap[ProID];
            vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
            for(auto it=SCre[ProH].begin();it!=SCre[ProH].end();it++){
                int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
                int cidH=Tree[rank[IDMap[Cid]]].height-1;
                int CidM=IDMap[Cid];

                map<int,int> Hnei; //Hnei.clear();
                vector<pair<int,int>> Lnei; //Lnei.clear();
                for(int j=0;j<Vert.size();j++){
                    if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                        Hnei[Vert[j].first]=Vert[j].second.first;
                    }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                        Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                    }
                }
                //check the affected shortcuts
                int hid,lid;
                for(int j=0;j<Tree[rank[CidM]].vert.size();j++){// for higher-order vertices
                    hid=Tree[rank[CidM]].vert[j].first;
                    if(Hnei.find(hid)!=Hnei.end()){
                        if(Cw+Hnei[hid]==Tree[rank[CidM]].vert[j].second.first){
                            Tree[rank[CidM]].vert[j].second.second-=1;
                            if(Tree[rank[CidM]].vert[j].second.second<1){
                                SCre[Tree[rank[CidM]].height].insert(hid);
                                if(Tree[rank[CidM]].height<MinH) MinH=Tree[rank[CidM]].height;
                                OCdis[make_pair(Cid,hid)]=Cw+Hnei[hid];
                            }
                        }
                    }
                }
                for(int j=0;j<Lnei.size();j++){//for lower-order vertices
                    lid=Lnei[j].first; lidM=IDMap[lid];
                    for(int k=0;k<Tree[rank[lidM]].vert.size();k++){
                        if(Tree[rank[lidM]].vert[k].first==Cid){
                            if(Tree[rank[lidM]].vert[k].second.first==Cw+Lnei[j].second){
                                Tree[rank[lidM]].vert[k].second.second-=1;
                                if(Tree[rank[lidM]].vert[k].second.second<1){
                                    SCre[Tree[rank[lidM]].height].insert(Cid);
                                    if(Tree[rank[lidM]].height<MinH) MinH=Tree[rank[lidM]].height;
                                    OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                                }
                            }
                            break;
                        }
                    }
                }

                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProIDM]].FN[cidH]){//if the label is obtained from shortcut
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[Cid]]].dis[i]){
                            Tree[rank[ProIDM]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProIDM]].FN[cidH]=false;
                    Tree[rank[ProIDM]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                        if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                            Tree[rank[ProIDM]].cnt[i]-=1;
                        }
                    }
                }

                //get the new value of shortcut
                //	cout<<Cw<<" increase to ";
                Cw=INF; int countwt=0;

                for(int i=0;i<Neighbors[ProID].size();i++){
                    if(Neighbors[ProID][i].first==Cid){
                        Cw=Neighbors[ProID][i].second;//the weight value in the original graph
                        countwt=1;
                        break;
                    }
                }

                int ssw,wtt,wid;
                vector<pair<int,int>> Wnodes; //Wnodes.clear();
                /*if(ProID<Cid)
                    Wnodes=SCconNodes[make_pair(ProID,Cid)]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodes[make_pair(Cid,ProID)];*/

                if(ProID<Cid)
                    Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
                else
                    Wnodes=SCconNodesMT[Cid][ProID];
                for(int i=0;i<Wnodes.size();i++){//for each supportive vertex of this shortcut
                    wid=Wnodes[i].first;
                    int widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<Cw){
                        Cw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==Cw){
                        countwt+=1;
                    }
                }

                //cout<<Cw<<endl;
                //refresh the shortcut to the new value
                for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                    if(Tree[rank[ProIDM]].vert[i].first==Cid){
                        Tree[rank[ProIDM]].vert[i].second.first=Cw;
                        Tree[rank[ProIDM]].vert[i].second.second=countwt;
                        break;
                    }
                }
            }

            if(influence){
                ProBeginID=ProID;
            }

            ProH-=1;
            ProID=Tree[Tree[rank[ProIDM]].pa].uniqueVertex;
        }
    }

    vector<int> line1; //line1.clear();
    line1.reserve(heightMax);
    pachid=Tree[Tree[rank[IDMap[ProBeginID]]].pa].uniqueVertex;
    while(Tree[rank[IDMap[pachid]]].height>1){
        line1.insert(line1.begin(),pachid);
        pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
    }
    line1.insert(line1.begin(),pachid);

    eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginID]],line1,ChangeNum,Tree,rank,VidtoTNid);

    //return ChangeNum;
}

void Graph::eachNodeProcessIncrease1Parti(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid){
    int childID=Tree[children].uniqueVertex;
    int childH=Tree[children].height-1;
//    int childM=IDMap[children];


    for(int i=0;i<Tree[children].dis.size();i++){
//        if(childID==84704){
//            cout<<childID<<": "<<i<<", "<<Tree[children].vAncestor[i]<< " "<<Tree[children].cnt[i]<< endl;
//        }
        if(Tree[children].cnt[i]==0){//the distance from child to line[i] should be updated
            if(i<BoundVertex[PartiTag[childID].first].size()){
                vUpdated[childID] = true;
            }

            changelabel+=1;
            //firstly, check which dis can be infected
            int disBF=Tree[children].dis[i];
            int tid;
            //chidlID
            for(int k=0;k<VidtoTNid[childID].size();k++){//check the tree nodes containing child
                tid=VidtoTNid[childID][k];
//                cout<<k<<" "<<tid<<" "<<IDMap[tid]<<" "<<childH<<" "<<Tree[tid].dis.size()<<" "<<Tree[tid].FN.size()<<endl;
                if(Tree[tid].FN[childH] && Tree[tid].dis[i]==disBF+Tree[tid].dis[childH]){//if the distance from tid to line[i] sources from child, valley path
                    Tree[tid].cnt[i]-=1;
                }
            }
//            cout<<"Flag 2"<<endl;
            //line[i]
            for(int k=0;k<VidtoTNid[line[i]].size();k++){
                tid=VidtoTNid[line[i]][k];
                if(Tree[tid].height>Tree[children].height && Tree[tid].vAncestor[childH] == childID){//children is the ancestor of tid
                    if(Tree[tid].FN[i] && Tree[tid].dis[childH]==disBF+Tree[tid].dis[i]){//if the distance from tid to child sources from line[i], peak path, out of scope
                        Tree[tid].cnt[childH]-=1;
                    }
                }
            }

            //secondly, calculate the actual distance
            int dis=INF; int count=0;
            int Dvb; int b,bH; int DDvb=INF;
            for(int j=0;j<Tree[children].vert.size();j++){
                Dvb=Tree[children].vert[j].second.first;
                b=Tree[children].vert[j].first;
                bH=Tree[rank[IDMap[b]]].height-1;
                if(bH<i){//if b is the ancestor of line[i]
                    if(Dvb+Tree[rank[IDMap[line[i]]]].dis[bH]<dis){
                        dis=Dvb+Tree[rank[IDMap[line[i]]]].dis[bH];
                        count=1;
                    }else if(Dvb+Tree[rank[IDMap[line[i]]]].dis[bH]==dis){
                        count+=1;
                    }
                }else if(bH==i){
                    DDvb=Dvb;//shortcut
                    if(Dvb<dis){
                        dis=Dvb;
                        count=1;
                    }else if(Dvb==dis){
                        count+=1;
                    }
                }else{//if line[i] is the ancestor of b
                    if(Dvb+Tree[rank[IDMap[b]]].dis[i]<dis){
                        dis=Dvb+Tree[rank[IDMap[b]]].dis[i];
                        count=1;
                    }else if(Dvb+Tree[rank[IDMap[b]]].dis[i]==dis){
                        count+=1;
                    }
                }
            }
            if(DDvb==dis) Tree[children].FN[i]=true;
            Tree[children].dis[i]=dis;
            Tree[children].cnt[i]=count;
        }
    }

    line.push_back(childID);
    for(int i=0;i<Tree[children].ch.size();i++){
        eachNodeProcessIncrease1Parti(Tree[children].ch[i],line,changelabel,Tree,rank,VidtoTNid);
    }
    line.pop_back();
}

void Graph::IncreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            if(Neighbor[a].find(b)!=Neighbor[a].end()){
                if(Neighbor[a][b]!=oldW){//only works for no-boundary
                    cout<<"Inconsistent! "<<Neighbor[a][b]<<" "<<oldW<<endl; exit(1);
                }
                Neighbor[a][b]=newW;
            }else{
                cout<<"Wrong for Neighbors!"<<endl; exit(1);
            }
            if(Neighbor[b].find(a)!=Neighbor[b].end()){
                if(Neighbor[b][a]!=oldW){
                    cout<<"Inconsistent! "<<Neighbor[b][a]<<" "<<oldW<<endl; exit(1);
                }
                Neighbor[b][a]=newW;
            }else{
                cout<<"Wrong for Neighbors!"<<endl; exit(1);
            }

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
                        if(Tree[rank[lid]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    ProBeginVertexSetOverlay.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        //each ProID corresponds to a line
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
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
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
                            OC.insert(OrderCompMin(Cid));
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
                                OC.insert(OrderCompMin(lid));
                                OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if(Tree[rank[ProID]].FN[cidH]){
                influence=true;
                //higher than Cid
                for(int i=0;i<cidH;i++){
                    if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                        Tree[rank[ProID]].cnt[i]-=1;
                    }
                }

                //equal to Cid
                Tree[rank[ProID]].FN[cidH]=false;
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

            for(auto it2=Neighbor[ProID].begin();it2!=Neighbor[ProID].end();++it2){
                if(it2->first==Cid){
                    Cw=it2->second;//the weight value in the original graph
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
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                r=rank[ProBeginVertexSetOverlay[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}

void Graph::IncreaseH2HBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW>=newW) {
            cout<<"Wrong update weight. "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl; exit(1);
        }
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
                    if(Tree[rank[lid]].vert[i].second.second<1){
                        OCdis[make_pair(lid,hid)]=oldW;
                        SCre[lid].insert(hid);
                        OC.insert(OrderCompMin(lid));
                    }
                }
                break;
            }
        }

    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
//    ProBeginVertexSetOverlay.clear();//has been cleared before
    vector<int> ProBeginVertexSetNew;
    bool influence,influenceExtend;
    int ProID; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        //each ProID corresponds to a line
        line=Tree[rank[ProID]].vAncestor;
//        line.clear(); line.reserve(heightMax);
//        int pachid=ProID;
//        while(Tree[rank[pachid]].height>1){
//            line.insert(line.begin(),pachid);
//            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
//        }
//        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
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
                            OC.insert(OrderCompMin(Cid));
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
                                OC.insert(OrderCompMin(lid));
                                OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }

            if(algoUpdate>PCH_No){
                if(PartiTags[ProID].first == -1){//if overlay vertex
                    //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                    if(Tree[rank[ProID]].FN[cidH]){
                        influence=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                                Tree[rank[ProID]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProID]].FN[cidH]=false;
                        Tree[rank[ProID]].cnt[cidH]-=1;

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                            if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                                Tree[rank[ProID]].cnt[i]-=1;
                            }
                        }
                    }
                }else{//if ProID is partition vertex
                    if(algoUpdate>=PH2H_Post){
                        if(PartiTags[Cid].first!=-1){//if Cid is partition vertex
                            if(Tree[rank[ProID]].FNPost[cidH]){
                                influence=true;
                                //higher than Cid
                                for(int i=0;i<cidH;i++){
                                    if(PartiTags[line[i]].first!=-1){//if partition vertex
                                        if(Tree[rank[ProID]].disPost[i]==Cw+Tree[rank[Cid]].disPost[i]){
                                            Tree[rank[ProID]].cntPost[i]-=1;
                                        }
                                    }
                                }

                                //equal to Cid
                                Tree[rank[ProID]].FNPost[cidH]=false;
                                Tree[rank[ProID]].cntPost[cidH]-=1;

                                //lower than Cid
                                for(int i=cidH+1;i<Tree[rank[ProID]].disPost.size();i++){
                                    if(Tree[rank[ProID]].disPost[i]==Cw+Tree[rank[line[i]]].disPost[cidH]){
                                        Tree[rank[ProID]].cntPost[i]-=1;
                                    }
                                }
                            }
                        }

                    }

                    if(algoUpdate==PH2H_Cross){
                        //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                        if(Tree[rank[ProID]].FN[cidH]){
//                            influence=true;
                            influenceExtend=true;
                            //higher than Cid
                            for(int i=0;i<cidH;i++){
                                if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                                    Tree[rank[ProID]].cnt[i]-=1;
                                }
                            }

                            //equal to Cid
                            Tree[rank[ProID]].FN[cidH]=false;
                            Tree[rank[ProID]].cnt[cidH]-=1;

                            //lower than Cid
                            for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                                if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                                    Tree[rank[ProID]].cnt[i]-=1;
                                }
                            }
                        }
                    }

                }

            }


            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            Cw=INF; int countwt=0;

            for(auto it2=Neighbor[ProID].begin();it2!=Neighbor[ProID].end();++it2){
                if(it2->first==Cid){
                    Cw=it2->second;//the weight value in the original graph
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
            if(PartiTags[ProID].first==-1) {
                ProBeginVertexSetNew.clear();
                ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size() + 1);
                ProBeginVertexSetNew.push_back(ProID);
                int rnew = rank[ProID], r;
                for (int i = 0; i < ProBeginVertexSetOverlay.size(); i++) {
                    r = rank[ProBeginVertexSetOverlay[i]];
                    if (LCAQueryOverlay(rnew, r) != rnew) {
                        ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                    }
                }
                ProBeginVertexSetOverlay = ProBeginVertexSetNew;
            }
            if(PartiTags[ProID].first!=-1){
                int pid=PartiTags[ProID].first;
                affectedParti.insert(pid);
                ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
                ProBeginVertexSetNew.push_back(ProID);
                int rnew=rank[ProID],r;
                for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                    r=rank[ProBeginVertexSetParti[pid][i]];
                    if(LCAQueryOverlay(rnew,r)!=rnew){
                        ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                    }
                }
                ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
            }
        }

        if(influenceExtend){
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size() + 1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew = rank[ProID], r;
            for (int i = 0; i < ProBeginVertexSetOverlay.size(); i++) {
                r = rank[ProBeginVertexSetOverlay[i]];
                if (LCAQueryOverlay(rnew, r) != rnew) {
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay = ProBeginVertexSetNew;
        }
    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

//            eachNodeProcessIncrease1PostMHLOverlay(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
            eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}

void Graph::IncreaseOverlayBatchPostMHL(vector<pair<pair<int, int>, pair<int, int>>> &wBatch,
                                     vector<Node> &Tree, vector<int> &rank,
                                     int heightMax, vector<map<int, vector<pair<int, int>>>> &SCconNodesMT,
                                     vector<vector<int>> &VidtoTNid, bool ifLabelU) {
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW>=newW) {
            cout<<"Wrong update weight. "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl; exit(1);
        }
//        if(Neighbor[a].find(b)!=Neighbor[a].end()){
//            Neighbor[a][b]=newW;
//        }else{
//            cout<<"Wrong for Neighbors!"<<endl; exit(1);
//        }
//        if(Neighbor[b].find(a)!=Neighbor[b].end()){
//            Neighbor[b][a]=newW;
//        }else{
//            cout<<"Wrong for Neighbors!"<<endl; exit(1);
//        }
//
//        for(int i=0;i<NeighborsOverlayV[a].size();i++){
//            if(NeighborsOverlayV[a][i].first==b){
////            cout<<Neighbor[a][i].second<<" "<<newW<<endl;
//                NeighborsOverlayV[a][i].second=newW;
//                break;
//            }
//        }
//        for(int i=0;i<NeighborsOverlayV[b].size();i++){
//            if(NeighborsOverlayV[b][i].first==a){
////            cout<<Neighbor[b][i].second<<" "<<newW<<endl;
//                NeighborsOverlayV[b][i].second=newW;
//                break;
//            }
//        }

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
                    if(Tree[rank[lid]].vert[i].second.second<1){
                        OCdis[make_pair(lid,hid)]=oldW;
                        SCre[lid].insert(hid);
                        OC.insert(OrderCompMin(lid));
                    }
                }
                break;
            }
        }

    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
//    ProBeginVertexSetOverlay.clear();//has been cleared before
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;

        //each ProID corresponds to a line
        line=Tree[rank[ProID]].vAncestor;
//        line.clear(); line.reserve(heightMax);
//        int pachid=ProID;
//        while(Tree[rank[pachid]].height>1){
//            line.insert(line.begin(),pachid);
//            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
//        }
//        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
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
                            OC.insert(OrderCompMin(Cid));
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
                                OC.insert(OrderCompMin(lid));
                                OCdis[make_pair(lid,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }

            if(algoUpdate>PCH_No){
                if(PartiTags[ProID].first != -1){
                    cout<<"Wrong for overlay shortcut update! "<<ProID<<" "<<PartiTags[ProID].first<<endl; exit(1);
                }
                //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
                if(Tree[rank[ProID]].FN[cidH]){
                    influence=true;
                    //higher than Cid
                    for(int i=0;i<cidH;i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }

                    //equal to Cid
                    Tree[rank[ProID]].FN[cidH]=false;
                    Tree[rank[ProID]].cnt[cidH]-=1;

                    //lower than Cid
                    for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                        if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                            Tree[rank[ProID]].cnt[i]-=1;
                        }
                    }
                }

            }


            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            Cw=INF; int countwt=0;

            for(auto it2=this->Neighbor[ProID].begin();it2!=this->Neighbor[ProID].end();++it2){
                if(it2->first==Cid){
                    Cw=it2->second;//the weight value in the original graph
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
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
                r=rank[ProBeginVertexSetOverlay[i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetOverlay.size();i++){
            ProBeginVertexID=ProBeginVertexSetOverlay[i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
            while(Tree[rank[pachidd]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1PostMHLOverlay(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}

void Graph::IncreaseOverlayBatchLabelPostMHL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid) {
    int ProBeginVertexID;
    int checknum=0;
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

        eachNodeProcessIncrease1PostMHLOverlay(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
    }
}

void Graph::IncreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid) {
    int ProBeginVertexID;
    int checknum=0;
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

        eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
    }
}

void Graph::IncreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            for(int i=0;i<Neighbors[a].size();i++){
                if(Neighbors[a][i].first==b){
                    Neighbors[a][i].second=newW;
                    break;
                }
            }
            for(int i=0;i<Neighbors[b].size();i++){
                if(Neighbors[b][i].first==a){
                    Neighbors[b][i].second=newW;
                    break;
                }
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    ProBeginVertexSetParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            Tree[rank[lid2M]].vert[k].second.second-=1;
                            if(Tree[rank[lid2M]].vert[k].second.second<1){
                                SCre[lid2].insert(Cid);
                                OC.insert(OrderCompMin(lid2));
                                OCdis[make_pair(lid2,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if(Tree[rank[ProIDM]].FN[cidH]){
                influence=true;
                //higher than Cid
                for(int i=0;i<cidH;i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }

                //equal to Cid
                Tree[rank[ProIDM]].FN[cidH]=false;
                Tree[rank[ProIDM]].cnt[cidH]-=1;

                //lower than Cid
                for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            Cw=INF; int countwt=0;

            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    Cw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
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

            if(PartiTag[ProID].second){//if boundary vertex
                updatedSC.emplace_back(make_pair(ProID,Cid),Cw);
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=Cw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }
        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}

void Graph::IncreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){//oldW may be incorrect!
            for(int i=0;i<Neighbors[a].size();i++){
                if(Neighbors[a][i].first==b){
                    if(oldW!=Neighbors[a][i].second){
                        cout<<"Edge weight inconsistent! "<<a<<" "<<b<<" "<<Neighbors[a][i].second<<" "<<oldW<<endl;
                        oldW=Neighbors[a][i].second;
                    }
                    Neighbors[a][i].second=newW;
                    break;
                }
            }
            for(int i=0;i<Neighbors[b].size();i++){
                if(Neighbors[b][i].first==a){
                    if(oldW!=Neighbors[b][i].second){
                        cout<<"Edge weight inconsistent! "<<b<<" "<<a<<" "<<Neighbors[b][i].second<<" "<<oldW<<endl;
                        oldW=Neighbors[b][i].second;
                    }
                    Neighbors[b][i].second=newW;
                    break;
                }
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    if(!ifLabelU){
        ProBeginVertexSetParti[pid].clear();
    }

    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            Tree[rank[lid2M]].vert[k].second.second-=1;
                            if(Tree[rank[lid2M]].vert[k].second.second<1){
                                SCre[lid2].insert(Cid);
                                OC.insert(OrderCompMin(lid2));
                                OCdis[make_pair(lid2,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }

//            if(ProID==84704){
//                cout<<"ProID "<<ProID<<endl;
//            }
            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            /*if(Tree[rank[ProIDM]].FN[cidH]){
                influence=true;
                //higher than Cid
                for(int i=0;i<cidH;i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }

                //equal to Cid
                Tree[rank[ProIDM]].FN[cidH]=false;
                Tree[rank[ProIDM]].cnt[cidH]-=1;

                //lower than Cid
                for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }
            }*/

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;

            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    newCw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
                        }
                    }

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }

            if(PartiTag[ProID].second){//if boundary vertex
//                cout<<"Boundary shortcut update: "<<ProID<<" "<<Cid<<" "<<Cw<<" "<<newCw<<endl;
                updatedSC.emplace_back(make_pair(ProID,Cid),newCw);
            }

            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=newCw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }
        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}

void Graph::IncreasePartiBatchPostMHLShortcut(int pid, vector<pair<pair<int, int>, pair<int, int>>> &wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int, int>, int>> &updatedSC, vector<int>& PropagateOverlay ) {
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW>=newW) {
            cout<<"Wrong update weight. "<<a<<" "<<b<<" "<<oldW<<" "<<newW<<endl; exit(1);
        }

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
                    if(Tree[rank[lid]].vert[i].second.second<1){
                        OCdis[make_pair(lid,hid)]=oldW;
                        SCre[lid].insert(hid);
                        OC.insert(OrderCompMin(lid));
                    }
                }
                break;
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    vector<int> ProBeginVertexSetNew;
    bool influence, influenceExtend;
    int ProID;
    vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x;
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProID]].vert;
        influence=false;
        influenceExtend=false;
        //each ProID corresponds to a line
        line=Tree[rank[ProID]].vAncestor;
//        line.clear(); line.reserve(heightMax);
//        int pachid=ProID;
//        while(Tree[rank[pachid]].height>1){
//            line.insert(line.begin(),pachid);
//            pachid=Tree[Tree[rank[pachid]].pa].uniqueVertex;
//        }
//        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[Cid]].height-1;

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }

            //check the affected shortcuts
            int hid2,lid2;
            if(PartiTags[Cid].first == -1){//if overlay vertex
//                if(Tree[rank[ProID]].height<MinHInf){
//                    MinHInf = Tree[rank[ProID]].height;
//                    ProBeginIDInf.clear();
//                    ProBeginIDInf.insert(ProID);
//                }else if(Tree[rank[ProID]].height==MinHInf){
//                    ProBeginIDInf.insert(ProID);
//                }
                int wsum,lidHeight;
                /// For Hnei, update the super edges between Cid and other interface vertex, i.e., the shortcuts in AdjaCore
                for(int i=0;i<Tree[rank[Cid]].vert.size();++i){
                    hid2=Tree[rank[Cid]].vert[i].first;
                    if(Hnei.find(hid2)!=Hnei.end()){
                        wsum=Cw+Hnei[hid2];
//                                if(hid == 48300 || hid == 27115)
//                                    cout<<"1 "<<Cid<<" "<<hid<<endl;
                        if(wsum==Tree[rank[Cid]].vert[i].second.first){
                            //since we do not know whether e(Cid,hid) is only supported by wsum, we will check it after all in-periphery shortcuts are processed
//                                SECoreForUpdate.insert(make_pair(Cid,hid));
                            sm->wait();
                            Tree[rank[Cid]].vert[i].second.second-=1;
//                                if(Tree[rank[Cid]].vert[i].second.second<1){
//                                    SCre[Cid].insert(hid2);
//                                    OC.insert(OrderCompMin(Cid));
//                                    OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
//                                }
                            sm->notify();
                            updatedSC.emplace_back(make_pair(Cid,hid2),wsum);//shortcut (Cid,hid2) may be affected

                        }
                    }
                }

                /// For Lnei, update the super edges of Lnei vertex
                for(int j=0;j<Lnei.size();j++){
                    lid2=Lnei[j].first;
                    if(PartiTags[lid2].first == -1){//if lid is interface vertex
                        for(int i=0;i<Tree[rank[lid2]].vert.size();++i) {
                            int vertid = Tree[rank[lid2]].vert[i].first;
                            if (vertid == Cid) {//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                                wsum=Cw+Lnei[j].second;
//                                    if(lid == 48300 || lid == 27115)
//                                        cout<<"2 "<<lid<<" "<<Cid<<endl;
                                if(Tree[rank[lid2]].vert[i].second.first==wsum){//update neighbor distance of Lnei
                                    sm->wait();
                                    Tree[rank[lid2]].vert[i].second.second-=1;
//                                    if(Tree[rank[lid2]].vert[i].second.second<1){
//                                        SCre[lid2].insert(Cid);
//                                        OC.insert(OrderCompMin(lid2));
//                                        OCdis[make_pair(lid2,Cid)]=Cw+Lnei[j].second;
//                                    }
                                    sm->notify();
                                    updatedSC.emplace_back(make_pair(lid2,Cid),wsum);//shortcut may be affected

                                }
                                break;
                            }
                        }
                    }
                    else{//if lid is periphery vertex
                        lidHeight=Tree[rank[lid2]].height-1;
                        for(int k=0;k<Tree[rank[lid2]].vert.size();k++){
                            int vertid=Tree[rank[lid2]].vert[k].first;
                            if(vertid==Cid){//Only deal with Cid for Lnei vertex; Tree[rank[lid]].vert[k].first
                                wsum=Cw+Lnei[j].second;
                                if(Tree[rank[lid2]].vert[k].second.first==wsum){///check 3: the interface entry of lid
                                    Tree[rank[lid2]].vert[k].second.second-=1;
                                    if(Tree[rank[lid2]].vert[k].second.second<1){
                                        SCre[lid2].insert(Cid);
                                        OC.insert(OrderCompMin(lid2));
                                        OCdis[make_pair(lid2,Cid)]=Cw+Lnei[j].second;

                                        /// interface
//                                        if(PartiTags[vertid].first == -1){//if the neighbor is interface vertex, never enter this branch
//                                            assert(Tree[rank[lid2]].disInf.find(vertid) != Tree[rank[lid2]].disInf.end());
//                                            if(Tree[rank[lid2]].disInf[vertid] == Cw+Lnei[j].second){
//                                                Tree[rank[lid2]].DisReInf.insert(vertid);//record the vertex id that the interface label should be updated
////                                        Tree[rank[lid]].disInf[vertid] = wsum;
//                                            }
//                                        }
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }
            else{//if partition vertex
                for(int j=0;j<Tree[rank[Cid]].vert.size();j++){
                    hid2=Tree[rank[Cid]].vert[j].first;
                    if(Hnei.find(hid2)!=Hnei.end()){
                        if(Cw+Hnei[hid2]==Tree[rank[Cid]].vert[j].second.first){
                            Tree[rank[Cid]].vert[j].second.second-=1;
                            if(Tree[rank[Cid]].vert[j].second.second<1){
                                SCre[Cid].insert(hid2);
                                OC.insert(OrderCompMin(Cid));
                                OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                            }
                        }
                    }
                }
                for(int j=0;j<Lnei.size();j++){
                    lid2=Lnei[j].first;
                    for(int k=0;k<Tree[rank[lid2]].vert.size();k++){
                        if(Tree[rank[lid2]].vert[k].first==Cid){
                            if(Tree[rank[lid2]].vert[k].second.first==Cw+Lnei[j].second){
                                Tree[rank[lid2]].vert[k].second.second-=1;
                                if(Tree[rank[lid2]].vert[k].second.second<1){
                                    SCre[lid2].insert(Cid);
                                    OC.insert(OrderCompMin(lid2));
                                    OCdis[make_pair(lid2,Cid)]=Cw+Lnei[j].second;
                                }
                            }
                            break;
                        }
                    }
                }
            }

            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if(algoUpdate>PCH_No){
                if(PartiTags[ProID].first == -1){
                    cout<<"Wrong for partition shortcut udpate! "<<pid<<" "<<ProID<<" "<<PartiTags[ProID].first<<endl; exit(1);
                }
                if(PartiTags[Cid].first != -1){//if not overlay index
                    if(Tree[rank[ProID]].FNPost[cidH]){
                        influence=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(PartiTags[line[i]].first==-1){//if overlay index
                                continue;
                            }
                            if(Tree[rank[ProID]].disPost[i]==Cw+Tree[rank[Cid]].disPost[i]){
                                Tree[rank[ProID]].cntPost[i]-=1;
//                                if(ProID==67975 && Tree[rank[ProID]].vAncestor[i]==67432){
//                                    cout<<"1 Find in SC! "<<ProID<<" "<<Tree[rank[ProID]].vAncestor[i]<<" "<<Tree[rank[ProID]].cntPost[i]<<endl;
//                                }
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProID]].FNPost[cidH]=false;
                        Tree[rank[ProID]].cntPost[cidH]-=1;
//                        if(ProID==67975 && Cid==67432){
//                            cout<<"1 Find in SC! "<<ProID<<" "<<Cid<<" "<<Tree[rank[ProID]].cntPost[cidH]<<endl;
//                        }

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProID]].disPost.size();i++){
                            if(Tree[rank[ProID]].disPost[i]==Cw+Tree[rank[line[i]]].disPost[cidH]){
                                Tree[rank[ProID]].cntPost[i]-=1;
//                                if(ProID==67975 && Tree[rank[ProID]].vAncestor[i]==67432){
//                                    cout<<"1 Find in SC! "<<ProID<<" "<<Tree[rank[ProID]].vAncestor[i]<<" "<<Tree[rank[ProID]].cntPost[i]<<endl;
//                                }
                            }
                        }
                    }
                }
                else{//if Cid is overlay vertex
                    int inf_i=BoundVertexMap[pid][Cid];
//                    if(Tree[rank[ProID]].FNInf[inf_i]){
                        influence=true;
                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProID]].disPost.size();i++){
                            if(PartiTags[line[i]].first==-1){//if overlay index
                                continue;
                            }
                            if(Tree[rank[ProID]].disPost[i]==Cw+Tree[rank[line[i]]].disInf[inf_i]){
                                Tree[rank[ProID]].cntPost[i]-=1;
//                                cout<<"Cid is overlay vertex. "<< ProID<<" ("<<Cid<<") "<<Tree[rank[ProID]].vAncestor[i]<<": "<<Tree[rank[ProID]].cntPost[i]<<endl;
//                                if(ProID==67975 && Tree[rank[ProID]].vAncestor[i]==67432){
//                                    cout<<"2 Find in SC! "<<ProID<<" "<<Tree[rank[ProID]].vAncestor[i]<<" "<<Tree[rank[ProID]].cntPost[i]<<endl;
//                                }
                            }
                        }
//                    }

                }

                if(algoUpdate==PH2H_Cross){
                    if(Tree[rank[ProID]].FN[cidH]){
                        influenceExtend=true;
                        //higher than Cid
                        for(int i=0;i<cidH;i++){
                            if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[Cid]].dis[i]){
                                Tree[rank[ProID]].cnt[i]-=1;
                            }
                        }

                        //equal to Cid
                        Tree[rank[ProID]].FN[cidH]=false;
                        Tree[rank[ProID]].cnt[cidH]-=1;

                        //lower than Cid
                        for(int i=cidH+1;i<Tree[rank[ProID]].dis.size();i++){
                            if(Tree[rank[ProID]].dis[i]==Cw+Tree[rank[line[i]]].dis[cidH]){
                                Tree[rank[ProID]].cnt[i]-=1;
                            }
                        }
                    }
                }
            }

//            if(ProID==84704){
//                cout<<"ProID "<<ProID<<endl;
//            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;
            for(auto it2=Neighbor[ProID].begin();it2!=Neighbor[ProID].end();++it2){
                if(it2->first==Cid){
                    newCw=it2->second;//the weight value in the original graph
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
            if(Wnodes.size()>0){
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

                    if(ssw+wtt<newCw){
                        newCw=ssw+wtt;
                        countwt=1;
                    }else if(ssw+wtt==newCw){
                        countwt+=1;
                    }
                }
            }


            //cout<<Cw<<endl;
            //refresh the shortcut to the new value
            for(int i=0;i<Tree[rank[ProID]].vert.size();i++){
                if(Tree[rank[ProID]].vert[i].first==Cid){
                    Tree[rank[ProID]].vert[i].second.first=newCw;
                    Tree[rank[ProID]].vert[i].second.second=countwt;
                    break;
                }
            }
        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[ProID],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[ProBeginVertexSetParti[pid][i]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }
        if(influenceExtend){//if the distance labeling is dectected changed
            sm->wait();
//            PropagateOverlay.push_back(ProID);
            ProBeginVertexSetNew.clear();
            ProBeginVertexSetNew.reserve(ProBeginVertexSetOverlay.size() + 1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew = rank[ProID], r;
            for (int i = 0; i < ProBeginVertexSetOverlay.size(); i++) {
                r = rank[ProBeginVertexSetOverlay[i]];
                if (LCAQueryOverlay(rnew, r) != rnew) {
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetOverlay[i]);
                }
            }
            ProBeginVertexSetOverlay = ProBeginVertexSetNew;
            sm->notify();
        }
    }

    sm->wait();
    affectedParti.insert(pid);
    sm->notify();



    //return checknum;
}

void Graph::IncreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU){
    int checknum=0;
    map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]) maintain the old distance before refreshed and avoid search in the adjacent list
    OCdis.clear();

    //NodeOrderss.clear();
//    NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
    vector<set<int>> SCre; //SCre.clear();
    set<int> ss; //ss.clear();
    SCre.assign(node_num,ss);//{vertexID, set<int>}
    set<OrderCompMin> OC; OC.clear();//vertexID in decreasing node order

    for(int k=0;k<wBatch.size();k++){
        int a=wBatch[k].first.first;
        int b=wBatch[k].first.second;
        int oldW=wBatch[k].second.first;
        int newW=wBatch[k].second.second;

        if(oldW!=newW){
            if(Neighbors[a].find(b)!=Neighbors[a].end()){
                Neighbors[a][b]=newW;
            }else{
                cout<<"Not found edge! "<<endl; exit(1);
            }

            if(Neighbors[b].find(a)!=Neighbors[b].end()){
                Neighbors[b][a]=newW;
            }else{
                cout<<"Not found edge! "<<endl; exit(1);
            }

            int lid,hid;
            if(NodeOrder[a]<NodeOrder[b]){
                lid=a;hid=b;
            }else{
                lid=b;hid=a;
            }
            int lidM=IDMap[lid];

            for(int i=0;i<Tree[rank[lidM]].vert.size();i++){
                if(Tree[rank[lidM]].vert[i].first==hid){
                    if(Tree[rank[lidM]].vert[i].second.first==oldW){
                        Tree[rank[lidM]].vert[i].second.second-=1;
                        if(Tree[rank[lidM]].vert[i].second.second<1){
                            OCdis[make_pair(lid,hid)]=oldW;
                            SCre[lid].insert(hid);
                            OC.insert(OrderCompMin(lid));
                        }
                    }
                    break;
                }
            }
        }
    }

//    vector<int> ProBeginVertexSet; ProBeginVertexSet.clear();
    ProBeginVertexSetParti[pid].clear();
    vector<int> ProBeginVertexSetNew;
    bool influence;
    int ProID, ProIDM; vector<int> line;
    while(!OC.empty()){
        ProID=(*OC.begin()).x; ProIDM=IDMap[ProID];
        OC.erase(OC.begin());
        vector<pair<int,pair<int,int>>> Vert=Tree[rank[ProIDM]].vert;
        influence=false;

        //each ProID corresponds to a line
        line.clear(); line.reserve(heightMax);
        int pachid=ProID;
        while(Tree[rank[IDMap[pachid]]].height>1){
            line.insert(line.begin(),pachid);
            pachid=Tree[Tree[rank[IDMap[pachid]]].pa].uniqueVertex;
        }
        line.insert(line.begin(),pachid);

        for(auto it=SCre[ProID].begin();it!=SCre[ProID].end();it++){
            int Cid=*it; int Cw=OCdis[make_pair(ProID,Cid)];
            int cidH=Tree[rank[IDMap[Cid]]].height-1;
            int CidM=IDMap[Cid];

            map<int,int> Hnei; //Hnei.clear();
            vector<pair<int,int>> Lnei; //Lnei.clear();
            for(int j=0;j<Vert.size();j++){
                if(NodeOrder[Vert[j].first]>NodeOrder[Cid]){
                    Hnei[Vert[j].first]=Vert[j].second.first;
                }else if(NodeOrder[Vert[j].first]<NodeOrder[Cid]){
                    Lnei.push_back(make_pair(Vert[j].first,Vert[j].second.first));
                }
            }
            //check the affected shortcuts
            int hid2,lid2,lid2M;
            for(int j=0;j<Tree[rank[CidM]].vert.size();j++){
                hid2=Tree[rank[CidM]].vert[j].first;
                if(Hnei.find(hid2)!=Hnei.end()){
                    if(Cw+Hnei[hid2]==Tree[rank[CidM]].vert[j].second.first){
                        Tree[rank[CidM]].vert[j].second.second-=1;
                        if(Tree[rank[CidM]].vert[j].second.second<1){
                            SCre[Cid].insert(hid2);
                            OC.insert(OrderCompMin(Cid));
                            OCdis[make_pair(Cid,hid2)]=Cw+Hnei[hid2];
                        }
                    }
                }
            }
            for(int j=0;j<Lnei.size();j++){
                lid2=Lnei[j].first; lid2M=IDMap[lid2];
                for(int k=0;k<Tree[rank[lid2M]].vert.size();k++){
                    if(Tree[rank[lid2M]].vert[k].first==Cid){
                        if(Tree[rank[lid2M]].vert[k].second.first==Cw+Lnei[j].second){
                            Tree[rank[lid2M]].vert[k].second.second-=1;
                            if(Tree[rank[lid2M]].vert[k].second.second<1){
                                SCre[lid2].insert(Cid);
                                OC.insert(OrderCompMin(lid2));
                                OCdis[make_pair(lid2,Cid)]=Cw+Lnei[j].second;
                            }
                        }
                        break;
                    }
                }
            }


            //before Cw=d(ProID,Cid) gets its new value, we first check which dis it will influence
            if(Tree[rank[ProIDM]].FN[cidH]){
                influence=true;
                //higher than Cid
                for(int i=0;i<cidH;i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[CidM]].dis[i]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }

                //equal to Cid
                Tree[rank[ProIDM]].FN[cidH]=false;
                Tree[rank[ProIDM]].cnt[cidH]-=1;

                //lower than Cid
                for(int i=cidH+1;i<Tree[rank[ProIDM]].dis.size();i++){
                    if(Tree[rank[ProIDM]].dis[i]==Cw+Tree[rank[IDMap[line[i]]]].dis[cidH]){
                        Tree[rank[ProIDM]].cnt[i]-=1;
                    }
                }
            }

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            Cw=INF; int countwt=0;

            for(auto it2=Neighbors[ProID].begin();it2!=Neighbors[ProID].end();++it2){
                if(it2->first==Cid){
                    Cw=it2->second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid,widM;
            vector<pair<int,int>> Wnodes;
            Wnodes.clear();
            /*if(SCconNodes.find(make_pair(ProID,Cid))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(ProID,Cid)];
            else if(SCconNodes.find(make_pair(Cid,ProID))!=SCconNodes.end())
                Wnodes=SCconNodes[make_pair(Cid,ProID)];*/
            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(Wnodes.size()>0){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i].first; widM=IDMap[wid];
                    for(int j=0;j<Tree[rank[widM]].vert.size();j++){
                        if(Tree[rank[widM]].vert[j].first==ProID){
                            ssw=Tree[rank[widM]].vert[j].second.first;
                        }
                        if(Tree[rank[widM]].vert[j].first==Cid){
                            wtt=Tree[rank[widM]].vert[j].second.first;
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
            for(int i=0;i<Tree[rank[ProIDM]].vert.size();i++){
                if(Tree[rank[ProIDM]].vert[i].first==Cid){
                    Tree[rank[ProIDM]].vert[i].second.first=Cw;
                    Tree[rank[ProIDM]].vert[i].second.second=countwt;
                    break;
                }
            }
        }

        if(influence){
            ProBeginVertexSetNew.clear(); ProBeginVertexSetNew.reserve(ProBeginVertexSetParti[pid].size()+1);
            ProBeginVertexSetNew.push_back(ProID);
            int rnew=rank[IDMap[ProID]],r;
            for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
                r=rank[IDMap[ProBeginVertexSetParti[pid][i]]];
                if(LCAQueryOverlay(rnew,r)!=rnew){
                    ProBeginVertexSetNew.push_back(ProBeginVertexSetParti[pid][i]);
                }
            }
            ProBeginVertexSetParti[pid]=ProBeginVertexSetNew;
        }

    }

    if(ifLabelU){
        int ProBeginVertexID;
        for(int i=0;i<ProBeginVertexSetParti[pid].size();i++){
            ProBeginVertexID=ProBeginVertexSetParti[pid][i];
            vector<int> linee; //linee.clear();
            linee.reserve(heightMax);
            int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
            while(Tree[rank[IDMap[pachidd]]].height>1){
                linee.insert(linee.begin(),pachidd);
                pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
            }
            linee.insert(linee.begin(),pachidd);

            eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
        }
    }

    //return checknum;
}

void Graph::IncreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax,  vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid) {
    int ProBeginVertexID;
    int checknum=0;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[IDMap[ProBeginVertexID]]].pa].uniqueVertex;
        while(Tree[rank[IDMap[pachidd]]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[IDMap[pachidd]]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);

        eachNodeProcessIncrease1Parti(rank[IDMap[ProBeginVertexID]], linee,checknum,Tree,rank,VidtoTNid);
    }
}

void Graph::IncreasePartiBatchLabelPostMHLExtend(int pid, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int> &ProBeginVertexSet, vector<vector<int>> &VidtoTNid) {
    int ProBeginVertexID=partiRoots[pid];
    int rootRank=rank[ProBeginVertexID];
    int checknum=0;
//    vector<int> linee; //linee.clear();
//    linee=Tree[Tree[rootRank].pa].vAncestor;
//    EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,Tree,rank);
//    cout<<"ProBeginVertexSetPartiExtend of Parti "<<pid<<": "<<ProBeginVertexSetPartiExtend[pid].size()<<endl;
    for(int i=0;i<ProBeginVertexSet.size();i++){
        ProBeginVertexID=ProBeginVertexSet[i];
//        cout<<"ProBeginVertexID: "<<ProBeginVertexID<<"("<<PartiTags[ProBeginVertexID].first<<","<<Tree[rank[ProBeginVertexID]].height<<")"<<endl;
        vector<int> linee; //linee.clear();
        linee.reserve(heightMax);
        int pachidd=Tree[Tree[rank[ProBeginVertexID]].pa].uniqueVertex;
        while(Tree[rank[pachidd]].height>1){
            linee.insert(linee.begin(),pachidd);
            pachidd=Tree[Tree[rank[pachidd]].pa].uniqueVertex;
        }
        linee.insert(linee.begin(),pachidd);
        eachNodeProcessIncrease1PostMHL(rank[ProBeginVertexID], linee,checknum,Tree,rank,VidtoTNid);
    }

}

void Graph::IncreasePartiBatchLabelPostMHLExtendV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>> &ProBeginVertexSetV, vector<vector<int>> &VidtoTNid) {
    int pid;
    for(int i=0;i<p.size();++i){
        pid=p[i];
        IncreasePartiBatchLabelPostMHLExtend(pid, Tree, rank, heightMax, ProBeginVertexSetV[pid], VidtoTNid);
    }
}

int Graph::ShortcutDisCheck(int ID1, int ID2){
    int d=INF;
    int ID;
    if(PartiTag[ID1].first!=PartiTag[ID2].first){
        cout<<ID1<<" and "<<ID2<<" are not in the same partition!"<<endl; exit(1);
    }
    int PID=PartiTag[ID1].first;
    int lid,hid;
    if(NodeOrder[ID1]>NodeOrder[ID2]){
        hid=ID1, lid=ID2;
    }else{
        hid=ID2, lid=ID1;
    }
    bool flag=false;
/*    for(int i=0;i<Trees[PID][ranks[PID][IDMap[lid]]].vert.size();i++){
        if(Trees[PID][ranks[PID][IDMap[lid]]].vert[i].first==hid){
            d=Trees[PID][ranks[PID][IDMap[lid]]].vert[i].second.first;
            flag=true;
            break;
        }
    }*/

    int Cw=INF;
    for(auto it2=NeighborsParti[lid].begin();it2!=NeighborsParti[lid].end();++it2){
        if(it2->first==hid){
            Cw=it2->second;//the weight value in the original graph
//            cout<<lid<<" "<<hid<<", Cw1: "<<Cw<<endl;
            flag=true;
            break;
        }
    }

    if(SCconNodesMTP[lid].find(hid) != SCconNodesMTP[lid].end()){//if found
        int widM;
        for(auto it=SCconNodesMTP[lid][hid].begin();it!=SCconNodesMTP[lid][hid].end();++it){
            ID=it->first;
            if(PartiTag[ID].second){//if boundary vertex
//                cout<<"Continue "<<ID<<endl;
                continue;
            }
            widM=IDMap[ID];
            int ssw=-1,wtt=-1;
            for(int j=0;j<Trees[PID][ranks[PID][widM]].vert.size();j++){
                if(Trees[PID][ranks[PID][widM]].vert[j].first==lid){
                    ssw=Trees[PID][ranks[PID][widM]].vert[j].second.first;
                    break;
                }

            }
            for(int j=0;j<Trees[PID][ranks[PID][widM]].vert.size();j++) {
                if (Trees[PID][ranks[PID][widM]].vert[j].first == hid) {
                    wtt = Trees[PID][ranks[PID][widM]].vert[j].second.first;
                    break;
                }
            }
            if(ssw==-1 || wtt==-1){
                cout<<"Wrong! "<<ssw<<" "<<wtt<<endl; exit(1);
            }

            if(ssw+wtt<Cw){
                Cw=ssw+wtt;
//                cout<<lid<<" "<<hid<<", Cw2: "<<Cw<<endl;
            }
        }
        //cout<<Cw<<endl;
        flag=true;
    }
//    else{
//        cout<<"Shortcut sc("<<ID1<<","<<ID2<<") does not exist!"<<endl;
//    }
    d=Cw;
    if(!flag){
        cout<<"Not found! "<<endl; exit(1);
    }
    return d;
}
