/*
 * Labelcon.cpp
 *
 *  Created on: 22 Dec 2020
 *      Author: zhangmengxuan
 */
#include "headH2H.h"

vector<int> NodeOrder_;//nodeID order
vector<int> _DD_;//true degree, temporal degree ,_DD2_
vector<int> NodeOrders;

///////////////////////////// Index Construction ////////////////////////////////

void Graph::IndexConstruction(int algo){
    algoIndex = algo;
    switch (algoIndex) {
        case 0:{
            cout<<"Dijkstra's algorithm."<<endl;
            break;
        }
        case 1:{
            cout<<"CH algorithm."<<endl;
            CHIndexConstruct();
            break;
        }
        case 2:{
            cout<<"H2H algorithm."<<endl;
            H2HIndexConstruct();
            break;
        }
        case 3:{
            cout<<"BiDijkstra's algorithm."<<endl;
            break;
        }
        default:{
            cout<<"Wrong SP index! "<<algoIndex<<endl; exit(1);
        }
    }
}

//// For CH index construction
void Graph::CHIndexConstruct(){
    string orderfile=sourcePath+dataset+".order";
    double runT1=0, runT2=0, runT3=0;
    Timer tt;

    tt.start();
    MDEContraction(orderfile);
    tt.stop();
    runT1=tt.GetRuntime();
    cout<<"Time for MDE contraction: "<<runT1<<" s."<<endl;

    IndexsizeCHWP();
}

//function for computing the index size
void Graph::IndexsizeCHWP(){
    unsigned long long m=0,m1=0,m2=0,m3=0,m4=0;
    //Overlay index
    for(int i=0;i<NeighborCon.size();i++){
        m1+=NeighborCon[i].size()*2*sizeof(int);//dis
    }

    for(int i=0;i< SCconNodesMT.size();i++){
        for(auto it=SCconNodesMT[i].begin(); it!=SCconNodesMT[i].end(); it++){
            m4+=sizeof(int)+(*it).second.size()*2*sizeof(int);
        }
    }

    //cout<<"Index size "<<(double)m1/1024/1024<<", Pruning point size "<<(double)m2/1024/1024<<endl;
    m=m1+m2+m3+m4;
    cout<<"CH label size: "<<(double)m1/1024/1024<<" MB"<<endl;
    cout<<"CH Update information size: "<<(double)m4/1024/1024<<" MB"<<endl;
    cout<<"Overall index size "<<(double)m/1024/1024<<" MB"<<endl;
}

void Graph::H2HIndexConstruct() {
    string orderfile=sourcePath+dataset+".order";
    double runT1=0, runT2=0, runT3=0;
    Timer tt;

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
        map<int, vector<int>> mi;
        SCconNodesMT.assign(node_num, mi);

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
        cout<<"Reading vertex ordering... "<<orderfile<<endl;
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

        map<int, vector<int>> mi;
        SCconNodesMT.assign(node_num, mi);//record the supportive vertices of a shortcut, only record edge once by leveraging the ID positions of endpoints

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
                            insertEorder(ID1,ID2,Neigh[i].second.first+Neigh[j].second.first);
                            if(ID1<ID2){
                                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                                    SCconNodesMT[ID1].insert({ID2,vector<int>()});
                                }
                                SCconNodesMT[ID1][ID2].push_back(x);//only record onece
                            }
                            else if(ID2<ID1){
                                if(SCconNodesMT[ID2].find(ID1)==SCconNodesMT[ID2].end()){//if not found
                                    SCconNodesMT[ID2].insert({ID1,vector<int>()});
                                }
                                SCconNodesMT[ID2][ID1].push_back(x);
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
                            thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
                        }
                        thread.join_all();
                    }else{
                        boost::thread_group thread;
                        for(int i=0;i<Neigh.size();i++){
                            pair<int,int> p;
                            p.first=i; p.second=(i+1);
                            thread.add_thread(new boost::thread(&Graph::NeighborComorder, this, boost::ref(Neigh), p, x));
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

    int treewidth=0;
	int nn;
	for(;len>=0;len--){
		int x=vNodeOrder[len];
		Node nod;
		nod.vert=NeighborCon[x];
		nod.uniqueVertex=x;
        if(treewidth<NeighborCon[x].size()+1){
            treewidth=NeighborCon[x].size()+1;
        }
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

    cout<<"Tree height: "<<heightMax<<" ; treewidth: "<<treewidth<<endl;
}
//function of H2H index construction
void Graph::makeIndex(){
    cout<<"Building H2H index..."<<endl;
	makeRMQ();

	//initialize
	vector<int> list; //list.clear();
	list.push_back(Tree[0].uniqueVertex);
	Tree[0].pos.clear();
	Tree[0].pos.push_back(0);
    Tree[0].vAncestor=list;

	for(int i=0;i<Tree[0].ch.size();i++){
		makeIndexDFS(Tree[0].ch[i],list);
	}

}
/// Functions for MDE contraction
void Graph::deleteEOrderGenerate(int u,int v){
    if(E[u].find(v)!=E[u].end()){
        E[u].erase(E[u].find(v));
        DD[u]--;
    }

    if(E[v].find(u)!=E[v].end()){
        E[v].erase(E[v].find(u));
        DD[v]--;
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
            insertEMTorder(ID1, ID2, w1+w2);
            if(ID1<ID2){
                if(SCconNodesMT[ID1].find(ID2)==SCconNodesMT[ID1].end()){//if not found
                    SCconNodesMT[ID1].insert({ID2,vector<int>()});
                }
                SCconNodesMT[ID1][ID2].push_back(x);
            }

        }
    }
//    sm->notify();
}
/// Functions for Tree contraction
int Graph::match(int x,vector<pair<int,pair<int,int>>> &vert){
    int nearest=vert[0].first;
    for(int i=1;i<vert.size();i++){
        if(rank[vert[i].first]>rank[nearest])
            nearest=vert[i].first;
    }
    int p=rank[nearest];
    return p;
}

/// Functions for Tree index contraction
void Graph::makeRMQ(){
    //EulerSeq.clear();
    toRMQ.assign(node_num,0);
    //RMQIndex.clear();
    makeRMQDFS(0, 1);
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
void Graph::makeRMQDFS(int p, int height){
    toRMQ[p] = EulerSeq.size();
    EulerSeq.push_back(p);
    for (int i = 0; i < Tree[p].ch.size(); i++){
        makeRMQDFS(Tree[p].ch[i], height + 1);
        EulerSeq.push_back(p);
    }
}

void Graph::makeIndexDFS(int p, vector<int>& list){
    //initialize
    int NeiNum=Tree[p].vert.size();
    Tree[p].pos.assign(NeiNum+1,0);
    Tree[p].dis.assign(list.size(),INF);
    Tree[p].cnt.assign(list.size(),0);
    Tree[p].FN.assign(list.size(),true);

    //pos
    //map<int,Nei> Nmap; Nmap.clear();//shortcut infor ordered by the pos ID
    for(int i=0;i<NeiNum;i++){
        for(int j=0;j<list.size();j++){
            if(Tree[p].vert[i].first==list[j]){//get the original distance information by shortcuts
                Tree[p].pos[i]=j;
                Tree[p].dis[j]=Tree[p].vert[i].second.first;
                Tree[p].cnt[j]=1;
                break;
            }
        }
    }
    Tree[p].pos[NeiNum]=list.size();
//    Tree[p].dis.push_back(0);//distance to itself
    Tree[p].vAncestor=list;
    Tree[p].vAncestor.push_back(Tree[p].uniqueVertex);//the last vertex is the tree node

    //dis
    for(int i=0;i<NeiNum;i++){
        int x=Tree[p].vert[i].first;
        int disvb=Tree[p].vert[i].second.first;
        int k=Tree[p].pos[i];//the k-th ancestor is x, the i-th neighbor is also the k-th ancesotr

        for(int j=0;j<list.size();j++){
            int y=list[j];//the j-th ancestor is y

            int z;//the distance from x to y
            if(k!=j){
                if(k<j)//if x is the ancestor of y
                    z=Tree[rank[y]].dis[k];
                else if(k>j)//if y is the ancestor of x
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


///////////////////////////// Query Processing ////////////////////////////////

//function for correctness check
void Graph::CorrectnessCheck(int runtimes){
    if(algoIndex==0||algoIndex==3){
        return;
    }
    Timer tt;
    double runT=0;
    srand (time(NULL));
    int s, t, d1=0, d2=0, d3;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ... ";
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
//        s=115113,t=178597;//
//        cout<<"Query "<<i<<": "<<s<<" "<<t<<endl;

//        if(runtimes == 1){
//            cout<<"s: "<<s<<" ; t: "<<t<<endl;
//        }
        tt.start();
        d1=Dijkstra(s,t,Neighbor);
        tt.stop();
        if(algoIndex==0){
            d2=d1;
        }else if(algoIndex==1){
            tt.start();
            d2=QueryCHWP(s,t);
            tt.stop();
        }else if(algoIndex==2){
            tt.start();
            d2=QueryH2H(s,t);
            tt.stop();
        }else if(algoIndex==3){
            tt.start();
            d2=BiDijkstra(s,t,Neighbor);
            tt.stop();
        }

        runT+=tt.GetRuntime();
//        cout<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<") "<<d2<<" "<<d1<<endl;
        if(d1!=d2){
            cout<<"InCorrect! "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<endl;
            exit(1);
        }
    }
    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms."<<endl;
}

void Graph::EffiCheck(int runtimes){
    string ODfile=sourcePath+dataset+".query";

    if(algoIndex==1){
        EffiCheckCH(ODfile,runtimes);//query efficiency test
//        EffiCheckCH(ODfile+"Parti",runtimes);//query efficiency test
//        EffiCheckCH(ODfile+"SameParti",runtimes);//query efficiency test
//        EffiCheckCH(ODfile+"CrossParti",runtimes);
    }else if(algoIndex==2){
        EffiCheckH2H(ODfile,runtimes);//query efficiency test
//        EffiCheckH2H(ODfile+"Parti",runtimes);//query efficiency test
//        EffiCheckH2H(ODfile+"SameParti",runtimes);//query efficiency test
//        EffiCheckH2H(ODfile+"CrossParti",runtimes);
    }else if(algoIndex==3){
        runtimes/=100;
        ifstream IF(ODfile);
        if(!IF){
            cout<<"Cannot open Map "<<ODfile<<endl;
            exit(1);
        }
        cout<<"Query file: "<<ODfile<<endl;
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
        double runT=0;
        int d1, d2;
        Timer tt;
        clock_t start = clock();

        vector<int> results(runtimes,-1);
        for(int i=0;i<runtimes;i++){
            s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
            tt.start();
            d2=BiDijkstra(s,t,Neighbor);
            tt.stop();
            runT += tt.GetRuntime();
            results[i]=d2;
//        cout<<d2<<" ";//<<endl;
//        if(d1!=d2){
//            cout<<"Wrong! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
//        }
        }

//    cout<<endl;

        cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms. "<<1000*(double)(clock() - start) / (CLOCKS_PER_SEC*runtimes)<<" ms."<<endl;
    }
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

            for(auto out=NeighborCon[topNodeIDForward].begin();out!=NeighborCon[topNodeIDForward].end();out++){
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

            for(auto in=NeighborCon[topNodeIDBackward].begin();in!=NeighborCon[topNodeIDBackward].end();in++){
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

//function for efficiency test
void Graph::EffiCheckH2H(string filename,int runtimes){
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
    double runT=0;
    int d1, d2;
    Timer tt;
    clock_t start = clock();

    vector<int> results(runtimes,-1);
    for(int i=0;i<runtimes;i++){
        s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
        tt.start();
        d2=QueryH2H(s,t);
        tt.stop();
        runT += tt.GetRuntime();
        results[i]=d2;
//        cout<<d2<<" ";//<<endl;
//        if(d1!=d2){
//            cout<<"Wrong! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
//        }
    }

//    cout<<endl;

    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms. "<<1000*(double)(clock() - start) / (CLOCKS_PER_SEC*runtimes)<<" ms."<<endl;
}

//function for correctness check
void Graph::CorrectnessCheckH2H(int runtimes){
    Timer tt;
    double runT=0;
    srand (time(NULL));
    int s, t, d1, d2, d3;
//    runtimes = 1;
    cout<<"Correctness check ("<<runtimes<<" rounds) ... ";
    for(int i=0;i<runtimes;i++){
//        if(i%100==0) cout<<i<<endl;
        s=rand()%node_num;
        t=rand()%node_num;
//        s=PartiVertex[9][rand()%PartiVertex[9].size()];
//        t=PartiVertex[9][rand()%PartiVertex[9].size()];
//        s=214318,t=159753;//NY
//        cout<<"Query "<<i<<": "<<s<<" "<<t<<endl;

        if(runtimes == 1){
//            cout<<"s: "<<s<<" ; t: "<<t<<endl;
        }
        d1=Dijkstra(s,t,Neighbor);

        tt.start();
        d2=QueryH2H(s,t);
        tt.stop();
        runT+=tt.GetRuntime();
//        cout<<s<<"("<<CoreTag[s]<<") "<<t<<"("<<CoreTag[t]<<") "<<d2<<" "<<d1<<endl;
        if(d1!=d2){
            cout<<"InCorrect! "<<i<<": "<<s<<"("<<NodeOrder[s]<<") "<<t<<"("<<NodeOrder[t]<<") "<<d2<<" "<<d1<<endl;
            exit(1);
        }
    }
    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms."<<endl;
}



//function for efficiency test
void Graph::EffiCheckCH(string filename,int runtimes){
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
        ODpair.emplace_back(ID1, ID2);
    }
    if(runtimes > num){
        runtimes = num;
    }
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    double runT=0;
    int d1, d2;
    Timer tt;
    clock_t start = clock();

    vector<int> results(runtimes,-1);
    for(int i=0;i<runtimes;i++){
        s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
        tt.start();
        d2=QueryCHWP(s,t);
        tt.stop();
        runT += tt.GetRuntime();
        results[i]=d2;
//        cout<<d2<<" ";//<<endl;
//        if(d1!=d2){
//            cout<<"Wrong! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
//        }
    }

//    cout<<endl;

    cout<<"Average Query Time: "<<1000*runT/runtimes<<" ms. "<<1000*(double)(clock() - start) / (CLOCKS_PER_SEC*runtimes)<<" ms."<<endl;
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

///////////////////////////// Index Maintenance ////////////////////////////////

//function for efficiency test
unsigned long long Graph::EffiCheckThroughput(vector<pair<int,int>>& ODpair, Timer& tRecord, int batchInterval){

    int s, t;
    double runT=0;
    int d1, d2;
    Timer tt;
    unsigned long long runtimes = 0;
    vector<int> results(ODpair.size(),-1);
    int i=0;
    double tNow=0;
    if(algoIndex==0){
        while(true){
            tRecord.stop();
            tNow=tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if(tNow>batchInterval){
                break;
            }

            s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
            tt.start();
            d2=Dijkstra(s,t,Neighbor);
            tt.stop();
            runT += tt.GetRuntime();
            results[i]=d2;
            ++runtimes;
            ++i;
            if(i==ODpair.size()){
                i=0;
            }
//        cout<<d2<<" ";//<<endl;
//        if(d1!=d2){
//            cout<<"Wrong! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
//        }
        }
    }else if(algoIndex==1){
        while(true){
            tRecord.stop();
            tNow=tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if(tNow>batchInterval){
                break;
            }

            s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
            tt.start();
            d2=QueryCHWP(s,t);
            tt.stop();
            runT += tt.GetRuntime();
            results[i]=d2;
            ++runtimes;
            ++i;
            if(i==ODpair.size()){
                i=0;
            }
//        cout<<d2<<" ";//<<endl;
//        if(d1!=d2){
//            cout<<"Wrong! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
//        }
        }
    }else if(algoIndex==2){
        while(true){
            tRecord.stop();
            tNow=tRecord.GetRuntime();
//        cout<<"tNow: "<<tNow<<endl;
            if(tNow>batchInterval){
                break;
            }

            s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
            tt.start();
            d2=QueryH2H(s,t);
            tt.stop();
            runT += tt.GetRuntime();
            results[i]=d2;
            ++runtimes;
            ++i;
            if(i==ODpair.size()){
                i=0;
            }
//        cout<<d2<<" ";//<<endl;
//        if(d1!=d2){
//            cout<<"Wrong! "<<s<<" "<<t<<" "<<d2<<" "<<d1<<endl; exit(1);
//        }
        }
    }


    cout<<"Throughput number: "<<runtimes;
    if(runtimes>0) {
        cout << " ; Average Query Time: " << 1000 * runT / runtimes << " ms.";
    }
    cout<<endl;
    return runtimes;
}

//function for efficiency test
unsigned long long Graph::EffiCheckThroughput(vector<pair<int,int>>& ODpair, int runtimes, int batchInterval, double& updateT, double& queryT){
    int s, t;
    double runT=0;
    int d1, d2;
    Timer tt;
    vector<int> results(ODpair.size(),-1);
    unsigned long long throughput=0;
    if(runtimes>ODpair.size()){
        runtimes=ODpair.size();
        cout<<"Actual runtimes: "<<runtimes<<endl;
    }

    if(algoIndex==1){
        runtimes=runtimes/10;
        for(int i=0;i<runtimes;++i){
            s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
            tt.start();
            d2=QueryCHWP(s,t);
            tt.stop();
            runT += tt.GetRuntime();
            results[i]=d2;
        }

    }else if(algoIndex==2){
        for(int i=0;i<runtimes;++i){
            s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
            tt.start();
            d2=QueryH2H(s,t);
            tt.stop();
            runT += tt.GetRuntime();
            results[i]=d2;
        }

    }
    double effiT=runT/runtimes;
    queryT+=effiT;
    if(updateT<batchInterval){
        throughput=(batchInterval-updateT)/effiT;
    }

    cout<<"Throughput number: "<<throughput << " ; Average Query Time: " << 1000 * effiT << " ms."<<endl;
    return throughput;
}


//function of testing the throughput of path-finding system, batchInterval is the time interval between two adjacent update batch (in seconds)
void Graph::SPThroughputTest(int updateType, bool ifBatch, int batchNum, int batchSize, int batchInterval, int runtimes) {
    cout<<"Shortest path query throughput test..."<<endl;
    // read updates
    string file = sourcePath+dataset + ".update";
    bool ifDebug=false;
//    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand (0);
    vector<pair<pair<int,int>,int>> updateData;
    ReadUpdate(file, updateData);

    if(batchNum*batchSize>updateData.size()){
        batchNum=floor(updateData.size()/batchSize);
    }
    cout<<"Update batch number: "<<batchNum<<" ; Batch size: "<<batchSize<<" ; Batch interval: "<< batchInterval<<endl;


    string filename = sourcePath+dataset + ".query";
    if(queryFName!=""){
        filename=queryFName;
    }
    ifstream IF(filename);
    if(!IF){
        cout<<"Cannot open Map "<<filename<<endl;
        exit(1);
    }
    cout<<"Query file: "<<filename<<endl;

    int num;
    vector<pair<int,int>> ODpair;
    IF>>num;
    for(int k=0;k<num;k++){
        IF>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    unsigned long long throughputNum=0;

    Timer tt;
    Timer tRecord;
    double runT1=0, runT2 = 0;
    double runT=0;
    double queryT=0;
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

                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;
                        oldW = updateData[update_i].second;
                        newW=oldW*0.5;
                        if(newW < 1) {
                            cout<<"New edge weight is not positive! "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                            exit(1);
                        }
                        if(ifDebug){
//                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                        ++update_i;
                    }
                    tRecord.start();
                    if(algoIndex==0){//Dijkstra
                        tt.start();
                        int a,b;
                        for(int k=0;k<wBatch.size();k++) {
                            a = wBatch[k].first.first;
                            b = wBatch[k].first.second;
                            newW = wBatch[k].second.second;

                            //modify the information in original graph
                            for (int i = 0; i < Neighbor[a].size(); i++) {
                                if (Neighbor[a][i].first == b) {
                                    Neighbor[a][i].second = newW;
                                    break;
                                }
                            }
                            for (int i = 0; i < Neighbor[b].size(); i++) {
                                if (Neighbor[b][i].first == a) {
                                    Neighbor[b][i].second = newW;
                                    break;
                                }
                            }
                        }
                        tt.stop();
                    }else if(algoIndex==1){//CH update
                        tt.start();
                        CHdecBat(wBatch);
                        tt.stop();
                    } else if(algoIndex==2){//H2H update
                        tt.start();
                        H2HdecBat(wBatch);
                        tt.stop();
                    }
                    runT=tt.GetRuntime();
                    runT1 += runT;
                    cout<<"Batch "<<u<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
//                    throughputNum += EffiCheckThroughput(ODpair,tRecord,batchInterval);//query efficiency test
                    throughputNum += EffiCheckThroughput(ODpair,runtimes,batchInterval,runT,queryT);

                    if(ifDebug){
//                        g2.CorrectnessCheckH2H(100);
                        CorrectnessCheck(100);
                    }
                }
                cout<<"\nOverall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average Decrease batch update Time: "<<runT1/batchNum<<" s; Average query time: "<<1000*queryT/batchNum<<" ms."<<endl;
            }
            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            if(ifBatch){//for batch update
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<batchNum;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;
                        oldW = updateData[update_i].second;
                        newW = oldW*1.5;
                        if(ifDebug){
//                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                        ++update_i;
                    }
                    tRecord.start();
                    if(algoIndex==0){
                        tt.start();
                        for(int wb=0;wb<wBatch.size();wb++) {
                            int a = wBatch[wb].first.first;
                            int b = wBatch[wb].first.second;
                            int oldW = wBatch[wb].second.first;
                            int newW = wBatch[wb].second.second;

                            //modify the original graph information
                            for (int i = 0; i < Neighbor[a].size(); i++) {
                                if (Neighbor[a][i].first == b) {
                                    Neighbor[a][i].second = newW;
                                    break;
                                }
                            }
                            for (int i = 0; i < Neighbor[b].size(); i++) {
                                if (Neighbor[b][i].first == a) {
                                    Neighbor[b][i].second = newW;
                                    break;
                                }
                            }
                        }
                        tt.stop();
                    }else if(algoIndex==1){//CH update
                        tt.start();
                        CHincBatMT(wBatch);
                        tt.stop();
                    } else if(algoIndex==2){//H2H update
                        tt.start();
                        H2HincBatMT(wBatch);
                        tt.stop();
                    }
                    runT=tt.GetRuntime();
                    runT2 += runT;
                    cout<<"Batch "<<u<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
//                    throughputNum += EffiCheckThroughput(ODpair,tRecord,batchInterval);//query efficiency test
                    throughputNum += EffiCheckThroughput(ODpair,runtimes,batchInterval,runT,queryT);
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }
                cout<<"\nOverall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average Increase batch update Time: "<<runT2/batchNum<<" s; Average query time: "<<1000*queryT/batchNum<<" ms."<<endl;
            }


            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

//function of testing the throughput on real-life updates
void Graph::RealUpdateThroughputTest(string updateFile){
    bool ifDebug=false;
    ifDebug=true;
    int runtimes=10000;
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    string line;
    vector<string> vs;
    int ID1,ID2,oldW,newW,weight;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==2);
    int batchNum=stoi(vs[0]);
    int batchInterval=stoi(vs[1]);
    int batchSize;
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNum);
    long long int aveBatchSize=0;
    int maxBatchSize=0, minBatchSize=INT32_MAX;
    for(int i=0;i<batchNum;++i){
        getline(IF,line);
        vs.clear();
        boost::split(vs,line,boost::is_any_of(" "));
        batchSize=stoi(vs[0]);
        assert(vs.size()==3*batchSize+1);
        aveBatchSize+=batchSize;
        if(maxBatchSize<batchSize) maxBatchSize=batchSize;
        if(minBatchSize>batchSize) minBatchSize=batchSize;
        for(int j=0;j<batchSize;++j){
            ID1=stoi(vs[3*j+1]), ID2=stoi(vs[3*j+2]), weight=stoi(vs[3*j+3]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    IF.close();
    aveBatchSize/=batchNum;
    cout<<"Update batch: "<<batchNum<<" ; Average batch size: "<<aveBatchSize<<" ; Maximal batch size: "<<maxBatchSize<<" ; Minimal batch size: "<<minBatchSize<<" ; Batch interval: "<< batchInterval<<endl;

    string queryF = sourcePath+dataset + ".query";

    ifstream IF2(queryF);
    if(!IF2){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF2>>num;
    for(int k=0;k<num;k++){
        IF2>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF2.close();

    //index maintenance
    Timer tRecord;
    double runT1=0, runT2 = 0;
    unsigned long long throughputNum=0;
    double runT=0;
    double queryT=0;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;
    Timer tt;
//    map<pair<int,int>,int> uEdges;
    for(int i=0;i<batchNum;++i){
        wBatchDec.clear(); wBatchInc.clear();
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[i].size();++j){
            ID1=batchUpdates[i][j].first.first, ID2=batchUpdates[i][j].first.second, weight=batchUpdates[i][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                if(it->first==ID2){
                    ifFind=true;
                    oldW=it->second;
                    if(oldW>weight){
                        wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                    }else if(it->second<weight){
                        wBatchInc.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                    }
                    break;
                }
            }
            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),weight});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<weight<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
//                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<i<<" . Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;

        //Step 1: Decrease updates
        /*if(!wBatchDec.empty()){
            cout<<"Decrease update."<<endl;
            tRecord.start();
            if(algoIndex==0){//Dijkstra
                tt.start();
                int a,b;
                for(int k=0;k<wBatchDec.size();k++) {
                    a = wBatchDec[k].first.first;
                    b = wBatchDec[k].first.second;
                    newW = wBatchDec[k].second.second;

                    //modify the information in original graph
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
                            Neighbor[a][i].second = newW;
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
                            Neighbor[b][i].second = newW;
                            break;
                        }
                    }
                }
                tt.stop();
            }else if(algoIndex==1){//CH update
                tt.start();
                CHdecBat(wBatchDec);
                tt.stop();
            } else if(algoIndex==2){//H2H update
                tt.start();
                H2HdecBat(wBatchDec);
                tt.stop();
            }
            runT=tt.GetRuntime();
            runT1 += runT;
            cout<<"Batch "<<i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
//                    throughputNum += EffiCheckThroughput(ODpair,tRecord,batchInterval);//query efficiency test
            throughputNum += EffiCheckThroughput(ODpair,runtimes,batchInterval,runT,queryT);
            if(ifDebug){
                CorrectnessCheck(100);
            }
        }*/

        //Step 2: Increase updates
        if(!wBatchInc.empty()){
            cout<<"Increase update."<<endl;
            tRecord.start();
            if(algoIndex==0){
                tt.start();
                for(int wb=0;wb<wBatchInc.size();wb++) {
                    int a = wBatchInc[wb].first.first;
                    int b = wBatchInc[wb].first.second;
                    int oldW = wBatchInc[wb].second.first;
                    int newW = wBatchInc[wb].second.second;

                    //modify the original graph information
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
                            Neighbor[a][i].second = newW;
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
                            Neighbor[b][i].second = newW;
                            break;
                        }
                    }
                }
                tt.stop();
            }else if(algoIndex==1){//CH update
                tt.start();
                CHincBatMT(wBatchInc);
                tt.stop();
            } else if(algoIndex==2){//H2H update
                tt.start();
                H2HincBatMT(wBatchInc);
                tt.stop();
            }
            runT=tt.GetRuntime();
            runT2 += runT;
            cout<<"Batch "<<i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
            throughputNum += EffiCheckThroughput(ODpair,runtimes,batchInterval,runT,queryT);
            if(ifDebug){
                CorrectnessCheck(100);
            }
        }



    }
    cout<<"\nOverall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average batch update Time: "<<runT2/batchNum<<" s; Average query time: "<<1000*queryT/batchNum<<" ms."<<endl;

}
//function of testing the throughput on real-life updates
void Graph::RandomUpdateThroughputTest(string updateFile, int batchNum, int batchSize, int batchInterval) {
    bool ifDebug=false;
//    ifDebug=true;
    int runtimes=10000;
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    string line;
    vector<string> vs;
    int ID1,ID2,oldW,newW,weight;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==1);
    int eNum=stoi(vs[0]);
    if(batchNum*batchSize>eNum){
        batchNum=eNum/batchSize;
        cout<<"Actual batch number: "<<batchNum<<endl;
    }
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNum);
    for(int i=0;i<batchNum;++i){
        for(int j=0;j<batchSize;++j){
            getline(IF,line);
            vs.clear();
            boost::split(vs,line,boost::is_any_of(" "));
            ID1=stoi(vs[0]), ID2=stoi(vs[1]), weight=stoi(vs[2]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    IF.close();
    cout<<"Update batch number: "<<batchNum<<" ; batch size: "<<batchSize<<" ; Batch interval: "<< batchInterval<<endl;
    string queryF = sourcePath+dataset + ".query";

    ifstream IF2(queryF);
    if(!IF2){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF2>>num;
    for(int k=0;k<num;k++){
        IF2>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF2.close();

    //index maintenance
    Timer tRecord;
    double runT1=0, runT2 = 0;
    unsigned long long throughputNum=0;
    double runT=0;
    double queryT=0;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;
    Timer tt;
//    map<pair<int,int>,int> uEdges;
    for(int i=0;i<batchNum;++i){
        wBatchDec.clear(); wBatchInc.clear(); runT=0;
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[i].size();++j){
            ID1=batchUpdates[i][j].first.first, ID2=batchUpdates[i][j].first.second, weight=batchUpdates[i][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
//                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            if(j<batchSize/2){//decrease update
                for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                    if(it->first==ID2){
                        ifFind=true;
                        oldW=it->second;
//                        weight=(0.5+0.5*rand()/(RAND_MAX+1.0))*oldW;
//                        cout<<weight<<endl;
                        weight=0.5*oldW;
                        if(weight>0 && weight<oldW){
//                        cout<<"Dec "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                            wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                        }

                        break;
                    }
                }
            }
            else{//increase update
                for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                    if(it->first==ID2){
                        ifFind=true;
                        oldW=it->second;
//                        weight=(1+1*rand()/(RAND_MAX+1.0))*oldW;
                        weight=1.5*oldW;
                        if(weight>0 && weight>oldW) {
//                        cout<<"Inc "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                            wBatchInc.emplace_back(make_pair(ID1, ID2), make_pair(oldW, weight));
                        }
                        break;
                    }
                }
            }

            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),weight});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<weight<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<i<<" . Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;

        //Step 1: Decrease updates
        if(!wBatchDec.empty()){
            cout<<"Decrease update. "<<wBatchDec.size()<<endl;
            tRecord.start();
            if(algoIndex==0){//Dijkstra
                tt.start();
                int a,b;
                for(int k=0;k<wBatchDec.size();k++) {
                    a = wBatchDec[k].first.first;
                    b = wBatchDec[k].first.second;
                    newW = wBatchDec[k].second.second;

                    //modify the information in original graph
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
                            Neighbor[a][i].second = newW;
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
                            Neighbor[b][i].second = newW;
                            break;
                        }
                    }
                }
                tt.stop();
            }else if(algoIndex==1){//CH update
                tt.start();
                CHdecBat(wBatchDec);
                tt.stop();
            } else if(algoIndex==2){//H2H update
                tt.start();
                H2HdecBat(wBatchDec);
                tt.stop();
            }
            runT+=tt.GetRuntime();
//            runT1 += runT;
//            cout<<"Batch "<<i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
////                    throughputNum += EffiCheckThroughput(ODpair,tRecord,batchInterval);//query efficiency test
//            throughputNum += EffiCheckThroughput(ODpair,runtimes,batchInterval,runT,queryT);
        }

        //Step 2: Increase updates
        if(!wBatchInc.empty()){
            cout<<"Increase update. "<<wBatchInc.size()<<endl;
            tRecord.start();
            if(algoIndex==0){
                tt.start();
                for(int wb=0;wb<wBatchInc.size();wb++) {
                    int a = wBatchInc[wb].first.first;
                    int b = wBatchInc[wb].first.second;
                    int oldW = wBatchInc[wb].second.first;
                    int newW = wBatchInc[wb].second.second;

                    //modify the original graph information
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
                            Neighbor[a][i].second = newW;
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
                            Neighbor[b][i].second = newW;
                            break;
                        }
                    }
                }
                tt.stop();
            }else if(algoIndex==1){//CH update
                tt.start();
                CHincBatMT(wBatchInc);
                tt.stop();
            } else if(algoIndex==2){//H2H update
                tt.start();
                H2HincBatMT(wBatchInc);
                tt.stop();
            }
            runT+=tt.GetRuntime();
//            runT2 += runT;
//            cout<<"Batch "<<i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
        }
        cout<<"Batch "<<i<<". Update time: "<<runT<<" s."<<endl;
        runT2 += runT;
        if(ifDebug){
            CorrectnessCheck(100);
        }
        throughputNum += EffiCheckThroughput(ODpair,runtimes/batchNum,batchInterval,runT,queryT);
    }
    cout<<"\nOverall throughput: "<<throughputNum<<" ; Average throughput: "<<throughputNum/batchNum<<" ; Average batch update Time: "<<runT2/batchNum<<" s; Average query time: "<<1000*queryT/batchNum<<" ms."<<endl;

}
//function of testing the throughput on random updates
void Graph::RandomUpdateThroughputTestQueueModel(int batchNum, int batchSize, int batchInterval, double T_r, int workerNum) {//T_r is in second
    bool ifDebug=false;
//    ifDebug=true;
    int runtimes=10000;
//    runtimes=1000;
    string updateFile=sourcePath+dataset+".update";
    ifstream IF(updateFile);
    if(!IF){
        cout<<"Cannot open file "<<updateFile<<endl;
        exit(1);
    }
    cout<<"update file: "<<updateFile<<endl;
    string line;
    vector<string> vs;
    int ID1,ID2,oldW,newW,weight;
    getline(IF,line);
    vs.clear();
    boost::split(vs,line,boost::is_any_of(" "));
    assert(vs.size()==1);
    int eNum=stoi(vs[0]);
    if(batchNum*batchSize>eNum){
        batchNum=eNum/batchSize;
        cout<<"Actual batch number: "<<batchNum<<endl;
    }
    vector<vector<pair<pair<int,int>,int>>> batchUpdates(batchNum);
    for(int i=0;i<batchNum;++i){
        for(int j=0;j<batchSize;++j){
            getline(IF,line);
            vs.clear();
            boost::split(vs,line,boost::is_any_of(" "));
            ID1=stoi(vs[0]), ID2=stoi(vs[1]), weight=stoi(vs[2]);
            batchUpdates[i].emplace_back(make_pair(ID1,ID2),weight);
        }
    }
    IF.close();
    cout<<"Update batch number: "<<batchNum<<" ; batch size: "<<batchSize<<" ; Batch interval: "<< batchInterval<<endl;

    string queryF = sourcePath+dataset + ".query";

    ifstream IF2(queryF);
    if(!IF2){
        cout<<"Cannot open file "<<queryF<<endl;
        exit(1);
    }
    cout<<"Query file: "<<queryF<<endl;
    int num;
    vector<pair<int,int>> ODpair;
    IF2>>num;
    for(int k=0;k<num;k++){
        IF2>>ID1>>ID2;
        ODpair.emplace_back(ID1, ID2);
    }
    IF2.close();

    //index maintenance
    Timer tt;
    double runT1=0, runT2 = 0;
    double throughputNum=0;
    double queryT=0;
    double updateTime=0;
    vector<vector<double>> stageUpdateT(2,vector<double>());//(stage, update time)
    vector<vector<double>> stageQueryT(2,vector<double>());//(stage, query times)
    vector<pair<pair<int,int>,pair<int,int>>> wBatchDec;
    vector<pair<pair<int,int>,pair<int,int>>> wBatchInc;

    double queryTime=EffiStageCheck(ODpair, runtimes, stageQueryT);

    for(int i=0;i<batchNum;++i){
        wBatchDec.clear(); wBatchInc.clear();
        map<pair<int,int>,int> uEdges;
        for(int j=0;j<batchUpdates[i].size();++j){
            ID1=batchUpdates[i][j].first.first, ID2=batchUpdates[i][j].first.second, weight=batchUpdates[i][j].second;
            bool ifFind=false;
            if(ID1>ID2){
                int temp=ID1;
                ID1=ID2, ID2=temp;
//                cout<<"ID2 is smaller!"<<ID1<<" "<<ID2<<endl;
            }
            if(j<batchSize/2){//decrease update
                for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                    if(it->first==ID2){
                        ifFind=true;
                        oldW=it->second;
//                        weight=(0.5+0.5*rand()/(RAND_MAX+1.0))*oldW;
//                        cout<<weight<<endl;
                        weight=oldW/2;
                        if(weight>0 && weight<oldW){
//                        cout<<"Dec "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                            wBatchDec.emplace_back(make_pair(ID1,ID2), make_pair(oldW,weight));
                        }

                        break;
                    }
                }
            }
            else{//increase update
                for(auto it=Neighbor[ID1].begin();it!=Neighbor[ID1].end();++it){
                    if(it->first==ID2){
                        ifFind=true;
                        oldW=it->second;
//                        weight=(1+1*rand()/(RAND_MAX+1.0))*oldW;
                        weight=2*oldW;
                        if(weight>0 && weight>oldW) {
//                        cout<<"Inc "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<weight<<endl;
                            wBatchInc.emplace_back(make_pair(ID1, ID2), make_pair(oldW, weight));
                        }
                        break;
                    }
                }
            }

            if(uEdges.find(make_pair(ID1,ID2))==uEdges.end()){//if not found
                uEdges.insert({make_pair(ID1,ID2),weight});
            }else{
                cout<<"Wrong. Find. "<<ID1<<" "<<ID2<<" "<<weight<<" "<<uEdges[make_pair(ID1,ID2)]<<" "<<oldW <<endl;
                exit(1);
            }

            if(!ifFind){
                cout<<"Wrong edge update. "<<ID1<<" "<<ID2<<" "<<endl; exit(1);
            }
        }
        cout<<"Batch "<<i<<" . Decrease update number: "<<wBatchDec.size()<<" ; Increase update number: "<<wBatchInc.size()<<endl;
        updateTime=0;
        //Step 1: Decrease updates
        if(!wBatchDec.empty()){
            cout<<"Decrease update. "<<wBatchDec.size()<<endl;
            if(algoIndex==0 || algoIndex==3){//Dijkstra
                tt.start();
                int a,b;
                for(int k=0;k<wBatchDec.size();k++) {
                    a = wBatchDec[k].first.first;
                    b = wBatchDec[k].first.second;
                    newW = wBatchDec[k].second.second;

                    //modify the information in original graph
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
                            Neighbor[a][i].second = newW;
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
                            Neighbor[b][i].second = newW;
                            break;
                        }
                    }
                }
                tt.stop();
            }else if(algoIndex==1){//CH update
                tt.start();
                CHdecBat(wBatchDec);
                tt.stop();
            } else if(algoIndex==2){//H2H update
                tt.start();
                H2HdecBat(wBatchDec);
                tt.stop();
            }
            updateTime+=tt.GetRuntime();
//            runT1 += runT;
//            cout<<"Batch "<<i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
////                    throughputNum += EffiCheckThroughput(ODpair,tRecord,batchInterval);//query efficiency test
//            throughputNum += EffiCheckThroughput(ODpair,runtimes,batchInterval,runT,queryT);
        }

        //Step 2: Increase updates
        if(!wBatchInc.empty()){
            cout<<"Increase update. "<<wBatchInc.size()<<endl;
            if(algoIndex==0 || algoIndex==3){
                tt.start();
                for(int wb=0;wb<wBatchInc.size();wb++) {
                    int a = wBatchInc[wb].first.first;
                    int b = wBatchInc[wb].first.second;
                    int oldW = wBatchInc[wb].second.first;
                    int newW = wBatchInc[wb].second.second;

                    //modify the original graph information
                    for (int i = 0; i < Neighbor[a].size(); i++) {
                        if (Neighbor[a][i].first == b) {
                            Neighbor[a][i].second = newW;
                            break;
                        }
                    }
                    for (int i = 0; i < Neighbor[b].size(); i++) {
                        if (Neighbor[b][i].first == a) {
                            Neighbor[b][i].second = newW;
                            break;
                        }
                    }
                }
                tt.stop();
            }else if(algoIndex==1){//CH update
                tt.start();
                CHincBatMT(wBatchInc);
                tt.stop();
            } else if(algoIndex==2){//H2H update
                tt.start();
                H2HincBatMT(wBatchInc);
                tt.stop();
            }
            updateTime+=tt.GetRuntime();
//            runT2 += runT;
//            cout<<"Batch "<<i<<". Update time: "<<tt.GetRuntime()<<" s."<<endl;
        }

        cout<<"Batch "<<i<<" . Update time: "<<updateTime<<" s."<<endl;

        if(ifDebug){
            CorrectnessCheck(100);
        }

//        for(int i=0;i<stageUpdateT.size();++i){
//            stageUpdateT[i].push_back(stageDurations[i]);
//        }


//        EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
        vector<double> duration = StageDurationCompute(batchInterval,updateTime);
        for(int j=0;j<duration.size();++j){
            stageUpdateT[j].push_back(duration[j]);
        }

    }
//    for(int i=0;i<stageDurations.size();++i){
//        stageDurations[i]/=batchNum;
//    }
//    EffiCheckStages(ODpair,runtimes,batchInterval,throughputNum,stageUpdateT,stageQueryT);
//    AverageStagePerformance(batchNum,stageUpdateT,stageQueryT);

//    double throughputNum2 = ThroughputEstimate(stageQueryT, stageUpdateT, T_r, batchInterval);
//    cout<<"Estimated throughput: "<<throughputNum2<<endl;
//    throughputNum = ThroughputSimulate(stageQueryT,stageUpdateT,batchInterval*batchNum,T_r,batchInterval,queryTime);

    pair<double,double> resultE = ThroughputEstimate(stageQueryT, stageUpdateT, T_r, batchInterval);
    cout<<"Estimated throughput: "<<resultE.first<<" ; fastest available query time: "<<resultE.second*1000<<" ms"<<endl;
//    throughputNum = ThroughputSimulate(stageQueryT,stageUpdateT,batchInterval*batchNum,T_r,batchInterval,resultE.second);
    throughputNum = ThroughputSimulate(stageQueryT,stageUpdateT,batchInterval*batchNum,T_r,batchInterval,resultE.second,workerNum);
//        vector<vector<double>> & simulated_query_costs, vector<double>& simulated_update_costs, int simulation_time, double threshold_time, double period_time, vector<vector<int> >& query_costs, vector<int>& batch_update_costs
    cout<<"\nPartiNum: "<<partiNum<<". Throughput: "<<throughputNum<<" ; Average batch update Time: "<<updateTime/batchNum<<" s."<<endl;
}

//function for efficiency test
double Graph::EffiStageCheck(vector<pair<int,int>> & ODpair, int runtimes, vector<vector<double>> & queryTimes){//return in seconds
    cout<<"Efficiency test. Run times: "<<runtimes<<endl;
    int s, t;
    Timer tt;
    double runT=0, runT1=0;
    int d1,d2;

    vector<int> results(runtimes,-1);
    /// For BiDijkstra
    for(int i=0;i<runtimes/100;i++){
        s=ODpair[i].first; t=ODpair[i].second;
//        d1=Dijkstra(s,t,Neighbor);
        tt.start();
        d2=BiDijkstra(s,t,Neighbor);
        tt.stop();
        runT1+=tt.GetRuntime();
        results[i]=d2;
        queryTimes[0].push_back(tt.GetRuntime());
    }
    runT1=runT1/(runtimes/100);
    cout<<"Average Efficiency of Q-Stage 1 : "<<1000*runT1<<" ms."<<endl;
    /// For SP index
    if(algoIndex==1) {
        for (int i = 0; i < runtimes / 10; i++) {
            s = ODpair[i].first;
            t = ODpair[i].second;
            tt.start();
            d2 = QueryCHWP(s, t);
            tt.stop();
            runT += tt.GetRuntime();
            results[i] = d2;
            queryTimes[1].push_back(tt.GetRuntime());
        }
        runT=runT/(runtimes/10);
        cout<<"Average Efficiency of Q-Stage 2 : "<<1000*runT<<" ms."<<endl;
    }
    else if(algoIndex==2){
        for (int i = 0; i < runtimes; i++) {
            s = ODpair[i].first;
            t = ODpair[i].second;
            tt.start();
            d2 = QueryH2H(s, t);
            tt.stop();
            runT += tt.GetRuntime();
            results[i] = d2;
            queryTimes[1].push_back(tt.GetRuntime());
        }
        runT=runT/runtimes;
        cout<<"Average Efficiency of Q-Stage 2 : "<<1000*runT<<" ms."<<endl;
    }else if(algoIndex==3 || algoIndex==0){
        queryTimes[1]=queryTimes[0];
        runT=runT1;
        cout<<"Average Efficiency of Q-Stage 2 : "<<1000*runT<<" ms."<<endl;
    }

    tt.stop();
    cout<<"Time for efficiency test: "<<tt.GetRuntime()<<" s."<<endl;
    return runT;
}

vector<double> Graph::StageDurationCompute(int intervalT, double updateTime){
//    double dt1=0,dt2=0,dt3=0,dt4=0,dt5=0;//actual duration for query stages
    vector<double> durations(2,0);
    double updateT=0;//overall update time

    //Stage 1: Dijkstra
    if(updateTime<intervalT){
        durations[0]=updateTime;
        //Stage 2: Index
        durations[1]=intervalT-updateTime;

    }else{
        durations[0]=intervalT;
    }

    return durations;
}

double Graph::analytical_update_first(double T_q, double T_u, double T_r, double T, double V_q) {//T is duration, T_u is the average update time for one update, T_r is in seconds, return throughput
//	M/M/1 return min(1.0 / T_q - 1.0 / T_r, 1.0 / T_q - tau * T_u / (T * T_q));
//	M/G/1
    if (T_r < T_q)
        return 0;
    double a = 2 * (T_r - T_q) / (V_q + 2 * T_r * T_q - T_q * T_q);
    double lambda_u = T_u / T;
    double b = (1 - lambda_u) / T_q;
    double c = a < b ? a : b;
    if (c < 0) return 0;
    return c;
}

// function for estimating the system throughput, old
/*double Graph::ThroughputEstimate(vector<vector<double>> &query_costs, vector<vector<double>> &update_costs, double threshold_time, double T) {
    double queryT, duration, queryT_var;
    double update_time;
    double throughput=0.0;

    for(int i=0;i<query_costs.size();++i){
        duration = get_mean(update_costs[i]);
        if(duration<=0)
            continue;
        queryT = get_mean(query_costs[i]);
        queryT_var = get_var(query_costs[i]);
        double tau = duration/T;
        double thr = tau*analytical_update_first(queryT, 0, threshold_time/1000,  duration, queryT_var);
        cout<<"Throughput of Q-stage "<<i+1<<" : "<<thr<<endl;
        throughput+=thr;
    }
    return throughput;
}*/

// function for estimating the system throughput
pair<double,double> Graph::ThroughputEstimate(vector<vector<double>> &query_costs, vector<vector<double>> &update_costs, double threshold_time, double T) {
    double queryT, duration, queryT_var;
    double update_time=0;
    double durationT=0;
    double throughput=0.0;
    double tau, thr;
    double query_time=INF;
//    vector<double> updateT;
//    updateT.push_back(0);

    for(int i=0;i<query_costs.size();++i){
        duration = get_mean(update_costs[i]);
//        updateT.push_back(update_time);
        if(duration<=0)
            continue;
        durationT+=duration;
        queryT = get_mean(query_costs[i]);
        queryT_var = get_var(query_costs[i]);
        tau = duration/T;
        if(queryT<query_time){
            query_time=queryT;
        }

        thr = analytical_update_first(queryT, update_time, threshold_time,  durationT, queryT_var);
//        cout<<i<<" . update time: "<<update_time<<" s ; estimate duration: "<<durationT<<" s"<<endl;
        update_time+=duration;
        cout<<"Throughput of Q-stage "<<i+1<<" : "<<thr<<" ; duration: "<< duration<<" s ; query time: "<<queryT*1000<<" ms"<<endl;
        throughput+=thr*tau;
    }
//    cout<<"Weighted average of throughput: "<<throughput<<endl;
    return make_pair(throughput,query_time);
}

// function for simulating the system throughput, old
//double Graph::ThroughputSimulate(vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime){//period_time is update interval time, simulation_time is the overall time for simulation
//    Timer tt;
//    tt.start();
//    double l = 0, r = 5000000; // be careful of 'r'
//    double throughput = 0.0;
//    int lambda;
//    int gap=10;
////    double qTime=INF;
////    for(int i=0;i<query_costs.size();++i){
////        if(qTime>query_costs[i][0] && query_costs[i][0]>0){
////            qTime=query_costs[i][0];
////        }
////    }
//    if(r>100/queryTime){
//        r=100/queryTime;
//        gap=max((int)r/1000,1);
//        gap=min(gap,10);
//    }
//
//    cout<<"Simulate the throughput of system... r: "<<r<<" ; gap: "<<gap<<endl;
//    while(l <= r){//why use different lambda? to get the suitable lambda
//        lambda = (l+r)/2;
//        auto queryList = generate_queryList(lambda, simulation_time);
//        cout<<"lambda: "<<lambda<<" ; the size of querylist: "<<queryList.size()<<endl;
//
//        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time);
//
////        pair<double,double> result = simulator_UpdateFirst(simulated_query_costs, simulated_update_costs, queryList, simulation_time, period_time);
//        if(result.second*1000 <= threshold_time){
//            if(throughput < result.first) {
////                cout<<"here:" << result.first<<endl;
//                throughput=result.first;
//
//                cout<<"lambda: "<<lambda<<" ; throughput: "<<result.first<< " ; response query time: "<< result.second*1000 <<" ms"<<endl;
//            }
//            l = lambda + gap;
////            l=l*1.1;
//        }else{
//            r = lambda - gap;
////            r=r/1.1;
//        }
//
//    }
//
//    tt.stop();
//    cout<<"Time for throughput simulation: "<<tt.GetRuntime()<<" s"<<endl;
//    return throughput;
//}

// function for simulating the system throughput
double Graph::ThroughputSimulate(vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime){//period_time is update interval time, simulation_time is the overall time for simulation
    Timer tt;
    tt.start();
    double l = 0, r = 2000000; // be careful of 'r'
    double throughput = 0.0;
    int lambda;
    int gap=1;

    if(r>10/queryTime){
        r=10/queryTime;
    }
    int num1=0;
    cout<<"Simulate the throughput of system... r: "<<r<<" ; gap: "<<gap<<endl;
    while(l <= r){//why use different lambda?
        lambda = (l+r)/2;
        if(lambda==0){
            lambda=1; num1++;
            if(num1>20){
                break;
            }
        }
        auto queryList = generate_queryList(lambda, simulation_time);
        cout<<"lambda: "<<lambda<<" ; the size of querylist: "<<queryList.size();

        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time);
//        pair<double,double> result = simulator_UpdateFirst(simulated_query_costs, simulated_update_costs, queryList, simulation_time, period_time);
        if(result.second <= threshold_time){
            if(throughput < result.first) {
//                cout<<"here:" << result.first<<endl;
                throughput=result.first;

                cout<<" ; throughput: "<<result.first<< " ; response query time: "<< result.second*1000 <<" ms";
            }
            l = lambda + gap;
//            l=l*1.1;
        }else{
            r = lambda - gap;
//            r=r/1.1;
        }
        cout<<endl;
    }

    tt.stop();
    cout<<"Time for throughput simulation: "<<tt.GetRuntime()<<" s"<<endl;
    return throughput;

}

double Graph::ThroughputSimulate(vector<vector<double>> & query_costs, vector<vector<double>>& update_costs, int simulation_time, double threshold_time, double period_time, double queryTime, int workerNum){//period_time is update interval time, simulation_time is the overall time for simulation
    Timer tt;
    tt.start();
    double l = 0, r = 2000000; // be careful of 'r'
    double throughput = 0.0;
    int lambda;
    int gap=1;

    r*=workerNum;

    if(r>workerNum*10/queryTime){
        r=workerNum*10/queryTime;
    }

    cout<<"Simulate the throughput of system... r: "<<r<<" ; gap: "<<gap<<" ; worker number: "<< workerNum<<endl;
    while(l < r){//why use different lambda?
        lambda = (l+r)/2;
        if(lambda==0){
            lambda=1;
        }
        auto queryList = generate_queryList(lambda, simulation_time);
        cout<<"lambda: "<<lambda<<" ; the size of querylist: "<<queryList.size();

//        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time);
        pair<double,double> result = simulator_UpdateFirst(query_costs, update_costs, queryList, simulation_time, period_time, workerNum);

//        pair<double,double> result = simulator_UpdateFirst(simulated_query_costs, simulated_update_costs, queryList, simulation_time, period_time);
        if(result.second <= threshold_time){
            if(throughput < result.first) {
//                cout<<"here:" << result.first<<endl;
                throughput=result.first;

                cout<<" ; throughput: "<<result.first<< " ; response query time: "<< result.second*1000 <<" ms";
            }
            l = lambda + gap;
//            l=l*1.1;
        }else{
            r = lambda - gap;
//            r=r/1.1;
        }
        cout<<endl;
    }

    tt.stop();
    cout<<"Time for throughput simulation: "<<tt.GetRuntime()<<" s"<<endl;
    return throughput;

}

//update first model
pair<double, double> Graph::simulator_UpdateFirst(vector<vector<double>> & query_cost, vector<vector<double>> & update_cost, vector<query> &queryList, int T, double period_time){
    int total_time; // terminal time for each period
    unsigned long long int count = 0; // the number of queries processed.

    double avg_response_time = 0.0;
    query* current_query;

    double c_time=0;  // current time for processing queries

//    T = T*microsecs_per_sec; // transform T to microseconds
    int num_updates = T/period_time;//number of batch updates
//    cout<<"batch number: "<<num_updates<<endl;
    int current_update_index=0;


    bool finished = false;

    int query_processed_in_one_update_slot = 0;

    double update_time_rest;


    while(current_update_index * period_time < T){

//        cout<<"batch "<<current_update_index<<endl;

        c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];
        double durationPoint=c_time;
        // start from c_time, end at total_time
        total_time = min((int)period_time * (current_update_index+1), T);
        query_processed_in_one_update_slot = 0;
        for(int i=0;i<query_cost.size();++i){
            update_time_rest = update_cost[i][current_update_index % update_cost.size()];//update time
//                cout<<"Duration of Q-Stage "<<i<<" : "<<update_time_rest<< " s"<<endl;
            c_time=durationPoint;
            if(update_time_rest==0)
                continue;

            durationPoint+=update_time_rest;
            while(c_time <= durationPoint) {//if there is remained time for querying
                if( count < queryList.size()) {
                    current_query = &queryList[count];
                }else {//if all queries are processed
                    finished = true;
//                    cout<<"!!! Finished the processing of all queries."<<endl;
                    break;
                }
//            cout<<"query "<<count<<": "<<current_query->init_time<<" "<<current_query->process_time<<" ";
                current_query->process_time = query_cost[i][query_processed_in_one_update_slot % query_cost[i].size()];//obtain the true query processing time
//            cout<<current_query->process_time<<endl;

                if (current_query->init_time> c_time){
                    c_time = current_query->init_time;
                }

                //c_time = max(c_time, current_query->init_time*1000000.0);
                if(c_time + current_query->process_time <= durationPoint) {
                    count++;
                    query_processed_in_one_update_slot++;
                    avg_response_time +=  (c_time -  current_query->init_time + current_query->process_time);//obtain the query response time
//                simulated_query_costs.push_back(current_query->process_time);
                    c_time += current_query->process_time;
                } else {//if the remained time is not enough for query processing
                    current_query->process_time = current_query->process_time - (durationPoint - c_time);
                    break;
                }
            }

        }
        if(finished) break;
//        if(c_time+1<total_time){
//            cout<<"Seems wrong. "<<c_time<<" "<<total_time<<" s"<<endl; exit(1);
//        }
        current_update_index++;
    }


    //cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
    return make_pair(count * 1.0 / T, avg_response_time / count);//return throughput (query per second) and average response time (in us)
}


pair<double, double> Graph::simulator_UpdateFirst(vector<vector<double>> & query_cost, vector<vector<double>> & update_cost, vector<query> &queryList, int T, double period_time, int workerNum){
    int total_time; // terminal time for each period
    unsigned long long int count = 0; // the number of queries processed.

    double avg_response_time = 0.0;
    query* current_query;

    double c_time=0;  // current time for processing queries

//    T = T*microsecs_per_sec; // transform T to microseconds
    int num_updates = T/period_time;//number of batch updates
//    cout<<"batch number: "<<num_updates<<endl;
    int current_update_index=0;


    bool finished = false;

    int query_processed_in_one_update_slot = 0;

    double update_time_rest;

    if(workerNum>1){//with multiple workers
        benchmark::heap<2, int, long long int> pqueue(workerNum);
        int topID; long long int topValue;
        vector<double> workers(workerNum,0.0);
        vector<double> c_times(workerNum,0.0);
        pair<double, int> minT;
        while(current_update_index * period_time < T){
//        cout<<"batch "<<current_update_index<<endl;
            c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];
            double durationPoint=c_time;
            // start from c_time, end at total_time
            total_time = min((int)period_time * (current_update_index+1), T);
            query_processed_in_one_update_slot = 0;
            for(int i=0;i<query_cost.size();++i){
                update_time_rest = update_cost[i][current_update_index % update_cost.size()];//update time
//                cout<<"Duration of Q-Stage "<<i<<" : "<<update_time_rest<< " s"<<endl;
                c_time=durationPoint;
                for(int wi=0;wi<c_times.size();++wi){
                    c_times[wi]=c_time;
                    pqueue.update(wi,c_time*Resolution);
                }
                if(update_time_rest==0)
                    continue;

                durationPoint+=update_time_rest;
                while(c_time <= durationPoint && !pqueue.empty()) {//if there is remained time for querying

                    pqueue.extract_min(topID, topValue);
                    c_time=c_times[topID];
                    if( count < queryList.size()) {
                        current_query = &queryList[count];
                    }else {//if all queries are processed
                        finished = true;
//                    cout<<"!!! Finished the processing of all queries."<<endl;
                        break;
                    }

                    current_query->process_time = query_cost[i][query_processed_in_one_update_slot % query_cost[i].size()];//obtain the true query processing time
                    if (current_query->init_time> c_times[topID]){
                        c_times[topID] = current_query->init_time;
                    }

                    //c_time = max(c_time, current_query->init_time*1000000.0);
                    if(c_times[topID] + current_query->process_time <= durationPoint) {
                        count++;
                        query_processed_in_one_update_slot++;
                        avg_response_time +=  (c_times[topID] -  current_query->init_time + current_query->process_time);//obtain the query response time
//                simulated_query_costs.push_back(current_query->process_time);
                        c_times[topID] += current_query->process_time;
                        pqueue.update(topID,c_times[topID]*Resolution);
                    } else {//if the remained time is not enough for query processing
                        current_query->process_time = current_query->process_time - (durationPoint - c_time);
                        c_times[topID] += current_query->process_time;
                        break;
                    }
                    if(!pqueue.empty()){
                        c_time=c_times[pqueue.top_id()];
                    }else{
                        break;
                    }


                }

            }
            if(finished) break;
//        if(c_time+1<total_time){
//            cout<<"Seems wrong. "<<c_time<<" "<<total_time<<" s"<<endl; exit(1);
//        }
            current_update_index++;
        }
    }
    else{//single worker
        while(current_update_index * period_time < T){
//        cout<<"batch "<<current_update_index<<endl;
            c_time = period_time * current_update_index;  // + batch_update_costs[current_update_index%batch_update_costs.size()];
            double durationPoint=c_time;
            // start from c_time, end at total_time
            total_time = min((int)period_time * (current_update_index+1), T);
            query_processed_in_one_update_slot = 0;
            for(int i=0;i<query_cost.size();++i){
                update_time_rest = update_cost[i][current_update_index % update_cost.size()];//update time
//                cout<<"Duration of Q-Stage "<<i<<" : "<<update_time_rest<< " s"<<endl;
                c_time=durationPoint;
                if(update_time_rest==0)
                    continue;

                durationPoint+=update_time_rest;
                while(c_time <= durationPoint) {//if there is remained time for querying
                    if( count < queryList.size()) {
                        current_query = &queryList[count];
                    }else {//if all queries are processed
                        finished = true;
//                    cout<<"!!! Finished the processing of all queries."<<endl;
                        break;
                    }
//            cout<<"query "<<count<<": "<<current_query->init_time<<" "<<current_query->process_time<<" ";
                    current_query->process_time = query_cost[i][query_processed_in_one_update_slot % query_cost[i].size()];//obtain the true query processing time
//            cout<<current_query->process_time<<endl;

                    if (current_query->init_time> c_time){
                        c_time = current_query->init_time;
                    }

                    //c_time = max(c_time, current_query->init_time*1000000.0);
                    if(c_time + current_query->process_time <= durationPoint) {
                        count++;
                        query_processed_in_one_update_slot++;
                        avg_response_time +=  (c_time -  current_query->init_time + current_query->process_time);//obtain the query response time
//                simulated_query_costs.push_back(current_query->process_time);
                        c_time += current_query->process_time;
                    } else {//if the remained time is not enough for query processing
                        current_query->process_time = current_query->process_time - (durationPoint - c_time);
                        break;
                    }
                }

            }
            if(finished) break;
//        if(c_time+1<total_time){
//            cout<<"Seems wrong. "<<c_time<<" "<<total_time<<" s"<<endl; exit(1);
//        }
            current_update_index++;
        }
    }




    //cout<<"Avg_response_time " << avg_response_time / count / 1000000 <<endl;
    return make_pair(count * 1.0 / T, avg_response_time / count);//return throughput (query per second) and average response time (in us)
}

/// CHWP algorithm
void Graph::IndexMaintenanceCHWP(int updateType, int updateSize, bool ifBatch, int batchSize) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = sourcePath+dataset + ".update";
    bool ifDebug=false;
    ifDebug=true;

    int ID1, ID2, oldW, newW;
    srand (0);
    int updateBatch=1;
    updateBatch=max(updateBatch,updateSize/batchSize);
    cout<<"Update size: "<<updateSize<<". Update batch: "<<updateBatch<<" ; Batch size: "<<batchSize<<endl;
    vector<pair<pair<int,int>,int>> updateData;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    ReadUpdate(file, updateData);
    Timer tt;
    double runT1=0, runT2 = 0;
    switch (updateType) {
        case 1:{
            //Decrease update
            cout<<"Update type: Decrease"<<endl;

            if(ifBatch){//for batch update
                if(updateBatch*batchSize>updateData.size()){
                    updateBatch=floor(updateData.size()/batchSize);
                }
                int update_i=0;

//                Graph g2=*this;
                for(int u=0;u<updateBatch;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
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
                        ++update_i;
                    }

                    tt.start();
                    CHdecBat(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }
                }
                cout<<"Average Decrease batch update Time: "<<runT1/(updateBatch*batchSize)<<" s.\n"<<endl;

            }

            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
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
                        newW=oldW*1.5;
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
//                        ++update_i;
                    }

                    tt.start();
                    CHincBatMT(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheck(100);
                    }

                }
                cout<<"Average Increase batch update Time: "<<runT2/(updateBatch*batchSize)<<" s.\n"<<endl;
            }

            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}


void Graph::CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    //maintain the index caused by the weight change
    //NodeOrders.clear();
    NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
    set<OrderCompCH> OC;
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the fresh distance and avoid search in the adjacent list
    //OC.clear(); OCdis.clear();

    int a,b,newW;//the weight of (a,b) decrease to newW
    for(int k=0;k<wBatch.size();k++){
        a=wBatch[k].first.first;
        b=wBatch[k].first.second;
        newW=wBatch[k].second.second;

        //modify the information in original graph
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

        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborCon[a].size();i++){
                if(NeighborCon[a][i].first==b){
                    if(NeighborCon[a][i].second.first>newW){
                        //cout<<OutNeighborCon[a][i].second.first<<"..........."<<newW<<endl;
                        NeighborCon[a][i].second.first=newW;
                        NeighborCon[a][i].second.second=1;

                        OCdis[make_pair(a,b)]=newW;
                        OC.insert(OrderCompCH(a,b));
                    }else if(NeighborCon[a][i].second.first==newW)
                        NeighborCon[a][i].second.second+=1;
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborCon[b].size();i++){
                if(NeighborCon[b][i].first==a){
                    if(NeighborCon[b][i].second.first>newW){
                        NeighborCon[b][i].second.first=newW;
                        NeighborCon[b][i].second.second=1;

                        OCdis[make_pair(b,a)]=newW;
                        OC.insert(OrderCompCH(b,a));
                    }else if(NeighborCon[b][i].second.first==newW)
                        NeighborCon[b][i].second.second+=1;
                    break;
                }
            }
        }
    }


    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];
        map<int,int> InM2t; //InM2t.clear();
        vector<pair<int,int>> InMLower; //InMLower.clear();
        for(int i=0;i<NeighborCon[s].size();i++){
            if(NodeOrder[NeighborCon[s][i].first]>NodeOrder[t])
                InM2t.insert(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
            else if(NodeOrder[NeighborCon[s][i].first]<NodeOrder[t])
                InMLower.push_back(make_pair(NeighborCon[s][i].first,NeighborCon[s][i].second.first));
        }
        int inID,inW,inWt;
        for(int i=0;i<NeighborCon[t].size();i++){
            inID=NeighborCon[t][i].first;
            if(InM2t.find(inID)!=InM2t.end()){
                inW=InM2t[inID];
                inWt=NeighborCon[t][i].second.first;
                if(inWt>inW+wt){
                    NeighborCon[t][i].second.first=inW+wt;
                    NeighborCon[t][i].second.second=1;
                    OCdis[make_pair(t,inID)]=inW+wt;
                    OrderCompCH oc={t,inID};
                    OC.insert(oc);
                }else if(inWt==inW+wt){
                    NeighborCon[t][i].second.second+=1;
                }
            }
        }

        for(int i=0;i<InMLower.size();i++){
            inID=InMLower[i].first; inW=InMLower[i].second;
            for(int j=0;j<NeighborCon[inID].size();j++){
                if(NeighborCon[inID][j].first==t){
                    inWt=NeighborCon[inID][j].second.first;
                    if(inWt>inW+wt){
                        NeighborCon[inID][j].second.first=inW+wt;
                        NeighborCon[inID][j].second.second=1;

                        OCdis[make_pair(inID,t)]=inW+wt;
                        OrderCompCH oc={inID,t};
                        OC.insert(oc);
                    }else if(inWt==inW+wt){
                        if(OCdis.find(make_pair(s,inID))==OCdis.end()) {//if not found
                            NeighborCon[inID][j].second.second += 1;
                        }
                    }
                    break;
                }
            }
        }
    }//finish change index
}

void Graph::CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
    //NodeOrders.clear();
    NodeOrders.assign(NodeOrder.begin(),NodeOrder.end());
    set<OrderCompCH> OC; //OC.clear();
    map<pair<int,int>,int> OCdis;//{(s,t),d} maintain the old distance before refreshed and avoid search in the adjacent list
    //OCdis.clear();

    for(int wb=0;wb<wBatch.size();wb++){
        int a=wBatch[wb].first.first;
        int b=wBatch[wb].first.second;
        int oldW=wBatch[wb].second.first;
        int newW=wBatch[wb].second.second;

        //modify the original graph information
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


        if(NodeOrder[a]<NodeOrder[b]){
            for(int i=0;i<NeighborCon[a].size();i++){
                if(NeighborCon[a][i].first==b){
                    if(NeighborCon[a][i].second.first==oldW){
                        NeighborCon[a][i].second.second-=1;
                        if(NeighborCon[a][i].second.second<1){
                            OrderCompCH oc={a,b};
                            OC.insert(oc);
                            OCdis[make_pair(a,b)]=oldW;
                        }
                    }
                    break;
                }
            }
        }else{
            for(int i=0;i<NeighborCon[b].size();i++){
                if(NeighborCon[b][i].first==a){
                    if(NeighborCon[b][i].second.first==oldW){
                        NeighborCon[b][i].second.second-=1;
                        if(NeighborCon[b][i].second.second<1){
                            OrderCompCH oc={b,a};
                            OC.insert(oc);
                            OCdis[make_pair(b,a)]=oldW;
                        }
                    }
                    break;
                }
            }
        }
    }

    while(!OC.empty()){
        int s=(*OC.begin()).x; int t=(*OC.begin()).y;
        int wt;
        OC.erase(OC.begin());
        wt=OCdis[make_pair(s,t)];//distance of s--->t before change
        int inID,inW;
        map<int,int> HigherIn; vector<pair<int,int>> LowerIn;
        //HigherIn.clear(); LowerIn.clear();
        //the shortcuts infected by s-->t
        for(int i=0;i<NeighborCon[s].size();i++){
            inID=NeighborCon[s][i].first;
            inW=NeighborCon[s][i].second.first;
            if(NodeOrder[inID]<NodeOrder[t]){
                LowerIn.emplace_back(inID,inW);
            }else if(NodeOrder[inID]>NodeOrder[t]){
                HigherIn.insert(make_pair(inID,inW));
            }
        }
        for(int i=0;i<NeighborCon[t].size();i++){
            inID=NeighborCon[t][i].first;
            if(HigherIn.find(inID)!=HigherIn.end()){
                inW=HigherIn[inID];
                if(NeighborCon[t][i].second.first==wt+inW){
                    NeighborCon[t][i].second.second-=1;
                    if(NeighborCon[t][i].second.second<1){
                        OrderCompCH oc={t,inID};
                        OC.insert(oc);
                        OCdis[make_pair(t,inID)]=wt+inW;
                    }
                }
            }
        }
        for(int i=0;i<LowerIn.size();i++){
            inID=LowerIn[i].first; inW=LowerIn[i].second;
            for(int j=0;j<NeighborCon[inID].size();j++){
                if(NeighborCon[inID][j].first==t){
                    if(NeighborCon[inID][j].second.first==inW+wt){
                        if(OCdis.find(make_pair(s,inID))==OCdis.end()) {//if not found, new
                            NeighborCon[inID][j].second.second -= 1;
                            if (NeighborCon[inID][j].second.second < 1) {
                                OrderCompCH oc = {inID, t};
                                OC.insert(oc);
                                OCdis[make_pair(inID, t)] = wt + inW;
                            }
                        }
                    }
                    break;
                }
            }
        }

        //get the new weight value of s-->t
        wt=INF; int countwt=0;
        for(int i=0;i<Neighbor[s].size();i++){
            if(Neighbor[s][i].first==t){
                wt=Neighbor[s][i].second;//the weight value in the original graph
                countwt=1;
                break;
            }
        }
        int ssw,wtt,wid;
        vector<int> Wnodes; //Wnodes.clear();
        if(s<t){
            //Wnodes=SCconNodes[make_pair(s,t)]; //cout<<"wid num "<<Wnodes.size()<<endl;
            Wnodes=SCconNodesMT[s][t];
        }else{
            //Wnodes=SCconNodes[make_pair(t,s)];
            Wnodes=SCconNodesMT[t][s];
        }

        for(int i=0;i<Wnodes.size();i++){
            wid=Wnodes[i];
            for(int j=0;j<NeighborCon[wid].size();j++){
                if(NeighborCon[wid][j].first==s){
                    ssw=NeighborCon[wid][j].second.first;
                }
                if(NeighborCon[wid][j].first==t){
                    wtt=NeighborCon[wid][j].second.first;
                }
            }

            if(ssw+wtt<wt){
                wt=ssw+wtt;
                countwt=1;
            }else if(ssw+wtt==wt){
                countwt+=1;
            }
        }

        //refresh the weight value of s--t in the index
        for(int i=0;i<NeighborCon[s].size();i++){
            if(NeighborCon[s][i].first==t){
//                cout<<"Refresh shortcut: "<<s<<" "<<t<<" "<<NeighborCon[s][i].second.first<<" "<<wt<<endl;
                NeighborCon[s][i].second.first=wt;
                NeighborCon[s][i].second.second=countwt;
                break;
            }
        }
    }
}

void Graph::IndexMaintenanceH2H(int updateType, int updateSize, bool ifBatch, int batchSize) {
    cout<<"Index update test..."<<endl;
    // read updates
    string file = sourcePath+dataset + ".update";
    bool ifDebug=false;
//    ifDebug=true;
    vector<pair<pair<int,int>,pair<int,int>>> wBatch;
    int ID1, ID2, oldW, newW;
    srand (0);
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

            if(ifBatch){//for batch update
                if(updateBatch*batchSize>updateData.size()){
                    updateBatch=floor(updateData.size()/batchSize);
                }
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                Graph g2=*this;
                for(int u=0;u<updateBatch;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
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
                        ++update_i;
                    }

                    tt.start();
                    g2.H2HdecBat(wBatch);
//                    H2HdecBat(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheckH2H(100);
                    }
                }
                cout<<"Average Decrease batch update Time: "<<runT1/(updateBatch*batchSize)<<" s."<<endl;
            }
            else{//for single-edge update
                for(int u=0;u<updateBatch;u++){
                    wBatch.clear();
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
                    wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    tt.start();
                    H2HdecBat(wBatch);
                    tt.stop();
                    runT1 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheckH2H(100);
                    }
                }
                cout<<"Average Decrease single-edge update Time: "<<runT1/updateBatch<<" s."<<endl;
            }

//            break;
        }
        case 2:{
            //Increase update
            cout<<"Update type: Increase"<<endl;
            if(ifBatch){//for batch update
                if(updateBatch*batchSize>updateData.size()){
                    updateBatch=floor(updateData.size()/batchSize);
                }
                int update_i=0;
                vector<pair<pair<int,int>,pair<int,int>>> wBatch;
                for(int u=0;u<updateBatch;u++){
                    wBatch.clear();
                    for(int i=0;i<batchSize;++i){
                        ID1 = updateData[update_i].first.first;
                        ID2 = updateData[update_i].first.second;
                        oldW = updateData[update_i].second;
                        newW=oldW*1.5;
                        if(ifDebug){
                            cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<oldW<<" "<<newW<<endl;
                        }
                        wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                        ++update_i;
                    }

                    tt.start();
                    H2HincBatMT(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
//                        CorrectnessCheckH2H(100);
                    }

                }
                cout<<"Average Increase batch update Time: "<<runT2/(updateBatch*batchSize)<<" s."<<endl;
            }
            else {//for single-edge update
                for(int u=0;u<updateBatch;u++){
                    wBatch.clear();
                    ID1 = updateData[u].first.first;
                    ID2 = updateData[u].first.second;
                    oldW = updateData[u].second;
                    newW=oldW*2;
                    if(ifDebug){
                        cout<<"Batch "<<u<<": "<<ID1<<" "<<ID2<<" "<<" "<<oldW<<" "<<newW<<endl;
                    }
                    wBatch.emplace_back(make_pair(ID1,ID2), make_pair(oldW,newW));
                    tt.start();
                    H2HincBatMT(wBatch);
                    tt.stop();
                    runT2 += tt.GetRuntime();
                    if(ifDebug){
                        CorrectnessCheckH2H(100);
                    }
                }
                cout<<"Average Increase single-edge update Time: "<<runT2/updateBatch<<" s."<<endl;
            }


            break;
        }
        default:{
            cout<<"Wrong update type!"<<endl;
            break;
        }

    }
}

void Graph::H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	map<int,int> checkedDis;

	for(int i=0;i<Tree.size();i++){
		Tree[i].DisRe.clear();//record the star weight change (causing the distance change)
	}

	//NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
	vector<set<OrderCompp>> SCre; //SCre.clear();
	SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
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
			int Cid=it->x; int Cw;
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
                Tree[rank[ProID]].cnt[cidH]=1;//new
			}else if(Tree[rank[ProID]].dis[cidH]==Cw){
				Tree[rank[ProID]].FN[cidH]=true;
                Tree[rank[ProID]].cnt[cidH]+=1;//new
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
                            if(SCre[ProID].find(lid)==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid]].vert[k].second.second += 1;
                            }
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
		EachNodeProBDis5(rank[ProBeginVertexID], linee, vertexIDChL,checkedDis);
	}
	//return checkedDis.size();
}

void Graph::EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis){
	bool ProIDdisCha=false;
    vector<int> cntNew(line.size(),0);
    vector<bool> flagUpdate(line.size(),false);
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
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
						}
                        else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];//use cntNew to redress the cnt value since the edge decrease may lead to more ways for dis (i.e., increase the cnt)
                            }
                        }
					}
					for(int i=bH+1;i<line.size();i++){
						checkedDis.insert(make_pair(child,i));
						if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
							Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
							Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
						}
                        else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
					}

				}else{//partial ancestor check

					if(vertexIDChL.find(b)!=vertexIDChL.end()){
						for(int i=0;i<bH;i++){
							checkedDis.insert(make_pair(child,i));
							if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
								Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
								Tree[child].FN[i]=false;
                                Tree[child].cnt[i]=1;//new
                                ProIDdisCha=true;
                                flagUpdate[i]=true;
                                cntNew[i]=1;
							}
                            else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                                cntNew[i]++;
                                if(flagUpdate[i]) {
                                    Tree[child].cnt[i]+=1;//new
                                }
                                else if(cntNew[i]>Tree[child].cnt[i]){
                                    Tree[child].cnt[i]=cntNew[i];
                                }
                            }
						}
					}
					for(int i=bH+1;i<line.size();i++){
						checkedDis.insert(make_pair(child,i));
						if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
							Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
							Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
						}
                        else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
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
						checkedDis.insert(make_pair(child,i));
						if(Tree[child].dis[i]>vbW+Tree[rank[b]].dis[i]){
							Tree[child].dis[i]=vbW+Tree[rank[b]].dis[i];
							Tree[child].FN[i]=false;
                            Tree[child].cnt[i]=1;//new
                            ProIDdisCha=true;
                            flagUpdate[i]=true;
                            cntNew[i]=1;
						}
                        else if(Tree[child].dis[i]==vbW+Tree[rank[b]].dis[i]){
                            cntNew[i]++;
                            if(flagUpdate[i]) {
                                Tree[child].cnt[i]+=1;//new
                            }
                            else if(cntNew[i]>Tree[child].cnt[i]){
                                Tree[child].cnt[i]=cntNew[i];
                            }
                        }
					}
				}
				for(int i=bH+1;i<line.size();i++){
					checkedDis.insert(make_pair(child,i));
					if(Tree[child].dis[i]>vbW+Tree[rank[line[i]]].dis[bH]){
						Tree[child].dis[i]=vbW+Tree[rank[line[i]]].dis[bH];
						Tree[child].FN[i]=false;
                        Tree[child].cnt[i]=1;//new
                        ProIDdisCha=true;
                        flagUpdate[i]=true;
                        cntNew[i]=1;
					}
                    else if(Tree[child].dis[i]==vbW+Tree[rank[line[i]]].dis[bH]){
                        cntNew[i]++;
                        if(flagUpdate[i]) {
                            Tree[child].cnt[i]+=1;//new
                        }
                        else if(cntNew[i]>Tree[child].cnt[i]){
                            Tree[child].cnt[i]=cntNew[i];
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
		EachNodeProBDis5(Tree[child].ch[i], line, vertexIDChL,checkedDis);
	}
	line.pop_back();

}

void Graph::H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch){
	int checknum=0;
	map<pair<int,int>,int> OCdis;//{(s,t),d} (order[s]<order[t]), original weight of the affected shortcut, maintain the old distance before refreshed and avoid search in the adjacent list
	//OCdis.clear();

	//NodeOrderss.clear();
//	NodeOrderss.assign(NodeOrder.begin(),NodeOrder.end());
	vector<set<OrderCompp>> SCre; //SCre.clear(); the affected shortcut pair
	SCre.assign(node_num,set<OrderCompp>());//{vertexID, set<int>}
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
			int Cid=it->x; int Cw=OCdis[make_pair(ProID,Cid)];
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
                            if(SCre[ProID].find(lid)==SCre[ProID].end()) {//if not found, avoid repeated count
                                Tree[rank[lid]].vert[k].second.second -= 1;
                                if (Tree[rank[lid]].vert[k].second.second < 1) {
                                    SCre[lid].insert(Cid);
                                    OC.insert(OrderCompp(lid));
                                    OCdis[make_pair(lid, Cid)] = Cw + Lnei[j].second;
                                }
                            }
						}
						break;
					}
				}
			}

            //get the new value of shortcut
            //	cout<<Cw<<" increase to ";
            int newCw=INF; int countwt=0;

            for(int i=0;i<Neighbor[ProID].size();i++){
                if(Neighbor[ProID][i].first==Cid){
                    newCw=Neighbor[ProID][i].second;//the weight value in the original graph
                    countwt=1;
                    break;
                }
            }

            int ssw,wtt,wid;
            vector<int> Wnodes;
            Wnodes.clear();

            if(ProID<Cid)
                Wnodes=SCconNodesMT[ProID][Cid]; //cout<<"wid num "<<Wnodes.size()<<endl;
            else
                Wnodes=SCconNodesMT[Cid][ProID];
            if(!Wnodes.empty()){
                for(int i=0;i<Wnodes.size();i++){
                    wid=Wnodes[i];
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

            if(newCw>Cw){
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

		eachNodeProcessIncrease1(rank[ProBeginVertexID], linee,checknum);
	}
	//return checknum;
}

void Graph::eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel){
	int childID=Tree[children].uniqueVertex;
	int childH=Tree[children].height-1;
	for(int i=0;i<Tree[children].dis.size();i++){
		if(Tree[children].cnt[i]<=0){//if the distance label to i-th ancestor should be maintained, use <= since cnt may not be accurate after update
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
//				if(Tree[PID].height>Tree[children].height){///modified for correctness
                if(Tree[PID].height>Tree[children].height && Tree[PID].vAncestor[childH] == childID){
//                    if(PID>Tree.size()){
//                        cout<<"PID error! "<<PID<<" "<<Tree.size()<<endl; exit(1);
//                    }
//                    if(childH>Tree[PID].dis.size()){
//                        cout<<"childH error! "<<childH<<" "<<Tree[PID].dis.size()<<": "<<children<<"("<<Tree[children].height<<") "<<PID<<"("<<Tree[PID].height<<")"<<endl; exit(1);
//                    }
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
		eachNodeProcessIncrease1(Tree[children].ch[i],line,changelabel);
	}
	line.pop_back();
}
