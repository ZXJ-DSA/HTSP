/*
 * head.h
 *
 *  Created on: 24 August 2022
 *      Author: Xinjie ZHOU
 */

#ifndef HEAD_H_
#define HEAD_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <boost/thread/thread.hpp>
#include <chrono>
#include <string>
#include "Heap.h"
//#include "labeling.hpp"
#include <omp.h>

#define INF 99999999
#define Dijk 0
#define PCH_No 1
#define PH2H_No 2
#define PH2H_Post 3
#define PH2H_Cross 4
#define CH 1
#define H2H 4


//typedef unsigned int vertex;
typedef int vertex;

using namespace std;
//using namespace boost;


extern vector<int> NodeOrder_;//nodeID order
extern vector<int> _DD_;
extern vector<int> _DD2_;

struct Nei{
	int nid;
	int w;
	int c;
};

struct OrderCompMax{// maximum-first, Higher-order first
    int ID;
    OrderCompMax(){ID=0;}
    OrderCompMax(int _ID){
        ID=_ID;
    }
    bool operator< (const OrderCompMax d) const{
        if(NodeOrder_[ID]!=NodeOrder_[d.ID])
            return NodeOrder_[ID]>NodeOrder_[d.ID];
        return ID>d.ID;
    }
};

struct OrderCompMin{//prior to reture the vertex with smaller order
    int x;
    OrderCompMin(int _x){
        x=_x;
    }
    bool operator< (const OrderCompMin& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
        }
    }
};

struct OrderCompp{//prior to return the vertex with smaller order
    int x;
    OrderCompp(int _x){
        x=_x;
    }
    bool operator< (const OrderCompp& d) const{
        if(x==d.x){//avoid the redundant
            return false;
        }else{
            if(x!=d.x)
                return NodeOrder_[x]<NodeOrder_[d.x];
        }
    }
};


//return smallest Dis with largest-order vertex
struct MinComp{
    int s, t, Dis;//s is the source vertex while t is the target vertex
    MinComp(){
        s=0, t=0, Dis=0;
    }
    MinComp(int _ID1, int _ID2, int _Dis){
        s=_ID1; t=_ID2; Dis=_Dis;
    }
    bool operator< (const MinComp d) const{
        if(Dis != d.Dis){
            return Dis<d.Dis;
        }else{
            if(NodeOrder_[s]!=NodeOrder_[d.s]){
                return NodeOrder_[s]>NodeOrder_[d.s];
            }
            return NodeOrder_[t]>NodeOrder_[d.t];
        }
    }
};


struct DegComp{//min-first
    int x;
    DegComp(int _x){
        x=_x;
    }
    bool operator < (const DegComp d) const{
        if(_DD_[x]!=_DD_[d.x])
            return _DD_[x]<_DD_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct DegComp2{//min-first
    int x;
    DegComp2(int _x){
        x=_x;
    }
    bool operator < (const DegComp2 d) const{
        if(_DD2_[x]!=_DD2_[d.x])
            return _DD2_[x]<_DD2_[d.x];
        return x<d.x;
    }
//    bool operator () (const DegComp1& a, const DegComp1& b) const{
//        return (_DD_[a.x] < _DD_[b.x]) || (_DD_[b.x] >= _DD_[a.x] && (a.x < b.x));
//    }

};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count(how many ways can lead to this super edge weight)
	vector<int> pos;
	vector<int> dis, cnt;//the distance value and corresponding count number (i.e., record how many path has the shortest distance)
    vector<int> disPost, cntPost;//the distance value of post-boundary strategy and corresponding count number
    vector<int> disInf;//the distance from this vertex to the boundary vertices and corresponding count number, cntInf
    vector<int> vAncestor;//the ancestors, which is corresponding to dis
    int rootPos;//index position of the root vertex of this partition in vAncestor
	vector<bool> FN;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert)
    vector<bool> FNPost;//another succinct way of FromNode, whether this distance label is directly obtained from shortcuts (vert)
    vector<bool> FNInf;//whether the interface distance is obtained from shortcuts (vert)
	set<int> DisRe;//record the vertex id that the distance label should be updated
    set<int> DisRePost;//record the vertex id that the interface label should be updated
	vector<int> ch;//position index of children
	int height, hdepth;//hdepty is the deepest node that a vertex still exists
	int pa;//parent, the pa of root vertex is 0, position index
	int uniqueVertex;//vertex id of this tree node
//	vector<int> piv;//pivot vetex, used in path retrieval
//    int treeroot;//the tree id of subtree root, i.e., rank[x]
	Node(){
		vert.clear();
		pos.clear();
		dis.clear(); cnt.clear();
        disPost.clear(); cntPost.clear();
        disInf.clear();
        vAncestor.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
		FN.clear(); FNInf.clear();
		DisRe.clear();
        DisRePost.clear();
//		piv.clear();
//        treeroot=-1;
	}
};

class Graph{
public:
    string graphfile;
    string dataset;
	int node_num=0;    //vertex number
	unsigned long long edge_num=0;    //edge number
	vector<vector<pair<vertex,int>>> Neighbor;//original graph
    vector<vector<pair<vertex,int>>> NeighborsParti;//<node_number,<in-partition adjacency lists>>
    vector<unordered_map<vertex,int>> NeighborsPartiPost;//adjacency lists of post-boundary partitions
    vector<unordered_map<vertex,int>> NeighborsOverlay;//<node_number,<adjacency lists of overlay graph>>
    vector<vector<pair<vertex,int>>> NeighborsOverlayV;//<node_number,<adjacency lists of overlay graph>>
//    vector<unordered_map<int,int>> BoundaryShortcuts;
    vector<pair<int,bool>> PartiTag;//<node_number,<partition_id,if_boundary>>
    vector<pair<int,set<int>>> PartiTags;//<node_number,<flag,set<partition_id>>>, flag=-1 means boundary vertex, flag>=0 means partition id
    vector<vector<vertex>> PartiVertex;//<partition_number,<in-partition vertices>>, in increasing vertex order, higher-rank vertex first
    vector<vector<vertex>> BoundVertex;//boundary vertices of each partition
//    vector<unordered_set<int>> BoundVertexSet;//boundary vertices of each partition
    vector<unordered_map<int,int>> BoundVertexMap;//map from boundary vertex id to its position id in BoundVertex
    vector<vertex> OverlayVertex;//overlay vertex in decreasing vertex order
//    vector<unordered_map<vertex,pair<int,int>>> repairShortcuts;//<ID1,ID2,<overlay weight, in-partition weight>>
    vector<map<int,int>> BoundShortcuts;//<ID1,ID2,weight>, source=-1 means stemming from core while source>=0 means stemming from specific partition
    int partiNum;   //partition number
    bool ifParallel = true;
    bool ifIncrease = false;
    int nodeNumOverlay = 0;//vertex number of overlay graph
    vector<double> stageDurations;//durations of querying stages
    vector<pair<int,int>> GraphLocation;//Longitude and Latitude, already times 1000,000

    //PostMHL
    double bRatioLower=0.2;//imbalance ratio of partitions
    double bRatioUpper=2;
    int bandWidth=100;//50;//bandwidth of PostMHL
    int HighestOrder;
    vector<bool> existOverlay;
    vector<int> partiRoots;//vertex id of partition root
    vector<map<int,map<int,int>>> SuppPartiID;//<ID1,<ID2,<pid,weight>>>, record the partition and its supportive value for a given interface edge
    vector<map<int,pair<int,set<int>>>> SuppPartiIDReal;//ID1,<ID2,<weight,set<pid>>>>, //record the partitions that really support a given interface edge
    set<int> affectedParti;//record the affected partitions
    vector<int> ProBeginIDV;
    vector<int> childNums;//children number of each vertex

	vector<int> DD; //intermediate variable in Contraction, DD2
	int threadnum=15;  //thread number
    int algoQuery=0;//algorithm for querying. 0:Dijkstra; 1: No-boundary; 2: Post-boundary; 3: Extended Label. 0: Dijkstra; 1: CH; 2:H2H. 0: Dijkstra; 1: PCH-No; 2: PH2H-No; 3: PCH-Post; 4: PH2H-Post; 5: PH2H-Extend
    int algoUpdate=0;//algorithm for core construction, (0: BPCL; 1: PCL; 2: PLL; 3: WPSL; 4: GLL; 5: Read)
    string algoParti="NC";
    int algoChoice=1;//HTSP system index. 1: Non-partition; 2: Partition

	//vertex order
	vector<int> NodeOrder;//nodeID order
	vector<int> vNodeOrder;//order nodeID
    int treewidth=0;//treewidth

    double tAncestor=0,tBoundary=0;

    /// Index Construction
//    vector<omp_lock_t> oml;
    unordered_map<int, Semaphore*> mSm;
    vector<Semaphore*> vSm;
    Semaphore* sm = new Semaphore(1);;// = new Semaphore(threadnum);

    //H2H index construction
    //intermediate variable and function used in the H2H index construction
    vector<map<int,pair<int,int>>> E;//ID1,ID2,(weight,count)
//    vector<map<OrderCompMin, pair<int,int>>> EOrder;
    vector<vector<pair<int,pair<int,int>>>> NeighborCon;//ID1,ID2,(weight,count)
    //for overlay graph
    vector<Node> Tree;
    vector<int> toRMQ;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndex;//?
    vector<int> rank;//rank[v]>0 indicates non-core vertex
    int heightMax;
//    vector<int> EulerSeq;//prepare for the LCA calculation, EulerSeq is the Euler tour, i.e., E[1,...,2n-1]
    //    vector<map<int, vector<int>>> SCconNodesMT;//supportive vertex, multiple thread of SCconNodes
    vector<map<int, vector<pair<int,int>>>> SCconNodesMT;//<ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNid;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for overlay

    //for no-boundary partitions
    vector<vertex> IDMap;//map the old id to new id in partitions
    vector<vector<Node>> Trees;//Trees for no-boundary
    vector<vector<int>> toRMQs;
    vector<vector<vector<int>>> RMQIndexs;
    vector<vector<int>> ranks;
    vector<int> heightMaxs;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMTP;//for partitions. <ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNidP;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for partition

    //for post-boundary partitions
//    vector<vertex> IDMapPost;//map the old id to new id in partitions
    vector<vector<Node>> TreesPost;//Trees for post-boundary
    vector<vector<int>> toRMQsPost;
    vector<vector<vector<int>>> RMQIndexsPost;
    vector<vector<int>> ranksPost;
    vector<int> heightMaxsPost;
    vector<map<int, vector<pair<int,int>>>> SCconNodesMTPost;//for partitions. <ID1,<ID2,<x,weight>>> where ID1<ID2, record the supportive vertices of a shortcut, only record edge once
    vector<vector<int>> VidtoTNidPost;//record the child tree nodes whose vert neighbors contain this tree node (nodeID--->tree node rank), for partition

    vector<bool> ifRepaired;

    vector<int> ProBeginVertexSetOverlay;
    set<int> vertexIDChLOverlay;
    vector<vector<int>> ProBeginVertexSetParti;
    vector<set<int>> vertexIDChLParti;
    vector<vector<int>> ProBeginVertexSetPartiExtend;//for post-boundary update
//    vector<set<int>> vertexIDChLPartiExtend;

    //extension label
    vector<Node> TreeExt;
    vector<int> toRMQExt;//toRMQ[p] records the first occurrence of node p in the Euler tour, i.e., R[p]
    vector<vector<int>> RMQIndexExt;//?
    vector<int> rankExt;//rank[v]>0 indicates non-core vertex
    int heightMaxExt;

    vector<bool> vUpdated;// flag of whether the label of vertex has been updated
    vector<int> bHeights;

    ~Graph(){
        clear();
    }
    void clear(){
        Neighbor.clear();
        Tree.clear();
        vSm.clear();
    }

    /// For non-partition SP index
    void HybridSPIndexConstruct();
    void MDEContraction(string orderfile);
//    void deleteEOrderGenerate(int u,int v);
    void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p);
    void insertEMTOrderGenerate(int u,int v,int w);
    //for contraction
    void deleteEorder(int u,int v);
    void insertEorder(int u,int v,int w);
    int match(int x,vector<pair<int,pair<int,int>>> &vert);
    void makeTree();
    void makeIndex();
    void IndexsizeH2H();  //Index size computation

    /// For partitioned SP index
    void HybridPSPIndexConstruct();
    void PMHLIndexConstruct();
    void PMHLIndexConstructOpt();
    void PostMHLIndexConstruct();

    /// For PH2H
    void PH2HIndexConstruct(); //PH2H index construction
    void ConstructBoundaryShortcutV(vector<int> & p, bool ifAllPair);
    void ConstructBoundaryShortcut(int pid);
    void ConstructBoundaryShortcutNoAllPair(int pid);
    void Construct_PartiIndex(bool ifParallel, bool ifLabelC);
    void Construct_OverlayGraph(bool ifParallel);
    void Construct_OverlayGraphNoAllPair(bool ifParallel);
    void Construct_OverlayIndex(bool ifLabelC);

    void PreConstructAllPairs(bool ifParallel);
    void PreConstructAllPairsPartiV(vector<int> & p);
    void PreConstructAllPairsParti(int pid);

    /// Post-boundary index
    void RefreshBoundaryEdgesAndLabelingPartiV(vector<int>& p);
    void RefreshBoundaryEdgesAndLabelingParti(int pid);
    void ConstructPartitionPost(bool ifParallel);
    void ConstructPostParti(int pid);
    void ConstructPostPartiV(vector<int>& p);
    void ConstructPartitionPostIndex(bool ifParallel, bool ifLabelU);
    void ConstructPartitionPostIndexOpt(bool ifParallel);
    void Repair_PartiIndex(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void Repair_PartiIndexPostMHLPost(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, double & runT);
    void Repair_PartiIndexForOpt(bool ifParallel, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexV(vector<int>& p, bool ifIncrease, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch, bool ifOpt);
    void RepairPartitionIndexDecrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexDecreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexIncrease(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RepairPartitionIndexIncreaseForOpt(int pid, map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);

    void PostMHLIndexUpdatePostPartV(vector<int>& p, bool ifIncrease);
    void PostMHLIndexUpdatePostPart(int pid, bool ifIncrease);
    void PostMHLIndexUpdatePostPartDFS(int pid, int p, vector<int>& ancestor, vector<int> & interface, vector<map<int,int>>& disInfs, int rootHeight, bool ifIncrease);
    void PostMHLIndexUpdatePostPartDFS(int pid, int p, vector<int>& ancestor, vector<int> & interface, vector<map<int,int>>& disInfs, set<int>& vertexIDChL, int rootHeight, bool ifIncrease, Timer& tt);

    /// Extension index
    void makeTreeExtDFS(int ID, vector<int> &list);
    void makeTreeIndexExtDFS(int ID, vector<int> &list);
    void ConstructExtensionLabelsNoAllPair();
    void ConstructExtensionLabelsNoAllPairTopDown();
    void ConstructExtensionLabels();
    void ConstructExtensionLabelPartiV(vector<int>& p, bool ifAllPair);
    void ConstructExtensionLabelParti(int pid);
    void ConstructExtensionLabelPartiNoAllPair(int pid);
    void RefreshExtensionLabelsNoAllPair(map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RefreshExtensionLabelsPostMHL(bool ifParallel, bool ifIncrease, double & runT);
    void RefreshExtensionLabels(map<int, vector<pair<pair<int, int>, pair<int, int>>>> & partiBatch);
    void RefreshExtensionLabelPartiV(vector<int>& p, bool ifTopDown);
    void RefreshExtensionLabelParti(int pid);
    void RefreshExtensionLabelPartiTopDown(int pid);
    void RefreshExtensionLabelTopDownV(vector<int>& p);
    void RefreshExtensionLabelTopDown(int ID);
    pair<int,int> ExtensionLabelCompute(int pid, int ID, int ancestor);

    /// PostMHL
    void TreeDecompositionPartitioningNaive();
    void TreeDecompositionPartitioning(int pNum, double bRatioUpper, double bRatioLower);
    void TDPContract();
    void TDPCreateTreeAndPartiNaive();
    void TDPCreateTreeAndParti(int pNum, double bRatioUpper, double bRatioLower);
    void MinimizeBoundaryNum(int upperB, int lowerB);
    void GetPartiRootCandidates(int p, int upperB, int lowerB, vector<pair<int,int>>& candidates, map<int,vector<int>>& candidatesMap);
    void GetChildrenNumDFS(int p, vector<int>& ancestors);
    void makePartitionDFS(vector<set<OrderCompMax>>& orderPartiV, int pid, int vid);
    void PostMHLIndexConstructOverlay();
    void PostMHLCreateIndex_Overlay();
    void PostMHLmakeIndexDFSOverlay(int p, vector<int>& ancestors);
    void PostMHLIndexConstructPost(bool ifParallel, double &t);
    void PostMHLPartitionIndexConstructV(vector<int>& p, bool ifPost);
    void PostMHLPartitionIndexConstruct(int pid, bool ifPost);
    void PostMHLmakeIndexDFSPartiPost(int p, vector<int>& ancestor, vector<int> & interface, map<int, unordered_map<int,int>>& disInfs);
    void PostMHLmakeIndexDFSParti(int p, vector<int>& ancestor, vector<int> & interface);
    void PostMHLIndexConstructExtend(bool ifParallel, double &t);
    void PostMHLPartitionIndexConstructExtend(int pid);
    void PostMHLPartitionIndexConstructExtendV(vector<int>& p);
    void PostMHLmakeIndexDFSPartiExtend(int p, vector<int>& ancestors);
    void IndexSizePostMHL();


    void BoundShortcutsCheck(bool ifParallel, bool ifIncrease);
    void BoundShortcutsCheckPartiV(pair<int,int> pRange, vector<int>& pids, vector<vector<pair<pair<int,int>,int>>>& affectedShortcuts, bool ifIncrease);
    void BoundShortcutsCheckParti(int pid, vector<pair<pair<int,int>,int>>& affectedShortcut, bool ifIncrease);
    void InterfacePropagate(int child, vector<int>& interfaces, vector<Node> &Tree, bool ifIncrease);
    void InterfacePropagateParallel(pair<int,int> pRange, vector<int>& pids, bool ifIncrease);
    void InterfaceEntryDecreaseUpdate(bool ifParallel);
    void AncestorEntryDecreaseUpdateParti(int pid);
    void AncestorEntryDecreaseUpdate(int pid, int child,vector<int>& line, vector<int>& interfaces, set<int>& vertexIDChL, map<int,int>& checkedDis, vector<Node> &Tree, vector<int> &rank, int rootHeight);

    void ThreadDistribute(vector<int>& vertices, vector<vector<int>>& processID);
    void ConstructPH2H_PartiV(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_PartiVCH(vector<int> & P, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_Parti(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
//    void ConstructPH2H_PartiNoLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void ConstructPH2H_PartiLabel(int pid, vector<vector<Node>>& Trees, vector<vector<int>>& ranks, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs);
    void H2HCreateTree_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP, vector<map<int, vector<pair<int,int>>>>& SCconNodesMTP, vector<vector<int>>& VidtoTNidP, vector<int>& heightMaxs);//Create tree for partition
    void H2HCreateIndex_Parti(int pid, vector<Node>& TreeP, vector<int>& rankP);//Create labels for partition

    void H2HCreateTree_Overlay();
    void H2HCreateIndex_Overlay();
    int ShortcutDisCheck(int ID1, int ID2);

    void makeTreeIndexDFSP(int p, vector<int>& list,  vector<Node>& TreeP, vector<int>& rankP);

    void makeRMQCoreP(int pid, vector<vector<int>>& toRMQs, vector<vector<vector<int>>>& RMQIndexs, vector<vector<Node>>& Trees);
    void makeRMQDFSCoreP(int pid, int p, int height, vector<int>& EulerSeqP, vector<vector<int>>& toRMQs, vector<vector<Node>>& Trees);
    void makeRMQCore();
    void makeRMQDFSCore(int p, int height, vector<int>& EulerSeq);
    void makeIndexDFS(int p, vector<int> &list);

    void makeRMQ(vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree);
    void makeRMQDFS(int p, int height, vector<int>& EulerSeq, vector<int>& toRMQ, vector<Node>& Tree);
//    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
//    void insertEMTorder(int u,int v,int w);
    void deleteECore(int u,int v);
    void insertECore(int u,int v,int w);
    void insertECoreMT(int u,int v,int w);
    int matchCore(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank);//vector<pair<int,int>> &vert
    int matchCoreParti(int x,vector<pair<int,pair<int,int>>> &vert, vector<int>& rank);//vector<pair<int,int>> &vert
    void IndexSizePH2H();  //Core-tree index size computation
    void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void NeighborComorderH2H(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
    void insertEMTorder(int u,int v,int w);
    void PH2HVertexOrdering(int type);
    void SketchGraphBuild();
    void OverlayOrderingBuild();
    void OverlayOrderingBuildBoundaryFirst(int nodenum, vector<vector<int>> Neighbor);
    void PartitionOrderingBuildMDE(bool ifParallel);
    void OrderingAssemblyMDEBoundaryFirst(int pNum);
    void OrderingAssemblyMDE(int pNum);
    void OrderingAssemblyBoundaryFirst(int pNum);
    void SketchOrder(vector<vector<pair<int,int>>> Neighbor, vector<int> &vNodeOrderSketch);

    void PartitionOrderingV(vector<int>& p);
    void PartitionOrdering(int pid);
    void deleteEOrderGenerate(int u,int v);
    void insertEOrderGenerate(int u,int v,int w);

//    vector<pair<pair<int,int>,int>> CutEdges;//the cut edges
    vector<vector<int>> NeighborSketch;
    vector<set<int>> NeighborSketchS;
//    vector<map<int,int>> vNodeOrderParti;
    vector<vector<int>> vNodeOrderParti;
    vector<int> vNodeOrderOverlay;

	///Query processing
    int QueryNP(int ID1, int ID2);
    int QueryCHWP(int ID1, int ID2);
    int QueryH2H(int ID1, int ID2);

	//PH2H
    int QueryCore(int ID1, int ID2);
    int QueryH2HPartition(int ID1, int ID2, int PID);
    int QueryH2HPartitionPost(int ID1, int ID2, int PID);
    int QueryPartiCore(int ID1, int ID2);
    int QueryPartiCoreDebug(int ID1, int ID2);
    int QueryPartiCoreExt(int ID1, int ID2);
    int QueryPartiCoreExtLCA(int ID1, int ID2);
    int QuerySameParti(int ID1, int ID2);
    int QuerySamePartiPost(int ID1, int ID2);
    int QuerySamePartiPostOpt(int ID1, int ID2);
    int QueryPartiParti(int ID1, int ID2);
    int QueryPartiPartiExt(int ID1, int ID2);
    int QueryPartiPartiExtLCA(int ID1, int ID2);
    int QueryPartiPartiExtLCADebug(int ID1, int ID2);



    //PCH
    int QueryCoreCH(int ID1, int ID2);
    int QueryCHPartition(int ID1, int ID2, int PID);
    int QueryPartiCoreCH(int ID1, int ID2);
//    int QuerySamePartiCH(int ID1, int ID2);
    int QueryPartiPartiCH(int ID1, int ID2);


    //PostMHL
    int QueryOverlayCH(int ID1, int ID2);
    int QueryOverlay(int ID1, int ID2);
    int QueryOverlayDebug(int ID1, int ID2);
    int QueryPostMHLSamePartiPost(int ID1, int ID2);
    int QueryPostMHLSamePartiPostDebug(int ID1, int ID2);
    int QueryPostMHLPartiParti(int ID1, int ID2);
    int QueryPostMHLPartiPartiDebug(int ID1, int ID2);
    int QueryPostMHLPartiPartiExt(int ID1, int ID2);
    int QueryPostMHLPartiPartiExtDebug(int ID1, int ID2);
    int QueryPostMHLPartiOverlay(int ID1, int ID2);
    int QueryPostMHLPartiOverlayDebug(int ID1, int ID2);


    void EffiCheck(string filename,int runtimes);
    void EffiCheckStages(vector<pair<int,int>> & ODpair, int runtimes, int intervalT, unsigned long long & throughputNum, vector<double>& stageUpdateT, vector<double>& stageQueryT);
    void StagePerformanceShow(int batchNum, vector<double>& stageUpdateT, vector<double>& stageQueryT);
    double EffiMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType);
    double EffiPMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType);
    double EffiPostMHLStage(vector<pair<int,int>> & ODpair, int runtimes, int queryType);
    int QueryPostMHL(int ID1, int ID2);
    int QueryPostMHLDebug(int ID1, int ID2);
    int QueryPMHLOpt(int ID1, int ID2);//algoQuery: 0: Dijkstra; 1: PCH-No; 4: PH2H-Post; 5: PH2H-Extend
    int QueryPMHL(int ID1, int ID2);//algoQuery: 0: Dijkstra; 1: PCH-No; 2: PH2H-No; 3: PCH-Post; 4: PH2H-Post; 5: PH2H-Extend
	int Query(int ID1, int ID2);//algoQuery: 0: Dijkstra; 2: PH2H-No; 4: PH2H-Post; 5: PH2H-Extend
    int QueryDebug(int ID1, int ID2);

    int QueryCoreDebug(int ID1, int ID2);
    int LCAQuery(int _p, int _q);
    int LCAQueryPartition(int _p, int _q, int PID);// query within partition
    int LCAQueryPartitionPost(int _p, int _q, int PID);// query within partition
    int LCAQueryOverlay(int _p, int _q);
    int LCAQuery(int _p, int _q, vector<int>& toRMQ, vector<vector<int>>& RMQIndex, vector<Node>& Tree);

    //Correctness Check
    void CorrectnessCheck(int runtimes);
    void CorrectnessCheckCore(int runtimes);
    void DFSTree(vector<int>& tNodes, int id);


    //Dijkstra
    int Dijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int BiDijkstra(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int Astar(int ID1, int ID2,vector<vector<pair<vertex,int>>> &Neighbor);
    int EuclideanDis(int s, int t);
    void RetrievePath(int ID1, int ID2, vector<int> & prece);
    int DijkstraCore(int ID1, int ID2);

    //// For throughput test
    void SPThroughputTest(int updateType, bool ifBatch, int batchNum, int batchSize, int batchInterval, int runtimes);
    void DecBatchThroughput(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);//process batch update
    void IncBatchThroughput(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);//process batch update
    void DecBatchThroughputNP(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);//process batch update
    void EachNodeProBDis5H2H(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis);
    void IncBatchThroughputNP(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);//process batch update
    void eachNodeProcessIncrease1H2H(int children, vector<int>& line, int& changelabel);

    void PMHLBatchUpdateDec(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PMHLBatchUpdateDecOpt(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PostMHLBatchUpdateDec(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PMHLBatchUpdateInc(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PMHLBatchUpdateIncOpt(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    void PostMHLBatchUpdateInc(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, int batch_i, double& runT1);
    unsigned long long EffiCheckThroughput(vector<pair<int,int>>& ODpair, Timer& tRecord, int batchInterval, unsigned long long& throughputNum);

    void IndexConstruction();

	/// Index update
    void IndexMaintenance(int updateType, int updateSize, bool ifBatch, int batchSize);
    void DecreaseSingle(int a, int b, int oldW, int newW);//process one update
    void IncreaseSingle(int a, int b, int oldW, int newW);//process one update
    void DecreaseBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//process batch update
    void IncreaseBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//process batch update
    void DecreaseOverlay(int a,int b, int newW, vector<unordered_map<vertex,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void DecreaseParti(int a,int b, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax);
    void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis,
    void EachNodeProBDis5PostMHLOverlay(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis,
    void EachNodeProBDis5Parti(int child,vector<int>& line,set<int>& vertexIDChL, vector<Node> &Tree, vector<int> &rank);//map<int,int>& checkedDis, map<int,int>& checkedDis,
    void IncreaseOverlay(int a, int b, int oldW, int newW, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1PostMHL(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1PostMHLOverlay(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void IncreaseParti(int a, int b, int oldW, int newW, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid);
    void eachNodeProcessIncrease1Parti(int children, vector<int>& line, int& changelabel, vector<Node> &Tree, vector<int> &rank, vector<vector<int>> &VidtoTNid);
    void DecreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);//batch decrease for overlay graph
    void DecreaseH2HBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void DecreaseOverlayBatchPostMHL(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);//batch decrease for overlay graph
    void DecreaseOverlayBatchLabelPostMHL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet,set<int>& vertexIDChL);
    void DecreasePartiBatchLabelPostMHLExtendV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>>& ProBeginVertexSetV, set<int> &vertexIDChL);
    void DecreasePartiBatchLabelPostMHLExtend(int pid, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, set<int> &vertexIDChL);
    void DecreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet,set<int>& vertexIDChL);
    void DecreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay);
    void DecreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt, vector<pair<pair<int,int>,int>>& updatedSC);
    void DecreasePartiBatchUpdateCheckPostMHL(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int,int>,pair<int,int>>>& overlayBatch, bool ifParallel);
    void DecreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void DecreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU, bool ifConstruct);
    void DecreasePartiBatchPostMHLShortcut(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch,  vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC);
    void DecreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>>& Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax, bool ifLabelU);
    void DecreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet,set<int>& vertexIDChL);
    void IncreaseOverlayBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreaseH2HBatch(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbor, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreaseOverlayBatchPostMHL(vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreaseOverlayBatchLabelPostMHL(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreaseOverlayBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchUpdateCheck(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay);
    void IncreasePartiBatchUpdateCheckCH(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<pair<pair<int,int>,pair<int,int>>>& weightOverlay, bool ifOpt,vector<pair<pair<int,int>,int>>& updatedSC);
    void IncreasePartiBatchUpdateCheckPostMHL(map<int, vector<pair<pair<int, int>, pair<int, int>>>>& partiBatch, vector<pair<pair<int, int>, pair<int, int>>> &overlayBatch, bool ifParallel);
    void IncreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void IncreasePartiBatchForOpt(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<vector<pair<int,int>>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, vector<pair<pair<int,int>,int>>& updatedSC, bool ifLabelU);
    void IncreasePartiBatchPostMHLShortcut(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch,  vector<Node> &Tree, vector<int> &rank, int heightMax, vector<pair<pair<int,int>,int>>& updatedSC, vector<int>& PropagateOverlay);
    void IncreasePartiBatch(int pid, vector<pair<pair<int,int>,pair<int,int>>>& wBatch, vector<unordered_map<int,int>> &Neighbors, vector<Node> &Tree, vector<int> &rank, int heightMax,vector<map<int, vector<pair<int,int>>>> &SCconNodesMT, vector<vector<int>> &VidtoTNid, bool ifLabelU);
    void IncreasePartiBatchLabel(vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchLabelPostMHLExtendV(vector<int>& p, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<vector<int>>& ProBeginVertexSetV, vector<vector<int>> &VidtoTNid);
    void IncreasePartiBatchLabelPostMHLExtend(int pid, vector<Node> &Tree, vector<int> &rank, int heightMax, vector<int>& ProBeginVertexSet, vector<vector<int>> &VidtoTNid);

	/// Graph Preprocessing
    void ReadGraph(string filename);
    void ReadCoordinate(string filename);
    void ReadUpdate(string filename,vector<pair<pair<int,int>,int>>& TestData);
    void ReadUpdate2(string filename,vector<pair<pair<int,int>,pair<int,int>>>& TestData);
    void ReadUpdate3(string filename,vector<pair<pair<int,int>,tuple<int,int,int>>>& TestData);
	void StainingMethod(int ID);
	void ODGene(int num, string filename);
    void ODGeneParti(int num, string filename);
    void ODGeneSameParti(int num, string filename);
    void ODGeneCrossParti(int num, string filename);
	void UpdateGene(int num, string filename);
    void QueryGenerationParti(bool ifSame);


    void WriteTreeIndexOverlay(string filename);
    void ReadTreeIndex(string file);
    void WriteTreeIndexParti(string filename);
    void WriteGraph(string graphfile);
    void WriteOrder(string filename);
    void ReadOrder(string filename);
    void CompareOrder(string filename1, string filename2);
    void GraphPartitionRead(string filename);

    vector<int> DFS_CC(vector<map<int,int>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);
    vector<int> DFS_CC(vector<vector<pair<int,int>>> & Edges, set<int> set_A, set<int> & set_B, int nodenum);

};

#endif // HEAD_H_
