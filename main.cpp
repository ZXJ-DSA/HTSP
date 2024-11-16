/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "head.h"

int main(int argc, char** argv){

    if( argc < 4 || argc > 19){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
//        printf("<arg3> HTSP system index, 1: CH+H2H; 2: PH2H; 3: PCH+PH2H; 4: Optimized PCH+PH2H; 5: PostMHL. default: 2\n");
        printf("<arg3> HTSP system index, 1: CH+H2H; 3: PCH+PH2H; 5: PostMHL. default: 5\n");
        printf("<arg4> (optional) average query response time requirement, in seconds. default: 0.5\n");
        printf("<arg5> (optional) partition number, e.g. 64\n");
        printf("<arg6> (optional) partition method, (NC: PUNCH; MT: METIS), default: NC\n");
        printf("<arg7> (optional) query strategy, (0:A*; 1: PCH (CH); 2: No-boundary; 3: Post-boundary; 4: Cross-boundary (H2H)), default: 4\n");
        printf("<arg8> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase), default: 0\n");
        printf("<arg9> (optional) batch number, default: 10\n");
        printf("<arg10> (optional) batch size, default: 10\n");
        printf("<arg11> (optional) batch interval (in seconds), default: 10\n");
        printf("<arg12> (optional) thread number, default: 15\n");
        printf("<arg13> (optional) query worker number, default: 15\n");
        printf("<arg14> (optional) same-partition query proportion, -1: random; 0-100: 0-100%,  default: -1 (random)\n");
        printf("<arg15> (optional) bandwidth, default: 50\n");
        printf("<arg16> (optional) lower bound ratio, default: 0.1\n");
        printf("<arg17> (optional) upper bound ratio, default: 2\n");
        printf("<arg18> (optional) preprocessing task (1: Partitioned MDE Ordering; 2: Partitioned Query Generation)\n");

        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";
    int algoChoice = 5;
    int partitionNum = 20;
    int algoQuery = 4;
    string algoParti = "NC";
    int updateType = 0;
    int runtimes = 10000;
//    runtimes = 1000;
    int batchNum = 10;
    batchNum = 2;
    int batchSize = 10;
    int batchInterval = 60;
//    bool ifBatch = false;
    int threadNum = 15;
    int preTask=0;
    int bandwidth=50;
    double lowerB=0.1;
    double upperB=2;
    int portion=-1;
    double T_r=0.5;//average query response time, s
    int workerNum = 15;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1] (Source Path): " << argv[1] << endl;//source path
        DesFile = argv[1];

        cout << "argv[2] (Dataset): " << argv[2] << endl;//dataset
        dataset = argv[2];

        cout << "argv[3] (System Index): " << argv[3] << endl;//system index
        algoChoice = stoi(argv[3]);

        if(argc > 4){
            cout << "argv[4] (Query Response Time, s): " << argv[4] << endl;//query response time
            T_r = stod(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5] (Partition Number): " << argv[5] << endl;//partition number
            partitionNum = stoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6] (Partition Method): " << argv[6] << endl;//partition method
            algoParti = argv[6];
//            if(algoParti != "NC" && algoParti != "MT"){
//                cout<<"Wrong partition method! "<<algoParti<<endl; exit(1);
//            }
        }
        if(argc > 7){
            cout << "argv[7] (Query Strategy): " << argv[7] << endl;//algorithm for query
            algoQuery = stoi(argv[7]);
        }

        if(argc > 8){
            cout << "argv[8] (Update Type): " << argv[8] << endl;//update type
            updateType = stoi(argv[8]);
        }

        if(argc > 9){
            cout << "argv[9] (Batch Number): " << argv[9] << endl;//batch number
            batchNum = stoi(argv[9]);
        }
        if(argc > 10){
            cout << "argv[10] (Batch Size): " << argv[10] << endl;//batch size
            batchSize = stoi(argv[10]);
        }
        if(argc > 11){
            cout << "argv[11] (Batch Interval): " << argv[11] << endl;//batch interval
            batchInterval = stoi(argv[11]);
        }
        if(argc > 12){
            cout << "argv[12] (Thread Number): " << argv[12] << endl;//thread number
            threadNum = stoi(argv[12]);
        }
        if(argc > 13){
            cout << "argv[13] (Query worker number): " << argv[13] << endl;
            workerNum = stoi(argv[13]);
        }
        if(argc > 14){
            cout << "argv[14] (Same-parti Proportion): " << argv[14] << endl;
            portion = stoi(argv[14]);
        }
        if(argc > 15){
            cout << "argv[15] (Bandwidth): " << argv[15] << endl;//bandwidth
            bandwidth = stoi(argv[15]);
        }
        if(argc > 16){
            cout << "argv[16] (Lower bound ratio): " << argv[16] << endl;//lower bound ratio
            lowerB = stod(argv[16]);
        }
        if(argc > 17){
            cout << "argv[17] (Upper bound ratio): " << argv[17] << endl;//upper bound ratio
            upperB = stod(argv[17]);
        }

        if(argc > 18){
            cout << "argv[18] (Preprocessing Task): " << argv[18] << endl;//preprocessing task
            preTask = stoi(argv[18]);
        }

    }

	//used for running time calculation
    Timer tt0;
    tt0.start();

    string sourcePath=DesFile+"/"+dataset+"/";
    string ODfile=sourcePath+dataset+".query";
    string updateFile=sourcePath+dataset+".update";



    Graph g;
    g.threadnum=threadNum;//thread number of parallel computation (can be changed)
    g.sourcePath=sourcePath;
    g.ifParallel = true;
    g.dataset=dataset;
    g.algoChoice=algoChoice;
    g.algoQuery=algoQuery;
    g.algoUpdate=algoQuery;
    g.algoParti=algoParti;
    g.partiNum=partitionNum;
    g.bandWidth=bandwidth;
    g.bRatioLower=lowerB;
    g.bRatioUpper=upperB;
    g.samePartiPortion=portion;
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"System Index: "<<algoChoice<<endl;
    if(g.algoQuery==Dijk){
        cout<<"Dijkstra's test !!!!!!!"<<endl;
    }else if(g.algoQuery==PH2H_No){
        cout<<"This is test for No-boundary strategy !!!!!!!"<<endl;
    }else if(g.algoQuery==PH2H_Post){
        cout<<"This is test for Post-boundary strategy !!!!!!!"<<endl;
    }else if(g.algoQuery==PH2H_Cross){
        cout<<"This is test for Cross-boundary strategy !!!!!!!"<<endl;
    }else if(g.algoQuery==PCH_No){
        if(g.algoChoice==2){
            cout<<"Wrong index choice."<<endl; exit(1);
        }
        cout<<"This is test for PCH-No !!!!!!!"<<endl;
    }else{
        cout<<"Wrong query strategy! "<<g.algoQuery<<endl; exit(1);
    }
    cout<<"Partition method: "<<g.algoParti<<endl;
    cout<<"Partition number: "<<g.partiNum<<endl;
    cout<<"Thread number: "<<g.threadnum<<endl;
    cout<<"Batch size: "<<batchSize<<endl;
    if(workerNum>threadNum){
        workerNum=threadNum;
    }
    cout<<"Query worker number: "<<workerNum<<endl;

    if(preTask==1){
//        g.PH2HVertexOrdering(0);//MDE ordering
//        g.PH2HVertexOrdering(1);//Boundary-first ordering
        g.PH2HVertexOrdering(2);//Boundary-first MDE ordering
    }
    else if(preTask==2){
        g.QueryGenerationParti(true);//same partition and real-world simulation
    }

//    g.ReadGraph(graphfile);//
//    g.StainingMethod(0);

    ///Task 1: Index construction
    g.IndexConstruction();
//    g.PH2HIndexConstruct();
//    g.WriteCTIndex(graphfile);

    ///Task 2: Query processing
//    g.CorrectnessCheck(100);
//    g.EffiCheck(ODfile,runtimes);//query efficiency test
//    g.EffiCheck(ODfile+"Parti",runtimes);//query efficiency test
//    g.EffiCheck(ODfile+"SameParti",runtimes);
//    g.EffiCheck(ODfile+"CrossParti",runtimes);
//    exit(0);
//    int id1=205055, id2=179964;
////        cout<<"Before update "<<i<<". "<<id1<<"("<<NodeOrder[id1]<<","<<PartiTags[id1].first<<") "<<id2<<"("<<NodeOrder[id2]<<","<<PartiTags[id2].first<<") "<<QueryPostMHL(id1,id2)<<"("<<Dijkstra(id1,id2,Neighbor)<<")"<<endl;
//    cout<<"Before all. "<<id1<<"("<<g.NodeOrder[id1]<<") "<<id2<<"("<<g.NodeOrder[id2]<<") "<<endl;
//    for(auto it=g.Tree[g.rank[id1]].vert.begin();it!=g.Tree[g.rank[id1]].vert.end();++it){
//        if(it->first==id2){
//            cout<<"CH index. "<<id1<<" "<<id2<<" "<<it->second.first<<" "<<it->second.second<<endl;
//            break;
//        }
//    }
//    g.PostMHLIndexCompareCH(g.sourcePath+dataset+".CHIndex19");
//    g.MHLIndexCompareCH(g.sourcePath+dataset+".CHIndex0");
//    exit(0);
    ///Task 3: Index update
//    if(dataset=="beijing" || dataset=="Guangdong" || dataset=="Test2"){//real-life updates
//        g.RealUpdateThroughputTest(sourcePath+dataset+"_20160105.updates");
//    }else{
        g.RandomUpdateThroughputTestQueueModel(batchNum, batchSize, batchInterval, T_r, workerNum);
//        g.RandomUpdateThroughputTest(sourcePath+dataset+".update", batchNum, batchSize, batchInterval);
//        g.SPThroughputTest(updateType, ifBatch, batchNum, batchSize, batchInterval, runtimes);
//    g.IndexMaintenance(updateType,updateSize, ifBatch, batchSize);//index maintenance
//    g.IndexMaiEachNntenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance
//    }


    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
//    exit(0);
    g.clear();
	return 0;
}
