/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "head.h"

int main(int argc, char** argv){

    if( argc < 4 || argc > 17){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> HTSP system index, 1: CH+H2H; 2: PH2H; 3: PCH+PH2H; 4: Optimized PCH+PH2H; 5: PostMHL. default: 2\n");
        printf("<arg4> (optional) partition number, e.g. 64\n");
        printf("<arg5> (optional) partition method, (NC: PUNCH; MT: METIS), default: NC\n");
        printf("<arg6> (optional) query strategy, (0:A*; 1: PCH (CH); 2: No-boundary; 3: Post-boundary; 4: Cross-boundary (H2H)), default: 4\n");
        printf("<arg7> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase), default: 0\n");
        printf("<arg8> (optional) whether batch update, (0: No; 1: Yes), default: 0\n");
        printf("<arg9> (optional) batch number, default: 10\n");
        printf("<arg10> (optional) batch size, default: 10\n");
        printf("<arg11> (optional) batch interval (in seconds), default: 10\n");
        printf("<arg12> (optional) thread number, default: 15\n");
        printf("<arg13> (optional) bandwidth, default: 50\n");
        printf("<arg14> (optional) lower bound ratio, default: 0.1\n");
        printf("<arg15> (optional) upper bound ratio, default: 2\n");
        printf("<arg16> (optional) preprocessing task (1: Partitioned MDE Ordering; 2: Partitioned Query Generation)\n");
        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";
    int algoChoice = 1;
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
    bool ifBatch = false;
    int threadNum = 15;
    int preTask=0;
    int bandwidth=50;
    double lowerB=0.1;
    double upperB=2;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1] (Source Path): " << argv[1] << endl;//source path
        DesFile = argv[1];

        cout << "argv[2] (Dataset): " << argv[2] << endl;//dataset
        dataset = argv[2];

        cout << "argv[3] (System Index): " << argv[3] << endl;//system index
        algoChoice = stoi(argv[3]);

        if(argc > 4){
            cout << "argv[4] (Partition Number): " << argv[4] << endl;//partition number
            partitionNum = stoi(argv[4]);
        }

        if(argc > 5){
            cout << "argv[5] (Partition Method): " << argv[5] << endl;//partition method
            algoParti = argv[5];
            if(algoParti != "NC" && algoParti != "MT"){
                cout<<"Wrong partition method! "<<algoParti<<endl; exit(1);
            }
        }
        if(argc > 6){
            cout << "argv[6] (Query Strategy): " << argv[6] << endl;//algorithm for query
            algoQuery = stoi(argv[6]);
        }

        if(argc > 7){
            cout << "argv[7] (Update Type): " << argv[7] << endl;//update type
            updateType = stoi(argv[7]);
        }
        if(argc > 8){
            cout << "argv[8] (Whether Batch Update): " << argv[8] << endl;//algorithm for update
            ifBatch = stoi(argv[8]);
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
            cout << "argv[13] (Bandwidth): " << argv[13] << endl;//bandwidth
            bandwidth = stoi(argv[13]);
        }
        if(argc > 14){
            cout << "argv[14] (Lower bound ratio): " << argv[14] << endl;//lower bound ratio
            lowerB = stod(argv[14]);
        }
        if(argc > 15){
            cout << "argv[15] (Upper bound ratio): " << argv[15] << endl;//upper bound ratio
            upperB = stod(argv[15]);
        }
        if(argc > 16){
            cout << "argv[16] (Preprocessing Task): " << argv[16] << endl;//preprocessing task
            preTask = stoi(argv[16]);
        }

    }

	//used for running time calculation
    Timer tt0;
    tt0.start();

//    string graphfile="/media/TraminerData/mengxuan/MengxuanGraphWPSL/Cond/CondWeighted";
    string graphfile=DesFile+"/"+dataset+"/"+dataset;
    string ODfile=graphfile+".query";
    string updateFile=graphfile+".update";



    Graph g;
    g.threadnum=threadNum;//thread number of parallel computation (can be changed)
    g.graphfile=graphfile;
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
    if(ifBatch){
        cout<<"Test for batch update! Batch size: "<<batchSize<<endl;

    }else{
        cout<<"Test for single-edge update!"<<endl;
        batchSize=1;
    }

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
    g.CorrectnessCheck(100);
    g.EffiCheck(ODfile,runtimes);//query efficiency test
//    g.EffiCheck(ODfile+"Parti",runtimes);//query efficiency test
//    g.EffiCheck(ODfile+"SameParti",runtimes);
//    g.EffiCheck(ODfile+"CrossParti",runtimes);
//    exit(0);
    ///Task 3: Index update
    g.SPThroughputTest(updateType, ifBatch, batchNum, batchSize, batchInterval, runtimes);
//    g.IndexMaintenance(updateType,updateSize, ifBatch, batchSize);//index maintenance
//    g.IndexMaiEachNntenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance

    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
//    exit(0);
    g.clear();
	return 0;
}
