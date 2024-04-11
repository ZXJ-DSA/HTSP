/*
 * main.cpp
 *
 *  Created on: 14 Oct 2022
 *      Author: Xinjie ZHOU
 */
#include "headH2H.h"

int main(int argc, char** argv){

    if( argc < 3 || argc > 11){//
        printf("usage:\n<arg1> source path, e.g. /export/project/xzhouby\n");
        printf("<arg2> name of dataset, e.g. NY\n");
        printf("<arg3> algorithm, 0:Dijkstra; 1: CH; 2:H2H. e.g. 1\n");
        printf("<arg4> (optional) update type, (0: No Update Test; 1: Decrease; 2: Increase), default: 0\n");
        printf("<arg5> (optional) whether batch update, (0: No; 1: Yes), default: 0\n");
        printf("<arg6> (optional) batch number, default: 10\n");
        printf("<arg7> (optional) batch size, default: 10\n");
        printf("<arg8> (optional) batch interval (in seconds), default: 10\n");
        printf("<arg9> (optional) thread number, default: 15\n");
        printf("<arg10> (optional) query file name\n");
        exit(0);
    }

    string DesFile="./data/";
    string dataset = "NY";

    int algorithm = 2;
    int updateType = 0;
    bool ifBatch = false;
    int runtimes = 10000;
//    runtimes=100;
    int batchNum = 10;
    int batchSize = 10;
    int batchInterval = 60;
    int threadNum = 15;
    string queryFName;

    if(argc > 1) {
        cout << "argc: " << argc << endl;
        cout << "argv[1] (Source Path): " << argv[1] << endl;//source path
        DesFile = argv[1];

        cout << "argv[2] (Dataset): " << argv[2] << endl;//dataset
        dataset = argv[2];

        cout << "argv[3] (Algorithm): " << argv[3] << endl;//algorithm
        algorithm = stoi(argv[3]);

        if(argc > 4){
            cout << "argv[4] (Update Type): " << argv[4] << endl;//update type
            updateType = stoi(argv[4]);
        }
        if(argc > 5){
            cout << "argv[5] (Batch Update): " << argv[5] << endl;//batch update
            ifBatch = stoi(argv[5]);
        }
        if(argc > 6){
            cout << "argv[6] (Batch Number): " << argv[6] << endl;//batch Number
            batchNum = stoi(argv[6]);
        }
        if(argc > 7){
            cout << "argv[7] (Batch Size): " << argv[7] << endl;//batch size
            batchSize = stoi(argv[7]);
        }
        if(argc > 8){
            cout << "argv[8] (Batch Interval): " << argv[8] << endl;//batch interval
            batchInterval = stoi(argv[8]);
        }
        if(argc > 9){
            cout << "argv[9] (Thread Number): " << argv[9] << endl;//thread number
            threadNum = stoi(argv[9]);
        }
        if(argc > 10){
            cout << "argv[10] (Same-parti Proportion): " << argv[10] << endl;
            queryFName = argv[10];
        }
    }


	//used for running time calculation
    Timer tt0;
    tt0.start();

    string sourcePath=DesFile+"/"+dataset+"/";
    string ODfile=sourcePath+dataset+".query";
    string updateFile=sourcePath+dataset+".update";

    Graph g;
    g.dataset=dataset;
    g.threadnum=threadNum;
    g.sourcePath=sourcePath;
    g.algoIndex=algorithm;
    g.queryFName=queryFName;
    cout<<"Dataset: "<<dataset<<endl;
    cout<<"Thread number: "<<threadNum<<endl;
    if(ifBatch){
        cout<<"Test for batch update! Batch size: "<<batchSize<<endl;
    }else{
        cout<<"Test for single edge update!"<<endl;
        batchSize=1;
    }


    g.ReadGraph(sourcePath+dataset);//
//    g.StainingMethod(0);

    ///Task 1: Index construction
    g.IndexConstruction(algorithm);
//    g.H2HIndexConstruct();
//    g.WriteCTIndex(graphfile);

    ///Task 2: Query processing
//    g.CorrectnessCheckCore(100);
    g.CorrectnessCheck(100);
    g.EffiCheck(runtimes);

//    g.SameTreeQueryTest(ODfile,runtimes);
//    exit(0);

    ///Task 3: Index update
    g.SPThroughputTest(updateType, ifBatch, batchNum, batchSize, batchInterval, runtimes);
//    g.IndexMaintenanceCHWP(updateType, updateSize, ifBatch, batchSize);//index maintenance
//    g.IndexMaintenanceH2H(updateType, updateSize, ifBatch, batchSize);//index maintenance
//    g.IndexMaintenance(updateFile+"ST",updateType,updateBatch);//same-tree index maintenance

    tt0.stop();
    cout<<"\nOverall runtime: "<<tt0.GetRuntime()<<" s."<<endl;
    cout<<"------------------\n"<<endl;
//    exit(0);
    g.clear();
	return 0;
}
