# HTSP
The implementation code of paper "High Throughput Shortest Path Query Processing for Large Dynamic Road Networks" (submitted to SIGMOD25). Please refer to the paper for the algorithm details. 

The full-version paper with an appendix section is also provided for your reference.

## Algorithms

The implementation code includes the index construction, query processing, and index update of our *MHL*, *PMHL*, and *PostMHL* algorithms. 


## Data
The datasets of this paper are sourced from [The 9th DIMACS Implementation Challenge - Shortest Paths](http://users.diag.uniroma1.it/challenge9/download.shtml). 

An example graph *NY* is provided in the directory *NY* for your reference. You can run our algorithms on the example graph by using the source path `./`. 


## Dependency

1. `g++` and `boost`

All the codes are runnable after `cmake` and `make`: go to the corresponding directory, `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.

