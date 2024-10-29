# China Datasets
We make three datasets (Guangdong, SouthChina, and EastChina) public, which are obtained from [NaviInfo](https://www.navinfo.com/en).

Each dataset contains two compressed files about the edges and nodes. For example, Guandong.zip stores the edge information while Guangdong.co.zip stores the node information. 


## Edge File
* For edge file, the first line shows the node number and edge number of this graph. 
* From the second line, each line stores one edge with the format `vertex_ID1, vertex_ID2, edge_weight`, where the edge weight is the travel time of this edge.

## Node File
* For node file, the first line shows the node number of this graph. 
* From the second line, each line stores the information of one node with the format `vertex_ID, longitudinal*10000000, latitudinal*10000000`.
