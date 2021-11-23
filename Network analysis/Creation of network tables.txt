Creation of network tables:

Step 1: Use Qiime2 to convert qza to json-type biom file. Details can be found in “spieceasi pipeline”

Step 2. Import json files into R and use SPIECEASI package. The txt file “SPIECEASI_workflow” shows how I did single matrices in R. This starts with a simple unrarefied biom file. You then go through some filtering steps that removes low abundance OUT and creates a “dummy” row so that sums are still kept. The spieceasi function then performs a log transformation and creates a matrix based on inverse covariance. Igraph then turns this network into a “node” and “edge” graph. We remove negative interactions and the final result is a write.graph file that can be imported into Cytoscape. The R file “combining matrices” does the same thing but includes creating two filtered phyloseq objects and a listed spieceasi function.

Step 3. Importing into Cytoscape. Before importing you must take the write.graph files and remove all of the “dummy” rows (should have values of “0” instead of OUT ids) and replace any space with a “,”. When importing the network in Cytoscape go to File>Import Network>from file and choose your files. The first column should be source node, the next connection node, and the final one edge attribute. Visit the Raes Lab for more information on using Cytoscape.

Step 4. Network Analysis. Simply chose the network, went to Tools>Network Analysis and selected it as an undirected network. 
