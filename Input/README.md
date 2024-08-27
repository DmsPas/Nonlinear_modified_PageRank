**NPR_datasets**: adjacency matrices from synthetic and real-world instances, intended for 
                  clustering and classification tasks. All matrices are in Matlab (.mat)
                  format.

The adjacency matrices in this folder were used in the numerical experiments in:
 
Costy Kodsi and Dimosthenis  Pasadakis,  "Nonlinear modified PageRank problem
for local graph partitioning," 2024.  

The full list of datasets used can be found at [permalink](https://drive.switch.ch/index.php/s/PEnKOcOYEWUILap).

### Graph creation

For all the graphs in the folder besides the LFR instances, for $n$ data points,
the connectivity matrix $G \in \mathbb{R}^{n\times n}$ is created from a k nearest
neighbors routine, with k set such that the resulting graph is connected. The
similarity matrix $S \in \mathbb{R}^{n\times n}$ between the data points is defined
as

$$
    s_{ij} = \mathrm{max} \{s_i(j), s_j(i)\} \;\; \text{with}\;
    s_i(j) = \mathrm{exp} (-4 \frac{\|x_i - x_j \|^2}{\sigma_i^2} )
$$

with $\sigma_i$ standing for the Euclidean distance between the $i$th data point
and its nearest k-nearest neighbor. The adjacency matrix $W$ is then created
as

$$
    W = G \odot S.
$$

Besides the adjacency matrices $W$, the node labels for each graph are part of
the submission.  If the graph has $c$ classes, the node labels are integers in
the range $0$ to $c-1$.

For a more detailed description of the datasets, and references of the sources
the datasets were obtained from, see the paper cited in the beginning of this
description.

### Synthetic datasets

| Name             | Nodes | Edges  | Classes |
|------------------|-------|--------|---------|
| LFR_10           | 1000  | 5062   | 19      |
| LFR_12           | 1000  | 5162   | 18      |
| LFR_14           | 1000  | 5323   | 21      |
| LFR_16           | 1000  | 5095   | 20      |
| LFR_18           | 1000  | 4796   | 18      |
| LFR_20           | 1000  | 4883   | 18      |
| LFR_22           | 1000  | 4960   | 18      |
| LFR_24           | 1000  | 5243   | 21      |
| LFR_26           | 1000  | 4887   | 17      |
| LFR_28           | 1000  | 5236   | 20      |
| LFR_30           | 1000  | 5229   | 20      |
| LFR_32           | 1000  | 4999   | 17      |
| LFR_34           | 1000  | 5002   | 17      |
| LFR_36           | 1000  | 5201   | 20      |
| LFR_38           | 1000  | 5179   | 19      |
| LFR_40           | 1000  | 4882   | 20      |
|------------------|-------|--------|---------|
| Gauss2_55_10NN   | 800   | 4902   | 2       |
| Gauss5_55_10NN   | 2000  | 12153  | 5       |
| Gauss8_55_10NN   | 3200  | 19319  | 8       |
| Gauss13_55_10NN  | 5200  | 31204  | 13      |
| Gauss18_55_10NN  | 7200  | 43013  | 18      |
| Gauss25_55_10NN  | 10000 | 59644  | 25      |
| Gauss32_55_10NN  | 12800 | 76044  | 32      |
| Gauss41_55_10NN  | 16400 | 97112  | 41      |
| Gauss50_55_10NN  | 20000 | 118348 | 50      |
| Gauss61_55_10NN  | 24400 | 144052 | 61      |


-------------------------------------------------------------------------------
Real-world datasets
-------------------------------------------------------------------------------
    
- Graphs from image classification problems. The RGB values are normalised 
   in the interval [0,1].

| Name              | Nodes | Edges  | Classes |
|-------------------|-------|--------|---------|
| MNIST_test_10NN   | 10000 | 72800  | 10      |
| USPS_10NN         | 11000 | 156900 | 10      |
| FMNIST_test_10NN  | 10000 | 79152  | 10      |


- Graphs from ORBIS: The Stanford geospatial network model of the Roman world, 
   [ORBIS](https://doi.org/10.2139/ssrn.2609654).

| Name        | Nodes | Edges | Classes |
|-------------|-------|-------|---------|
| Orbis_Km    | 677   | 1104  | 47      |
| Orbis_Days  | 677   | 1104  | 47      |
| Orbis_Cost  | 677   | 1104  | 47      |
