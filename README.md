# HOLaGraF

Name stands for **H**igher **O**rder **La**placian: **Gra**dient **F**low.

### TODO

- [ ] change eigenvector algorithm in `np.linal.eigh` to match the accuracy and correct exit in the matlab version
- [ ] fix `hom2.ipynb` random weights tester (and nice plotter) with the theoretical gradient instead of the numerical
- [ ] add some kind of the `gridSearch` for the `thrs["mu"]` parameter
- [ ] write an adaptive strategy for the `thrs["aplha"]`
- [ ] add the free-gradient-flow transition on the epsilon change
- [x] ~~write theoretical gradient~~
  - [x]  ~~test with random weights and perturbations~~
- [x] ~~step increase only after the second consequent "accept"~~
  

### Project structure for the matlab version

```
.
+-- eps_flow4.m (main file, epsilon flow)
|   +-- alpha_flow2.m (subflow by alpha)
|   |   +-- single_flow.m (flow for a given params)
|   +-- getAnswerFromBack.m (extracts minimum and distance)
+-- Graph building
|   +-- readEdges.m (reads list of edges from file)
|   +-- readTrigs.m (reads list of trigs from file)
|   +-- B1fromEdges.m (extract B1 from edges)
|   +-- B2fromTrig.m (extract B2 from trigs and edges)
+-- Graph's matrices
|   +-- getAdjB1.m (adjacency matrix of the graph)
|   +-- getAdjWB1.m (weighted adjacency matrix of the graph)
|   +-- getL0.m (normalized classical Laplacian)
|   +-- HodgeLW_fr.m (1-order Hodge Laplacian)
|   |   +-- getDt.m (triangle weights extractor)
|   |   +-- getM.m (positions of minimal weights for triangles)
+-- Functional
|   +-- getFk_l2.m (main file for the functional)
|   +-- getFk1.m (first term, eigs of L1)
|   +-- getFk2.m (second term, connectedness)
|   +-- getFk3.m (third term, versor penalisation, deprecated)
+-- Gradient
|   +-- getGradNum.m (numerical gradient)
|   +-- getGrad.m (theoretical gradient)
+-- Plotters
|   +-- simpleDrawB1.m (simple graph plotter)
+-- Testers
|   +-- testL0.m (checks adequecy of the normalized L0)
|   +-- test_Grad.m (compares numertical and theoretical gradient)
+-- Support
|   +-- Sym.m (the usual sym function)
|   +-- getEigSort.m (get eigenvector and eigenvalues and sort them)

```

### Project structure for the python version
- `hom.py` main file with all the functions
  - `Graph` class implemented!
- `hom2.ipynb` random flow tester
- `test.grad.py` tester for the theoretical gradient