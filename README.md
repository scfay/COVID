# COVID
This is the code that supports the manuscript "Simple control for complex pandemics" by Sarah C. Fay, Dalton J. Jones, Munther A. Dahleh, and A.E. (Peko) Hosoi.

This code to simulate outbreaks of disease to illustrate the trade-offs of social distancing, random testing, and contact tracing. 
Two models of the underlying network structure are used: a branching model and an individual-based model.

The main simulation occurs in main_SimpleControlForComplexPandemics.m.

It makes use of the two simulators for the different types of network models: infectionBranching.m and infectionIBM.m.

infectionIBM.m relies on generating an graph that describes people as nodes and their interactions as edges. The functions that generate these graphs are:
  generateGraphBinom.m  (The number of edges that each node has follows a binomial distribution)
  generateGraphBinom.m  (The number of edges that each node has follows a uniform distribution)
  
The anaylictal results in the manuscript are plotted with literature data in tradeoffs.m.
