# Dijkstra’s Algorithm with Predictions

This repository contains the code which was used for the experiments described in the paper:

> “Dijkstra’s Algorithm with Predictions to Solve the Single-Source Many-Targets Shortest-Path Problem”   
> by Willem Feijen and Guido Schäfer  
> arXiv link: [https://arxiv.org/abs/2112.11927](https://arxiv.org/abs/2112.11927)

We are grateful to Ruben Brokkelkamp for creating an initial implementation of our algorithm in C++, which we used as a starting point to build the rest of our code around.

## How to use

Download the code from github by cloning the repository:

    git clone https://github.com/w-feijen/dijkstra-predictions-for-SSMTSP.git

Below, we use `<local-path>` to refer to the (absolute) path of the directory containing the code.

### Docker

We use [Docker](https://www.docker.com) to package up our code and all its dependencies to make it self-contained and platform independent. Please make sure that you have Docker installed (e.g., Docker Desktop) and that it is running in the background.

1.  In the directory `docker`, edit the script `create_docker.sh` to adjust the path to your `<local-path>` accordingly:
    
        docker run -i -t -v <local-path>:/home/learning-in-dijkstra dijkstra /bin/bash

2.  Make sure that Docker is installed and the Docker daemon is running (e.g., by running Docker Desktop in the background).

3.  Build the Docker image by running the `build_docker.sh` script. This may take some time.

4.  Run the Docker image by executing the script `create_docker.sh`.

This will create a new environment in your terminal containing the following sub-directories:

    root@50cbb3aeb5ff:/home# ls
    FunctionalPlus  eigen  frugally-deep  json  learning-in-dijkstra

Now, change to directory containing the c++ source code:

    root@50cbb3aeb5ff:/home# cd learning-in-dijkstra/learning-in-dijkstra/src/

### Source Code

There are three main programs that can be used for different testing purposes.

- `random`
- `testset`
- `fortunate`

You can compile them all at once by

    make all
    
Or compile them separately by one of the following

    make random
    make testset
    make fortunate

The main files differ in the graph instances they test the shortest path algorithms on.
Each main file runs four different algorithms to solve the SSMTSP problem:

1.  Dijkstra: plain Dijkstra algorithm adapted to the many-targets setting
2.  Dijkstra pruning: adapted Dijksta algorithm pruning all edges exceeding bound B
3.  Dijkstra oracle: same as in 2, but pruning all edges exceeding the exact shortest path distance D
4.  Dijkstra prediction: adapted Dijkstra combining pruning with prediction of shortest path distance
    For Dijkstra prediction, different methods are implemented to obtain a prediction
     4.a The oracle prediction (ORPR), an oracle gives the actual shortest path value at start of algorithm
     4.b The ML prediction (PRED), obtained by a trained neural network in 10th iteration
     4.c The BFS prediction (BFS)
     4.d The weighted BFS prediction (BFS-w)

A summary of the number of queue operations and running times is given at the end. Note that the report of CPU running time requires to set the `-DOUTPUT_TIMING` flag when compiling (uncomment corresponding line in `Makefile`).

In the Makefile, on can also choose between a priority queue based on a 'Fibonacci' or a 'Binomial' heap, by commenting out the correct line and set the HEAP flag either to `-DFIBONACCI` or `-DBINOMIAL`.

#### Program `random`

    make random

`random` runs the above algorithms (1,2,3,4.a) on random instances (as defined in the paper). The program can be called with the following arguments: 

    ./random n t [c] [q] [alpha] [beta]

Here n and t specify the number of nodes and trials, respectively. Optional the average degree c, target-node probability q, and parameters alpha and beta (as defined in the paper) can be specified.

#### Program `testset`

    make testset 
    
`testset` runs the above algorithms (1,2,3,4.a,4.b,4.c,4.d) on the set of test instances (as defined in the paper). The program can be called with the following arguments:

    ./random t [alpha] [beta]
    
Here t specifies the number of trials of the test set are evalueated. t should not exceed 100, which is the number of instances in the test set which are shared in this repository. The other 9.90O graphs in the test set can be shared upon request, please send an e-mail to Willem. Optional parameters alpha and beta (as defined in the paper) can be specified. 
    
#### Program `fortunate`

    make fortunate

`fortunate` runs the above algorithms (1,2,3,4.a) on so-called fortunate instances. We created these instances to test how (and how fast) our algorithm works with many Insert and Decrease Priority operations. The program can be called with the following arguments:

    ./random n t [r] [alpha] [beta]
    
Here n and t specify the number of nodes and trials, respectively.  Optional parameters alpha and beta (as defined in the paper) can be specified. The other optional parameter is r, which is an input parameter for creation of the so-called fortunate instance. 

A fortunate graph instance is defined as follows. Two input parameters must be given: the number of nodes, $n\in \mathbb{N}$, and the fraction of nodes which is on the shortest path, $r\in \mathbb{R}\_+ $. From these parameters, the number of nodes on the shortest path, $x$, can be derived: $x:= \lfloor r\cdot n\rfloor$. We denote these nodes on the shortest path as $u_i$, for $i=0,\ldots, x-1$. Furthermore, we label $u_0$ as the start $s$ and label $u_{x-1}$ as the target $t$. The rest of the $n-x$ nodes in the graph are denoted by $v_j$, for $j=1, 2, \ldots, n-x$. Between these $n$ nodes, there are $(x-1)(1+n-x)$ edges in the fortunate graph. There is an edge $(u_i, u_{i+1})$ for $i=0,1, \dots, x-2$, each with weight $\frac{1}{x-1}$. Together, these edges form the shortest path from $s$ to $t$, counting up to a shortest path length of exactly 1. Furthermore, there are edges $(u_i, v_j)$ for $i=0,1, \ldots, x-2$ and $j=1,2,\dots,n-x$. Each edge from $u_i$ to $v_j$ has length $2- \frac{2i}{x-1}$.

#### Clean-up

You can remove all object and/or executable files by calling: 

    make clean
    make clean-all
    
### Create NN Model Notebook
The traces which were used to build the neural network model are in the `\b_files` folder. The notebook `\python\create_nn_model.ipynb` shows an example of how these traces can be used to train a neural network model. It outputs a `model.h5` file. This should be converted into a file which frugally-deep can read. To do so, run the `convert_model.py` model with the `model.h5` file as input. The `convert_model.py` method can be found in the fdeep folder after building the docker container. See https://github.com/Dobiasd/frugally-deep for more explanation. The `convert_model.py` outputs a `fdeep_model.json` file which can be read by fdeep in our C++ code.

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
