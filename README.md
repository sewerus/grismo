# Grismo

Graph Isomorphism Problem Solver by Seweryn Panek Wrocław University of Science and Technology.

## Table of contents
* [Gismo](#grismo)
	* [Implemented algorithms](#implemented-algorithms)
	* [Types of tests](#types-of-tests)
	* [Setup](#setup)
	* [Usage](#usage)
* [Project](#project)
	* [Introduction](#introduction)
	* [Research problem](#research-problem)
	* [Description of the methods](#description-of-the-methods)
	* [Research plan and simulation program](#research-plan-and-simulation-program)
	* [Conducted research](#conducted-research)
	* [The results of conducted research](#the-results-of-conducted-research)
	* [Summary](#summary)
	* [References](#references)

To read this description with math expressions, please, see README.ipynb. Github doesn't render LATEX code in README.md on the projects' main pages.

## Implemented algorithms
-   Brute method,
-   Method based on BFS algorithm,
-   Method based on DFS algorithm,
-   Weisfeiler-Lehman's method,
-   Method based on graph spectrum comparison.

## Types of tests
-   single test for given graphs,
-   many tests for an increasing number of vertices,
-   many tests for an increasing number of edges.

## Setup
Clone this repo to your desktop and install all the dependencies.

## Usage
Just run command `python main.py` to run an application. Then choose in the right column type of test you want to conduct and it's parameters. If you want start a symulation click button "Przeprowadź badanie".

# Project

The aim of the project is to study already known and new created methods to solving problem of graph isomorphism. For this purpose algorithms have been implemented and a simulation tool has been prepared to compare these methods. The simulation program allows to study the impact of parameters such as the number of vertices and the number of edges of the graphs. The methods were compared in terms of their reliability and calculation time. Based on the collected results, it was possible to draw conclusions and make further hypotheses.

## Introduction

A graph is a mathematical structure consisting of a set of vertices and a set of edges having ends at two selected vertices. It is used to present and describe relationships between objects.

A graph can be presented by set of vertices *V* and set of edges *E*, where *E = {{x, y}: x, y in V}*. Each edge is a pair of vertices which are neighbours. In labeled graph all vertices have unique number.

In this article all discussed graphs are simple, connected and symmetric. It means, that there are no multiple edges and no own loops, between any two vertices there is a path and if vertex *v_1* is a neighbour of vertex *v_2*, then vertex *v_2* is the neighbour of vertex *v_1*.

A graph can be represented by an adjacency matrix and an adjacency list. Adjacency matrices for all discussed graphs are square (0, 1)-matrices with zeros on diagonals. If in *i*-row and *j* column of an adjacency matrix there is 1, then vertices *v_i* and *v_j* are neighbours. If there is 0, then these vertices are not neighbours. An adjacency list represent a graph by a list of neighbours' labels for each vertex.

![Example of a graph with 6 vertices.](./readme_images/graph_1.png)

![Adjacency matrix of the graph](./readme_images/matrix_1.png)

![Adjacency list of the graph](./readme_images/list_1.png)

## Research problem

There is an isomorphism between two graphs when the vertices of one of them can be relabeled in such a way that the vertices in both graphs have exactly the same neighbors. An example of two isomorphic graphs is shown in Fig. below, where each vertex is given a corresponding label and in both graphs each vertex has exactly the same neighbors.

Research whether two graphs are isomorphic is the problem of isomorphism resolution. Graph isomorphism preserves all graph properties, for example: number of vertices, number of edges, and consistency. Therefore, isomorphic graphs are usually identified. In this paper two examined graphs will be called "graph *A*" and "graph *B*".

![Example of two isomorphic graphs.](./readme_images/2_isomorphism.png)

Solving this problem has been used in chemistry, comparing the structure of molecules and atoms and the bonds between them. Also in some anti-plagiarism systems, two texts are presented using graphs and their similarity is checked, e.g. by testing isomorphism.

The problem of isomorphism resolution of two graphs belongs to the NP class, but has not yet been shown to be NP-complete. On the other hand, there are no known deterministic, probabilistic or quantum polynomial algorithms solving this problem. It is also not known whether the problem belongs to the co-NP class (complementary complexity class for NP decision problems).

## Description of the methods

This chapter describes 5 methods for solving the problem of graph isomorphism, which will be investigated later in this work.

### Brutal Method

In this method, for the graph *B*, all possible other graphs resulting from the change of the order of its vertices are generated.

For *n* vertices there is *n!* possible ways to arrange them in order. For each graph resulting from the rearrangements of the vertices of the graph *B*, a comparison is made with the graph *A* by comparing all neighbor relations in the adjacency matrix. The adjacency matrix is a square matrix with dimensions *n \times n*, so comparing two matrices has a square complexity.
Hence the computational complexity of the whole method is *O (n! \cdot n^2)*.

The creation of each vertex permutation is done vertex by vertex, i.e. the first vertex from the *n* unused vertices is drawn first, then the second from the *n-1* unused vertices until all vertices are used.

After each drawn vertex, a subgraph of the original graph is obtained. For example: after drawing only 5 vertices out of 10, we get a subgraph consisting of only 5 vertices.

To improve this method, algorithm compares the subgraphs at each stage. Thanks to this, it is possible to make an earlier decision about the lack of isomorphism for this draw - before the full permutation of vertices is drawn, the checked permutation will be interrupted and the next one starts.

### Method that uses BFS algorithm

Breadth-first search (BFS) is one of the graph searching algorithms. The graph is traversed from a given vertex *s* and consists in visiting all vertices (neighbors) reachable from it. The result of the algorithm is a search tree rooted in *s*, containing all the nodes reachable from *s*.

This algorithm was used to create the isomorphism check method, introducing the concept of "graph transition stage". In the first step, a subgraph is created from the base vertex *s* and all its neighbors. In the second stage, the subgraph is supplemented with vertices' neighbors from the previous stage. The successive stages of the subgraph transition enlarge the subgraph until the original graph is obtained (visiting all vertices). An example of traversing the graph using the BFS algorithm and new vertices at each stage is presented in Fig. below, where the vertex *s* is the vertex labeled as 1.

![Example of graph transition stages.](./readme_images/3_bfs_ex.png)

While the brute force method examined all possible vertex permutations of the graph *B* to see if any of them would produce an identical graph to the graph *A*, the BFS graph search method limits the number of permutations. If the subgraph obtained at the *i*-th transition stage of the graph *B* is isomorphic to the graph obtained at the *i*-th transition stage of the graph *A*, in order to prove the isomorphism in the next stage, it is enough to check all possible permutations of new vertices in the next stage.

This method is not effective for inconsistent graphs because then the BFS search method does not find all vertices.

### Method that uses DFS algorithm

A depth-first search (DFS) is another graph searching algorithm. The traversing of the graph starts from the given vertex s and consists in examining all edges coming from the vertex s. It is a recursive algorithm.

Using this algorithm you can get a set of all possible paths coming from any vertex. For example, Fig. below presents two graphs with all their paths.

![Representation of sample graphs using a list of edges.](./readme_images/3_dfv_example.png)
Research isomorphism can be carried out by trying to assign the edges of the graph *B* representing the edges of the graph representing the graph *A*. At the beginning of this method, an edges list is created for the graph *A* from vertex number 1. Then, for graph *B*, lists of edges coming from each vertex are created. Edges are grouped according to their length. If for a given initial vertex of graph *B* the number of groups and the number of their edges do not correspond to the representation of graph *A*, then such a vertex is immediately rejected.

The graphs in Fig. are represented by 3 groups of edges with the lengths of 5 vertices, 4 vertices and 3 vertices consecutively. The first group has two edges, the others have one edge.

### Weisfeiler-Lehman's algorithm

First, all the vertices of both graphs are colored with the same color (a color can be represented by unique integer number). Then, for each vertex, a multi-set (set with repetitions) of the colors of its neighbors is created. Such pairs are saved for each graph, and if the graphs are isomorphic, identical sets of color pairs and the aforementioned multisets should be produced.

The next step is to recolor the vertices so that the two vertices have the same color only if they previously created the same color and multiset pair.

If the multisets of colors of both graphs are different, then the method determines the lack of isomorphism. Otherwise, the isomorphism has not been resolved.

Therefore, coloring and creating multi-color sets can be repeated in circles. If the obtained multisets of colors are identical after the duplications, then the result confirming the isomorphism of the examined graphs is returned.

It is then the Weisfeiler-Lehman's  algorithm of the *k* dimension. Note that it can state isomorphism when in fact the two graphs are not isomorphic. It has been proved, however, that this method correctly distinguishes all planar graphs for *k = 3*.

### The method of comparing the spectrum

Let *X* be a graph. The eigenvalues and eigenvectors of the adjacency matrix of the graph *X* are called the eigenvalues and eigenvectors of the graph *X*, respectively. By the *Spec X*, spectrum of the graph *X*, we mean the set of its eigenvalues together with information about their multiplicity. The spectrum of graph *X* is written in the form of a matrix, the first row of which contains the roots of the  characteristic polynomial of the graph *X*, and the second row of their multiplicities, respectively.

In this method, the spectrums of both examined graphs are calculated. The eigenvalues in the spectrum are sorted. If the obtained matrices representing the spectrums are identical, the result of this method will be confirmation of the isomorphism of the examined graphs.

## Research plan and simulation program

The chosen programming language is Python 3.7. The development environment used for this project is PyCharm. The program interface was created using the PyQt5 library. Other libraries that provided ready-made mechanisms, including for time measurement and calculation of eigenvalues of matrices are time, math, bisec, numpy. The pyqtgraph library was used to generate the charts.

The simulation program provides for 3 types of tests.

### Single test

In order to check single results for previously determined graphs, the simulation program should allow the insertion of such graphs, e.g. through an adjacency matrix. It is also necessary to be able to draw graphs.

After selecting the methods and validating the test, the results for each method should be summarized in the form of a table, which will include information about the verdict whether the graphs are isomorphic and the calculation time.

For the Weifeiler-Lehman method, it is also necessary to select the dimension of this method.

This test will not check the effect of the number of vertices or edges on the performance of the proposed methods, but it will be helpful when examining graphs entered by the user.

### Investigation of the effect of the number of vertices

In order to determine the dependence, which is influenced by the number of graph vertices on the selected methods, it will be necessary to perform repeated tests for different graphs for the these methods.

To do this, it will be needed to establish a range for the number of vertices and the number of simple tests to run for each number in that range. Moreover, it should be possible to choose whether each of the drawn graphs has a fixed number of edges, or whether the number of edges should change with the number of vertices. It should also be possible to choose whether all graphs drawn for testing are isomorphic or random.

Due to the methods based on searching the graph with the BFS and DFS algorithms, all drawn graphs must be consistent.

An example test should have the following parameters:
* Number of vertices from 5 to 15
* Number of edges variable: 50\% of all possible edges
* All graphs are isomorphic
* 50 simple tests were performed for each number of vertices.

In a case that a given method shows too long computation times compared to others, it will be tested only to as many vertices as it is possible to perform the tests. For a larger number of vertices, those methods will be tested that allow obtaining results in a reasonable time.

### Investigation of the effect of the number of edges

In order to determine the dependence of the selected methods on the number of graph edges, it will be necessary to perform repeated tests for different graphs for these methods.

To do this, it will be necessary to establish the number of vertices and the range of the number of edges for the simple tests that will be run.

An example test should have the following parameters:
* Number of vertices equal to 10
* The number of edges ranges from 11 to 45 (i.e. to the full graph)
* All random graphs (isomorphic or not)
* 100 simple tests were performed for each number of vertices.

## Conducted research

### Investigation of the effect of the number of vertices

The first test was performed for all tested methods for a constant number of edges equal to 16. The range of vertices is from 7 to 15 and for each number of vertices 20 graphs (isomorphic or not) were drawn. The results of this test are presented in Fig. below using a chart of the average execution time of the algorithm on a vertical axis on a logarithmic scale depending on the number of vertices.

In the next test the number of edges will be variable, but the ratio of the graph edges used to all their edges will be constant and equal to 50\%. The selected range of the number of vertices is from 5 to 14. As in the previous test, the drawn graphs will be random (isomorphic or not) and 20 simple tests will be performed for each number of vertices. The results of this test are presented in Fig. below.

In order to better investigate the spectrum comparison method, the Weisfeiler-Lehman algorithm and the BFS algorithm, a larger number of vertices were tested. This test has the same settings as the previous test, but the vertex range is from 5 to 30. The results of this test are presented in Fig. below.

![Results of the first test of the effect of the number of vertices](./readme_images/4_vert_1.png)
![Results of the second test of the effect of the number of vertices](./readme_images/4_vert_2.png)
![Results of the third test of the effect of the number of vertices](./readme_images/4_vert_3.png)

### Investigation of the effect of the number of edges

The first test was performed for the number of vertices equal to 7. The range of the tested edges was set from 6 to 21. The test ran as expected up to an edge number of 17. The eighteenth drawn 17-edge graph was checked by the brutal method and the BFS method, but the single test using the DFS method took more than 85 hours, after that time the test was discontinued. Hence the conclusion that for some graphs the method using the DFS algorithm is at risk of combinatorial explosion resulting from too many possible paths coming from one vertex.

It was decided to perform the test and discuss the results for all methods for the number of edges from 6 to 16, the results are shown in Fig. below. The second test was carried out for the number of edges from 6 to 21, but with the exception of the method using the DFS algorithm, the results are shown in  Fig. below.

It was also decided to test all methods for the number of edges from 9 to 37 and all methods except the method using the DFS algorithm for the number of edges from 9 to 45. The results of the first test are presented in Fig. below. The results of the second test are presented in Fig. below.

In order to better investigate the influence of the number of edges on the graph spectrum comparison method and the Weisfeiler-Lehman's algorithm, one more test was performed only for these methods. Graphs with 20 vertices and the number of edges from 50 to 190 were selected. The results of this test are presented in Fig. below.

![Results of the first test of the effect of the number of edges](./readme_images/4_edg_1.png)
![Results of the second test of the effect of the number of edges](./readme_images/4_edg_2.png)
![Results of the third test of the effect of the number of edges](./readme_images/4_edg_3.png)
![Results of the fourth test of the effect of the number of edges](./readme_images/4_edg_4.png)
![Results of the fifth test of the effect of the number of edges](./readme_images/4_edg_5.png)

## The results of conducted research

### Reliability of the methods

The least reliable method is the method that uses the DFS algorithm; the method that returned an incorrect result only once in several thousand tests performed is the method that compares spectra graphs, the methods that have always returned correct results are:
* brutal method,
* method using the BFS algorithm,
* Weisfeiler-Lehman's algorithm (for dimension 3).

Weisfeiler-Lehman's algorithm of dimension 3 certainly returns correct results for planar graphs [2]. In this project, all examined graphs were consistent, which allows for the hypothesis that dimension 3 of this algorithm allows for obtaining certain results for connected graphs.

Another fact that should be compared with information found in publications [2] is that not all graphs have identical spectrums. So far, it has only been proven that all tree algorithms can be distinguished in terms of isomorphism by comparing their spectra. This project shows that for connected graphs, this method rarely returns an erroneous result.

### Calculation times

Brutal method:
* the computation time of the brutal method depends on the number of vertices in an exponential way,
* for a fixed number of vertices, this method also depends on the number of edges,
* for a small number of edges, the calculation time of this method decreases with the increasing number of edges,
* for almost full graphs, the computation time of this method increases with the increasing number of edges,
* the exception are full graphs, for which this method returns the correct result very quickly, what is a result of the algorithm y and the fact that two complete graphs with the same number of vertices are isomorphic.

Method using the BFS algorithm:
* the calculation time of the method using the BFS algorithm increases with the increase of the number of vertices, but without a clear exponential or linear tendency,
* has very short computation times for graphs with few edges (e.g. up to 16 edges for graphs with 10 vertices), but for a larger number of edges the computation time for this method is extended. The longest computation times are observed for full graphs.

Method using the DFS algorithm:
* the calculation time of this method depends mainly on the number of edges - it grows exponentially,
* for a constant number of edges, the time of this method decreases with an increase in the number of vertices - this is because the ratio of the number of edges of the graph decreases compared to the number of all possible edges,
* for a given number of vertices there is a number of edges for which a combinatorial explosion can be observed - too many possible paths resulting from the DFS algorithm itself causes a very rapid increase in computation time.

Weisfeiler-Lehman's algorithm:
* the time of this method has always been the shortest of all methods,
* its time grows very slowly both with the increase in the number of graph edges and with the increase in the number of their vertices.

Method for comparing graph spectrums:
* the computation time of this method is similar to the computation time of the Weisfeiler-Lehman's algorithm,
* its time increases very slowly as the number of graph vertices increases, but it is independent of the number of edges.

## Summary

Both in terms of reliability and calculation times, the method using the Weisfeiler-Lehman's algorithm turned out to be the best one.

The times of the remaining reliable methods (brutal and using the BFS algorithm) increase too rapidly for large graphs, making them impractical to use for problems where there is a need of resolving large graph isomorphism problem.

An alternative to the Weisfeiler-Lehman's algorithm is the method of comparing graph spectrumsAn alternative to the Weisfeiler-Lehman's algorithm is the method of comparing graph spectrums. For several thousand tests, it only returned an incorrect result once.

For several thousand tests, it only returned an incorrect result once.

The worst method both in terms of computation time and reliability turned out to be the method using the DFS algorithm.

## References

[1]: Brendan L. Douglas, "The Weisfeiler-Lehman Method and Graph Isomorphism Testing"

[2]: Zając K., "Algorytmy kombinatoryczne i graficzne w spektralnej klasyfikacji skończonych bigrafów oraz sieciowych systemów pierwiastków" (doctoral dissertation)

[3]: Bartoszuk M., "System do oceny podobieństwa kodów źródłowych w językach funkcyjnych oparty na metodach uczenia maszynowego i agregacji danych" (doctoral dissertation)

[4]: Kowalik Ł., "Problem izomorfizmu grafów", http://www.deltami.edu.pl/temat/informatyka/algorytmy/2018/11/26/Problem_izomorfizmu_grafow/, accessed at 2020-03-20.

[5]: Institut Teknologi Bandung, Lecture recording with calculating the graph spectrum, https://www.youtube.com/watch?v=BjoPgQPMisQ, accessed at 2020-04-20
