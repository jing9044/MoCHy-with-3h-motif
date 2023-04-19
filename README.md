# MoCHy-with-3h-motif
Source code for the paper Hypergraph motifs and their extensions beyond binary, Geon Lee, Seokbum Yoon, Jihoon Ko, Hyunju Kim, Kijung Shin, VLDBJ.

This work is an extended version of [Hypergraph Motifs: Concepts, Algorithms, and Discoveries](https://arxiv.org/abs/2003.01853), [VLDB 2020](https://vldb2020.org/).

We additionaly propose **Ternary Hypergraph Motifs (3h-motifs)**, which are extensions of h-motifs beyond binary.  

While h-motifs consider only on the emptiness of seven subsets derived from intersections among three hyperedges, 3h-motifs further differntiate patterns based on the cardinality of these subsets. 

We add only source codes for counting 3h-motifs to https://github.com/geon0325/MoCHy. Please refer this link for information of h-motifs and MoCHy.

## Datasets
* The sample dataset is available [here](https://gist.github.com/pszufe/02666497d2c138d1b2de5b7f67784d2b#sec_dblp).
* The real-world datasets used in the paper are available [here](https://www.cs.cornell.edu/~arb/data/) or [here](http://dmlab.kaist.ac.kr/hmotif/).
* In the paper, we used datasets with unique hyperedges, where duplicated hyperedges are removed. 

## Input & Output Format
* The input format should be lines of hyperedges, where each line represents the nodes contained in each hyperedge.
* The index of the nodes should start from 0.
* For example, with 3 hyperedges: {0, 1, 2}, {2, 3}, and {1, 3, 4, 5}, the input file should be:
```
0,1,2
2,3
1,3,4,5
```
* If you type "none" as argument then the output of the code will be counts of h-motifs:
```
h-motif 1: 123
h-motif 2: 22
...
h-motif 26: 31
```
* If you type "ab1" as argument then the output of the code will be counts of 3h-motifs:
```
3h-motif 1: 12
3h-motif 2: 7
...
3h-motif 431: 2
```

## Running Demo
You can run demo with the sample dataset (dblp_graph.txt).
1. To run **MoCHy-E**, type 'run_exact.sh'.
2. To run *parallelized* **MoCHy-E**, type 'run_exact_par.sh'.
3. To run **MoCHy-A**, type 'run_approx_ver1.sh'.
4. To run **MoCHy-A+**, type 'run_approx_ver2.sh'.
5. To run *parallelized* **MoCHy-A+**, type 'run_approx_ver2_par.sh'.
6. To run *memory-bounded* **MoCHy-A+**, type 'run_approx_ver2_memory.sh'.

## Terms and Conditions
If you use this code as part of any published research, please acknowledge our VLDBJ paper.
```
@article{,
  title={Hypergraph motifs and their extensions beyond binary},
  author={Lee, Geon and Yoon, Seokbum and Ko, Jihoon and Kim, Hyunju and Shin, Kijung},
  journal={},
  year={2023},
  publisher={}
}
```

## Contact Information
If you have any questions, please contact [Seokbum Yoon](jing9044@kaist.ac.kr).
