# Dual-Channel Hybrid Community Detection in Attributed Networks

This repository provides a reference implementation of *DHCD* introduced in the paper "Dual-Channel Hybrid Community Detection in Attributed networks" (https://www.sciencedirect.com/science/article/pii/S0020025520310963), which has been accepted by Information Sciences.

### Abstract
This study considers the problem of hybrid community detection in attributed networks based on the information of network topology and attributes with the aim to address the following two shortcomings of existing hybrid community detection methods. First, many of these methods are based on the assumption that network topology and attributes carry consistent information but ignore the intrinsic mismatch correlation between them. Second, network topology is typically treated as the dominant source of information, with attributes employed as the auxiliary source; the dominant effect of attributes is seldom explored or indeed considered. To address these limitations, this paper presents a novel Dual-channel Hybrid Community Detection (DHCD) method that considers the dominant effects of topology and attributes separately. The concept of transition relation between the topology and attribute clusters is introduced to explore the mismatch correlation between the two sources and learn the behavioral and content diversity of nodes. An extended overlapping community detection algorithm is introduced based on the two types of diversity. By utilizing network attributes, DHCD can simultaneously derive the community partitioning membership and corresponding semantic descriptions. The superiority of DHCD over state-of-the-art community detection methods is demonstrated on a set of synthetic and real-world networks.

### Usage

Note that the initliazation strategy may significantly influence the performance of DHCD, i.e., accuracy between the community detection result and ground-truth. Including the strategy in the paper, this repository demonstrates three initialization strategies:
1. Random initialization for **X**, **Y**, **Z**, **U**, and **V**;
2. Hybird initialization, i.e., NNDSVD initialization for **X**, **Y**, and **Z** as well as random initialization for **U** and **V**;
3. NNDSVD initialization for **X**, **Y**, **Z**, **U**, and **V**.

Since the ranomd initialization of NMF has the inherent limitation of getting local optimal solution but not the global optimum, we independently run DHCD T-A and DHCD A-T 10 times in the demonstration of strategy 1 and 2. The community detection result with the minimum objective function value in the 10 independent runs is reported. The random initialization may poentially result in better community detection results (with high quality metrics), but the reported result may not be unique for multiple runs.

In contrast, the community detection result w.r.t. the NNDSVD initialization strategy is unique for a specific settings of the hyper-parameters. The total optimization time of DHCD w.r.t. the NNDSVD initialization is also much shorter than that with ranomd initialization, as we only need to run the DHCD algorithm once for the NNDSVD initialization. However, the corresponding performance may not be as good as the random initialization.

Descriptions of the demostration code are as follow:
1. **demo_dis_rand.m**: disjoint community detection with random initialization, in which the quality metrics w.r.t. different parameter settings will be saved in 'DHCD_T-A(dis_rand).txt' and 'DHCD_A-T(dis_rand).txt';
2. **demo_dis_hyd.m**: disjoint community detection with hybrid initiliazation, in which the quality metrics w.r.t. different parameter settings will be saved in 'DHCD_T-A(dis_hyd).txt' and 'DHCD_A-T(dis_hyd).txt';
3. **demo_dist_SVD.m**: disjoint community detection with NNDSVD initialization, in which the quality metrics w.r.t. different parameter settings will be saved in 'DHCD_T-A(dis_SVD).txt' and 'DHCD_A-T(dis_SVD).txt';
4. **demo_over_rand.m**: overlapping community detecton with random initialization, in which the quality metrics w.r.t. different parameter settings will be saved in 'DHCD_T-A(over_rand).txt' and 'DHCD_A-T(over_rand).txt';
5. **demo_over_hyd.m**: overlapping community detection with hybrid initialization, in which the quality metrics w.r.t. different parameter settings will be saved in 'DHCD_T-A(over_hyd).txt' and 'DHCD_A-T(over_hyd).txt';
6. **demo_over_SVD.m**: overlapping community detection with NNDSVD initialization, in which the quality metrics w.r.t. different parameter settings will be saved in 'DHCD_T-A(over_SVD).txt' and 'DHCD_A-T(over_SVD).txt';
7. **res/**: some example results of Washington, Texas, Cora, and Citeseer are given in the 'res' directory.

### Citing
Please cite the following paper if you use *DHCD* in your research:
```
@article{qin2021dual,
  title={Dual-Channel Hybrid Community Detection in Attributed Networks},
  author={Qin, Meng and Lei, Kai},
  journal={Information Sciences},
  volume={551},
  pages={146--167},
  year={2021},
  publisher={Elsevier}
}
```
If you have any questions, you can contact the author via [mengqin_az@foxmail.com].
