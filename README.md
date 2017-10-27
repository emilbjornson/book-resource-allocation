Optimal Resource Allocation in Coordinated Multi-Cell Systems
==================

This is a code package is related to the follow research book:

Emil Björnson, Eduard Jorswieck, “[Optimal Resource Allocation in Coordinated Multi-Cell Systems](http://kth.diva-portal.org/smash/get/diva2:608533/FULLTEXT01),” Foundations and Trends in Communications and Information Theory, vol. 9, no. 2-3, pp. 113-381, 2013

The package contains Matlab implementations of many of the algorithms described in the book. The use of these algorithms is exemplified by Matlab scripts that generate some of the figures shown in the book. *We encourage you to also perform reproducible research!*

## Abstract of Book

The use of multiple antennas at base stations is a key component in the design of cellular communication systems that can meet high-capacity demands in the downlink. Under ideal conditions, the gain of employing multiple antennas is well-recognized: the data throughput increases linearly with the number of transmit antennas if the spatial dimension is utilized to serve
many users in parallel. The practical performance of multi-cell systems is, however, limited by a variety of nonidealities, such as insufficient channel knowledge, high computational complexity, heterogeneous user conditions, limited backhaul capacity, transceiver impairments, and the constrained level of coordination between base stations.

This tutorial presents a general framework for modeling different multi-cell scenarios, including clustered joint transmission, coordinated beamforming, interference channels, cognitive radio, and spectrum sharing between operators. The framework enables joint analysis and insights that are both scenario independent and dependent.

The performance of multi-cell systems depends on the resource allocation; that is, how the time, power, frequency, and spatial resources are divided among users. A comprehensive characterization of resource allocation problem categories is provided, along with the signal processing algorithms that solve them. The inherent difficulties are revealed: (a) the overwhelming spatial degrees-of-freedom created by the multitude of transmit antennas; and (b) the fundamental tradeoff between maximizing aggregate system throughput and maintaining user fairness. The tutorial provides a pragmatic foundation for resource allocation where the system utility metric can be selected to achieve practical feasibility. The structure of optimal resource allocation is also derived, in terms of beamforming parameterizations and optimal operating points.

This tutorial provides a solid ground and understanding for optimization of practical multi-cell systems, including the impact of the nonidealities mentioned above. The Matlab code is available online for some of the examples and algorithms in this tutorial.


## Content of Code Package

The code package contains implementations of 10 algorithms and precoding schemes. The package also contains the code that reproduces 7 figures in the book: 2.8, 2.10, 3.1, 3.2, 3.6, 4.5 and 4.9. See the file "documentation.pdf" and each file for further documentation. 

The convex optimization problems are implemented using the modeling language [CVX](http://cvxr.com/cvx/).


## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
