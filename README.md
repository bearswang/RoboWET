# ROBOWET: A Simulation Platform For Robotic Wireless Energy Transfer

Wireless energy transfer (WET) is a ground-breaking technology for cutting the last wire between mobile sensors and power grids in smart cities. Yet, WET only offers effective transmission of energy over a short distance. Robotic WET is an emerging paradigm that mounts the energy transmitter on a mobile robot and navigates the robot through different regions in a large area to charge remote energy harvesters.

The platform consists of the following elements:

1) Mathematical models such as channel fading, nonlinear energy harvesting, robot motion time;

2) Planning algorithms such as anchor selection, route planning, resource allocation, collision avoidance.

NOTE: Different solvers (including SDPT3, Sedumi, Mosek) could be tried to address the numerical issues arising in solving convex optimization problems. 

The overall architecture and experimental validation of robotic WET are provided in the following paper:

[RoboWET] S. Wang, R. Han, Y. Hong, Q. Hao, M. Wen, L. Musavian, S. Mumtaz, and D. W. K. Ng, 
``Robotic wireless energy transfer in dynamic environments: System design and experimental validation,''
IEEE Communications Magazine, Mar. 2022.


### Energy Harvesting Model

The "model" folder contains the Matlab code for the nonlinear energy harvesting model in the paper:

[TWC17] S Wang, M Xia, K Huang, and Y.-C. Wu, 
``Wirelessly powered two-way communication with nonlinear energy harvesting model: Rate regions under fixed and mobile relay,'' 
IEEE Trans. Wireless Commun., vol. 16, no. 12, pp. 8190-8204, Dec. 2017.

### Wirelessly Powered Communication

The "fixed" folder contains the Matlab code for the wirelessly powered communication in the paper:

[JSTSP16] S. Wang, M. Xia, and Y.-C. Wu, 
``Multi-pair two-way relay network with harvest-then-transmit users: resolving pairwise uplink-downlink coupling,'' 
IEEE Journal of Selected Topics in Signal Processing, vol. 10, no. 8, pp. 1506-1521, Dec. 2016.

### Robotic Wirelessly Powered Two-Way Relay

The "mobile" folder contains the Marlab code for the robotic two way relay in the paper:

[TWC17] S Wang, M Xia, K Huang, and Y.-C. Wu, 
``Wirelessly powered two-way communication with nonlinear energy harvesting model: Rate regions under fixed and mobile relay,'' 
IEEE Trans. Wireless Commun., vol. 16, no. 12, pp. 8190-8204, Dec. 2017.

### Robotic Wirelessly Powered Backscatter Communication

The "global_planning" folder contains the Matlab code for the robotic IoT in the paper:

[TWC20] S Wang, M Xia, and Y.-C. Wu,
``Backscatter data collection with unmanned ground vehicle: Mobility management and power allocation,''
IEEE Trans. Wireless Commun., vol. 18, no. 4, pp. 2314-2328, 2020.

### Collision Avoidance Local Planning

The "local_planning" folder contains the associated Python/C code for the collision avoidance problem in the paper:

[IROS20] R Han, S Chen, and Q Hao, ``A Distributed Range-Only Collision Avoidance Approach for Low-cost Large-scale Multi-Robot Systems,'' in Proc. IEEE IROS, Las Vegas, NV, USA, 2020.

