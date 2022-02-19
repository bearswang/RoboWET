# ROBOWET: A Simulation Platform For Robotic Wireless Energy Transfer

Wireless energy transfer (WET) is a ground-breaking technology for cutting the last wire between mobile sensors and power grids in smart cities. Yet, WET only offers effective transmission of energy over a short distance. Robotic WET is an emerging paradigm that mounts the energy transmitter on a mobile robot and navigates the robot through different regions in a large area to charge remote energy harvesters.

The platform consists of

1) mathematical models such as channel fading, nonlinear energy harvesting, robot motion time.

2) planning algorithms such as anchor selection, route planning, resource allocation, collision avoidance.

NOTE: Different solvers (including SDPT3, Sedumi, Mosek) could be tried to address the numerical issues arising in solving convex optimization problems. 

[RoboWET] S. Wang, R. Han, Y. Hong, Q. Hao, M. Wen, L. Musavian, S. Mumtaz, D. W. K. Ng, 
``Robotic wireless energy transfer in dynamic environments: System design and experimental validation,''
IEEE Communications Magazine, Mar. 2022.


# Energy Harvesting Model

The "model" folder is the Matlab code for the nonlinear energy harvesting model in the paper:

[TWC17] S Wang, M Xia, K Huang, YC Wu, 
``Wirelessly powered two-way communication with nonlinear energy harvesting model: Rate regions under fixed and mobile relay,'' 
IEEE Trans. Wireless Commun., vol. 16, no. 12, pp. 8190-8204, Dec. 2017.

# Wirelessly Powered Communication

The "fixed" folder is the Matlab code for the wirelessly powered communication in the paper:

[JSTSP16] S. Wang, M. Xia and Y.-C. Wu, 
``Multi-pair two-way relay network with harvest-then-transmit users: resolving pairwise uplink-downlink coupling,'' 
IEEE Journal of Selected Topics in Signal Processing, vol. 10, no. 8, pp. 1506-1521, Dec. 2016.

# Robotic Wirelessly Powered Two-Way Relay

The "mobile" fplder is the Marlab code for the mobile relay case in the paper:

[TWC17] S Wang, M Xia, K Huang, YC Wu, 
``Wirelessly powered two-way communication with nonlinear energy harvesting model: Rate regions under fixed and mobile relay,'' 
IEEE Trans. Wireless Commun., vol. 16, no. 12, pp. 8190-8204, Dec. 2017.

# Robotic Wirelessly Powered Backscatter Communication

The "global_planning" folder is the Matlab code for the robotic IoT in the paper:

[TWC20] S Wang, M Xia, YC Wu,
``Backscatter data collection with unmanned ground vehicle: Mobility management and power allocation,''
IEEE Trans. Wireless Commun., vol. 18, no. 4, pp. 2314-2328, 2020.

# Collision Avoidance Local Planning

The "local_planning" folder is the code for the collision avoidance problem in the paper:

[IROS20] R Han, S Chen, Q Hao, ``A Distributed Range-Only Collision Avoidance Approach for Low-cost Large-scale Multi-Robot Systems,'' in Proc. IEEE IROS, 2020.

