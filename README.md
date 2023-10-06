# Optimal Energy-Efficient Transmit Beamforming for Multi-User MISO Downlink

This repo is cloned from [this CodeOcean capsule](https://codeocean.com/capsule/7272942/tree/v1), which contains the code for the following scientific paper:

Oskari Tervo, Le-Nam Tran, and Markku Juntti, "Optimal Energy-Efficient Transmit Beamforming for Multi-User MISO Downlink," IEEE Trans. Signal Process., vol. 63, no. 20, pp. 5574-5587, Oct. 2015.

## Instructions
The "**main.m**" script plots the convergence of Algorithms 1 and 3 in the paper. Note that, for Algorithm 1, we avail of the **optimizer** function in YALMIP that creates an object with a pre-compiled low-level numerical format which can be used to efficiently solve a series of the feasibility problem. Further information on this interesting functionality can be found [here](https://yalmip.github.io/command/optimizer/).
