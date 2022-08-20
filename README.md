# Particle-Swarm-Optimization-Based-First-Order-Plus-Dead-Time-System-Identification
This project applies one of the most used and efficient optimization algorithms namely Particle Swarm Optimization for system identification purpose.
PSO is a stochastic optimization technique inspired by swarm behaviour observed in nature such as bird and fish schooling,these particles(fish/birds) are flown through hyperdimensional search space. Each particle has a position and velocity and is evaluated through an objective function that we want to maximize or minimize. Once all particles are assessed, we select the best one (solution that maximize or minimize our objective function) and we tune the other particles'position and velocities so to track it, this particle is usually called pbest(best particle). As PSO is an iterative process,particles continue to move around the search space and look for more best solutions. After each iteration, a global best solution is selected based upon the performance of each particle. The Dynamic of particles is adjusted by the mean of a set of coeffcients to speed up or slow down some particles'position and velocities according to their vicinity to potential optimal solution or simply reverse their direction to track the good path led by the particle holding best knowledge.

Before designing our control law for our process (Temperature lab) the first thing to start with is to build up as accurate as possible a model that describes the best our process's dynamics. Our process is an arduino-based lab made up of a heater (actuator) and a sensor to measure the temperature. We started by acquiring data response (temperature) of the process subjected to step signal. The step response was close to the one of a first order system with time delay. Therefore, we chose an empirical description of it.As yet, our model consists of process gain,process time constant and process dead time.
Our objective is to find these three parameters to describe the process. First method is graphic-based, from data response, we try to find the process gain at steady state. the same can be done by evaluating process time constant (time at which the response reaches 63% of its maximum) and we measure also process dead time which is the shift existing between the beginning of input and output.
The second method is by using regression , which means trying to minimize by the mean of an optimization algorithm a performance indicator called ISE(Integral Squared Error).
We have used the second method as it is more efficient and an optimal solution can be obtained.
