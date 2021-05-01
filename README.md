# Abstract

In this project, we explore the basis for an evolutionary hybrid controller to improve bipedal locomotion in a NAO robot.  The controller primarily controls the upper body while the lower body is moved through arbitrary open loop motions.  We split the ZMP space between two controllers.  While in the ‘safe’ ZMP region, the first controller operates with the goal of reducing COM velocity / ZMP distance.  When the ZMP falls out of the safe region and the robot is at risk of falling, the robot switches to the second controller and attempts to stabilize.  Each controller is separately simulated in the Webots environment.

# Background

The locomotion of bipedal robots has consistently been an area of great interest. Many potential operating environments are unfriendly to wheeled robots, from man-made features such as stairs to natural obstacles like roots, rocks, and uneven ground.  Such elements make the strong ground mobility of bipedal robots highly desirable. However, bipedal locomotion is difficult to solve analytically.  Several methods have been widely used in the last decade: primarily the zero-moment-pole (ZMP) [1] and the spring-loaded inverted pendulum (SLIP) [2] methods.  Unfortunately, both methods assume there is no slip at the foot-ground interface and can thus encounter problems when operating in a setting which allows slip.
  
The continual increase of computing power has led to the application of machine learning to many fields, including bipedal locomotion [3].  Machine learning is attractive due to its ease of application to systems where not all variables are explicitly known.  In the context of bipedal locomotion, this means that machine learning based controllers can adapt to differing robot and environmental dynamics.  Machine-learning based controllers can also be used to supplement existing controllers through hybrid methods.  Hybrid methods are similar to 1D piecewise functions.  By dividing the state space into subregions and allocating each subregion to a controller, the system can switch between controllers depending on the current state of the robot [4].  One such implementation may include having a high-performance controller operate in a stable operating region, and in the case of a disturbance that would destabilize the robot the system may switch to a robust controller which brings the robot back into the stable operating region.
     
Both traditional and machine-learning based controllers often simplify the dynamics by removing the arms (and occasionally torso) from the model [5] [6].  While this simplifies the initial controller development, controllers often face difficulty in moving from operating in simulations to controlling physical robots. This problem is particularly relevant to machine-learning based controllers [4] [7], where even slightly inaccurate simulation assumptions such as initial state can cause the controller to fail.

In this paper, we develop the basis for a hybrid controller to improve bipedal locomotion through upper body motion.  We split the ZMP space of the supporting leg between these two controllers: controller “A” and controller “B”. While in the safe ZMP region, the first controller (A) operates with the goal of reducing COM velocity and ZMP distance.  When the ZMP falls out of the safe region and the robot is at risk of falling, the robot switches to the second controller (B) and attempts to stabilize.  For the evolution of both controllers, open loop lower body motions are used.  Controller A receives a forward gait motion, and controller B receives a destabilizing movement of a single leg.  Each controller is seperately simulated and evolved in Webots, an open-source robotics simulator.

# Methodology

## Enviroment / NAO

All evolutions and simulations were done in Webots (C++).  Webots is an open-source 3D robotics simulator widely used in robotics research [8] [9] [10].  Webots basic libraries include various premade 3D models of existing robots.  In this paper, we use the NAO model (Figure 1).

![fig1](https://user-images.githubusercontent.com/16945020/116784728-ab1c6a00-aa4a-11eb-8572-eeff7b0e2855.png)

NAO is a bipedal robot with 25 degrees of freedom.  Each foot has four force sensors located on the vertices, which allow for ZMP calculations.  The controller aiming to improve the standard gait cycle (controller A) controls the shoulder and elbow motors.  The controller aiming to maintain one-legged stability controls the same upper body motors as controller A, and additionally the supporting / grounded leg.  The genetic controllers output target positions for the motors, which are then controlled internally by the respective motor PID controllers (Figure 2).  The optional supervisor node is enabled to record the dynamics and quickly reset the simulation to iterate through each generation.

![fig2](https://user-images.githubusercontent.com/16945020/116784744-be2f3a00-aa4a-11eb-9b97-088d5b09b8a4.png)

## Genetic Algorithm

A genetic algorithm is used to evolve the controllers.  First, a population of n organisms (controllers) is created.  Each controller in the population is simulated (with controller A receiving an open loop forward gait input to the lower body, and controller B receiving an open loop destabilizing motion to the right leg), and the specific simulation results are stored in the organism (controller) object.  The simulation results of each organism are then used to calculate the organism’s fitness score (Eq. 1), which is a measure of how well the organism performed in the simulation.  The population is then sorted by descending fitness score, and the lowest scoring half of the population is removed, leaving n/2 organisms in the population.  The surviving half of the population is then randomly bred to produce n/4 child organisms.  Next, n/4 organisms are randomly selected from the surviving half of the population and then copied and mutated.  These n/4 children and n/4 mutated organisms are then added to the surviving population of n/2 organisms, bringing the total population size back to n organisms.  The population now begins the next generation and repeats the above process.

The fitness score calculation is the most important part of the genetic algorithm as the fitness score is used to select which organisms survive and are subsequently bred / mutated.  The fitness score is calculated as follows:

![eq1](https://user-images.githubusercontent.com/16945020/116784774-e4ed7080-aa4a-11eb-98d2-e8da96b19641.png)

where Wt, Wzmp, Wx, and Wcom are constant positive scalars used to weight the fitness components.  t is the total stable time, ZMPtotal distance is the total distance of the supporting leg from the ZMP across the simulation (see Eq. 2), x is the total distance translated in x, and vcom y,z is the total COM velocity in the y, z axes (see Eq. 5).  Note that since all weight coefficients are positive scalars, each controller is rewarded for its stable time and distance travelled in x (since these products are added) and punished for increased ZMP distance and COM velocity (since these products are subtracted).  Note that when characteristics are irrelevant or simplified evolution is desired, weights can be set to zero to remove their impact on the genetic algorithm.  For example, in the evolution of the controller B organisms which is merely trying to stabilize, the weight of the x translation (Wx) is set to zero.  It is also typically desirable to use the weighting coefficients to keep components in the same order of magnitude.

The ZMP is useful for determining bipedal stability in a LIPM sense.  In addition to the fitness function, the ZMP values are also used to determine controller switching (see Figure 4) and, in conjunction with the total foot force values, if the robot has fallen over so that a given organism’s simulation can be ended.  The total ZMP distance is calculated as follows:

![eq2](https://user-images.githubusercontent.com/16945020/116785076-7dd0bb80-aa4c-11eb-87ca-32340b3a1130.png)

where the x and y ZMP coordinates for a given foot are calculated respectively in Eq. 3 and Eq. 4 [12] [13].  The two-norm of the ZMP based on these coordinates is found at each control step in the simulation, multiplied by the length of the control step, and summed over all steps.  The sum is then divided by the total time of the simulation, t.

![eq3,4](https://user-images.githubusercontent.com/16945020/116785065-74475380-aa4c-11eb-85e3-7ad65d9e49cd.png)


W and L are respectively the width and length of the foot (assumed to be rectangular).  f1, f2, f3 and f4 are the forces in each corner of the foot.

The final fitness component is the total COM velocity.  The total COM velocity is calculated as:

![eq5](https://user-images.githubusercontent.com/16945020/116784808-0e0e0100-aa4b-11eb-9136-69e159ecfbf1.png)

## Class Structure

The class hierarchy is shown below in Figure 3.  The highest-level class is the Population class.  A population contains a vector of organisms, and also stores population-wide information such as the number of input and output variables, the mutation chance, the total runtime, and the generation number.

![fig3](https://user-images.githubusercontent.com/16945020/116784828-254cee80-aa4b-11eb-8965-093d27fa1ed7.png)

The Population class also contains the member functions to breed and mutate members of the population, and the read/write functions to save and load population files.

Each population contains a number of Organisms (controllers).  The organism objects store the simulation data, such as the total stable time, the total ZMP distance, the total x translation, the total COM velocity, and the number of simulations.  The Organism class also contains a function to get the fitness score of the organism, as per Eq. 1.  Each organism has some defining genetics, stored as a member vector of Gene objects.

Each organism has one gene per output variable, such that the genes can be iterated through to calculate each output variable.  The Gene class has the function ‘calculateValue()’ for this purpose (see Eq. 6).  Each output variable is a nonlinear function of all of the input variables.  To facilitate this each gene object has a vector of Expressions, one Expression object per input variable.

The Expression class contains several vectors of scalars.  These vectors are initialized randomly and are of varying length.  Likewise, the scalars are also initialized randomly.  These vectors allow each Organism to be uniquely represented as a series of scalars.

![eq6](https://user-images.githubusercontent.com/16945020/116784843-3695fb00-aa4b-11eb-9716-cb79781af2de.png)

The above equation calculates an output variable yi.  Each value of n is a random positive integer (signifying the random lengths of the vectors in each Expression object).  m is the number of input variables, and i ranges from zero to the number of output variables.  The variables A through H are the scalars contained in the Expression class’s vectors.  For the purposes of this paper, the coefficient scalars were initialized in [-2:2], and the exponent scalars were initialized to a random integer in [0:5].  This scaling is important to speed up evolution, both because large coefficients can slow convergence due to motor positional bounds and because it is desirable to keep different coefficients within an order of magnitude of each other to reduce the chance of single dominating scalars.

## Switching in ZMP space

While time constraints prevented real-time simulations of the coupled controllers switching, the code structure is designed to be extended to switch between controller A (forward gait) and B (single-legged stability) based on the current position in the ZMP space.  Figure 1 shows the proposed allocation of the ZMP space.

![fig4](https://user-images.githubusercontent.com/16945020/116784860-49a8cb00-aa4b-11eb-803f-61bd4b768b73.png)

Note that the interior region is rectangular – this region is defined as:

![eq7,8](https://user-images.githubusercontent.com/16945020/116784866-55948d00-aa4b-11eb-9b0a-ec20a03d2ff4.png)

As defined in Eq. 3 and Eq. 4, W and L are respectively the width and length of NAO’s foot.  The outer region is unbounded.  While in the safe region of the ZMP space, the system uses controller A with the primary goal of reducing ZMP distance and COM velocity in the y and z axes.  When the system moves outside of the safe region of the ZMP space (such as in the case of an external disturbance), the system switches to controller B, which attempts to stabilize the robot and move the ZMP coordinates back into the safe region.   

# Complications and Considerations

## Emergent Behavior

Various complications were encountered.  The most common obstacle when using genetic algorithms is undesired emergent behaviour, stemming from improperly selected fitness functions.  Such undesired behaviour occurred in the evolution of controller A (Figure 5).  In this case, the fitness function used the unmodified x translation to calculate the fitness score.

![fig5](https://user-images.githubusercontent.com/16945020/116784878-74931f00-aa4b-11eb-9868-ac0c3c012ba4.png)

The robot is rewarded for achieving a greater total x translation, so the controller learns to throw its arms forward at the end of the gait cycle, causing NAO to fall forward and end the simulation with a greater total x translation.  This behaviour was prevented by capping the possible rewarded x translation to 0.5 m.

Simulation inaccuracies and constraints are another source of errors.  The simulations calculate forces through bounding boxes, meaning that robot-ground clipping is allowed.  When the position or rotation of the robot is manually reset while clipping is occurring, an accumulation of velocity and acceleration error can occur as the simulation continuously attempts to push the robot out of the ground.  One such symptom of this is the robot begins to glide without external stimulation.  Such simulation constraints also appear when attempting to simulate sagittal support.  Sagittal support is often used in developing controllers for bipedal robots, as it allows the simplification of a 3D system into a 2D plane.  A third controller to generate gait cycles and trajectories was originally planned, but the inability to well support the robot in the sagittal plane (see Figure 6) led to the use of open loop gait motion instead.

![fig6](https://user-images.githubusercontent.com/16945020/116784896-9096c080-aa4b-11eb-9afb-a764d3932030.png)

Another consideration was the selection of the mutation chance.  Traditionally, a popular choice has been to define the mutation chance as the inverse of the sum of the number of input and number of output variables, such that on average one Expression object will change in an organism per mutation call.  However, there is no single optimal mutation constant as the optimal mutation value varies per system and between generations [14] [15], and that dynamic mutation values provide superior space searching [16].  Mutation chances with magnitudes that are too small are slow to evolve and often do not efficiently converge.  Mutation chances with magnitudes that are too high do not sufficiently search the space around the current local maxima and may miss regions of superior performance.  With this in mind, a dynamic mutation chance dependent on the normalized population standard deviation of the population was introduced (Eq. 9).

![eq9](https://user-images.githubusercontent.com/16945020/116784912-9ee4dc80-aa4b-11eb-8107-19cc09d372f6.png)

ninput and noutput are respectively the number of input and output variables.  σ0 is a constant, in this case selected to be 0.785.  σi is the population standard deviation (see Eq. 10) of the fitness scores (a useful scalar which captures the quantifying components) for the current generation.  c is a constant integer, and the optimal value will vary per use-case.  In these simulations, the optimal value of c was found to be 1 (Figure 7).

![fig7](https://user-images.githubusercontent.com/16945020/116784929-b754f700-aa4b-11eb-9068-f52f02f5b943.png)

Note that the previously mentioned population standard deviation is calculated as follows:

![eq10](https://user-images.githubusercontent.com/16945020/116784937-c340b900-aa4b-11eb-9f34-6af7b6b14344.png)

Where the individual values (in this case the fitness scores) are denoted as sj, μ is the mean of the fitness components, and N is the number of organisms in the population.

One final consideration is the size of the search space.  Genetic algorithms search a given space for local maxima, and as such reducing the size of the search space has a direct impact on the speed and quality of evolution.  In this case we exploit symmetry to reduce the search space by 75% from the naïve case.  NAO is physically symmetrical about the XZ axis.  This means that given some left and right inputs states which produce some given output positions for the left and right side motors, swapping these left and right input states should produce swapped left and right output positions.

# Simulation Results and Conclusions

Using a mutation constant value of c = 1, controller A (handling the open loop gait motion) was evolved over 1000 generations.  The results of Eq. 9 are shown below in Figure 8.

![fig8](https://user-images.githubusercontent.com/16945020/116784953-db183d00-aa4b-11eb-9bdb-f5897139a8c6.png)

Note how as the population standard deviation continues to increase (due to the increasing fitness disparity between mutated copies and star organisms), the mutation chance continues to decrease in order to search the state space around the high performing organisms.  The generational fitness is shown below in Figure 9.
     
![fig9](https://user-images.githubusercontent.com/16945020/116784966-e9feef80-aa4b-11eb-936b-df096ef57777.png)

Note that the mean fitness lags the max fitness.  As new local maxima are found, the respective organism’s genes are spread into the population over subsequent generations through breeding and mutated copies.  However, the generated motion (Figure 10) after 1000 generations for controller A is not natural / human-like.

![fig10](https://user-images.githubusercontent.com/16945020/116785138-cc7e5580-aa4c-11eb-8c3a-adb24e1e99e9.png)

While there are semblances to natural gait, the robot favours holding its arms behind its body, as opposed to a human gait in which the arms would swing back-and-forth.  This is likely due to the punishing ZMP component in the fitness score.  The knees are bent forward during the gait cycle due to the LIPM nature of the open loop gait, which moves the COM and ZMP coordinates slightly forward.  It is also possible that the unnatural gait could be improved by reducing the search space, such as by solving the form of the governing equations to provide a better template for the output variable calculation (in comparison to Eq. 6).

![fig11](https://user-images.githubusercontent.com/16945020/116784984-08fd8180-aa4c-11eb-91c2-c230080df684.png)

![image](https://user-images.githubusercontent.com/16945020/116784998-16b30700-aa4c-11eb-9c08-b30041d1b755.png)

Due to time constraints, controller A and controller B were tested separately, and testing the hybrid switching proposed in Figure 4 will have to be saved for future work.  Performance for controller B (Figure 11) did converge more quickly but will require coupling with controller A to provide a larger variety of motion plants to evolve more generally, as well as to encounter aspects such as non-zero initial velocity which were not part of the simplified model.

# References

[1] 	R. Tedrake, S. Kuindersma, R. Deits and K. Miura, “A closed-form solution for real-time zmp gait generation and feedback stabilization,” IEEE-RAS 15th International Conference on Humanoid Robots, pp. 936-940, Nov 2015. 

[2] 	H. R. Vejdani, A. Wu, H. Geyer and J. W. Hurst, “Touch-down angle control for spring-mass walking,” IEEE International Conference on Robotics and Automation (ICRA), pp. 5101-5106, May 2015. 

[3] 	R. Wang, S. Hudson, Y. Li, H. Wu and C. Zhou, “Normalized Neural Network for Energy Efficient Bipedal Walking Using Nonlinear Inverted Pendulum Model,” IEEE International Conference on Robotics and Biomimetics, pp. 1400-1406, 2019. 

[4] 	B. Semwal, M. Raj and G. Nanda, “Hybrid Model for Passive Locomotion Control of a Biped Humanoid: The Artificial Neural Network Approach,” International Journal of Interactive Multimedia and Artificial Intelligence, vol. 5, pp. 40-46, 2018. 

[5] 	W. L. Ma, Y. Or and A. Ames, “Dynamic Walking on Slippery Surfaces: Demonstrating Stable Bipedal Gaits with Planned Ground Slippage,” IEEE International Conference on Robotics and Automation, May 2019. 

[6] 	X. Xiao and F. Asano, “Analytical Solution of Target Walking Speed Generation by Underactuated Compass-like Bipedal Walker,” IEEE International Conference on Robotics and Biomimetics, pp. 1017-1022. 

[7] 	Z. Li, X. Cheng, X. Peng, P. Abbeel, S. Levine, G. Berseth and K. Sreenath, “Reinforcement Learning for Robust Parameterized Locomotion Control of Bipedal Robots,” International Conference on Robotics and Automation, 2021. 

[8] 	V. Ganapathy, C. Soh and W. Lui, “Utilization of Webots and Khepera II as a platform for Neural Q-Learning controllers,” IEEE Symposium on Industrial Electronics Applications, vol. 2, pp. 783-788, 2009. 

[9] 	Y. Pan and X. Ma, “Design of Industrial Robot Sorting System with Visual Guidance Based on Webots,” 3rd International Conference on Computer and Communication Systems, pp. 516-521, 2018. 

[10] 	L. Hohl, R. Tellez, O. Michel and A. Janijspeert, “Aibo and Webots: Simulation, wireless remote control and controller transfer,” Robotics and Autonomous Systems, vol. 54, no. 6, pp. 472-485, 2006. 

[11] 	A. Massah, A. Sharifi, Y. Salehinia and F. Najafi, “An Open Loop Walking on Different Slopes for NAO Humanoid Robot,” IEEE International Symposium on Robotics and Intelligent Sensors, pp. 296-304, 2012. 

[12] 	C. Gil, H. Calvo and H. Sossa, “Learning an Efficient Gait Cycle of a Biped Robot Based on Reinforcement Learning and Artificial Neural Networks,” Applied Sciences - Advanced Mobile Robotics, pp. 502-526, 2019. 

[13] 	J.-L. Lin and W.-C. J. Y.-J. C. Kao-Shing Hwang, “Gait Balance and Acceleration of a Biped Robot,” IEEE Access, vol. 4, pp. 2439-2449, 2016. 

[14] 	M. Serpell and J. Smith, “Self-Adaptation of Mutation Operator and Self-Adaptation of Mutation Operator and in Genetic Algorithms,” Evolutionary Computation, vol. 18, no. 3, pp. 491-514, 2010. 

[15] 	R. Greenwell and J. Angus, “Optimal Mutation Probability for Genetic Algorithm,” Mathematical annd Computer Modeling, vol. 21, no. 8, pp. 1-11, 1995. 

[16] 	A. Hassanat, K. Almohammadi, E. Alkafaween, E. Abunawas, A. Hammouri and A. Hammouri, “Choosing Mutation and Crossover Ratios for Genetic Algorithms—A Review with a New Dynamic Approach,” Information, vol. 10, no. 390, 2019. 

[17] 	T. Li, H. Geyer, C. G. Atkeson and A. Rai, “Using Deep Reinforcement Learning to Learn High-Level Policies on the ATRIAS Biped,” IEEE International Conference on Robotics and Automation, pp. 263-269, 2019. 

[18] 	R. J. Griffin, G. Wiedebach, S. Bertand, A. Leonessa and J. Pratt, “Straight-Leg Walking Through Underconstrained Whole-Body Control,” IEEE International Conference on Robotics and Automation, 2018. 

[19] 	Y. Huang, Q. Wang, B. Chen and G. Xie, “Modeling and gait selection of passivity-based seven-link bipeds with dynamic series of walking phases,” Robotica, vol. 30, no. 1, pp. 39-51, 2012. 

[20] 	M. Susi, “Gait Analysis for Pedestrian Navigation Using MEMS Handheld Devices,” University of Calgary, 2015. 





