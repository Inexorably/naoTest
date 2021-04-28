# naoTest

(Currently outdated and does not cover gait evolution)

In creating controllers for bipedal locomotion, the upper body is often neglected and modeled as a rigid, armless torso.  In this repo, the population of organisms (controllers) learn to use the robot's upper body and left hip + left ankle roll to stabilize the robot while a motion plant moves all motors in the right leg and left ankle pitch to attempt to destabilize the system.  Note that the organism controllers currently take the left foot zmp and associated time derivatives as inputs, as it is assumed that the left foot is grounded.  An additional reversed copy of the controller can be added to extend the allowed cases to either feet, but restricting to a single foot speeds up learning.

<img src="https://user-images.githubusercontent.com/16945020/114114382-9a3a6780-9895-11eb-8176-22d95daae787.png" width="720">

A simple nonlinear structure is used for the organism controllers.  Each organism takes NUM_INPUT_VARS inputs and produces NUM_OUTPUT_VARS outputs (all doubles).  Weights can be examined by opening the .pop or .organism files generated from running runEvolutions().  Note that the value at the zero state of all control functions is zero.  For a given output variable y<sub>i</sub>, an organism will calculate that output variable per the following equation:

![gene equation](https://raw.githubusercontent.com/Inexorably/naoTest/4c677763f9eeb2a55ff75f23151b81af5fd559d8/media/NAO%20expression.png)

Where x is the input variable (in the current set up, the zmp x y coordinates of the left foot and the first + second time derivatives comprise the 6 input variables).  Note that n is a random constant scalar.  The Expression class can be examined for more detail.

Evolutions can currently be ran by calling runEvolutions() in main.  The destabilizing lower body motions are called from runEvolutions().  These may either be function calls (such as moveRightLeg()) or direct manipulations in runEvolutions() (such as calling LAnklePitch->setPosition(0.2\*sin(1.2\*robot->getTime()))).

![stability test gif](https://raw.githubusercontent.com/Inexorably/naoTest/media/media/naoTest3compressed.gif)

Generational data can be obtained by calling writePopulationInfo() in main.

<img src="https://raw.githubusercontent.com/Inexorably/naoTest/media/media/naoplot.jpg" width="480">
