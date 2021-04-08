# naoTest

In creating controllers for bipedal locomotion, the upper body is often neglected and modeled as a rigid, armless torso.  In this repo, the population of organisms (controllers) learn to use the robot's upper body and left hip + left ankle roll to stabilize the robot while a motion plant moves all motors in the right leg and left ankle pitch to attempt to destabilize the system.  Note that the organism controllers currently take the left foot zmp and associated time derivatives as inputs, as it is assumed that the left foot is grounded.  An additional reversed copy of the controller can be added to extend the allowed cases to either feet, but restricting to a single foot speeds up learning.

<img src="https://user-images.githubusercontent.com/16945020/114021636-c1f1e700-9825-11eb-84e6-cf89739f7717.png" width="720">

A simple nonlinear structure is used for the organism controllers.  Each organism takes NUM_STATE_VARS inputs and produces NUM_OUTPUT_VARS outputs (all doubles).  Each output from an organism is produced by the following formula, where the constants A:F are the genetic weights which allow for reproduction / mutation.  Note that for a given organism, the weights vary per output variable -- ie the A<sub>1</sub> value to produce the first output variable is a different value from the A<sub>1</sub> value used to produce the second output variable.  Weights can be examined by opening the .pop or .organism files generated from running runEvolutions().  Note that the value at the zero state of all control functions is zero.

![image](https://user-images.githubusercontent.com/16945020/114023883-43e30f80-9828-11eb-98c6-80a24b3b371b.png)


Evolutions can currently be ran by calling runEvolutions() in main.  The destabilizing lower body motions are called from runEvolutions().  These may either be function calls (such as moveRightLeg()) or direct manipulations in runEvolutions() (such as calling LAnklePitch->setPosition(0.2\*sin(1.2\*robot->getTime()))).

![Stability test gif](https://raw.githubusercontent.com/Inexorably/naoTest/media/media/naoTest3compressed.gif)

Generational data can be obtained by calling writePopulationInfo() in main.

<img src="https://raw.githubusercontent.com/Inexorably/naoTest/media/media/naoplot.jpg" width="480">
