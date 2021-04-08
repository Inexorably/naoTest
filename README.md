# naoTest

Evolutionary learning to stabilize a NAO robot in the Webots simulator.

In creating controllers for bipedal locomotion, the upper body is often neglected and modeled as a rigid, armless torso.  In this repo, the population of organisms (controllers) learn to use the robot's upper body and left hip + left ankle roll to stabilize the robot while a motion plant moves all motors in the right leg and left ankle pitch to attempt to destabilize the system.

![Stability test gif](https://raw.githubusercontent.com/Inexorably/naoTest/media/media/naoTest3compressed.gif)

Generational data can be obtained by calling writePopulationInfo() in main.

<img src="https://raw.githubusercontent.com/Inexorably/naoTest/media/media/naoplot.jpg" width="480">
