# naoTest

Evolutionary learning to stabilize a NAO robot in the Webots simulator.

In creating controllers for bipedal locomotion, the upper body is often neglected and modeled as a rigid, armless torso.  In this repo, the population of organisms (controllers) learn to use the robot's upper body and left hip + left ankle roll to stabilize the robot while a motion plant moves all motors in the right leg and left ankle pitch to attempt to destabilize the system.

![Stability test gif](https://raw.githubusercontent.com/Inexorably/naoTest/media/media/naoTest3compressed.gif)

Generational data can be obtained by calling writePopulationInfo() in main.

<img src="https://user-images.githubusercontent.com/16945020/113967561-0f4e6400-97e6-11eb-9068-59a5d3dd001b.png" width="480">
