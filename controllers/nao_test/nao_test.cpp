// File:          nao_test.cpp
// Date:
// Description:
// Author:
// Modifications:

// You may need to add webots include files such as
// <webots/DistanceSensor.hpp>, <webots/Motor.hpp>, etc.
// and/or to add some other includes
#include <webots/Robot.hpp>
#include <webots/Motor.hpp>
#include <webots/camera.hpp>

// For supervisor / COM etc
#include <webots/Node.hpp>

#include <webots/accelerometer.hpp>
#include <webots/camera.hpp>
#include <webots/DistanceSensor.hpp>
#include <webots/gps.hpp>
#include <webots/gyro.hpp>
#include <webots/InertialUnit.hpp>
#include <webots/keyboard.hpp>
#include <webots/led.hpp>
#include <webots/robot.hpp>
#include <webots/TouchSensor.hpp>
#include <webots/utils/motion.hpp>

#include <webots/Supervisor.hpp>

#include <iostream>
#include <math.h>

// All the webots classes are defined in the "webots" namespace
using namespace webots;

// Simulated devices
static WbDeviceTag CameraTop, CameraBottom;  // cameras

// This is the main program of your controller.
// It creates an instance of your Robot instance, launches its
// function(s) and destroys it at the end of the execution.
// Note that only one instance of Robot should be created in
// a controller program.
// The arguments of the main function can be specified by the
// "controllerArgs" field of the Robot node
int main(int argc, char **argv) {
  // create the Robot instance.
  Supervisor *robot = new Supervisor();

  // get the time step of the current world.
  int timeStep = (int)robot->getBasicTimeStep();
  
  // get handle to robot's translation field
  Node *robot_node = robot->getFromDef("NAO");
  if (robot_node == NULL) {
    std::cerr << "No DEF NAO node found in the current world file" << std::endl;
    exit(1);
  }
  Field *trans_field = robot_node->getField("translation");

  // You should insert a getDevice-like function in order to get the
  // instance of a device of the robot. Something like:
  //  Motor *motor = robot->getMotor("motorname");
  //  DistanceSensor *ds = robot->getDistanceSensor("dsname");
  //  ds->enable(timeStep);
  
  Motor *headYaw = robot->getMotor("HeadYaw");
  Motor *RElbowRoll = robot->getMotor("RElbowRoll");
  Motor *RHipPitch = robot->getMotor("RHipPitch");
  Motor *LShoulderPitch = robot->getMotor("LShoulderPitch");
  Motor *RShoulderPitch = robot->getMotor("RShoulderPitch");
  
  Camera *cameraTop = robot->getCamera("CameraTop");
  Camera *cameraBottom = robot->getCamera("CameraBottom");
  cameraTop->enable(4 * timeStep);
  cameraBottom->enable(4 * timeStep);


  // Main loop:
  // - perform simulation steps until Webots is stopping the controller
  while (robot->step(timeStep) != -1) {
    double t = robot->getTime();
    // Read the sensors:
    // Enter here functions to read sensor data, like:
    //  double val = ds->getValue();

    // Process sensor data here.

    // Enter here functions to send actuator commands, like:
    //  motor->setPosition(10.0);
    
    // Print some info:
    static int ticker = 0;
    ticker++;
    if (ticker > 60) {
      const double *pos = trans_field->getSFVec3f();
      const double *com = robot_node->getCenterOfMass();
      std::cout << "t = " << t << ", Position: " << pos[0] << ' ' << pos[1] << ' ' << pos[2];
      std::cout << ", COM: " << com[0] << ' ' << com[1] << ' ' << com[2] << std::endl;
      ticker = 0;
    }
    
    // Set arms down.
    LShoulderPitch->setPosition(2*sin(t/10));
    RShoulderPitch->setPosition(2*sin(t-M_PI)/10);
    
    // Do some linear interpolation with sin centered around sin(t-pi/2).  Then:
    // max at t%2pi == pi, min at t%2pi == 0.
    // interpolation: y = y1 + ((x – x1) / (x2 – x1)) * (y2 – y1)
    double interpolated = 0.7*(RHipPitch->getMinPosition() + ((sin(robot->getTime()/3.0) - -1) / (1 - -1) * (0 - RHipPitch->getMinPosition())));
    //RHipPitch->setPosition(interpolated);
/*
    interpolated = LShoulderPitch->getMinPosition() + ((sin(robot->getTime()/3.0-M_PI)+1 - 0) / (2 - 0) * (LShoulderPitch->getMaxPosition() - 0));
    LShoulderPitch->setPosition(LShoulderPitch->getMaxPosition());
    

    interpolated = RShoulderPitch->getMinPosition() + ((sin(robot->getTime()) - -1) / (1 - -1) * (0 - RShoulderPitch->getMinPosition()));
    RShoulderPitch->setPosition(1.5);
  */  

  };

  // Enter here exit cleanup code.

  delete robot;
  return 0;
}
