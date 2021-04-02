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
#include <assert.h>
#include <vector>

// All the webots classes are defined in the "webots" namespace
using namespace webots;

// globals - todo: make globals.h //

const double FOOT_WIDTH = 0.08;  // per http://simspark.sourceforge.net/wiki/index.php/Models
const double FOOT_LENGTH = 0.16;



////////////////////////////////

// STRUCTS - TODO: ADD TO SEPERATE SOURCE FILES ///

struct point {
  double x;
  double y;
};

///////////////////////////////////////////////////

// Utility functions /////////////////////////

double clamp(double value, double min, double max) {
  if (min > max) {
    assert(0);
    return value;
  }
  return value < min ? min : value > max ? max : value;
}

// Takes the left and right foot force sensors, and prints the values to console.
void printFootSensors (TouchSensor *fsrL, TouchSensor *fsrR) {
  const double *fsv[2] = {fsrL->getValues(), fsrR->getValues()};  // force sensor values

  double l[4], r[4];
  double newtonLeft = 0, newtonRight = 0;

  // The coefficients were calibrated against the real
  // robot so as to obtain realistic sensor values.
  l[0] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Front Left
  l[1] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Front Right
  l[2] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Rear Right
  l[3] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Rear Left

  r[0] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Front Left
  r[1] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Front Right
  r[2] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Rear Right
  r[3] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Rear Left

  int i;
  for (i = 0; i < 4; ++i) {
    l[i] = clamp(l[i], 0, 25);
    r[i] = clamp(r[i], 0, 25);
    newtonLeft += l[i];
    newtonRight += r[i];
  }

  printf("----------foot sensors----------\n");
  printf("   left       right\n");
  printf("+--------+ +--------+\n");
  printf("|%3.1f  %3.1f| |%3.1f  %3.1f|  front\n", l[0], l[1], r[0], r[1]);
  printf("|        | |        |\n");
  printf("|%3.1f  %3.1f| |%3.1f  %3.1f|  back\n", l[3], l[2], r[3], r[2]);
  printf("+--------+ +--------+\n");
  printf("total: %g Newtons, %g kilograms\n", newtonLeft + newtonRight, (newtonLeft + newtonRight) / 9.81);
}

// Returns a vector of 2 point elements consisting of the ZMP coordinates for each foot.
std::vector<point> getZMPCoordinates (TouchSensor *fsrL, TouchSensor *fsrR) {
  const double *fsv[2] = {fsrL->getValues(), fsrR->getValues()};  // force sensor values
  double l[4], r[4];
  double newtonLeft = 0, newtonRight = 0;
  
  l[0] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Front Left
  l[1] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Front Right
  l[2] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Rear Right
  l[3] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Rear Left

  r[0] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Front Left
  r[1] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Front Right
  r[2] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Rear Right
  r[3] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Rear Left

  int i;
  for (i = 0; i < 4; ++i) {
    l[i] = clamp(l[i], 0, 25);
    r[i] = clamp(r[i], 0, 25);
    newtonLeft += l[i];
    newtonRight += r[i];
  }
  
  // Create return object.
  std::vector<point> zmps;
  point temp;
  
  // Left foot.
  temp.x = FOOT_WIDTH * ((l[1]+l[3])-(l[0]+l[2]))/(2*(l[0]+l[1]+l[2]+l[3]));
  temp.y = FOOT_LENGTH * ((l[0]+l[1])-(l[2]+l[3]))/(2*(l[0]+l[1]+l[2]+l[3]));
  zmps.push_back(temp);
  
  // Right foot.
  temp.x = FOOT_WIDTH * ((r[1]+r[3])-(r[0]+r[2]))/(2*(r[0]+r[1]+r[2]+r[3]));
  temp.y = FOOT_LENGTH * ((r[0]+r[1])-(r[2]+r[3]))/(2*(r[0]+r[1]+r[2]+r[3]));
  zmps.push_back(temp);
  
  return zmps; 
}
/////////////////////////////////////////////

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

  //////////////////////////////////////////////////////////////////////////

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

  TouchSensor *fsrL = robot->getTouchSensor("LFsr");
  TouchSensor *fsrR = robot->getTouchSensor("RFsr");
  fsrL->enable(timeStep);
  fsrR->enable(timeStep);
  
  //////////////////////////////////////////////////////////////////////////

  // Main loop:
  // - perform simulation steps until Webots is stopping the controller
  while (robot->step(timeStep) != -1) {
    double t = robot->getTime();
    // Read the sensors:
    // Enter here functions to read sensor data, like:
    //  double val = ds->getValue();
    
    // Find the 

    // Process sensor data here.

    // Enter here functions to send actuator commands, like:
    //  motor->setPosition(10.0);
    
    // Print some info:
    static int ticker = 0;
    ticker++;
    if (ticker > 60) {
      /*
      const double *pos = trans_field->getSFVec3f();
      const double *com = robot_node->getCenterOfMass();
      std::cout << "t = " << t << ", Position: " << pos[0] << ' ' << pos[1] << ' ' << pos[2];
      std::cout << ", COM: " << com[0] << ' ' << com[1] << ' ' << com[2] << std::endl;
      */
      
      //printFootSensors(fsrL, fsrR);
      
      std::vector<point> zmps = getZMPCoordinates(fsrL, fsrR);
      
      std::cout << zmps[0].x << ", " << zmps[0].y << std::endl;
      
      ticker = 0;
    }
    
    //Move arms.
    LShoulderPitch->setPosition(2*sin(t/10));
    RShoulderPitch->setPosition(2*sin(t/10));
    

    
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
