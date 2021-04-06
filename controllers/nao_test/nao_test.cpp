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


#include <math.h>

#include <vector>

#include "utilities.h"
#include "globals.h"
#include "motions.h"
#include "genetics.h"

// All the webots classes are defined in the "webots" namespace
using namespace webots;

// Avoid memory leak.
static int cleanUp(Supervisor *supervisor) {
  delete supervisor;
  return 0;
}

// This is the main program of your controller.
// It creates an instance of your Robot instance, launches its
// function(s) and destroys it at the end of the execution.
// Note that only one instance of Robot should be created in
// a controller program.
// The arguments of the main function can be specified by the
// "controllerArgs" field of the Robot node
int main(int argc, char **argv) {
  //////////////////////////////////////////////////////////////////////////
  
  // Initializing robot.
  
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
  Field *rot_field = robot_node->getField("rotation");

  //////////////////////////////////////////////////////////////////////////

  // Get motors per NAO proto file.
  
  //Upper body.
  Motor *headYaw = robot->getMotor("HeadYaw");
  Motor *HeadPitch = robot->getMotor("HeadPitch");
  Motor *RShoulderPitch = robot->getMotor("RShoulderPitch");
  Motor *RShoulderRoll = robot->getMotor("RShoulderRoll");
  Motor *RElbowYaw = robot->getMotor("RElbowYaw");
  Motor *RElbowRoll = robot->getMotor("RElbowRoll");
  Motor *LShoulderPitch = robot->getMotor("LShoulderPitch");
  Motor *LShoulderRoll = robot->getMotor("LShoulderRoll");
  Motor *LElbowYaw = robot->getMotor("LElbowYaw");
  Motor *LElbowRoll = robot->getMotor("LElbowRoll");
  
  //Lower body.
  Motor *RHipYawPitch = robot->getMotor("RHipYawPitch");
  Motor *RHipRoll = robot->getMotor("RHipRoll");
  Motor *RHipPitch = robot->getMotor("RHipPitch");
  Motor *RKneePitch = robot->getMotor("RKneePitch");
  Motor *RAnklePitch = robot->getMotor("RAnklePitch");
  Motor *RAnkleRoll = robot->getMotor("RAnkleRoll");
  Motor *LHipYawPitch = robot->getMotor("LHipYawPitch");
  Motor *LHipRoll = robot->getMotor("LHipRoll");
  Motor *LHipPitch = robot->getMotor("LHipPitch");
  Motor *LKneePitch = robot->getMotor("LKneePitch");
  Motor *LAnklePitch = robot->getMotor("LAnklePitch");
  Motor *LAnkleRoll = robot->getMotor("LAnkleRoll");
  
  // For iterating through:
  Motor *controlMotors[10] = {RShoulderPitch, RShoulderRoll, RElbowYaw, RElbowRoll, LShoulderPitch, LShoulderRoll, LElbowYaw, LElbowRoll, LHipRoll, LAnkleRoll};

  //Misc sensors
  Camera *cameraTop = robot->getCamera("CameraTop");
  Camera *cameraBottom = robot->getCamera("CameraBottom");
  cameraTop->enable(4 * timeStep);
  cameraBottom->enable(4 * timeStep);
  TouchSensor *fsrL = robot->getTouchSensor("LFsr");
  TouchSensor *fsrR = robot->getTouchSensor("RFsr");
  fsrL->enable(timeStep);
  fsrR->enable(timeStep);
  
  //////////////////////////////////////////////////////////////////////////
  
  // Generate an organism population of controllers.
  Population p(POPULATION_SIZE);
  
  // For debugging purposes, output the best current stable time.
  double bestStableTime = 0;
  
  // Track the zmp derivatives every STEPS_PER_CONTROL steps to use as additional control inputs.
  // Assume starting from rest.
  int zmplx_prev = 0;
  int zmply_prev = 0;
  int zmplxdt = 0;
  int zmplydt = 0;
  int zmplxdt_prev = 0;
  int zmplydt_prev = 0;
  int zmplxd2t = 0;
  int zmplyd2t = 0;
  
  // Iterate through each generation and evolve.
  for (int i = 0; i < NUM_GENERATIONS; i++) {
    std::cout << "Entering generation " << i << std::endl;
    
    // For each organism in the population, run the simulation in order to generate fitness values.
    for (Organism& o : p.m_organisms) {
    
      // Get the start time of each organism's simulation.
      double t0 = robot->getTime();
    
      // Simulate robot.
      while (robot->step(timeStep) != -1) {
        // Move right leg randomly.
        moveRightLeg(robot->getTime(), RHipYawPitch, RHipRoll, RHipPitch, RKneePitch, RAnklePitch, RAnkleRoll);
        
        // Don't attempt to control every step.  Waiting more steps can reduce noise.
        static int ticker = 0;
        ticker++;
        if (ticker > STEPS_PER_CONTROL) {
          // Get the inputs (zmpx, zmpy, respective motor target position).
          // We are moving the right foot, so we are interested in the zmps of the left foot.
          // The left foot zmps is the first element returned by getZMPCoordinates.
          auto zmps = getZMPCoordinates(fsrL, fsrR);
          double zmplx = zmps[0].m_x;
          double zmply = zmps[0].m_y;
          
          // Check the derivatives of the zmps.
          zmplxdt = (zmplx - zmplx_prev) / (timeStep * STEPS_PER_CONTROL);
          zmplydt = (zmply - zmply_prev) / (timeStep * STEPS_PER_CONTROL);
          zmplxd2t = (zmplxdt - zmplxdt_prev) / (timeStep * STEPS_PER_CONTROL);
          zmplyd2t = (zmplydt - zmplydt_prev) / (timeStep * STEPS_PER_CONTROL);
          
          // Put the control inputs into a vector to pass to the control function.
          std::vector<double> x = {zmplx, zmply, zmplxdt, zmplydt, zmplxd2t, zmplyd2t};
            
          // Generate each output variable based on the input (ie, loop through the system of equations).
          for (int j = 0; j < NUM_OUTPUT_VARS; j++) {
            // Find the respective motor input position and clamp it to the min:max bounds of that motor.
            double input = o.m_genetics[j].calculateValue(x);
            input = clamp(input, controlMotors[j]->getMinPosition(), controlMotors[j]->getMaxPosition());   
            controlMotors[j]->setPosition(input);
            
            // Set the current zmp values as the previous values so we can calculate derivatives in the next step.
            zmplx_prev = zmplx;
            zmply_prev = zmply;
            zmplxdt_prev = zmplxdt;
            zmplydt_prev = zmplydt;
          }
          ticker = 0;
        }
      
        // If the robot falls, break.
        if (!isStable(fsrL, fsrR)) {
          break;
        }
      }
      
      // Get the time the robot was stable.
      double stableTime = robot->getTime() - t0;
      
      if (stableTime > bestStableTime) {
        bestStableTime = stableTime;
        std::cout << "New longest stable time: " << bestStableTime << std::endl;
      }
      
      // Increment the runs counter and stable time tracker variables.
      o.m_totalStableTime += stableTime;
      ++o.m_numSimulations;
      
      // Reset simulation.
      robot->simulationReset();
      robot->step(timeStep);
      robot->simulationReset();
    }
    
    // Each organism in this generation of the population has been simulated now.
    // We sort by the fitness score.
    p.sortOrganisms();
    
    // Save best performing half of population (POPULATION_SIZE/2).
    p.m_organisms.erase(p.m_organisms.begin() + POPULATION_SIZE/2, p.m_organisms.end());
    
    // Breed the best performing half of the population (POPULATION_SIZE/2/2 == POPULATION_SIZE/4).
    p.reproduceOrganisms();
    
    // Make copies of random members of the previous generation and children, 
    // mutate them, and add them to population.
    // Adds another POPULATION_SIZE/4, and restores our population to POPULATION_SIZE.
    p.mutateOrganisms();
    
    // Save the population to a file or update the active file.
  }
  
  // Enter here exit cleanup code.

  delete robot;
  return 1;
}
