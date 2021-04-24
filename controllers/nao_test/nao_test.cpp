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

#include <windows.h>

#include <math.h>
#include <vector>
#include <chrono>

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


// Creates a population and evolves the population of organisms through breeding / mutation.
// Constants are set in globals.h.
int runEvolutions(int argc, char **argv) {
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
  
  // If the pops/ folder where population files are stored does not exist, create it.
  CreateDirectory("pops", NULL);

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
  Population p(1000, 6, 10);
  
  // Load p from the default file.  If no such file exists or file is corrupted,
  // the random p created upon construction will not be changed.
  p.load(DEFAULT_POPULATION_FILENAME, false);
  
  // For debugging purposes, output the best current fitness score.
  double bestFitnessScore = -11111111;
  
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
  
  // We may not start at the first generation.
  int initialGeneration = p.m_generation;
  
  // Iterate through each generation and evolve.
  for (; p.m_generation - initialGeneration < NUM_GENERATIONS; p.m_generation++) {
    std::cout << "Entering generation " << p.m_generation << std::endl;
    
    // Track the runtime of each generation to add to m_runtime of the population object.
    auto start = std::chrono::system_clock::now();
    
    // For each organism in the population, run the simulation in order to generate fitness values.
    unsigned int progressTickerA = 0;   // For printing to console / debugging.
    unsigned int progressTickerB = 0;
    for (Organism& o : p.m_organisms) {
      // Print progress to console every 10% of the current population.
      progressTickerA++;
      if (progressTickerA > p.m_numOrganisms/10) {
        progressTickerB++;
        progressTickerA = 0;
        std::cout << "Current generation (" << p.m_generation << ") progress: " << progressTickerB*10 << "%\n";
      }
      
      // Get the start time of each organism's simulation.
      double t0 = robot->getTime();
    
      // Simulate robot.  If robot stays stable for more than SIMULATION_TIME_MAX seconds, break.
      while (robot->step(timeStep) != -1 && robot->getTime() < SIMULATION_TIME_MAX) {
        // Move right leg randomly.
        moveRightLeg(robot->getTime(), RHipYawPitch, RHipRoll, RHipPitch, RKneePitch, RAnklePitch, RAnkleRoll);
        LAnklePitch->setPosition(0.2*sin(1.2*robot->getTime()));
      }
      
      // Get the time the robot was stable.
      double stableTime = robot->getTime() - t0;
      
      // Increment the runs counter and stable time tracker variables.
      o.m_totalStableTime += stableTime;
      ++o.m_numSimulations;
      
      // Get the fitness score of this Organism.  This needs to be done after updating
      // m_totalStableTime and m_numSimulations.
      double fitnessScore = o.getFitness();
      
      // If this is a new best fitness scoring organism, we save it and print to console.
      if (fitnessScore > bestFitnessScore) {
        bestFitnessScore = fitnessScore;
        std::cout << "New best fitness score: " << bestFitnessScore << std::endl;
        std::string bestOrganismFilename = "pops/star.organism";
        std::cout << "Saving best organism to: " << bestOrganismFilename << std::endl;
        o.save(bestOrganismFilename);
      }
      
      // Reset simulation.
      robot->simulationReset();
      robot->step(timeStep);
      robot->simulationReset();
    }
    
    // Each organism in this generation of the population has been simulated now.
    // We sort by the fitness score.
    p.sortOrganisms();
    
    // Track the runtime of each generation to add to m_runtime of the population object.
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    p.m_runtime += diff.count();
    
    // Save this generation for possible plotting purposes.
    std::cout << "Saving to historic generation population file at pops/generation_" << p.m_generation << ".pop\n";
    std::string generationFilename = "pops/generation_" + std::to_string(p.m_generation) + ".pop";
    p.save(generationFilename);
    
    // Save best performing half of population (POPULATION_SIZE/2).
    std::cout << "Pruning weakest half of population.\n";
    p.m_organisms.erase(p.m_organisms.begin() + p.m_numOrganisms/2, p.m_organisms.end());
    
    // Breed the best performing half of the population (POPULATION_SIZE/2/2 == POPULATION_SIZE/4).
    std::cout << "Breeding population.\n";
    p.reproduceOrganisms();

    
    // Make copies of random members of the previous generation and children, 
    // mutate them, and add them to population.
    // Adds another POPULATION_SIZE/4, and restores our population to POPULATION_SIZE.
    std::cout << "Mutating population.\n";
    p.mutateOrganisms();
    
    // Update the default population file which holds the latest population.
    // Default file name is DEFAULT_POPULATION_FILENAME.
    std::cout << "Saving to latest population file at " << DEFAULT_POPULATION_FILENAME << '\n';
    p.save();
  }
  
  // Clean up and return.
  delete robot;
  return 1;
}

// Creates a population and evolves the population of organisms through breeding / mutation.
// Constants are set in globals.h.
// Uses Gait classes (GaitGene, GaitOrganism, GaitPopulation).
int runGaitEvolutions(int argc, char **argv) {
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
  
  // If the pops/ folder where population files are stored does not exist, create it.
  CreateDirectory("pops", NULL);

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
  std::vector<Motor*> controlMotors = {RAnklePitch, RKneePitch, RHipPitch, LHipPitch, LKneePitch, LAnklePitch};

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
  GaitPopulation p(1000);
  
  // Load p from the default file.  If no such file exists or file is corrupted,
  // the random p created upon construction will not be changed.
  p.load(DEFAULT_POPULATION_FILENAME, false);
  
  // For debugging purposes, output the best current fitness score.
  double bestFitnessScore = -11111111;
  
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
  
  // We may not start at the first generation.
  int initialGeneration = p.m_generation;
  
  // Iterate through each generation and evolve.
  for (; p.m_generation - initialGeneration < NUM_GENERATIONS; p.m_generation++) {
    std::cout << "Entering generation " << p.m_generation << std::endl;
    
    // Track the runtime of each generation to add to m_runtime of the population object.
    auto start = std::chrono::system_clock::now();
    
    // For each organism in the population, run the simulation in order to generate fitness values.
    unsigned int progressTickerA = 0;   // For printing to console / debugging.
    unsigned int progressTickerB = 0;
    for (GaitOrganism& o : p.m_organisms) {
      // Print progress to console every 10% of the current population.
      progressTickerA++;
      if (progressTickerA > p.m_numOrganisms/10) {
        progressTickerB++;
        progressTickerA = 0;
        std::cout << "Current generation (" << p.m_generation << ") progress: " << progressTickerB*10 << "%\n";
      }
      
      // Get the start time of each organism's simulation.
      double t0 = robot->getTime();
      
      // The last time we swapped the gait cycle.
      double tg = 0;
      
      // False means that right leg is support, true means left leg is support.
      bool swingPhase = false;
    
      // Simulate robot.  If robot stays stable for more than SIMULATION_TIME_MAX seconds, break.
      while (robot->step(timeStep) != -1 && robot->getTime() < SIMULATION_TIME_MAX) {/*
        // Hold the robot in the air while we test.
        const double pos[3] = {0.0, 0.6, 0.0};
        trans_field->setSFVec3f(pos);
        
        // Reset the velocities.
        const double vel[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        robot_node->setVelocity(vel);
        */
        
        if (robot->getTime() < 0.6) {
          // Lower the arms while testing.
          RShoulderPitch->setPosition(1.5708);
          LShoulderPitch->setPosition(1.5708);
          continue;
        }
        
        /*
        const double rot[4] = {1.0, 0.0, 0.0, -1.5708};
        rot_field->setSFRotation(rot);
        */
        
        // Don't attempt to control every step.  Waiting more steps can reduce noise.
        static int ticker = 0;
        ticker++;
        if (ticker > STEPS_PER_CONTROL) {
          // Get the inputs (zmpx, zmpy, respective motor target position).
          // Record the dominant food zmps.
          auto zmps = getZMPCoordinates(fsrL, fsrR);
          double zmplx = zmps[0].m_x;
          double zmply = zmps[0].m_y;
          if (zmps[1].m_isSupporting) {
            zmplx = zmps[1].m_x;
            zmply = zmps[1].m_y;
            // TODO: This will cause a discontinuity on the step that we switch supporting feet
            // in the derivatives.  Handle this.
          }
          
          // Add to the organism's m_totalZMPDistance member so that we can reward
          // keeping the zmp x y closer to zero state.
          // We weight by how much time has been spent at this zmp coordinate.
          o.m_totalZMPDistance += sqrt(pow(zmplx, 2)+pow(zmply, 2))/(STEPS_PER_CONTROL*timeStep);

          // Check the derivatives of the zmps.
          zmplxdt = (zmplx - zmplx_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmplydt = (zmply - zmply_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmplxd2t = (zmplxdt - zmplxdt_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmplyd2t = (zmplydt - zmplydt_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          
          // Put the control inputs into a vector to pass to the control function.
          double alpha = 1.0;
          // Check if we need to swap the support / swing legs.
          double tChopped = robot->getTime() - tg;
          if (tChopped > 1/alpha) {
            swingPhase = !swingPhase;
          }
          
          std::vector<double> x = {alpha, tChopped};
          
          // Generate the output vector of the relative angles q.
          std::vector<double> q = o.m_gaitGene.calculateValue(x);
          
          std::vector<double> test = q;   // We will overwrite values, but make the same size.  Is cheap.
          
          // We re-order the motors based on swingPhase.
          // If the right leg is the support leg.
          if (swingPhase) { 
            controlMotors = {RAnklePitch, RKneePitch, RHipPitch, LHipPitch, LKneePitch, LAnklePitch};
          }
          else {
            controlMotors = {LAnklePitch, LKneePitch, LHipPitch, RHipPitch, RKneePitch, RAnklePitch};
          }
          
          test[3] = -1*q[3]/2.0;          // q4 is split between the hip joints, so we half.
          test[2] = q[3]/2.0 + q[2];      // add q3 to the supporting hip joint.
          test[1] = 0;                    // We set the knee joint of the support leg to 0.
            
          // Loop through the control motors and do the inputs.
          for (int i = 0; i < 5; ++i) {
            double input = test[i];
            input = clamp(input, controlMotors[i]->getMinPosition(), controlMotors[i]->getMaxPosition());
            controlMotors[i]->setPosition(input);
          }
          
          // Set the current zmp values as the previous values so we can calculate derivatives in the next step.
          zmplx_prev = zmplx;
          zmply_prev = zmply;
          zmplxdt_prev = zmplxdt;
          zmplydt_prev = zmplydt;
          
          ticker = 0;
        }
      
        // If the robot falls, break.
        if (!isStable(fsrL, fsrR)) {
          break;
        }
      }
      
      // Get the time the robot was stable.
      double stableTime = robot->getTime() - t0;
      
      // Get the total distance in x travelled.
      const double* trans = trans_field->getSFVec3f();
      o.m_totalTranslationX += trans[0];
      
      // Increment the runs counter and stable time tracker variables.
      o.m_totalStableTime += stableTime;
      ++o.m_numSimulations;
      
      // Get the fitness score of this Organism.  This needs to be done after updating
      // m_totalStableTime and m_numSimulations.
      double fitnessScore = o.getFitness();
      
      // If this is a new best fitness scoring organism, we save it and print to console.
      if (fitnessScore > bestFitnessScore) {
        bestFitnessScore = fitnessScore;
        std::cout << "New best fitness score: " << bestFitnessScore << std::endl;
        std::string bestOrganismFilename = "pops/star.organism";
        std::cout << "Saving best organism to: " << bestOrganismFilename << std::endl;
        o.save(bestOrganismFilename);
      }
      
      // Reset simulation.
      robot->simulationReset();
      robot->step(timeStep);
      robot->simulationReset();
    }
    
    // Each organism in this generation of the population has been simulated now.
    // We sort by the fitness score.
    p.sortOrganisms();
    
    // Track the runtime of each generation to add to m_runtime of the population object.
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    p.m_runtime += diff.count();
    
    // Save this generation for possible plotting purposes.
    std::cout << "Saving to historic generation population file at pops/generation_" << p.m_generation << ".pop\n";
    std::string generationFilename = "pops/generation_" + std::to_string(p.m_generation) + ".pop";
    p.save(generationFilename);
    
    // Save best performing half of population (POPULATION_SIZE/2).
    std::cout << "Pruning weakest half of population.\n";
    p.m_organisms.erase(p.m_organisms.begin() + p.m_numOrganisms/2, p.m_organisms.end());
    
    // Breed the best performing half of the population (POPULATION_SIZE/2/2 == POPULATION_SIZE/4).
    std::cout << "Breeding population.\n";
    p.reproduceOrganisms();

    
    // Make copies of random members of the previous generation and children, 
    // mutate them, and add them to population.
    // Adds another POPULATION_SIZE/4, and restores our population to POPULATION_SIZE.
    std::cout << "Mutating population.\n";
    p.mutateOrganisms();
    
    // Update the default population file which holds the latest population.
    // Default file name is DEFAULT_POPULATION_FILENAME.
    std::cout << "Saving to latest population file at " << DEFAULT_POPULATION_FILENAME << '\n';
    p.save();
  }
  
  // Clean up and return.
  delete robot;
  return 1;
}

int test(int argc, char **argv) {
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
  
  // If the pops/ folder where population files are stored does not exist, create it.
  CreateDirectory("pops", NULL);

  //////////////////////////////////////////////////////////////////////////

  // Get motors per NAO proto file.
  
  //Upper body.
  Motor *HeadYaw = robot->getMotor("HeadYaw");
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
  
  // Load motion files.
  Motion hand_wave("../../motions/HandWave.motion");
  Motion forwards("../../motions/Forwards50.motion");
  Motion forwardsTest("../../motions/Forwards50_1.2x_v_p.motion");   // Modified walking motion with 2x speed and no hip / ankle roll
  Motion backwards("../../motions/Backwards.motion");
  Motion side_step_left("../../motions/SideStepLeft.motion");
  Motion side_step_right("../../motions/SideStepRight.motion");
  Motion turn_left_60("../../motions/TurnLeft60.motion");
  Motion turn_right_60("../../motions/TurnRight60.motion");
  
  // For iterating through:
  // The motors we use as inputs for the controller (in addition to the zmp coords).
  // We get these with getTargetPosition().
  // In order to force symmetry, we use mirrored input vectors where the order of L and R are swapped.
  std::vector<Motor*> inputMotorsR = {RHipYawPitch,RHipRoll,RHipPitch,RKneePitch,RAnklePitch,RAnkleRoll};
  std::vector<Motor*> inputMotorsL = {LHipYawPitch,LHipRoll,LHipPitch,LKneePitch,LAnklePitch,LAnkleRoll};
  
  // The motors we use as outputs that the controller provides target positions to.
  // We use these with setPosition().
  // In order to enforce symmetry, we split the output motors into the left and right side of the body.
  // This means that for a given organism, we can call calculateValues and provide mirrored inputs
  // and get mirrored outputs.  This should reduce the search space by 75% (n^2, and reducing n by half).
  std::vector<Motor*> outputMotorsR = {RShoulderPitch, RShoulderRoll, RElbowYaw, RElbowRoll};
  std::vector<Motor*> outputMotorsL = {LShoulderPitch, LShoulderRoll, LElbowYaw, LElbowRoll};

  //Misc sensors
  Camera *cameraTop = robot->getCamera("CameraTop");
  Camera *cameraBottom = robot->getCamera("CameraBottom");
  //cameraTop->enable(4 * timeStep);
  //cameraBottom->enable(4 * timeStep);
  TouchSensor *fsrL = robot->getTouchSensor("LFsr");
  TouchSensor *fsrR = robot->getTouchSensor("RFsr");
  fsrL->enable(timeStep);
  fsrR->enable(timeStep);
  
  //////////////////////////////////////////////////////////////////////////

  // Generate an organism population of controllers.
  Population p(100, inputMotorsR.size(), outputMotorsR.size());
  
  // Load p from the default file.  If no such file exists or file is corrupted,
  // the random p created upon construction will not be changed.
  p.load(DEFAULT_POPULATION_FILENAME, false);
  
  // For debugging purposes, output the best current fitness score.
  // Set bestFitnessScore to the best fitness score in the current population.
  // Iterate through the organisms of the population and set the fitness score.
  // If none of the organisms have been simulated (ie all return the default fitness score),
  // bestFitnessScore will be set to FITNESS_FLOOR.
  double bestFitnessScore = FITNESS_FLOOR;
  for (auto o : p.m_organisms) {
    if (o.getFitness() > bestFitnessScore) {
      bestFitnessScore = o.getFitness();
    }
  }
  std::cout << "Initial bestFitnessScore for star organism set to: " << bestFitnessScore << std::endl;
  
  // We may not start at the first generation.
  int initialGeneration = p.m_generation;
  
  // Iterate through each generation and evolve.
  for (; p.m_generation - initialGeneration < NUM_GENERATIONS; p.m_generation++) {
    std::cout << "Entering generation " << p.m_generation << std::endl;
    
    // Track the runtime of each generation to add to m_runtime of the population object.
    auto start = std::chrono::system_clock::now();
    
    // Track the organism properties of each generation, so that we can calculate the
    // population standard deviations in order to increase the mutation rate when stagnation is
    // detected.  Note that we use the fitness 'components' of each property, which have been weighted,
    // as these are all of the same order of magnitude.
    std::vector<double> timeComponents;
    std::vector<double> zmpComponents;
    std::vector<double> translationXComponent;
    std::vector<double> comVelocityComponent;

    // For each organism in the population, run the simulation in order to generate fitness values.
    unsigned int progressTickerA = 0;   // For printing to console / debugging.
    unsigned int progressTickerB = 0;
    for (Organism& o : p.m_organisms) {
      // Print progress to console every 10% of the current population.
      progressTickerA++;
      if (progressTickerA > p.m_numOrganisms/10) {
        progressTickerB++;
        progressTickerA = 0;
        std::cout << "Current generation (" << p.m_generation << ") progress: " << progressTickerB*10 << "%\n";
      }
      
      // We skip the simulation step if we have already calculated the organism's values.
      if (o.m_numSimulations < 1) {
        // Get the start time of each organism's simulation.
        double t0 = robot->getTime();
        
        // The last time we swapped the gait cycle.
        double tg = 0;
        
        // False means that right leg is support, true means left leg is support.
        bool swingPhase = false;
        
        // We want to minimize the kinetic energy by minimizing COM movement, so we track the
        // positions of the com in order to track the velocity.  We initialize it to the inital value.
        const double* com_0 = robot_node->getCenterOfMass();
        std::vector<double> com_prev;
        
        // Store the previous values in a vector, as the array pointed to is deallocated each time step.
        for (int i = 0; i < 3; i++) {
          com_prev.push_back(com_0[i]);
        }
        
        // Move into the initial position.
        double delay = 1.2;      
        while (robot->step(timeStep) != -1 && robot->getTime() < delay) {
          LHipYawPitch->setPosition(0);
          LHipRoll->setPosition(0);
          LHipPitch->setPosition(-0.7332);
          LKneePitch->setPosition(1.4664);
          LAnklePitch->setPosition(-0.7332);
          LAnkleRoll->setPosition(0);
          RHipYawPitch->setPosition(0);
          RHipRoll->setPosition(0);
          RHipPitch->setPosition(-0.7332);
          RKneePitch->setPosition(1.4664);
          RAnklePitch->setPosition(-0.7332);
          RAnkleRoll->setPosition(0);                
        }
        
        // Play the motion.
        forwardsTest.play();
        
        // Set to true if robot fell and broke loop early.
        bool fell = false;
        
        // Strangely, neither simulationReset does nor simulationResetPhysics work in setting
        // resetting the velocities and accelerations of the robot.  Thus, we bring give time to slow NAO
        // to prevent one organism from messing up the next.
        double postDelay = 1;
      
        // Simulate robot.  If robot is still stable at the end of the motion file, break.
        while (robot->step(timeStep) != -1 && robot->getTime() < forwardsTest.getDuration()/1000.0+delay+postDelay) {
          // Strangely, neither simulationReset does nor simulationResetPhysics work in setting
          // resetting the velocities and accelerations of the robot.  Thus, we bring give time to slow NAO
          // to prevent one organism from messing up the next.
          if (forwardsTest.isOver()) {
            //std::cout << "Motion playback finished at " << robot->getTime() << " s.\n";
            continue;
          }
        
          // Don't attempt to control every step.  Waiting more steps can reduce noise.
          static int ticker = 0;
          ticker++;
          if (ticker > STEPS_PER_CONTROL) {
            // Get the center of mass.
            const double* com = robot_node->getCenterOfMass();
            
            // Find the magnitude of the velocity vector of the COM.
            double comV = sqrt(0*pow(com[0]-com_prev[0],2)+pow(com[1]-com_prev[1],2)+pow(com[2]-com_prev[2],2));
            
            // Add the velocity magnitude to the organism's member variable.
            o.m_totalCOMVelocity += comV / (timeStep * STEPS_PER_CONTROL);
            
            // Set the current COM to be the previous COM.
            for (int i = 0; i < 3; i++) {
              com_prev.push_back(com[i]);
            }
          
            // Get the inputs (zmpx, zmpy, respective motor target position).
            auto zmps = getZMPCoordinates(fsrL, fsrR);
            double zmplx = zmps[0].m_x;
            double zmply = zmps[0].m_y;
            double zmprx = zmps[1].m_x;
            double zmpry = zmps[1].m_y;
            
            // Add to the organism's m_totalZMPDistance member so that we can reward
            // keeping the zmp x y closer to zero state.
            // We weight by how much time has been spent at this zmp coordinate.
            // Add the dominant foot zmps.
            if (zmps[0].m_isSupporting) {
              o.m_totalZMPDistance += sqrt(pow(zmplx, 2)+pow(zmply, 2))/(STEPS_PER_CONTROL*timeStep);          
            }
            else {
              o.m_totalZMPDistance += sqrt(pow(zmprx, 2)+pow(zmpry, 2))/(STEPS_PER_CONTROL*timeStep);          
            }
            
            // Put the control inputs into a vector to pass to the control function.
            std::vector<double> xR;
            std::vector<double> xL;
            
            // Loop through the control / input motors and put their target positions into the x vector.
            for (auto m : inputMotorsR) {
              xR.push_back(m->getTargetPosition());
            }
            for (auto m : inputMotorsL) {
              xL.push_back(m->getTargetPosition());
            }
              
            // Generate each output variable based on the input (ie, loop through the system of equations).
            // Note that we need to multiply some inputs by negative to correctly mirror the motions.
            // The ShoulderRoll, ElbowYaw, and ElbowRoll should be multiplied by -1.  The ShoulderPitch should not.
            std::vector<double> mirror = {1, -1, -1, -1};
            for (int j = 0; j < p.m_numOutputVars; j++) {
              // Find the respective motor input position and clamp it to the min:max bounds of that motor.
              double inputR = o.m_genetics[j].calculateValue(xL);
              double inputL = mirror[j]*o.m_genetics[j].calculateValue(xR);
              
              // Handle NaN results.  If is NaN, set to zero.  Else, don't modify.            
              inputR = std::isnan(inputR) ? 0 : inputR;
              inputL = std::isnan(inputL) ? 0 : inputL;
              
              inputR = clamp(inputR, outputMotorsR[j]->getMinPosition(), outputMotorsR[j]->getMaxPosition());   
              inputL = clamp(inputL, outputMotorsL[j]->getMinPosition(), outputMotorsL[j]->getMaxPosition());   
              outputMotorsR[j]->setPosition(inputR);
              outputMotorsL[j]->setPosition(inputL);
            }
            
            ticker = 0;
          }
        
          // If the robot falls, break.
          if (!isStable(fsrL, fsrR)) {
            fell = true;
            break;
          }
        }
        
        // If we did not travel at least 0.5m in x and stay stable the whole way, we don't reward the organism
        // by increasing its characteristics.
        // Get the total distance in x travelled.
        const double* trans = trans_field->getSFVec3f();
        if (!fell && trans[0] > 0.5) {
          // Get the time the robot was stable.
          //double stableTime = robot->getTime() - t0 - delay - postDelay;
          o.m_totalStableTime += forwardsTest.getDuration()/1000.0; //stableTime;
          o.m_totalTranslationX += 0.5;    // We limit benefit to 0.5 m as we want to mainly reward reducing comv.
        }
        
        // Increment runs counter.
        ++o.m_numSimulations;
      
        // Get the fitness score of this Organism.  This needs to be done after updating
        // m_totalStableTime and m_numSimulations.
        double fitnessScore = o.getFitness();
        
        // If this is a new best fitness scoring organism, we save it and print to console.
        if (fitnessScore > bestFitnessScore) {
          bestFitnessScore = fitnessScore;
          std::cout << "-------------------------------\n";
          std::cout << "New best fitness score: " << bestFitnessScore << std::endl;
          std::string bestOrganismFilename = "pops/star.organism";
          std::cout << "Saving best organism to: " << bestOrganismFilename << "\n";
          o.printFitnessComponents();
          std::cout << "-------------------------------\n";
          o.save(bestOrganismFilename);
        }
      }
      
      // Even if we skip the simulation run due to have already calcuating it, we still
      // will check the stdev of the generation for adjusting the mutation chance.
      // Push back the individual organism weighted components for each property onto their
      // respective vectors so that we can find the population std dev for each generation
      // to adjust the mutation rate going forward.
      std::vector<double> components = o.getFitnessComponents();
      timeComponents.push_back(components[0]);
      zmpComponents.push_back(components[1]);
      translationXComponent.push_back(components[2]);
      comVelocityComponent.push_back(components[3]);
      
      // Reset simulation.  There seems to be a bug where inertia is not reset, so we have to reset
      // twice to counter this (?).
      robot->simulationReset();
      robot->step(timeStep);
      robot->simulationReset();
      //robot->simulationResetPhysics();
    }
    
    // Each organism in this generation of the population has been simulated now.
    // We sort by the fitness score.
    p.sortOrganisms();
    
    // Track the runtime of each generation to add to m_runtime of the population object.
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    p.m_runtime += diff.count();
    
    // Save this generation for possible plotting purposes.
    std::cout << "Saving to historic generation population file at pops/generation_" << p.m_generation << ".pop\n";
    std::string generationFilename = "pops/generation_" + std::to_string(p.m_generation) + ".pop";
    p.save(generationFilename);
    
    // Adjust the mutation rate based on the standard deviation of the fitness components to prevent stagnation.
    // The mutation chance for each generation_i is of the form:
    // chance_i = (num_inputs*num_outputs)^-1 * (stddev_0 / s_i)^c
    // Where stddev_0 is a base approximation of the standard deviation for a completely new population
    // and c is a constant.  
    // Get the base mutation chance which is (num_inputs*num_outputs)^-1.
    double baseMutation = 1.0/static_cast<double>(p.m_numInputVars * p.m_numOutputVars);
    double totalStdDev = sqrt(pow(getStdDev(timeComponents),2)+pow(getStdDev(zmpComponents),2)+pow(getStdDev(translationXComponent),2)+pow(getStdDev(comVelocityComponent),2));
    std::cout << "Component standard deviation: " << totalStdDev << std::endl;
    p.m_chanceMutation = baseMutation*pow(MUTATION_CHANCE_STD_DEV_BASE/totalStdDev,MUTATION_CHANCE_C);
    std::cout << "Adjusting mutation chance to: " << p.m_chanceMutation << std::endl;
    
    // Save best performing half of population (POPULATION_SIZE/2).
    std::cout << "Pruning weakest half of population.\n";
    p.m_organisms.erase(p.m_organisms.begin() + p.m_numOrganisms/2, p.m_organisms.end());
    
    // Breed the best performing half of the population (POPULATION_SIZE/2/2 == POPULATION_SIZE/4).
    std::cout << "Breeding population.\n";
    p.reproduceOrganisms();

    
    // Make copies of random members of the previous generation and children, 
    // mutate them, and add them to population.
    // Adds another POPULATION_SIZE/4, and restores our population to POPULATION_SIZE.
    std::cout << "Mutating population.\n";
    p.mutateOrganisms();
    
    // Update the default population file which holds the latest population.
    // Default file name is DEFAULT_POPULATION_FILENAME.
    std::cout << "Saving to latest population file at " << DEFAULT_POPULATION_FILENAME << '\n';
    p.save();
  }
  
  // Clean up and return.
  delete robot;
  return 1;
}

// Writes historic data in csv format to filename for plotting in an external program.
// Assumes files are in pops/generation_X.pop, where X varies from 0 to n.
// Column 1: Generation number
// Column 2: Mean stable time
// Column 3: Max stable time
// Column 4: Mean zmp distance
// Column 5: Min zmp distance
// Column 6: Mean trans x distance
// Column 7: Max trans x distance
// Column 8: Mean COM velocity
// Column 9: Min COM velocity
// Column 10: Generation stddev
// Column 11: Mutation chance
// Column 12: Fitness score mean
// Column 13: Fitness score max
// TODO: Make not require known population size n, input vars i, and output vars o.
void writePopulationInfo(const std::string& outfilename, const int& n, const int& i, const int&o) {
  std::cout << "Beginning historic data writing to: " << outfilename << std::endl;

  // Create a dummy population object.
  Population p(n, i, o);
  
  // Create the ofstream.
  std::ofstream outfile;
  outfile.open(outfilename);

  // Loop through the generations.
  for (int i = 0; true; i++) {
    // Create the generation filename we are checking.
    std::string generationFilename = "pops/generation_" + std::to_string(i) + ".pop";    
  
    // Create the ifstream so we can check file existance.
    std::ifstream infile(generationFilename);
    
    // If the file does not exist, break.
    std::string line;
    std::getline(infile, line);
    if (line.empty()) {
      break;
    }
    
    // Load the generation into p.
    p.load(generationFilename);
    
    // For each generation, we store the fitness components so we can calculate the generational stddev.
    std::vector<double> timeComponents;
    std::vector<double> zmpComponents;
    std::vector<double> translationXComponent;
    std::vector<double> comVelocityComponent;
    
    // Loop through p and calculate the values of interest.  Initialize the min/max variable trackers to the first
    // valid member of the generation's population's value.
    double stableTimeMean = 0;
    double stableTimeMax;
    double zmpDistanceMean = 0;
    double zmpDistanceMin;
    double transXMean = 0;
    double transXMax;
    double comVMean = 0;
    double comVMin;
    double fitnessMean = 0;
    double fitnessMax;
    
    // Find the first valid organism and break.  TODO: Can combine loops and increase efficiency slightly.
    for (Organism o : p.m_organisms) {   
      // Only record organisms with good runs, ie stable time > 5s and translation x > 0.5 m.
      if (o.m_totalStableTime/static_cast<double>(o.m_numSimulations) < 5.0 || o.m_totalTranslationX/static_cast<double>(o.m_numSimulations) <  0.5) {
        continue;
      }
      stableTimeMax = o.m_totalStableTime/static_cast<double>(o.m_numSimulations);
      zmpDistanceMin = o.m_totalZMPDistance/static_cast<double>(o.m_numSimulations);
      transXMax = o.m_totalTranslationX/static_cast<double>(o.m_numSimulations);
      comVMin = o.m_totalCOMVelocity/static_cast<double>(o.m_numSimulations);
      fitnessMax = o.getFitness();
      break;
    }
    
    // Record the number of non-skipped organisms in each generation so that we can divide the totals/n to get mean.
    int numOrganismsRecorded = 0;
    
    for (Organism o : p.m_organisms) { 
      // Get the fitness components needed for the generational stddev calculations.  The stdev data should
      // include all organisms, not just successful ones, so we push the fitness components onto their respective
      // vectors prior to checking if the organism is successful (and continuing to next iteration if it is not).
      std::vector<double> components = o.getFitnessComponents();
      timeComponents.push_back(components[0]);
      zmpComponents.push_back(components[1]);
      translationXComponent.push_back(components[2]);
      comVelocityComponent.push_back(components[3]);
      
      // Only record organisms with good runs, ie stable time > 5s and translation x > 0.5 m.
      if (o.m_totalStableTime/static_cast<double>(o.m_numSimulations) < 5.0 || o.m_totalTranslationX/static_cast<double>(o.m_numSimulations) <  0.5) {
        continue;
      }
      
      // Increment the recorded organisms ticker.
      numOrganismsRecorded++;
     
      // Add the mean organism values to the respective objects.  We will divide by the number of organisms
      // before writing to the csv.
      double organismMeanStableTime = o.m_totalStableTime/static_cast<double>(o.m_numSimulations);
      stableTimeMean += organismMeanStableTime;
      double organismZMPDistanceMean = o.m_totalZMPDistance/static_cast<double>(o.m_numSimulations);
      zmpDistanceMean += organismZMPDistanceMean;
      double organismTransXMean = o.m_totalTranslationX/static_cast<double>(o.m_numSimulations);
      transXMean += organismTransXMean;
      double organismCOMVMean = o.m_totalCOMVelocity/static_cast<double>(o.m_numSimulations);
      comVMean += organismCOMVMean;
      double organismFitness = o.getFitness();
      fitnessMean += organismFitness;
      
      // Check if the organism values presented thus far are the greatest magnitude so far.
      // If so, update the max value for each characteristic for the current generation.
      if (organismMeanStableTime > stableTimeMax) {
        stableTimeMax = organismMeanStableTime;
      }
      if (organismZMPDistanceMean < zmpDistanceMin) {
        zmpDistanceMin = organismZMPDistanceMean;
      }
      if (organismTransXMean > transXMax) {
        transXMax = organismTransXMean;
      }
      if (organismCOMVMean < comVMin) {
        comVMin = organismCOMVMean;
      }
      if (organismFitness > fitnessMax) {
        fitnessMax = organismFitness;
      }
    }
    // The mean variables are currently the sum of the mean values of the organisms in this generation.
    // Divide by the number of organisms to get mean for the generation.
    stableTimeMean = stableTimeMean/static_cast<double>(numOrganismsRecorded);
    zmpDistanceMean = zmpDistanceMean/static_cast<double>(numOrganismsRecorded);
    transXMean = transXMean/static_cast<double>(numOrganismsRecorded);
    comVMean = comVMean/static_cast<double>(numOrganismsRecorded);
    fitnessMean = fitnessMean/static_cast<double>(numOrganismsRecorded);
    
    // Get the generational stddev.
    double totalStdDev = sqrt(pow(getStdDev(timeComponents),2)+pow(getStdDev(zmpComponents),2)+pow(getStdDev(translationXComponent),2)+pow(getStdDev(comVelocityComponent),2));
    
    // Write our values for the generation to our output file.
    outfile << i;
    outfile << "," << stableTimeMean << "," << stableTimeMax;
    outfile << "," << zmpDistanceMean << "," << zmpDistanceMin;
    outfile << "," << transXMean << "," << transXMax;
    outfile << "," << comVMean << "," << comVMin;
    outfile << "," << totalStdDev << "," << p.m_chanceMutation;
    outfile << "," << fitnessMean << "," << fitnessMax;
    outfile << '\n';
  }

  // We have read all the successive generation files.
  std::cout << "Wrote " << p.m_generation - 1 << " generations of data to " << outfilename << '\n';
  return;
}

int main(int argc, char **argv) {    
  // Testing controller evolution.
  //test(argc, argv);

  // Evolve the controllers (single leg stability).
  //runEvolutions(argc, argv);
  
  // Print generation data to output csv file.
  writePopulationInfo("pops/historicalData.csv", 100, 6, 4);
  
  return 3;
}