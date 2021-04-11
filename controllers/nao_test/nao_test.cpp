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
  Motor *controlMotors[12] = {RShoulderPitch, RShoulderRoll, RElbowYaw, RElbowRoll, LShoulderPitch, LShoulderRoll, LElbowYaw, LElbowRoll, LHipRoll, LAnkleRoll, RHipRoll, RAnkleRoll};

  //Misc sensors
  Camera *cameraTop = robot->getCamera("CameraTop");
  Camera *cameraBottom = robot->getCamera("CameraBottom");
  //cameraTop->enable(4 * timeStep);
  //cameraBottom->enable(4 * timeStep);
  TouchSensor *fsrL = robot->getTouchSensor("LFsr");
  TouchSensor *fsrR = robot->getTouchSensor("RFsr");
  fsrL->enable(timeStep);
  fsrR->enable(timeStep);
  
  // Load motion files.
  Motion hand_wave("../../motions/HandWave.motion");
  Motion forwards("../../motions/Forwards50.motion");
  Motion forwardsTest("../../motions/Forwards50testNoRoll.motion");   // Modified walking motion with 2x speed and no hip / ankle roll
  Motion backwards("../../motions/Backwards.motion");
  Motion side_step_left("../../motions/SideStepLeft.motion");
  Motion side_step_right("../../motions/SideStepRight.motion");
  Motion turn_left_60("../../motions/TurnLeft60.motion");
  Motion turn_right_60("../../motions/TurnRight60.motion");
  
  //////////////////////////////////////////////////////////////////////////
  
  // Generate an organism population of controllers.
  Population p;
  
  // Load p from the default file.  If no such file exists or file is corrupted,
  // the random p created upon construction will not be changed.
  p.load(DEFAULT_POPULATION_FILENAME, false);
  
  // For debugging purposes, output the best current fitness score.
  double bestFitnessScore = -11111111;
  
  // Track the zmp derivatives every STEPS_PER_CONTROL steps to use as additional control inputs.
  // Assume starting from rest.
  
  // Left foot zmp data.
  int zmplx_prev = 0;
  int zmply_prev = 0;
  int zmplxdt = 0;
  int zmplydt = 0;
  int zmplxdt_prev = 0;
  int zmplydt_prev = 0;
  int zmplxd2t = 0;
  int zmplyd2t = 0;

  // Right foot zmp data.
  int zmprx_prev = 0;
  int zmpry_prev = 0;
  int zmprxdt = 0;
  int zmprydt = 0;
  int zmprxdt_prev = 0;
  int zmprydt_prev = 0;
  int zmprxd2t = 0;
  int zmpryd2t = 0;
  
  // We may not start at the first generation.
  int initialGeneration = p.m_generation;
  
  // Iterate through each generation and evolve.
  for (; p.m_generation - initialGeneration < NUM_GENERATIONS; p.m_generation++) {
    std::cout << "Entering generation " << p.m_generation << std::endl;
    
    // Track the runtime of each generation to add to m_runtime of the population object.
    auto start = std::chrono::system_clock::now();
    
    // For each organism in the population, run the simulation in order to generate fitness values.
    int progressTickerA = 0;   // For printing to console / debugging.
    int progressTickerB = 0;
    for (Organism& o : p.m_organisms) {
      // Print progress to console every 10% of the current population.
      progressTickerA++;
      if (progressTickerA > POPULATION_SIZE/10) {
        progressTickerB++;
        progressTickerA = 0;
        std::cout << "Current generation (" << p.m_generation << ") progress: " << progressTickerB*10 << "%\n";
      }
      
      // Get the start time of each organism's simulation.
      double t0 = robot->getTime();
      
      // Test the walk forward motion.
      LShoulderPitch->setPosition(1.57);
      RShoulderPitch->setPosition(1.57);
      bool played = false;
    
      // Simulate robot.  If robot stays stable for more than SIMULATION_TIME_MAX seconds, break.
      while (robot->step(timeStep) != -1 && robot->getTime() < SIMULATION_TIME_MAX) {
        if (robot->getTime() > 0.5 && !played) {
          played = true;
          forwardsTest.play();
          continue;
        }
        
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
          double zmprx = zmps[1].m_x;
          double zmpry = zmps[1].m_y;
          
          // Add to the organism's m_totalZMPDistance member so that we can reward
          // keeping the zmp x y closer to zero state.
          // We weight by how much time has been spent at this zmp coordinate.
          // We use the foot with more weight on it.
          // If the left foot has more weight on it:
          if (zmps[0].m_dominant) {
            o.m_totalZMPDistance += sqrt(pow(zmplx, 2)+pow(zmply, 2))/(STEPS_PER_CONTROL*timeStep);
          }
          // Else, the right foot is supporting more weight:
          else {
            o.m_totalZMPDistance += sqrt(pow(zmprx, 2)+pow(zmpry, 2))/(STEPS_PER_CONTROL*timeStep);
          }

          // Check the derivatives of the zmps.
          zmplxdt = (zmplx - zmplx_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmplydt = (zmply - zmply_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmplxd2t = (zmplxdt - zmplxdt_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmplyd2t = (zmplydt - zmplydt_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmprxdt = (zmprx - zmprx_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmprydt = (zmpry - zmpry_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmprxd2t = (zmprxdt - zmprxdt_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          zmpryd2t = (zmprydt - zmprydt_prev) / static_cast<double>(timeStep * STEPS_PER_CONTROL);
          
          // Put the control inputs into a vector to pass to the control function.
          std::vector<double> x = {zmplx, zmply, zmplxdt, zmplydt, zmplxd2t, zmplyd2t, zmprx, zmpry, zmprxdt, zmprydt, zmprxd2t, zmpryd2t};
            
          // Generate each output variable based on the input (ie, loop through the system of equations).
          for (int j = 0; j < NUM_OUTPUT_VARS; j++) {
            // Find the respective motor input position and clamp it to the min:max bounds of that motor.
            double input = controlMotors[j]->getTargetPosition()+o.m_genetics[j].calculateValue(x);
            input = clamp(input, controlMotors[j]->getMinPosition(), controlMotors[j]->getMaxPosition());   
            controlMotors[j]->setPosition(input);
            
            // Set the current zmp values as the previous values so we can calculate derivatives in the next step.
            zmplx_prev = zmplx;
            zmply_prev = zmply;
            zmplxdt_prev = zmplxdt;
            zmplydt_prev = zmplydt;
            zmprx_prev = zmprx;
            zmpry_prev = zmpry;
            zmprxdt_prev = zmprxdt;
            zmprydt_prev = zmprydt;
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
      
      // Get the x translation.
      const double *pos = trans_field->getSFVec3f();
      o.m_translationX += pos[0];
      
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
    p.m_organisms.erase(p.m_organisms.begin() + POPULATION_SIZE/2, p.m_organisms.end());
    
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
// Column 2: Min stable time
// Column 3: Mean stable time
// Column 4: Max stable time
void writePopulationInfo(const std::string& outfilename) {
  std::cout << "Beginning historic data writing to: " << outfilename << std::endl;

  // Create a dummy population object.
  Population p;
  
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
    
    // Loop through p and calculate the values of interest.
    double stableTimeMin = 0;
    double stableTimeMean = 0;
    double stableTimeMax = 0;
    for (Organism o : p.m_organisms) {
      double organismMeanStableTime = o.getFitness();
      stableTimeMean += organismMeanStableTime;
      if (organismMeanStableTime < stableTimeMin) {
        stableTimeMin = organismMeanStableTime;
        continue;
      }
      if (organismMeanStableTime > stableTimeMax) {
        stableTimeMax = organismMeanStableTime;
        continue;
      }
    }
    // stableTimeMean is currently the sum of the mean stable times of the organisms in that generation.
    // Divide by the number of organisms to get mean.
    stableTimeMean = stableTimeMean/p.m_organisms.size();
    
    // Write our values to our output file.
    outfile << i << "," << stableTimeMin << "," << stableTimeMean << "," << stableTimeMax << '\n';
  }

  // We have read all the successive generation files.
  std::cout << "Wrote " << p.m_generation - 1 << " generations of data to " << outfilename << '\n';
  return;
}

int main(int argc, char **argv) {
  // Evolve the controllers.  
  runEvolutions(argc, argv);
  
  // Print generation data to output csv file.
  //writePopulationInfo("pops/historicalData.csv");
  
  return 2;
}