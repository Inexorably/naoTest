#include "motions.h"

// Randomly move the right leg.
void moveRightLeg(const double t, Motor *RHipYawPitch, Motor *RHipRoll, Motor *RHipPitch, Motor *RKneePitch, Motor *RAnklePitch, Motor *RAnkleRoll) {
  // Begin the simulation by lifting the right hip / knee so that we are only supported on
  // the left foot, to simplify learning.  Wait 1 second for this to occur.
  if (t < 1){
    RHipPitch->setPosition(-t);
    RKneePitch->setPosition(t);
    return;
  }


  // Generate a random number to deci
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(0, 1);
  //std::uniform_int_distribution<int> dist(-1, 1);

  std::vector<Motor*> rightLegMotors;
  rightLegMotors.push_back(RHipYawPitch);
  rightLegMotors.push_back(RHipRoll);
  rightLegMotors.push_back(RHipPitch);
  rightLegMotors.push_back(RKneePitch);
  rightLegMotors.push_back(RAnklePitch);
  rightLegMotors.push_back(RAnkleRoll);
   
  // For each motor, choose some random position target that is within x +- x*LAMBDA of the current position.
  for (auto m : rightLegMotors) {
    // Get the current target position.
    double targ = m->getTargetPosition();
    
    // Randomly modify the target position to targ +- [-1:1]*targ*LAMBDA.
    targ += sin(t/3)*dist(mt)*LAMBDA;
    
    // Clamp the target to min:max positions.
    targ = clamp(targ, m->getMinPosition(), m->getMaxPosition());
    
    //Set the new target position.
    m->setPosition(targ);
  }
}