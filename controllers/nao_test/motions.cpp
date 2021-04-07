#include "motions.h"

// Randomly move the right leg.
void moveRightLeg(const double t, Motor *RHipYawPitch, Motor *RHipRoll, Motor *RHipPitch, Motor *RKneePitch, Motor *RAnklePitch, Motor *RAnkleRoll) {
  // Begin the simulation by lifting the right hip / knee so that we are only supported on
  // the left foot, to simplify learning.  Wait 1 second for this to occur.
  if (t < 1){
    RHipPitch->setPosition(t*RHipPitch->getMinPosition());
    RKneePitch->setPosition(t*RKneePitch->getMaxPosition());
    return;
  }

  std::vector<Motor*> rightLegMotors;
  rightLegMotors.push_back(RHipYawPitch);
  rightLegMotors.push_back(RHipRoll);
  //rightLegMotors.push_back(RHipPitch);
  rightLegMotors.push_back(RKneePitch);
  rightLegMotors.push_back(RAnklePitch);
  rightLegMotors.push_back(RAnkleRoll);
   
  // For each motor, choose some random position target that is within x +- x*LAMBDA of the current position.
  double d = 0;    // Index for frequency modifier.  We do this to increase motion variation.
  for (auto m : rightLegMotors) {
    // Get the current target position.
    double targ = m->getTargetPosition();
    
    // Randomly modify the target position to targ +- [-1:1]*targ*LAMBDA.
    double delta = sin((1.0+0.05*d)*2*t)*1.5*LAMBDA;
    d++;       // Increment the index for the frequency modifier.
    
    // Make it progressively harder as t increases.  At 20 seconds, target changes double.
    // By 40 seconds, target changes 5x in magnitude.
    targ += delta + delta * pow(t/20.0, 2);
    
    // Clamp the target to min:max positions.
    targ = clamp(targ, m->getMinPosition(), m->getMaxPosition());
    
    //Set the new target position.
    m->setPosition(targ);
  }
}