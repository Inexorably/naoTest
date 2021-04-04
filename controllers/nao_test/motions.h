#pragma once

// Execute various lower body motions in order to act as a plant for the controller
// to learn to stabilize through upper body motions.

// We only need to focus on learning to stabilize around a specific leg, because
// we can mirror the final controller for stabilizing the other leg as NAO is symmetrical.

#include <webots/Motor.hpp>

#include <random>
#include <iostream>
#include <math.h>

#include "globals.h"
#include "utilities.h"

using namespace webots;

// Randomly move the right leg.
void moveRightLeg(const double t, Motor *RHipYawPitch, Motor *RHipRoll, Motor *RHipPitch, Motor *RKneePitch, Motor *RAnklePitch, Motor *RAnkleRoll);