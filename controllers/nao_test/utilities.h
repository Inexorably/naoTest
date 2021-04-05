#pragma once

// Various utility functions such as printing variables to console or misc math utilities.

#include <webots/TouchSensor.hpp>

#include <vector>
#include <assert.h>
#include <iostream>
#include <random>

#include "globals.h"

// All the webots classes are defined in the "webots" namespace
using namespace webots;

// Simple cartesian point structure.
struct Point {
  double m_x;
  double m_y;
};

// Takes the left and right foot force sensors, and prints the values to console.
void printFootSensors(TouchSensor *fsrL, TouchSensor *fsrR);

// If value is between min:max, return value.  Else, return max if value > max or min if value < min.
double clamp(const double value, const double min, const double max);

// Returns a vector of 2 point elements consisting of the ZMP coordinates for each foot.
std::vector<Point> getZMPCoordinates(TouchSensor *fsrL, TouchSensor *fsrR);

// Returns true with probability p.
bool trueWithProbability(const double p);

// Returns false if the robot has fallen over, as defined per total grf on feet < FOOT_FORCE_MIN.
// Else returns true.
bool isStable(TouchSensor *fsrL, TouchSensor *fsrR);
