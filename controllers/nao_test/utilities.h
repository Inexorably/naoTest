#pragma once

// Various utility functions such as printing variables to console or misc math utilities.

#include <webots/TouchSensor.hpp>

#include <vector>
#include <assert.h>
#include <iostream>
#include <random>
#include <math.h>

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

// Take the pointer to the rotational field [rot_field->getSFRotation()], in axis-angle form.
// Axis-angle form is format of unit x y z vector + rotation angle.
// Return the euclidean rotational angles (in xz, xy, yz, meaning order is heading, attitude, bank).
// https://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToEuler/index.htm
/* Example usage in some function near top/main() level:
 *  const double* rot = rot_field->getSFRotation();
 *  std::vector<double> angles = getRotationalAngles(rot);
 */
std::vector<double> getRotationalAngles(const double* rot);

// Returns true with probability p.
bool trueWithProbability(const double p);

// Returns false if the robot has fallen over, as defined per total grf on feet < FOOT_FORCE_MIN.
// Else returns true.
bool isStable(TouchSensor *fsrL, TouchSensor *fsrR);
