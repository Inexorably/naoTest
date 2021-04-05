#pragma once

// Various utility functions such as printing variables to console or misc math utilities.

#include <webots/TouchSensor.hpp>

#include <vector>
#include <assert.h>
#include <iostream>

#include "globals.h"

// All the webots classes are defined in the "webots" namespace
using namespace webots;

// Simple cartesian point structure.
struct Point {
  double m_x;
  double m_y;
};

// Takes the left and right foot force sensors, and prints the values to console.
void printFootSensors (TouchSensor *fsrL, TouchSensor *fsrR);

// If value is between min:max, return value.  Else, return max if value > max or min if value < min.
double clamp(double value, double min, double max);

// Returns a vector of 2 point elements consisting of the ZMP coordinates for each foot.
std::vector<Point> getZMPCoordinates(TouchSensor *fsrL, TouchSensor *fsrR);