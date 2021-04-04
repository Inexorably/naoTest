#pragma once

// NAO physical properties
const double FOOT_WIDTH = 0.08;  // per http://simspark.sourceforge.net/wiki/index.php/Models
const double FOOT_LENGTH = 0.16;

// For lower body motion, new target positions will be within x +- [0:1]*x*LAMBDA of current position.
// Motors are servos, so target position is in rads.
const double LAMBDA = 0.01;