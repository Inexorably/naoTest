#pragma once

/////////////////////////////////////// NAO physical properties //////////////////////////////////////////

const double FOOT_WIDTH = 0.08;  // per http://simspark.sourceforge.net/wiki/index.php/Models
const double FOOT_LENGTH = 0.16;

/////////////////////////////////////// Random motion / plant ////////////////////////////////////////////

// For lower body motion, new target positions will be within x +- [0:1]*x*LAMBDA of current position.
// Motors are servos, so target position is in rads.
const double LAMBDA = 0.01;

/////////////////////////////////////// Genetics constants ///////////////////////////////////////////////

// Determines the max number of instances of a type of function on a state variable an expression can have.
// See Expression class in "genetics.h".
const int MAX_EXPRESSION_SUBLENGTHS = 4;

// Expression constant min / max values.
const double EXPRESSION_CONST_MIN = -20;
const double EXPRESSION_CONST_MAX = 20;

// Number of state variables.
const int NUM_STATE_VARS = 3;

