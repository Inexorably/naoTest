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
const int EXPRESSION_MAX_SUBLENGTHS = 4;

// Expression constant min / max values for non exponential values.
const double EXPRESSION_CONST_MIN = -0.2;
const double EXPRESSION_CONST_MAX = 0.2;

// Exponential constant bounds.
const double EXPRESSION_CONST_EXP_MAX = 2;
const double EXPRESSION_CONST_EXP_MIN = -3;

// Number of state variables.
const int NUM_STATE_VARS = 3;

// Number of output vars.
const int NUM_OUTPUT_VARS = 8;

// Mutation probability between 0 to 1.
const double MUTATION_CHANCE = 1/static_cast<double>(NUM_STATE_VARS * NUM_OUTPUT_VARS);