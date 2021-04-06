#pragma once

/////////////////////////////////////// Simulation constants /////////////////////////////////////////////

// If the total grf / force feet supporting drops below FOOT_FORCE_MIN, the robot has fallen over
// and we end the simulation.
const double FOOT_FORCE_MIN = 40;

// Max simulation time before breaking in seconds.
const int SIMULATION_TIME_MAX = 60;

// The number of generations of our population we iterate through.
const int NUM_GENERATIONS = 100;

// The number of steps between each control attempt.  Higher values allow less noise but slower response.
const int STEPS_PER_CONTROL = 4;

/////////////////////////////////////// NAO physical properties //////////////////////////////////////////

const double FOOT_WIDTH = 0.08;  // per http://simspark.sourceforge.net/wiki/index.php/Models
const double FOOT_LENGTH = 0.16;

/////////////////////////////////////// Random motion / plant ////////////////////////////////////////////

// Modifier for delta position target in motion plant.
// Motors are servos, so target position is in rads.
const double LAMBDA = 0.02;

// Have slightly different frequency coefficients for each motor to introduce more motion variance.
// Note: this must be of length NUM_OUTPUT_VARS.
// TODO: Add verification check on OMEGA.
const std::vector<double> OMEGA = {1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50};

/////////////////////////////////////// Genetics constants ///////////////////////////////////////////////

// Determines the max number of instances of a type of function on a state variable an expression can have.
// See Expression class in "genetics.h".
const int EXPRESSION_MAX_SUBLENGTHS = 4;

// Expression constant min / max values for non exponential values.
const double EXPRESSION_CONST_MIN = -2;
const double EXPRESSION_CONST_MAX = 2;

// Exponential constant bounds.  Using doubles for exponents can cause domain errors (ie -2^2.34) and NaN.
// Allowing negative exponents creates inf values easily (ie 0.000000001^-5).
const int EXPRESSION_CONST_EXP_MAX = 5;
const int EXPRESSION_CONST_EXP_MIN = 0;

// Number of state variables.
const int NUM_STATE_VARS = 6;

// Number of output vars.
const int NUM_OUTPUT_VARS = 10;

// Mutation probability between 0 to 1.
const double MUTATION_CHANCE = 1/static_cast<double>(NUM_STATE_VARS * NUM_OUTPUT_VARS);

// Number of organisms in population.
const int POPULATION_SIZE = 1000;

/////////////////////////////////////// File constants ///////////////////////////////////////////////

// Default population filename.
const std::string DEFAULT_POPULATION_FILENAME = "pops/population.pop";

// XML-ish blocks for parsing files.
const std::string FILE_BLOCK_POPULATION = "<population>\n";
const std::string FILE_BLOCK_POPULATION_SIZE = "\t<size>\n";
const std::string FILE_BLOCK_NUM_STATE_VARS = "\t<NUM_STATE_VARS>\n";
const std::string FILE_BLOCK_NUM_OUTPUT_VARS = "\t<NUM_OUTPUT_VARS>\n";
const std::string FILE_BLOCK_ORGANISM = "<organism>\n";
const std::string FILE_BLOCK_INDEX = "\t<index>\n";
const std::string FILE_BLOCK_TOTAL_STABLE_TIME = "\t<m_totalStableTime>\n";
const std::string FILE_BLOCK_NUM_SIMULATIONS = "\t<m_numSimulations>\n";
const std::string FILE_BLOCK_GENETICS = "\t<m_genetics>\n";
const std::string FILE_BLOCK_EXPRESSIONS = "\t\t<m_expressions>\n";
const std::string FILE_BLOCK_POLY = "\t\t\t<m_poly>\n";
const std::string FILE_BLOCK_LOG = "\t\t\t<m_log>\n";
const std::string FILE_BLOCK_SIN = "\t\t\t<m_sin>\n";
const std::string FILE_BLOCK_COS = "\t\t\t<m_cos>\n";
const std::string FILE_BLOCK_EXP = "\t\t\t<m_exp>\n";

// Delimitter stripped versions of the above due to getline etc stripping \n when we want to use for
// comparisons.
const std::string FILE_BLOCK_POPULATION_STRIPPPED = "<population>";
const std::string FILE_BLOCK_POPULATION_SIZE_STRIPPPED = "\t<size>";
const std::string FILE_BLOCK_NUM_STATE_VARS_STRIPPPED = "\t<NUM_STATE_VARS>";
const std::string FILE_BLOCK_NUM_OUTPUT_VARS_STRIPPPED = "\t<NUM_OUTPUT_VARS>";
const std::string FILE_BLOCK_ORGANISM_STRIPPPED = "<organism>";
const std::string FILE_BLOCK_INDEX_STRIPPPED = "\t<index>";
const std::string FILE_BLOCK_TOTAL_STABLE_TIME_STRIPPPED = "\t<m_totalStableTime>";
const std::string FILE_BLOCK_NUM_SIMULATIONS_STRIPPPED = "\t<m_numSimulations>";
const std::string FILE_BLOCK_GENETICS_STRIPPPED = "\t<m_genetics>";
const std::string FILE_BLOCK_EXPRESSIONS_STRIPPPED = "\t\t<m_expressions>";
const std::string FILE_BLOCK_POLY_STRIPPPED = "\t\t\t<m_poly>";
const std::string FILE_BLOCK_LOG_STRIPPPED = "\t\t\t<m_log>";
const std::string FILE_BLOCK_SIN_STRIPPPED = "\t\t\t<m_sin>";
const std::string FILE_BLOCK_COS_STRIPPPED = "\t\t\t<m_cos>";
const std::string FILE_BLOCK_EXP_STRIPPPED = "\t\t\t<m_exp>";