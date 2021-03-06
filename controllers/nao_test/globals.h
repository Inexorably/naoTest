#pragma once

#include <string>

/////////////////////////////////////// Simulation constants /////////////////////////////////////////////

// If the total grf / force feet supporting drops below FOOT_FORCE_MIN, the robot has fallen over
// and we end the simulation.
const double FOOT_FORCE_MIN = 30;

// Max simulation time before breaking in seconds.
const int SIMULATION_TIME_MAX = 60;

// The number of generations of our population we iterate through.
const int NUM_GENERATIONS = 10001;

// The number of steps between each control attempt.  Higher values allow less noise but slower response.
const int STEPS_PER_CONTROL = 4;

// Some very low value to act as a min bound for the fitness score, that will never actually be reached.
const int FITNESS_FLOOR = -9999999;

/////////////////////////////////////// NAO physical properties //////////////////////////////////////////

const double FOOT_WIDTH = 0.08;  // Per http://simspark.sourceforge.net/wiki/index.php/Models.
const double FOOT_LENGTH = 0.16; // In meters.

/////////////////////////////////////// Random motion / plant ////////////////////////////////////////////

// Modifier for delta position target in motion plant.
// Motors are servos, so target position is in rads.
const double LAMBDA = 0.02;

/////////////////////////////////////// Genetics constants ///////////////////////////////////////////////

// Determines the max number of instances of a type of function on a input variable an expression can have.
// See Expression class in "genetics.h".
const int EXPRESSION_MAX_SUBLENGTHS = 2;

// Expression constant min / max values for non exponential values.
const double EXPRESSION_CONST_MIN = -2;
const double EXPRESSION_CONST_MAX = 2;

// Exponential constant bounds.  Using doubles for exponents can cause domain errors (ie -2^2.34) and NaN.
// Allowing negative exponents creates inf values easily (ie 0.000000001^-5).
const int EXPRESSION_CONST_EXP_MAX = 5;
const int EXPRESSION_CONST_EXP_MIN = 0;

// GaitGene constant min and max values.
const double GAITGENE_CONST_MIN = -2.5;
const double GAITGENE_CONST_MAX = 2.5;

// FITNESS_WEIGHT_ZMP_TRANSITION_TIME: average stable time at which we transition
// to highly valuing zmp stability.  To have any affect, must be less than SIMULATION_TIME_MAX.
// FITNESS_WEIGHT_ZMP_TRANSITION_COEF: the multiplier by which we multiply the zmp 
// component to add more weighting to balancing.
const double FITNESS_WEIGHT_ZMP_TRANSITION_TIME = 30;
const double FITNESS_WEIGHT_ZMP_TRANSITION_COEF = 0;

// Weighting values for Organism::getFitness().
const double FITNESS_WEIGHT_ZMP_COEF = 200;
const double FITNESS_WEIGHT_TRANSLATION_X_COEF = 3;
const double FITNESS_WEIGHT_TIME_COEF = 0.2;
const double FITNESS_WEIGHT_COMV_COEF = 1300;

// The mutation chance for each generation_i is of the form:
// chance_i = (num_inputs*num_outputs)^-1 * (stddev_0 / s_i)^c
// Where stddev_0 is a base approximation of the standard deviation for a completely new population
// and c is a constant.  
const double MUTATION_CHANCE_STD_DEV_BASE = 0.785;
extern int MUTATION_CHANCE_C;

/////////////////////////////////////// File constants ///////////////////////////////////////////////

// Filename / paths for saving historical generation data.
// Not explictly constant so can change in main and do various conditions.
extern std::string DEFAULT_POPULATION_FILENAME;
extern std::string DEFAULT_GENERATION_FILE_PREFIX;

// XML-ish blocks for parsing files.
const std::string FILE_BLOCK_POPULATION = "<population>\n";
const std::string FILE_BLOCK_RUNTIME = "\t<runtime>\n";
const std::string FILE_BLOCK_GENERATION = "\t<generation>\n";
const std::string FILE_BLOCK_POPULATION_SIZE = "\t<m_numOrganisms>\n";
const std::string FILE_BLOCK_NUM_INPUT_VARS = "\t<m_numInputVars>\n";
const std::string FILE_BLOCK_NUM_OUTPUT_VARS = "\t<m_numOutputVars>\n";
const std::string FILE_BLOCK_CHANCE_MUTATION = "\t<m_chanceMutation>\n";
const std::string FILE_BLOCK_ORGANISM = "<organism>\n";
const std::string FILE_BLOCK_GAIT_ORGANISM = "<GaitOrganism>\n";
const std::string FILE_BLOCK_INDEX = "\t<index>\n";
const std::string FILE_BLOCK_TOTAL_STABLE_TIME = "\t<m_totalStableTime>\n";
const std::string FILE_BLOCK_NUM_SIMULATIONS = "\t<m_numSimulations>\n";
const std::string FILE_BLOCK_TOTAL_ZMP_DISTANCE = "\t<m_totalZMPDistance>\n";
const std::string FILE_BLOCK_TOTAL_TRANSLATION_X = "\t<m_totalTranslationX>\n";
const std::string FILE_BLOCK_TOTAL_COM_VELOCITY = "\t<m_totalCOMVelocity>\n";
const std::string FILE_BLOCK_GENETICS = "\t<m_genetics>\n";
const std::string FILE_BLOCK_EXPRESSIONS = "\t\t<m_expressions>\n";
const std::string FILE_BLOCK_CONSTANTS = "\t\t<m_constants>\n";
const std::string FILE_BLOCK_POLY = "\t\t\t<m_poly>\n";
const std::string FILE_BLOCK_LOG = "\t\t\t<m_log>\n";
const std::string FILE_BLOCK_SIN = "\t\t\t<m_sin>\n";
const std::string FILE_BLOCK_COS = "\t\t\t<m_cos>\n";
const std::string FILE_BLOCK_EXP = "\t\t\t<m_exp>\n";

// Delimitter stripped versions of the above due to getline etc stripping \n when we want to use for
// comparisons.
const std::string FILE_BLOCK_POPULATION_STRIPPED = "<population>";
const std::string FILE_BLOCK_RUNTIME_STRIPPED = "\t<runtime>";
const std::string FILE_BLOCK_POPULATION_SIZE_STRIPPED = "\t<m_numOrganisms>";
const std::string FILE_BLOCK_NUM_INPUT_VARS_STRIPPED = "\t<m_numInputVars>";
const std::string FILE_BLOCK_NUM_OUTPUT_VARS_STRIPPED = "\t<m_numOutputVars>";
const std::string FILE_BLOCK_ORGANISM_STRIPPED = "<organism>";
const std::string FILE_BLOCK_INDEX_STRIPPED = "\t<index>";
const std::string FILE_BLOCK_TOTAL_STABLE_TIME_STRIPPED = "\t<m_totalStableTime>";
const std::string FILE_BLOCK_NUM_SIMULATIONS_STRIPPED = "\t<m_numSimulations>";
const std::string FILE_BLOCK_TOTAL_ZMP_DISTANCE_STRIPPED = "\t<m_totalZMPDistance>";
const std::string FILE_BLOCK_TOTAL_TRANSLATION_X_STRIPPED = "\t<m_totalTranslationX>";
const std::string FILE_BLOCK_TOTAL_COM_VELOCITY_STRIPPED = "\t<m_totalCOMVelocity>";
const std::string FILE_BLOCK_GENETICS_STRIPPED = "\t<m_genetics>";
const std::string FILE_BLOCK_EXPRESSIONS_STRIPPED = "\t\t<m_expressions>";
const std::string FILE_BLOCK_CONSTANTS_STRIPPED = "\t\t<m_constants>";
const std::string FILE_BLOCK_POLY_STRIPPED = "\t\t\t<m_poly>";
const std::string FILE_BLOCK_LOG_STRIPPED = "\t\t\t<m_log>";
const std::string FILE_BLOCK_SIN_STRIPPED = "\t\t\t<m_sin>";
const std::string FILE_BLOCK_COS_STRIPPED = "\t\t\t<m_cos>";
const std::string FILE_BLOCK_EXP_STRIPPED = "\t\t\t<m_exp>";