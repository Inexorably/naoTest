#pragma once

// This header prototypes the classes needed for the genetics of the controller, i.e
// base controller structure, mutations, reproduction, populations, etc.

#include <vector>
#include <random>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <string>

#include "globals.h"
#include "utilities.h"

// Expressions are for each of the 3 state vars comprising an equation (Genetic).
// Member vectors are of size 2*EXPRESSION_MAX_SUBLENGTHS.
struct Expression {
  // Construct a random subexpression (single variable).
  Expression();

  // Coefficients for polynomial A1*x^B1 ...
  std::vector<double> m_poly;
  
  // Coefficients for log, A1*log(|B1*x|) ...
  // TEMPORARILY DISABLED FOR TESTING if non f(0) == 0 functions are appropriate.
  std::vector<double> m_log;

  // Coefficients for sin, A1*sin(B1*x) ...
  std::vector<double> m_sin;
  
  // Coefficients for cos, A1*cos(B1*x) ...
  // TEMPORARILY DISABLED FOR TESTING if non f(0) == 0 functions are appropriate.
  std::vector<double> m_cos;
  
  // Coefficients for exp, A1*exp(B1*x)-A1 ...
  std::vector<double> m_exp;
};

// The genetics of an organism.  Returns a value based on the NUM_STATE_VARS state variables.
struct Gene {
  // Construct a genome with NUM_STATE_VARS expression objects in m_expressions.
  Gene();

  // One Expression per variable, so size NUM_STATE_VARS vector.
  std::vector<Expression> m_expressions;

  // Takes the current values of the state variables, and returns the output per m_expressions.
  // Note: currently hard coded to 3 state vars.
  double calculateValue(const double& x1, const double& x2, const double& x3) const;
  
  // Takes the current values of the state variables, and returns the output per m_expressions.
  // Takes a vector of doubles of size NUM_STATE_VARS.
  double calculateValue(const std::vector<double>& x) const;
};


// The individual organisms of a population can mutate and reproduce.
struct Organism {
    // Construct an organism with NUM_OUTPUT_VARS Gene members in m_genetics.
    Organism();

    // Mutate the current organism.
    void mutate();
    
    // Creates a child organism with the genetics of *this and partner organism.
    Organism reproduce(const Organism& partner) const;
    
    // Holds the genetics of the organism.  We have NUM_STATE_VARS inputs of interest.
    // We have NUM_OUTPUT_VARS outputs of interest.
    // So we have NUM_OUTPUT_VARS varying equations (genetics) each with
    // NUM_STATE_VARS input variables, ie vector is size NUM_OUTPUT_VARS.
    std::vector<Gene> m_genetics;
    
    // The fitness of the organism (determined by average time before falling in simulation).
    // Returns m_totalStableTime/m_numSimulations.
    double getFitness() const;
    
    // Total time stable accross all simulations.
    double m_totalStableTime;

    // The number of times this controller has been simulated.
    int m_numSimulations;
   
    // Defining comparison operators of organism for sorting / pruning purposes.
    // Compares by getFitness().
    bool operator < (const Organism& rhs) const;
    bool operator > (const Organism& rhs) const;
};


// Holds a population / generation of controllers (organisms).
struct Population {
  // Initialize a population with n random organisms.
  Population(const int& n);
  
  // This will be kept in descending order.  > operator for organism will be defined.
  std::vector<Organism> m_organisms;
  
  // Sort the m_organisms vector using Organism::operator<.
  void sortOrganisms();
  
  // Breed the current population with 2 partners creating 1 child organism and add the children to
  // m_organisms.  Increases size of m_organisms by 50%.
  void reproduceOrganisms();
  
  // Copy and mutate random members until m_organsisms.size() + number of mutated members == POPULATION_SIZE.
  // Then, append the mutated members to m_organsisms.
  void mutateOrganisms();
  
  // Saves the population to a given filename.
  void save(const std::string& filename) const;
  
  // Loads the population from a given file.
  void load(const std::string& filename);
};