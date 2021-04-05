#pragma once

// This header prototypes the classes needed for the genetics of the controller, i.e
// base controller structure, mutations, reproduction, populations, etc.

#include <vector>
#include <random>
#include <math.h>

#include "globals.h"

// Expressions are for each of the 3 state vars comprising an equation (Genetic).
// Member vectors are of size 2*MAX_EXPRESSION_SUBLENGTHS.
struct Expression {
  // Construct a random subexpression (single variable).
  Expression();

  // Coefficients for polynomial A1*x^B1 ...
  std::vector<double> m_poly;
  
  // Coefficients for log, A1*log(B1*x) ...
  std::vector<double> m_log;

  // Coefficients for sin, A1*sin(B1*x) ...
  std::vector<double> m_sin;
  
  // Coefficients for cos, A1*cos(B1*x) ...
  std::vector<double> m_cos;
  
  // Coefficients for exp, A1*exp(B1*x) ...
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
};


// The individual organisms of a population can mutate and reproduce.
struct Organism {
    // Construct an organism with 8 Gene members in m_genetics.
    Organism();

    // Mutate the current organism.
    void mutate();
    
    // Creates a child organism with the genetics of *this and partner organism.
    Organism reproduce(const Organism& partner) const;
    
    // Holds the genetics of the organism.  We have NUM_STATE_VARS inputs of interest:
    // the x and y zmps of the foot holding the load (left), and the z com.
    // We have 8 outputs of interest (left and right shoulder pitch / roll, 
    // elbow yaw / roll.  So we have 8 varying equations (genetics) each with
    // NUM_STATE_VARS input variables, ie vector is size 8.
    std::vector<Gene> m_genetics;
    
    // The fitness of the organism (determined by average time before falling in simulation).
    // Returns m_totalStableTime/m_numSimulations.
    double getFitness();
    
    // Total time stable accross all simulations.
    double m_totalStableTime;

    // The number of times this controller has been simulated.
    int m_numSimulations;
   
    // Defining comparison operators of organism for sorting / pruning purposes.
    bool operator < (const Organism& rhs) const;
    bool operator > (const Organism& rhs) const;
};


// Holds a population / generation of controllers (organisms).
struct Population {
  // Initialize a population with n random organisms.
  Population(const int& n);
  
  // This will be kept in descending order.  > operator for organism will be defined.
  std::vector<Organism> m_organisms;





};