#pragma once

// This header prototypes the classes needed for the genetics of the controller, i.e
// base controller structure, mutations, reproduction, populations, etc.

#include <vector>
#include <random>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>

#include "globals.h"
#include "utilities.h"

// Expressions are for each of the 3 input vars comprising an equation (Genetic).
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
  
  // Coefficients for cos, A1*cos(B1*x)-A1 ...
  std::vector<double> m_cos;
  
  // Coefficients for exp, A1*exp(B1*x) ... can function as a constant with low B1.
  std::vector<double> m_exp;
};

// A generic nonlinear gene comprised of polynomial, sin, cos, and exp terms.
// The genetics of an organism.  Returns a value based on the m_numInputVars input variables.
struct Gene {  
  // Construct a Gene with i (m_numInputVars) expression objects in m_expressions.
  Gene(const int& i);
  
  // Gene needs to know how many input variables in order to push back the correct
  // number of Expression objects onto m_expressions.
  int m_numInputVars;
  
  // Takes the current values of the input variables, and returns the output per m_expressions.
  // Takes a vector of doubles of size m_numInputVars.
  double calculateValue(const std::vector<double>& x) const;
  
  // One Expression per variable, so size m_numInputVars vector.
  std::vector<Expression> m_expressions;
};

// The gait form of a gene, per equation 14 in https://journals.sagepub.com/doi/full/10.5772/58845.
// See equation image in op of https://github.com/Inexorably/naoTest/issues/12.
struct GaitGene {
  // Construct a Gene with i (m_numInputVars) expression objects in m_expressions.
  // This is a CPG / gait generator, so this does depends on input alpha and time t.
  GaitGene();
  
  // Note that we accept an x vector to match the style of the Gene class, but the vector
  // should be of size 2 as per our comments on the GaitGene constructor.
  // ie, the input should be alpha and time t.
  std::vector<double> calculateValue(const std::vector<double>& x) const;
  
  // The variables of interest are:
  /* A1
   * A2 = 0
   * A3
   * A4
   * A5  
   */
  // We store these constants in m_constants.
  // Note that omega = 2pi*alpha (the input), and that K_r = alpha.
  // Note that we define this in the top level Gene class so that we can push different types of Genes
  // on an organism level.
  std::vector<double> m_constants;
};

// The individual organisms of a population can mutate and reproduce.
struct Organism {
    // Construct an organism with o (m_numOutputVars) Gene members in m_genetics.
    Organism(const int& i, const int& o);

    // Mutate the current organism.
    void mutate();
    
    // Creates a child organism with the genetics of *this and partner organism.
    Organism reproduce(const Organism& partner) const;
    
    // Holds the genetics of the organism.  We have m_numInputVars inputs of interest.
    // We have m_numOutputVars outputs of interest.
    // So we have m_numOutputVars varying equations (genetics) each with
    // m_numInputVars input variables, ie vector is size m_numOutputVars.
    std::vector<Gene> m_genetics;
    
    // Same values as parent Population if organism is part of Population object.
    // Number of input/output variables.  Moved inside the Population struct so that we can have multiple objects
    // of differing input / outputs, which using globals for num_input/output vars prevented.
    int m_numInputVars;
    int m_numOutputVars;
    
    // Mutation probability between 0 to 1 for a given gene.
    double m_chanceMutation;
    
    // The fitness of the organism, determined by average time before falling in simulation and
    // the average distance of zmp coordinates from the origin.
    double getFitness() const;
    
    // Print the fitness components to console so we can see if we need to adjust the weights.
    // Returns the fitness components in a vector of doubles.
    // Order: time, zmp, trans_x, comv_yz.
    std::vector<double> getFitnessComponents() const;
    
    // Print the fitness components to console so we can see if we need to adjust the weights.
    void printFitnessComponents() const;
    
    // Total time stable accross all simulations, in seconds.
    double m_totalStableTime;

    // The number of times this controller has been simulated.
    int m_numSimulations;
    
    // The TOTAL zmp distance from 0, 0, ie if zmp is at 1, 1 for 2 seconds, the
    // m_totalZMPDistance value would be sqrt(2)*2.  Units are meters.
    double m_totalZMPDistance;
    
    // Total translated distance in x.
    double m_totalTranslationX;
    
    // The total 'velocity' of the COM.
    // At each control step, the current COM velocity * control time step is added to this
    // variable.  The units are m/s * s.  Ie, m_totalCOMVelocity/m_totalStableTime gives the
    // average COM velocity.
    double m_totalCOMVelocity;
   
    // Defining comparison operators of organism for sorting / pruning purposes.
    // Compares by getFitness().
    bool operator < (const Organism& rhs) const;
    bool operator > (const Organism& rhs) const;
    
    // Save the current organism to a file.
    void save(const std::string& filename) const;
};

// A version of the Organism struct with GaitGene members for gait controller evolution.
struct GaitOrganism {
    // Construct an organism with a single GaitGene member.
    GaitOrganism();

    // Mutate the current organism.
    void mutate();
    
    // Creates a child organism with the genetics of *this and partner GaitOrganism.
    GaitOrganism reproduce(const GaitOrganism& partner) const;
    
    // Holds the genetics of the organism.
    GaitGene m_gaitGene;

    // Same values as parent Population if organism is part of Population object.
    // Number of input/output variables.  Moved inside the Population struct so that we can have multiple objects
    // of differing input / outputs, which using globals for num_input/output vars prevented.
    int m_numInputVars;
    int m_numOutputVars;
    
    // Mutation probability between 0 to 1 for a given gene.
    double m_chanceMutation;
    
    // The fitness of the organism, determined by average time before falling in simulation and
    // the average distance of zmp coordinates from the origin.
    double getFitness() const;
    
    // Total time stable accross all simulations, in seconds.
    double m_totalStableTime;

    // The number of times this controller has been simulated.
    int m_numSimulations;
    
    // The TOTAL zmp distance from 0, 0, ie if zmp is at 1, 1 for 2 seconds, the
    // m_totalZMPDistance value would be sqrt(2)*2.  Units are meters.
    double m_totalZMPDistance;
    
    // Total translated distance in x.
    double m_totalTranslationX;
   
    // Defining comparison operators of GaitOrganism for sorting / pruning purposes.
    // Compares by getFitness().
    bool operator < (const GaitOrganism& rhs) const;
    bool operator > (const GaitOrganism& rhs) const;
    
    // Save the current organism to a file.
    void save(const std::string& filename) const;
};


// Holds a population / generation of controllers (organisms).
struct Population {
  // Initialize a population with n random organisms, i input vars, and o output vars.
  Population(const int& n, const int& i, const int& o);
  
  // This will be kept in descending order.  > operator for organism will be defined.
  std::vector<Organism> m_organisms;
  
  // The current generation of this population.
  int m_generation;
  
  // The total runtime in seconds of the population, including all preceeding generations.
  double m_runtime;
  
  // Number of input/output variables.  Moved inside the Population struct so that we can have multiple objects
  // of differing input / outputs, which using globals for num_input/output vars prevented.
  int m_numInputVars;
  int m_numOutputVars;
  
  // Mutation probability between 0 to 1 for a given gene.
  double m_chanceMutation;
  
  // The population size.  This is the size of m_organisms, and we store this for when we are manipulating
  // m_organisms.
  unsigned int m_numOrganisms;
  
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
  
  // Saves the population to the default filename DEFAULT_POPULATION_FILENAME.
  void save() const;
  
  // Loads the population from a given file.  If ignoreHistory is true, do not load the
  // m_totalStableTime and m_numSimulation members from the file. 
  void load(const std::string& filename, const bool& ignoreHistory);
  
  // Loads the population from a given file.
  // TODO: Load validation such as confirming right population size etc.
  void load(const std::string& filename);
  
  // Loads the population from the default filename DEFAULT_POPULATION_FILENAME.
  void load();
};

// A population struct holding n GaitOrganisms.
struct GaitPopulation {
  // Initialize a GaitPopulation with n random organisms.
  GaitPopulation(const int& n);
  
  // This will be kept in descending order.  > operator for organism will be defined.
  std::vector<GaitOrganism> m_organisms;
  
  // The current generation of this population.
  int m_generation;
  
  // The total runtime in seconds of the population, including all preceeding generations.
  double m_runtime;
  
  // Number of input/output variables.  Moved inside the Population struct so that we can have multiple objects
  // of differing input / outputs, which using globals for num_input/output vars prevented.
  int m_numInputVars;
  int m_numOutputVars;
  
  // Mutation probability between 0 to 1 for a given gene.
  double m_chanceMutation;
  
  // The population size.  This is the size of m_organisms, and we store this for when we are manipulating
  // m_organisms.
  unsigned int m_numOrganisms;
  
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
  
  // Saves the population to the default filename DEFAULT_POPULATION_FILENAME.
  void save() const;
  
  // Loads the population from a given file.  If ignoreHistory is true, do not load the
  // m_totalStableTime and m_numSimulation members from the file. 
  void load(const std::string& filename, const bool& ignoreHistory);
  
  // Loads the population from a given file.
  // TODO: Load validation such as confirming right population size etc.
  void load(const std::string& filename);
  
  // Loads the population from the default filename DEFAULT_POPULATION_FILENAME.
  void load();
};