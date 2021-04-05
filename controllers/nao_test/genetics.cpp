#include "genetics.h"

///////////////// Expression /////////////////////////////////

// Construct a random subexpression (single variable).
Expression::Expression() {
  // Upon construction, fill the member vectors to represent some random expressions.

  // Create a double and int rng device.
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> distI(0, EXPRESSION_MAX_SUBLENGTHS);
  std::uniform_int_distribution<int> distExp(EXPRESSION_CONST_EXP_MIN, EXPRESSION_CONST_EXP_MAX);
  std::uniform_real_distribution<double> distD(EXPRESSION_CONST_MIN, EXPRESSION_CONST_MAX);
  
  // For each of the member vectors, fill randomly.
  
  // m_poly
  int temp = distI(mt);             //The number of elements being added.
  for (int i = 0; i < temp; i++) {
    // For each iteration, we push back two random doubles, since the basic forms of each subexpression
    // require an A and B value.  Could simply loop i < 2*temp, but this better represents the idea.
    m_poly.push_back(distD(mt));
    m_poly.push_back(distExp(mt));
  }
  
  /*
  // m_log
  temp = distI(mt);             //The number of elements being added.
  for (int i = 0; i < temp; i++) {
    // For each iteration, we push back two random doubles, since the basic forms of each subexpression
    // require an A and B value.  Could simply loop i < 2*temp, but this better represents the idea.
    m_log.push_back(distD(mt));
    m_log.push_back(distD(mt));
  }
  */
  
  // m_sin
  temp = distI(mt);             //The number of elements being added.
  for (int i = 0; i < temp; i++) {
    // For each iteration, we push back two random doubles, since the basic forms of each subexpression
    // require an A and B value.  Could simply loop i < 2*temp, but this better represents the idea.
    m_sin.push_back(distD(mt));
    m_sin.push_back(distD(mt));
  }
  
  /*
  // m_cos
  temp = distI(mt);             //The number of elements being added.
  for (int i = 0; i < temp; i++) {
    // For each iteration, we push back two random doubles, since the basic forms of each subexpression
    // require an A and B value.  Could simply loop i < 2*temp, but this better represents the idea.
    m_cos.push_back(distD(mt));
    m_cos.push_back(distD(mt));
  }
  */
  
  // m_exp
  temp = distI(mt);             //The number of elements being added.
  for (int i = 0; i < temp; i++) {
    // For each iteration, we push back two random doubles, since the basic forms of each subexpression
    // require an A and B value.  Could simply loop i < 2*temp, but this better represents the idea.
    m_exp.push_back(distD(mt));
    m_exp.push_back(distExp(mt));
  }
}

///////////////// Gene /////////////////////////////////

// Construct a genome with 3 expression objects in m_expressions.
Gene::Gene() {
  // Push back NUM_STATE_VARS subexpressions - one per state variable.
  for (int i = 0; i < NUM_STATE_VARS; i++) { 
    Expression temp;         // Default constructor creates a random expression.
    m_expressions.push_back(temp);
  }
}

// Takes the current values of the state variables, and returns the output per m_expressions.
// Note: this is currently hard coded to 3 state vars.
double Gene::calculateValue(const double& x1, const double& x2, const double& x3) const {
  // Array so we can iterate through state variables.
  double x[3] = {x1, x2, x3};
  
  // The value to be returned.
  double result = 0;
  
  // Iterate through the state variables as per the subexpressions.
  for (int i = 0; i < NUM_STATE_VARS; i++) {
    // For each state variable, apply the various nonlinear expressions.
    
    // poly
    for (size_t j = 0; j < m_expressions[i].m_poly.size(); j += 2) {
      result += m_expressions[i].m_poly[j] * pow(x[i], static_cast<int>(m_expressions[i].m_poly[j+1]));
    }
    
    //std::cout << "poly: " << result << std::endl;
    
    // log
    for (size_t j = 0; j < m_expressions[i].m_log.size(); j += 2) {
      result += m_expressions[i].m_log[j] * log(abs(x[i] * m_expressions[i].m_log[j+1]));
    }
    
    //std::cout << "log: " << result << std::endl;
    
    // sin
    for (size_t j = 0; j < m_expressions[i].m_sin.size(); j += 2) {
      result += m_expressions[i].m_sin[j] * sin(x[i] * m_expressions[i].m_sin[j+1]);
    }
    
    //std::cout << "sin: " << result << std::endl;
    
    // cos
    for (size_t j = 0; j < m_expressions[i].m_cos.size(); j += 2) {
      result += m_expressions[i].m_cos[j] * cos(x[i] * m_expressions[i].m_cos[j+1]);
    }
    
    //std::cout << "cos: " << result << std::endl;
    
    // exp
    for (size_t j = 0; j < m_expressions[i].m_exp.size(); j += 2) {
      result += m_expressions[i].m_exp[j] * exp(x[i] * m_expressions[i].m_exp[j+1]) - m_expressions[i].m_exp[j];
    }
    
    //std::cout << "exp: " << result << std::endl;
  }

  return result;
}

///////////////// Organism /////////////////////////////////

// Construct an organism with NUM_OUTPUT_VARS Gene members in m_genetics.
Organism::Organism() {
  for (int i = 0; i < NUM_OUTPUT_VARS; i++) { 
    Gene temp;
    m_genetics.push_back(temp);
  }
  
  m_totalStableTime = 0;
  m_numSimulations = 0;
}

// Mutate the current organism.  Each expression has a MUTATION_CHANCE chance of changing.
void Organism::mutate() {
  //Loop through the m_genetics and m_expressions and randomly change the expressions.
  // m_genetics is of size NUM_OUTPUT_VARS, and each Gene g.m_expressions is of size NUM_STATE_VARS.
  for (Gene& g : m_genetics) {
    for (Expression& e : g.m_expressions) {
      // Randomly replace one of the expressions with a new one.
      if (trueWithProbability(MUTATION_CHANCE)) {
        e = Expression();
      }
    }
  }
}

// The fitness of the organism (determined by average time before falling in simulation).
// Returns m_totalStableTime/m_numSimulations (the average stable time).
double Organism::getFitness() const {
  return m_totalStableTime/m_numSimulations;
}

// Defining comparison operators of organism for sorting / pruning purposes.
// Compares by getFitness().
bool Organism::operator < (const Organism& rhs) const {
  return getFitness() < rhs.getFitness();
}

// Defining comparison operators of organism for sorting / pruning purposes.
// Compares by getFitness().
bool Organism::operator > (const Organism& rhs) const {
  return getFitness() > rhs.getFitness();
}

///////////////// Population /////////////////////////////////

// Initialize a population with n random organisms.
Population::Population(const int& n) {
  for (int i = 0; i < n; i++) {
    m_organisms.push_back(Organism());
  }
}

// Sort the m_organisms vector using Organism::operator<.
void Population::sortOrganisms() {
  std::sort(m_organisms.begin(), m_organisms.end());
}