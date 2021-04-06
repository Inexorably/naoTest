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

// TODO: Depreciate.
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

// Takes the current values of the state variables, and returns the output per m_expressions.
// Takes a vector of doubles of size NUM_STATE_VARS.
double Gene::calculateValue(const std::vector<double>& x) const {
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

// Creates a child organism with the genetics of *this and partner organism.
Organism Organism::reproduce(const Organism& partner) const {
  // Randomly mix the expressions of each partner.
  
  // Create the child object to be returned.
  // TODO: Can optimize by adding constructor that skips rng coefficient creation since
  // immediately overwritten here, but performance gains are likely negligible.
  Organism child;
  
  // Loop through the genetics (NUM_OUTPUT_VARS many) of the child object.
  for (int i = 0; i < NUM_OUTPUT_VARS; i++) {
    // Loop through the expressions (NUM_STATE_VARS many) in each gene.
    for (int j = 0; j < NUM_STATE_VARS; j++) {
      // For each expression, set the child expression to one of the parent's.
      // If true with 50% chance, set to parent one (*this).
      if (trueWithProbability(0.5)) {
        child.m_genetics[i].m_expressions[j] = m_genetics[i].m_expressions[j];
      }
      // else set to parent two (partner).
      else {
        child.m_genetics[i].m_expressions[j] = partner.m_genetics[i].m_expressions[j];
      }
    }
  }
  
  // Return the child.
  return child;
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

// Save the current organism to a file.
void Organism::save(const std::string& filename) const {
  // Open a file stream.
  std::ofstream outfile;
  outfile.open(filename);

  outfile << "<organism>\n";
  outfile << "\t<index>\n";
  outfile << "\t\t0\n";
  outfile << "\t<m_totalStableTime>\n";
  outfile << "\t\t" << std::to_string(m_totalStableTime) << '\n';
  outfile << "\t<m_numSimulations>\n";
  outfile << "\t\t" << std::to_string(m_numSimulations) << '\n';

  // Write each organisms m_genetics, of size NUM_OUTPUT_VARS.
  outfile << "\t<m_genetics>\n";
  
  // Loop through the genetics of each organism.
  for (Gene g : m_genetics) {
    // Write the expressions of each gene.
    outfile << "\t\t<m_expressions>\n";
    
    // Loop through the expressions and write the variable values.
    // m_expressions is of size NUM_STATE_VARS.
    for (Expression e : g.m_expressions) {
      // Loop through and write m_poly coeffs
      outfile << "\t\t\t<m_poly>\n";
      for (double d : e.m_poly) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
      
      // Loop through and write m_log coeffs
      outfile << "\t\t\t<m_log>\n";
      for (double d : e.m_log) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
      
      // Loop through and write m_sin coeffs
      outfile << "\t\t\t<m_sin>\n";
      for (double d : e.m_sin) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
      
      // Loop through and write m_cos coeffs
      outfile << "\t\t\t<m_cos>\n";
      for (double d : e.m_cos) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
      
      // Loop through and write m_exp coeffs
      outfile << "\t\t\t<m_exp>\n";
      for (double d : e.m_exp) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
    }
  }
}

///////////////// Population /////////////////////////////////

// Initialize a population with n random organisms.
Population::Population(const int& n) {
  for (int i = 0; i < n; i++) {
    m_organisms.push_back(Organism());
  }
}

// Sort the m_organisms vector in descending order.
void Population::sortOrganisms() {
  std::sort(m_organisms.begin(), m_organisms.end(), std::greater <>());
}

// Breed the current population with 2 partners creating 1 child organism.
void Population::reproduceOrganisms() {
  // Shuffle the population so that we can easily select random partners.
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(m_organisms.begin(), m_organisms.end(), g);

  // Create a vector to store the children.
  std::vector<Organism> children;

  // Loop through the shuffled organisms and breed them.
  for (size_t i = 0; i < m_organisms.size(); i += 2) {
    children.push_back(m_organisms[i].reproduce(m_organisms[i+1]));
  }
  
  // Append the children organism vector to the current population.
  m_organisms.insert(std::end(m_organisms), std::begin(children), std::end(children));
  
  // Sort the m_organisms vector.  While not useful now, can be beneficial in future
  // depending on design choice and is of negligible cpu cost.
  sortOrganisms();
}

// Copy and mutate random members until m_organsisms.size() + number of mutated members == POPULATION_SIZE.
// Then, append the mutated members to m_organsisms.
void Population::mutateOrganisms() {
  // Create a vector to store the mutated copies.
  std::vector<Organism> mutations;
  
  // Create an int rng device to select random members of m_organisms with.
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> distI(0, m_organisms.size());

  // While we do not have enough organisms for the population size ...
  while (m_organisms.size() + mutations.size() < POPULATION_SIZE) {
    // Create a copy of a random organism in the population.
    Organism copy = m_organisms[distI(mt)];
    
    // Mutate the copy organism.
    copy.mutate();
    
    // Push back the mutated copy.
    mutations.push_back(copy);
  }
  
  // Append the mutations vector to m_organisms.
  m_organisms.insert(std::end(m_organisms), std::begin(mutations), std::end(mutations));
  
  // Sort the m_organisms vector.  While not useful now, can be beneficial in future
  // depending on design choice and is of negligible cpu cost.
  sortOrganisms();
}

// Saves the population to a given filename.
void Population::save(const std::string& filename) const {
  // Open a file stream.
  std::ofstream outfile;
  outfile.open(filename);
  
  // Begin by writing population header information.
  outfile << "<population>\n";
  outfile << "\t<size>\n";
  outfile <<"\t\t" << std::to_string(m_organisms.size()) << '\n';
  outfile << "\t<NUM_STATE_VARS>\n";
  outfile << "\t\t" << std::to_string(NUM_STATE_VARS) << '\n';
  outfile << "\t<NUM_OUTPUT_VARS>\n";
  outfile << "\t\t" << std::to_string(NUM_OUTPUT_VARS) << '\n';
  
  // Output each organism to the file.
  for (unsigned int i = 0; i < m_organisms.size(); i++) {
    outfile << "<organism>\n";
    outfile << "\t<index>\n";
    outfile << "\t\t" << std::to_string(i) << '\n';
    outfile << "\t<m_totalStableTime>\n";
    outfile << "\t\t" << std::to_string(m_organisms[i].m_totalStableTime) << '\n';
    outfile << "\t<m_numSimulations>\n";
    outfile << "\t\t" << std::to_string(m_organisms[i].m_numSimulations) << '\n';

    // Write each organisms m_genetics, of size NUM_OUTPUT_VARS.
    outfile << "\t<m_genetics>\n";
    
    // Loop through the genetics of each organism.
    for (Gene g : m_organisms[i].m_genetics) {
      // Write the expressions of each gene.
      outfile << "\t\t<m_expressions>\n";
      
      // Loop through the expressions and write the variable values.
      // m_expressions is of size NUM_STATE_VARS.
      for (Expression e : g.m_expressions) {
        // Loop through and write m_poly coeffs
        outfile << "\t\t\t<m_poly>\n";
        for (double d : e.m_poly) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
        
        // Loop through and write m_log coeffs
        outfile << "\t\t\t<m_log>\n";
        for (double d : e.m_log) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
        
        // Loop through and write m_sin coeffs
        outfile << "\t\t\t<m_sin>\n";
        for (double d : e.m_sin) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
        
        // Loop through and write m_cos coeffs
        outfile << "\t\t\t<m_cos>\n";
        for (double d : e.m_cos) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
        
        // Loop through and write m_exp coeffs
        outfile << "\t\t\t<m_exp>\n";
        for (double d : e.m_exp) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
      }
    }
  }
  
  // Close the ofstream.
  outfile.close();
}

// Saves the population to the default filename DEFAULT_POPULATION_FILENAME.
void Population::save() const {
  save(DEFAULT_POPULATION_FILENAME);
}