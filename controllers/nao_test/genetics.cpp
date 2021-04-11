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
  
  // m_cos
  temp = distI(mt);             //The number of elements being added.
  for (int i = 0; i < temp; i++) {
    // For each iteration, we push back two random doubles, since the basic forms of each subexpression
    // require an A and B value.  Could simply loop i < 2*temp, but this better represents the idea.
    m_cos.push_back(distD(mt));
    m_cos.push_back(distD(mt));
  }
  
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
  // Push back NUM_INPUT_VARS subexpressions - one per input variable.
  for (int i = 0; i < NUM_INPUT_VARS; i++) { 
    Expression temp;         // Default constructor creates a random expression.
    m_expressions.push_back(temp);
  }
}

// Takes the current values of the input variables, and returns the output per m_expressions.
// Takes a vector of doubles of size NUM_INPUT_VARS.
double Gene::calculateValue(const std::vector<double>& x) const {
  // The value to be returned.
  double result = 0;
  
  // Iterate through the input variables as per the subexpressions.
  for (int i = 0; i < NUM_INPUT_VARS; i++) {
    // For each input variable, apply the various nonlinear expressions.
    
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
      result += m_expressions[i].m_cos[j] * cos(x[i] * m_expressions[i].m_cos[j+1]) - m_expressions[i].m_cos[j];
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
  m_totalZMPDistance = 0;
}

// Mutate the current organism.  Each expression has a MUTATION_CHANCE chance of changing.
void Organism::mutate() {
  //Loop through the m_genetics and m_expressions and randomly change the expressions.
  // m_genetics is of size NUM_OUTPUT_VARS, and each Gene g.m_expressions is of size NUM_INPUT_VARS.
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
    // Loop through the expressions (NUM_INPUT_VARS many) in each gene.
    for (int j = 0; j < NUM_INPUT_VARS; j++) {
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

// The fitness of the organism, determined by average time before falling in simulation and
// the average distance of zmp coordinates from the origin.
double Organism::getFitness() const {
  // Returns -1 if num_simulations == 0.
  if (m_numSimulations == 0) { // More representative of concept than implictly casting as bool.
    return -1;
  }
  
  // We base the fitness score on the following components.
  double timeComponent = m_totalStableTime/static_cast<double>(m_numSimulations);
  double zmpComponent = (m_totalStableTime*sqrt(pow(FOOT_WIDTH, 2)+pow(FOOT_LENGTH, 2))-m_totalZMPDistance)/static_cast<double>(m_numSimulations);
  double translationXComponent = 15*m_translationX/static_cast<double>(m_numSimulations);
  
  // If we exceed the transition time, greatly increase the zmp component weighting to
  // encourage the robot to minimize zmp.
  if (timeComponent > FITNESS_WEIGHT_ZMP_TRANSITION_TIME) {
    zmpComponent = zmpComponent * FITNESS_WEIGHT_ZMP_TRANSITION_COEF;
  }
  
  return FITNESS_WEIGHT_TIME_COEF*timeComponent + FITNESS_WEIGHT_TIME_ZMP*zmpComponent + FITNESS_WEIGHT_TRANSLATION_X_COEF*translationXComponent;
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

  outfile << FILE_BLOCK_ORGANISM;
  outfile << FILE_BLOCK_INDEX;
  outfile << "\t\t" << std::to_string(0) << '\n';
  outfile << FILE_BLOCK_TOTAL_STABLE_TIME;
  outfile << "\t\t" << std::to_string(m_totalStableTime) << '\n';
  outfile << FILE_BLOCK_NUM_SIMULATIONS;
  outfile << "\t\t" << std::to_string(m_numSimulations) << '\n';
  outfile << FILE_BLOCK_TOTAL_ZMP_DISTANCE;
  outfile << "\t\t" << std::to_string(m_totalZMPDistance) << '\n';
  outfile << FILE_BLOCK_TOTAL_TRANSLATION_X;
  outfile << "\t\t" << std::to_string(m_translationX) << '\n';

  // Write each organisms m_genetics, of size NUM_OUTPUT_VARS.
  outfile << FILE_BLOCK_GENETICS;
    
  // Loop through the genetics of each organism.
  for (Gene g : m_genetics) {
    // Write the expressions of each gene.
    outfile << FILE_BLOCK_EXPRESSIONS;
    
    // Loop through the expressions and write the variable values.
    // m_expressions is of size NUM_INPUT_VARS.
    for (Expression e : g.m_expressions) {
      // Loop through and write m_poly coeffs
      outfile << FILE_BLOCK_POLY;
      for (double d : e.m_poly) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
      
      // Loop through and write m_log coeffs
      outfile << FILE_BLOCK_LOG;
      for (double d : e.m_log) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
      
      // Loop through and write m_sin coeffs
      outfile << FILE_BLOCK_SIN;
      for (double d : e.m_sin) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
      
      // Loop through and write m_cos coeffs
      outfile << FILE_BLOCK_COS;
      for (double d : e.m_cos) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
      
      // Loop through and write m_exp coeffs
      outfile << FILE_BLOCK_EXP;
      for (double d : e.m_exp) {
        outfile << "\t\t\t\t" << std::to_string(d) << '\n';
      }
    }
  }
}

///////////////// Population /////////////////////////////////

// Initialize a population with n random organisms.
Population::Population(const int& n) {
  m_generation = 0;
  for (int i = 0; i < n; i++) {
    m_organisms.push_back(Organism());
  }
}

// Create a population with POPULATION_SIZE random organisms.
Population::Population() {
  m_generation = 0;
  m_runtime = 0;
  for (int i = 0; i < POPULATION_SIZE; i++) {
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
  outfile << FILE_BLOCK_POPULATION;
  outfile << FILE_BLOCK_RUNTIME;
  outfile <<"\t\t" << std::to_string(m_runtime) << '\n';
  outfile << FILE_BLOCK_GENERATION;
  outfile <<"\t\t" << std::to_string(m_generation) << '\n';
  outfile << FILE_BLOCK_POPULATION_SIZE;
  outfile <<"\t\t" << std::to_string(m_organisms.size()) << '\n';
  outfile << FILE_BLOCK_NUM_INPUT_VARS;
  outfile << "\t\t" << std::to_string(NUM_INPUT_VARS) << '\n';
  outfile << FILE_BLOCK_NUM_OUTPUT_VARS;
  outfile << "\t\t" << std::to_string(NUM_OUTPUT_VARS) << '\n';
  
  // Output each organism to the file.
  for (unsigned int i = 0; i < m_organisms.size(); i++) {
    outfile << FILE_BLOCK_ORGANISM;
    outfile << FILE_BLOCK_INDEX;
    outfile << "\t\t" << std::to_string(i) << '\n';
    outfile << FILE_BLOCK_TOTAL_STABLE_TIME;
    outfile << "\t\t" << std::to_string(m_organisms[i].m_totalStableTime) << '\n';
    outfile << FILE_BLOCK_NUM_SIMULATIONS;
    outfile << "\t\t" << std::to_string(m_organisms[i].m_numSimulations) << '\n';
    outfile << FILE_BLOCK_TOTAL_ZMP_DISTANCE;
    outfile << "\t\t" << std::to_string(m_organisms[i].m_totalZMPDistance) << '\n';
    outfile << FILE_BLOCK_TOTAL_TRANSLATION_X;
    outfile << "\t\t" << std::to_string(m_organisms[i].m_translationX) << '\n';

    // Write each organisms m_genetics, of size NUM_OUTPUT_VARS.
    outfile << FILE_BLOCK_GENETICS;
    
    // Loop through the genetics of each organism.
    for (Gene g : m_organisms[i].m_genetics) {
      // Write the expressions of each gene.
      outfile << FILE_BLOCK_EXPRESSIONS;
      
      // Loop through the expressions and write the variable values.
      // m_expressions is of size NUM_INPUT_VARS.
      for (Expression e : g.m_expressions) {
        // Loop through and write m_poly coeffs
        outfile << FILE_BLOCK_POLY;
        for (double d : e.m_poly) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
        
        // Loop through and write m_log coeffs
        outfile << FILE_BLOCK_LOG;
        for (double d : e.m_log) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
        
        // Loop through and write m_sin coeffs
        outfile << FILE_BLOCK_SIN;
        for (double d : e.m_sin) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
        
        // Loop through and write m_cos coeffs
        outfile << FILE_BLOCK_COS;
        for (double d : e.m_cos) {
          outfile << "\t\t\t\t" << std::to_string(d) << '\n';
        }
        
        // Loop through and write m_exp coeffs
        outfile << FILE_BLOCK_EXP;
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

// Loads the population from a given file.  If ignoreHistory is true, do not load the
// m_totalStableTime and m_numSimulation members from the file.
void Population::load(const std::string& filename, const bool& ignoreHistory) {
  // Create the ifstream.
  std::ifstream infile(filename);
  
  // Mostly skip the fixed header information.
  std::string line;                  // Some possible values:
  std::getline(infile, line);        // <population>

  //If the first line is empty, return.  
  if (line.empty()) {
    return;
  }

  // More initial header info:  
  std::getline(infile, line);        // <runtime>
  std::getline(infile, line);        // 632635481.079
  
  // Clean and set the m_generation value if we are not ignoring history.
  line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
  line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  if (!ignoreHistory) {
    m_runtime = std::stod(line);
  }
  
  // Skip more header information.
  std::getline(infile, line);        // <generation>
  std::getline(infile, line);        // 3
  
  // Clean and set the m_generation value if we are not ignoring history.
  line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
  line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
  if (!ignoreHistory) {
    m_generation = std::stod(line) + 1.0;   //We create the object as the next generation, to correctly resume progress.
  }
  
  //Continue handling the header and skipping lines.
  std::getline(infile, line);        // <size>
  std::getline(infile, line);        // 1000
  std::getline(infile, line);        // <NUM_INPUT_VARS>
  std::getline(infile, line);        // 6
  std::getline(infile, line);        // <NUM_OUTPUT_VARS>
  std::getline(infile, line);        // 10
  
  // organismTicker is the index of the current organism we are working on.  When we
  // encounter a <organism> line, we increment this ticker.  Thus, we start at -1 as
  // we will encounter such a line immediately for the organism at 0 index.
  int organismTicker = -1;
  
  // The index of the working gene.
  // When we encounter an <m_expressions> line, we increment this.  When we switch to
  // the next organism (ie organismTicker increments) we re-set this to -1.
  // We start at -1 as we increment when we first encounter m_expressions, and the first
  // index should be 0.
  int geneTicker = -1;
  
  // See organismTicker and geneTicker comments.
  int expressionTicker = -1;
  
  // Iterate through each organism and add.
  for (; std::getline(infile, line); ) {
    // If we are encountering a new organism object, increment the index of the organism
    // in *this.m_organisms that we are working on.  We also reset geneTicker and expressionTicker.
    if (line == FILE_BLOCK_ORGANISM_STRIPPED) {
      organismTicker++;
      //std::cout << line << ": incrementing organismTicker: " << organismTicker << std::endl;
      geneTicker = -1;
      expressionTicker = -1;
      continue;
    }
    
    // The index is for people manually reading the file and referencing organisms, so we don't care about it.
    // Index is already equal to organismTicker.
    if (line == FILE_BLOCK_INDEX_STRIPPED) {
      //std::cout << line << ": skipping next line" << std::endl;
      // We aren't interested in the index, so get it and then continue so that we skip the value.
      std::getline(infile, line);
      //std::cout << line << ": did not record this" << std::endl;
      continue;
    }
    
     // Get and set m_totalStableTime.
     if (line == FILE_BLOCK_TOTAL_STABLE_TIME_STRIPPED) {
       // Get the value contained in the next line, strip the \t and \n chars,
       // and set the m_totalStableTime value.
       std::getline(infile, line);
       line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
       line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
       
       // If ignoreHistory is true, we don't load the totalStableTime value.
       if (!ignoreHistory) {
         m_organisms[organismTicker].m_totalStableTime = std::stod(line);
       }
       //std::cout << line << ": recording total stable time" << std::endl;
       continue;
     }
     
     // Get and set m_numSimulations.
     if (line == FILE_BLOCK_NUM_SIMULATIONS_STRIPPED) {
       // Get the value contained in the next line, strip the \t and \n chars,
       // and set the m_totalStableTime value.
       std::getline(infile, line);
       line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
       line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
       // If ignoreHistory is true, we don't load the m_numSimulations value.
       if (!ignoreHistory) {
         m_organisms[organismTicker].m_numSimulations = std::stod(line);
       }
       //std::cout << line << ": recording num simulations" << std::endl;
       continue;
     }
     
     // Get and set m_totalZMPDistance.
     if (line == FILE_BLOCK_TOTAL_ZMP_DISTANCE_STRIPPED) {
       // Get the value contained in the next line, strip the \t and \n chars,
       // and set the m_totalStableTime value.
       std::getline(infile, line);
       line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
       line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
       // If ignoreHistory is true, we don't load the m_totalZMPDistance value.
       if (!ignoreHistory) {
         m_organisms[organismTicker].m_totalZMPDistance = std::stod(line);
       }
       //std::cout << line << ": recording total zmp distance" << std::endl;
       continue;
     }
     
      // Get and set m_translationX.
     if (line == FILE_BLOCK_TOTAL_TRANSLATION_X_STRIPPED) {
       // Get the value contained in the next line, strip the \t and \n chars,
       // and set the m_totalStableTime value.
       std::getline(infile, line);
       line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
       line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
       // If ignoreHistory is true, we don't load the m_translationX value.
       if (!ignoreHistory) {
         m_organisms[organismTicker].m_translationX = std::stod(line);
       }
       //std::cout << line << ": recording total translation x distance" << std::endl;
       continue;
     }
     
     // If detect m_genetics block we just continue to next line.
     if (line == FILE_BLOCK_GENETICS_STRIPPED) {
       //std::cout << line << ": continuing" << std::endl;
       continue;
     }
     
     // If we detect an m_expressions block, we increment the gene index and reset the expression index.
     if (line == FILE_BLOCK_EXPRESSIONS_STRIPPED) {
       geneTicker++;
       expressionTicker = -1;
       //std::cout << line << ": incrementing geneTicker: " << geneTicker;
       //std::cout << ", resetting expressionTicker: " << expressionTicker << std::endl;
       continue;
     }
     
     // Since the expression members have multiple members due to multiple input variables,
     // we need to handle an m_poly block being in the 'line' string at the bottom of the main loop.
     // Thus, we cannot simply use continue; as otherwise the getline at the top of the mainloop will
     // overwrite this.  So, we loop through the expressions until the final line is not a
     // m_poly block.
     while (line == FILE_BLOCK_POLY_STRIPPED) {
     
       // If we detect a m_poly block, we clear the respective vector and fill it with values
       // until we encounter a non number line.
       // The m_poly block is the first subexpression block that we will encounter.
       // Thus, we increment the expressionTicker when we do so.
       if (line == FILE_BLOCK_POLY_STRIPPED) {  // Redundant in while loop but kept to keep style consistent with other expression blocks.
         expressionTicker++;
         //std::cout << line << ": incrementing expressionTicker: " << expressionTicker << std::endl;
         
         // Clear the current m_poly vector so that we can fill it with the values being read.
         m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_poly.clear();
         
         // Loop through the lines until we encounter the m_log block which follows.       
         while (std::getline(infile, line) && line != FILE_BLOCK_LOG_STRIPPED) {
           // Clean line and push it back onto the expression vector.
           line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
           line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
           m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_poly.push_back(std::stod(line));
           //std::cout << line << ": pushed back onto m_poly vector." << std::endl;
         }
         // Note that we don't continue here, as we already have the next expression block in line.
       }
       
       // If we detect a m_log block, we clear the respective vector and fill it with values
       // until we encounter a non number line.
       if (line == FILE_BLOCK_LOG_STRIPPED) {
         //std::cout << line << std::endl;
         
         // Clear the current m_log vector so that we can fill it with the values being read.
         m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_log.clear();
         
         // Loop through the lines until we encounter the m_log block which follows.       
         while (std::getline(infile, line) && line != FILE_BLOCK_SIN_STRIPPED) {
           // Clean line and push it back onto the expression vector.
           line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
           line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
           m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_log.push_back(std::stod(line));
           //std::cout << line << ": pushed back onto m_log vector." << std::endl;
         }
         // Note that we don't continue here, as we already have the next expression block in line.
       }
       
       // If we detect a m_sin block, we clear the respective vector and fill it with values
       // until we encounter a non number line.
       // The m_sin block is the first subexpression block that we will encounter.
       // Thus, we increment the expressionTicker when we do so.
       if (line == FILE_BLOCK_SIN_STRIPPED) {
         //std::cout << line << std::endl;
         
         // Clear the current m_sin vector so that we can fill it with the values being read.
         m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_sin.clear();
         
         // Loop through the lines until we encounter the m_sin block which follows.       
         while (std::getline(infile, line) && line != FILE_BLOCK_COS_STRIPPED) {
           // Clean line and push it back onto the expression vector.
           line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
           line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
           m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_sin.push_back(std::stod(line));
           //std::cout << line << ": pushed back onto m_sin vector." << std::endl;
         }
         // Note that we don't continue here, as we already have the next expression block in line.
       }    
     
       // If we detect a m_cos block, we clear the respective vector and fill it with values
       // until we encounter a non number line.
       if (line == FILE_BLOCK_COS_STRIPPED) {
         //std::cout << line << std::endl;
         
         // Clear the current m_cos vector so that we can fill it with the values being read.
         m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_cos.clear();
         
         // Loop through the lines until we encounter the m_cos block which follows.       
         while (std::getline(infile, line) && line != FILE_BLOCK_EXP_STRIPPED) {
           // Clean line and push it back onto the expression vector.
           line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
           line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
           m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_cos.push_back(std::stod(line));
           //std::cout << line << ": pushed back onto m_cos vector." << std::endl;
         }
         // Note that we don't continue here, as we already have the next expression block in line.
       }
    
       // If we detect a m_exp block, we clear the respective vector and fill it with values
       // until we encounter a non number line.
       if (line == FILE_BLOCK_EXP_STRIPPED) {
         //std::cout << line << std::endl;
         
         // Clear the current m_exp vector so that we can fill it with the values being read.
         m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_exp.clear();
         
         // Unlike for the previous expressions, there are three blocks we can
         // encounter that we neeed to handle.  We can either encounter a new <m_expression> block,
         // or we can encounter a new <organism> block.  We handle those cases after this if statement,
         // once again without calling continue, as we cannot go backward when using std::getline().
         // We can also encounter a new <m_log> block.
         while (std::getline(infile, line) && line != FILE_BLOCK_EXPRESSIONS_STRIPPED && line != FILE_BLOCK_ORGANISM_STRIPPED && line != FILE_BLOCK_POLY_STRIPPED) {
           // Clean line and push it back onto the expression vector.
           line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
           line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
           m_organisms[organismTicker].m_genetics[geneTicker].m_expressions[expressionTicker].m_exp.push_back(std::stod(line));
           //std::cout << line << ": pushed back onto m_exp vector." << std::endl;
         }
         // Note that we don't continue here, as we already have the next expression block in line.
      }
      
      // If we get here, line must be the poly block, in which case the loop will increment
      // and we will begin storing the next expression members.  Thus, we don't need to write
      // any code here.
    }
    
    // As per our previous comments in the while loop of the exp block, we need to handle the
    // <organism> and <m_expression> block cases without going to the next loop iteration, as otherwise
    // getline will be called at the start of the loop iteration and will overwrite line before we can handle
    // it appropriately.  If line is neither, it should be a poly block, and we continue the subloop after
    // incrementing expressionTicker.
     
    // If line is an expression block, increment the geneTicker and continue.
    if (line == FILE_BLOCK_EXPRESSIONS_STRIPPED) {
      geneTicker++;
      expressionTicker = -1;
      //std::cout << line << ": incrementing geneTicker: " << geneTicker << std::endl;
      //std::cout << ", resetting expressionTicker: " << expressionTicker << std::endl;
      continue;
    }
     
    // Else if the line is an organism block:
    // If we are encountering a new organism object, increment the index of the organism
    // in *this.m_organisms that we are working on.  We also reset geneTicker and expressionTicker.
    if (line == FILE_BLOCK_ORGANISM_STRIPPED) {
      organismTicker++;
      //std::cout << line << ": incrementing organismTicker: " << organismTicker << std::endl;
      geneTicker = -1;
      expressionTicker = -1;
      continue;
    }
  }
}

// Loads the population from a given file.
// TODO: Load validation such as confirming right population size etc.
void Population::load(const std::string& filename) {
  load(filename, false);
}

// Loads the population from the default filename DEFAULT_POPULATION_FILENAME.
void Population::load() {
  load(DEFAULT_POPULATION_FILENAME);
}