#include "utilities.h"

// Takes the left and right foot force sensors, and prints the values to console.
void printFootSensors (TouchSensor *fsrL, TouchSensor *fsrR) {
  const double *fsv[2] = {fsrL->getValues(), fsrR->getValues()};  // force sensor values

  double l[4], r[4];
  double newtonLeft = 0, newtonRight = 0;

  // The coefficients were calibrated against the real
  // robot so as to obtain realistic sensor values.
  l[0] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Front Left
  l[1] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Front Right
  l[3] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Rear Left
  l[2] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Rear Right

  r[0] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Front Left
  r[1] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Front Right
  r[3] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Rear Left
  r[2] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Rear Right

  int i;
  for (i = 0; i < 4; ++i) {
    l[i] = clamp(l[i], 0, 25);
    r[i] = clamp(r[i], 0, 25);
    newtonLeft += l[i];
    newtonRight += r[i];
  }

  printf("----------foot sensors----------\n");
  printf("   left       right\n");
  printf("+--------+ +--------+\n");
  printf("|%3.1f  %3.1f| |%3.1f  %3.1f|  front\n", l[0], l[1], r[0], r[1]);
  printf("|        | |        |\n");
  printf("|%3.1f  %3.1f| |%3.1f  %3.1f|  back\n", l[2], l[3], r[2], r[3]);
  printf("+--------+ +--------+\n");
  printf("total: %g Newtons, %g kilograms\n", newtonLeft + newtonRight, (newtonLeft + newtonRight) / 9.81);
}

// If value is between min:max, return value.  Else, return max if value > max or min if value < min.
double clamp(const double value, const double min, const double max) {
  if (min > max) {
    assert(0);
    return value;
  }
  return value < min ? min : value > max ? max : value;
}

// Returns a vector of 2 point elements consisting of the ZMP coordinates for each foot.
std::vector<Point> getZMPCoordinates (TouchSensor *fsrL, TouchSensor *fsrR) {
  const double *fsv[2] = {fsrL->getValues(), fsrR->getValues()};  // force sensor values
  double l[4], r[4];
  double newtonLeft = 0, newtonRight = 0;
  
  l[0] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Front Left
  l[1] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Front Right
  l[3] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Rear Left
  l[2] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Rear Right

  r[0] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Front Left
  r[1] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Front Right
  r[3] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Rear Left
  r[2] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Rear Right

  int i;
  for (i = 0; i < 4; ++i) {
    l[i] = clamp(l[i], 0, 25);
    r[i] = clamp(r[i], 0, 25);
    newtonLeft += l[i];
    newtonRight += r[i];
  }
  
  // Create return object.
  std::vector<Point> zmps;
  Point temp;
  
  // Left foot.
  temp.m_x = FOOT_WIDTH * ((l[1]+l[3])-(l[0]+l[2]))/(2*(l[0]+l[1]+l[2]+l[3]));
  temp.m_y = FOOT_LENGTH * ((l[0]+l[1])-(l[2]+l[3]))/(2*(l[0]+l[1]+l[2]+l[3]));
  zmps.push_back(temp);
  
  // Right foot.
  temp.m_x = FOOT_WIDTH * ((r[1]+r[3])-(r[0]+r[2]))/(2*(r[0]+r[1]+r[2]+r[3]));
  temp.m_y = FOOT_LENGTH * ((r[0]+r[1])-(r[2]+r[3]))/(2*(r[0]+r[1]+r[2]+r[3]));
  zmps.push_back(temp);
  
  return zmps; 
}

// Returns true with probability p.
bool trueWithProbability(const double p) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::bernoulli_distribution d(p);

  return d(mt);
}

// Returns false if the robot has fallen over, as defined per total grf on feet < FOOT_FORCE_MIN.
// Else returns true.
bool isStable(TouchSensor *fsrL, TouchSensor *fsrR){
  const double *fsv[2] = {fsrL->getValues(), fsrR->getValues()};  // force sensor values

  double l[4], r[4];
  double newtonLeft = 0, newtonRight = 0;

  // The coefficients were calibrated against the real
  // robot so as to obtain realistic sensor values.
  l[0] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Front Left
  l[1] = fsv[0][2] / 3.4 + 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Front Right
  l[3] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Rear Left
  l[2] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Rear Right

  r[0] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Front Left
  r[1] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Front Right
  r[3] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Rear Left
  r[2] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Rear Right

  int i;
  for (i = 0; i < 4; ++i) {
    l[i] = clamp(l[i], 0, 25);
    r[i] = clamp(r[i], 0, 25);
    newtonLeft += l[i];
    newtonRight += r[i];
  }
  
  return newtonLeft + newtonRight > FOOT_FORCE_MIN;
}
