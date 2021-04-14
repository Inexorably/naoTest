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
  l[2] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] - 1.15 * fsv[0][1];  // Left Foot Rear Right
  l[3] = fsv[0][2] / 3.4 - 1.5 * fsv[0][0] + 1.15 * fsv[0][1];  // Left Foot Rear Left

  r[0] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Front Left
  r[1] = fsv[1][2] / 3.4 + 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Front Right
  r[2] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] - 1.15 * fsv[1][1];  // Right Foot Rear Right
  r[3] = fsv[1][2] / 3.4 - 1.5 * fsv[1][0] + 1.15 * fsv[1][1];  // Right Foot Rear Left

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
  temp.m_x = FOOT_WIDTH * ((l[1]+l[2])-(l[0]+l[3]))/(2*(l[0]+l[1]+l[3]+l[2]));
  temp.m_y = FOOT_LENGTH * ((l[0]+l[1])-(l[3]+l[2]))/(2*(l[0]+l[1]+l[3]+l[2]));
  
  // Check for NaN errors.
  if (std::isnan(temp.m_x))
    temp.m_x = 0;
  if (std::isnan(temp.m_y))
    temp.m_y = 0;
    
  zmps.push_back(temp);
  
  // Right foot.
  temp.m_x = FOOT_WIDTH * ((r[1]+r[3])-(r[0]+r[2]))/(2*(r[0]+r[1]+r[2]+r[3]));
  temp.m_y = FOOT_LENGTH * ((r[0]+r[1])-(r[2]+r[3]))/(2*(r[0]+r[1]+r[2]+r[3]));
  
  // Check for NaN errors.
  if (std::isnan(temp.m_x))
    temp.m_x = 0;
  if (std::isnan(temp.m_y))
    temp.m_y = 0;
    
  zmps.push_back(temp);
  
  return zmps; 
}

// Take the pointer to the rotational field [rot_field->getSFRotation()], in axis-angle form.
// Axis-angle form is format of unit x y z vector + rotation angle.
// Return the euclidean rotational angles (in xy, xz, yz).
// https://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToEuler/index.htm
std::vector<double> getRotationalAngles(const double* rot) {
  // Pull the values from rot.  Note that theta is in radians.
  double x = rot[0];
  double y = rot[1];
  double z = rot[2];
  double theta = rot[3];
  
  // heading rotation around y (in xz), attitude around z (in xy), bank around x (in yz).
  double heading = atan2(y * sin(theta)- x * z * (1 - cos(theta)) , 1 - (pow(y,2) + pow(z,2)) * (1 - cos(theta)));
  double attitude = asin(x * y * (1 - cos(theta)) + z * sin(theta));
  double bank = atan2(x * sin(theta)-y * z * (1 - cos(theta)) , 1 - (pow(x,2) + pow(z,2)) * (1 - cos(theta)));
  
  std::vector<double> angles = {heading, attitude, bank};
  
  return angles;
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
