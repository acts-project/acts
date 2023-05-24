#this code is trying to implement the alignment of superstructures

#include <iostream>
#include <vector>
#include <cmath>

// our class Surface that represents surface element and holds 6 parameters (T_x, T_y, T_z, R_x, R_y, R_z)
// Notation used for rotation parameters: R_z - descibes the rotation around z-axis, R_y - describes the rotation around y-axis, R_x - describes the rotation around x-axis
class Surface {
public:
  double x_position;  // x position
  double y_position;  // y position
  double z_position;  // z position
  double Phi;  // R_z
  double Theta;  // R_y
  double Psi;  // R_x

  // Constructor
  Surface(double posX, double posY, double posZ, double rotPhi, double rotTheta, double rotPsi)
      : x_position(posX), y_position(posY), z_position(posZ), Phi(rotPhi), Theta(rotTheta), Psi(rotPsi) {}

   // Defining the function rotate that will rotate surface element - around x,y,z axis
  void rotate(double rx, double ry, double rz) {
    Phi += rx;
    Theta += ry;
    Psi += rz;
  }


  // Defining the function translate that will translate the surface element around the x,y,z axis 
  void translate(double tx, double ty, double tz) {
    x_position += tx;
    y_position += ty;
    z_position += tz;
  }


  // Check out - information on the position and orientation of the surface: 
  void print() const {
    std::cout << "Position of a surface in x,y,z direction: (" << x_position << ", " << y_position << ", " << z_position << ")\n";
    std::cout << "Orientation of a surface in Phi, Theta, Psi direction: (" << Phi << ", " << Theta << ", " << Psi << ")\n";
  }
};

// Let's define the function that will align superstructures
void alignSuperstructure(std::vector<Surface>& superstructure) {
  const int maxIterations = 100; //iterativly aligns 
  const double Threshold = 0.002;

// For each surface in the superstructure, preform the alignment (position and orientation)
  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    for (auto& surface : superstructure) {
      
      // this is try: without predefined data set, I chose to test the code with this set of arbitary values (getRandomValueInRange is called for it - within a range)
      double dx_pos = getRandomValueInRange(-0.2, 0.2);
      double dy_pos = getRandomValueInRange(-0.2, 0.2);
      double dz_pos = getRandomValueInRange(-0.2, 0.2);
      double dPhi_pos = getRandomValueInRange(-0.02, 0.02);
      double dTheta_pos = getRandomValueInRange(-0.02, 0.02);
      double dPsi_pos = getRandomValueInRange(-0.02, 0.02);

      surface.translate(dx_pos, dy_pos, dz_pos);
      surface.rotate(dPhi_pos, dTheta_pos, dPsi_pos);
    }

    // desired level of accuracy shouldn't be >  specific value (= convergenceThreshold) then label alignment as not acceptable 
    double overallMisalignment = calculateOverallMisalignment(superstructure);
    if (overallMisalignment < Threshold) {
      std::cout << "Success! The superstructure is aligned! \n";
      break;
    }
  }
}

// Function to calculate the overall misalignment of the superstructure
double calculateOverallMisalignment(const std::vector<Surface>& superstructure) {
  // Example: Calculate the overall misalignment based on the positions and orientations of the surfaces
  // You can define your own criteria or use mathematical calculations to determine the overall misalignment

  double overallMisalignment = 0.0;
  for (const auto& surface : superstructure) {
    overallMisalignment += calculateSurfaceMisalignment(surface);
  }

  return overallMisalignment;
}

double getRandomValueInRange(double minValue, double maxValue) {
  double random = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  return minValue + random * (maxValue - minValue);
}

int main() {
  

  // Superstructures with initial values
  std::vector<Surface> superstructure;
  superstructure.emplace_back(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);  // Surface 1
  superstructure.emplace_back(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);  // Surface 2
  superstructure.emplace_back(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);  // Surface 3
  superstructure.emplace_back(1.0, 1.0, 0.0, 0.0, 0.0, 0.0);  // Surface 4

  


  // at the end --> call the alignSuperstructure that will align superstructure
  alignSuperstructure(superstructure);

