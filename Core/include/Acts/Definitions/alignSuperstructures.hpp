#ifndef SUPERSTRUCTURE_ALIGNMENT_H
#define SUPERSTRUCTURE_ALIGNMENT_H

#include <iostream>
#include <vector>
#include <cmath>

class Surface {
public:
    double x_position;
    double y_position;
    double z_position;
    double Phi;
    double Theta;
    double Psi;

    Surface(double posX, double posY, double posZ, double rotPhi, double rotTheta, double rotPsi)
        : x_position(posX), y_position(posY), z_position(posZ), Phi(rotPhi), Theta(rotTheta), Psi(rotPsi) {}

    void rotate(double rx, double ry, double rz);
    void translate(double tx, double ty, double tz);
    void print() const;
};

void alignSuperstructure(std::vector<Surface>& superstructure);
double calculateOverallMisalignment(const std::vector<Surface>& superstructure);
double getRandomValueInRange(double minValue, double maxValue);

#endif // SUPERSTRUCTURE_ALIGNMENT_H
