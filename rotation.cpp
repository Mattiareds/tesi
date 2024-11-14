
//ROTAZIONE DI 45 GRADI RISPETTO ALLA COSTRUZIONE INIZIALE PER AVERE I PIANI 111 IN POSIZIONE SENSATA
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

struct Atom {
    double x, y, z;
};

void rotateYZPlane(std::vector<Atom>& atoms, double angle_degrees) {
    double angle_radians = angle_degrees * M_PI / 180.0;
    double cos_angle = cos(angle_radians);
    double sin_angle = sin(angle_radians);

    for (auto& atom : atoms) {
        double y_new = cos_angle * atom.y - sin_angle * atom.z;
        double z_new = sin_angle * atom.y + cos_angle * atom.z;
        
        // Update atom's y and z coordinates
        atom.y = y_new;
        atom.z = z_new;
    }
}

int main() {
    std::ifstream infile("siti.xyz");
    std::ofstream outfile("rotated_siti.xyz");

    if (!infile || !outfile) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    std::vector<Atom> atoms;
    Atom atom;

    // Read atom positions from the file
    while (infile >> atom.x >> atom.y >> atom.z) {
        atoms.push_back(atom);
    }

    // Rotate atoms by 45 degrees in the yz-plane
    rotateYZPlane(atoms, 45.0);

    // Write rotated positions to the output file
    for (const auto& rotated_atom : atoms) {
        outfile << rotated_atom.x << " " << rotated_atom.y << " " << rotated_atom.z << std::endl;
    }

    std::cout << "Rotation complete. Check rotated_atoms.txt for results." << std::endl;

    return 0;
}
