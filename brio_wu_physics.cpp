//file for physics of brio wu shock tube

//what we're going to do here is use an HLL solver for this:
//instead of calculating every single wave speed,
//we're just going to use 1 wave going left, and another going right

struct Vector3D{
    double vector[3];
    double dot_product(const Vector3D& vector2) const {
        return vector[0]*vector2.vector[0] + vector[1]*vector2.vector[1] + vector[2]*vector2.vector[2];
    }
    double norm_squared() const {
        return dot_product(*this);
    }
};

struct GridCell{
    double rho;  double P;  Vector3D B;     Vector3D v;
    // density   pressure   magnetic field  velocity
    double epsilon;           double h;
    // energy per unit mass   specific enthalpy
};

class Brio_Wu_Physics{
    public:
        const double gamma = 2.0; //adiabatic index, indicates speed of sound & thermal compression temperature
        //for the brio wu shock tube, gamma is always 2, that's why I made it a constant

        double get_epsilon(double Pressure, double rho_density, double gamma){
            return Pressure/(rho_density*(gamma-1)); //calculation for energy per unit mass
        }
        double get_h(double epsilon_energy, double Pressure, double rho_density){
            return 1.0 + epsilon_energy + Pressure/rho_density; //calculation for specific enthalpy
        }


        //we'll start with a simple method for finding the relativistic sound speed:

};

//first, let's approximate the speeds of these 2 (left&right) waves