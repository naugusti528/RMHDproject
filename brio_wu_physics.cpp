//file for physics of brio wu shock tube

//what we're going to do here is use an HLL solver for this:
//instead of calculating every single wave speed,
//we're just going to use 1 wave going left, and another going right

#include <cmath>

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

        double get_epsilon(double Pressure, double rho_density){
            return Pressure/(rho_density*(gamma-1.0)); //calculation for energy per unit mass
        } //this varies in each grid cell, so it needs to be derived from this formula

        double get_h(double epsilon_energy, double Pressure, double rho_density){
            return 1.0 + epsilon_energy + Pressure/rho_density; //calculation for specific enthalpy
        } //like epsilon, this also varies in each grid cell

        double get_relativistic_sound_squared(double Pressure, double rho_density){
            double epsilon = get_epsilon(Pressure, rho_density);
            double h = get_h(epsilon, Pressure, rho_density);
            return gamma*Pressure/(rho_density*h);
        } //extreme conditions mean sound travels impossibly fast, close to relativistic speeds

        double get_lorentz_factor(const Vector3D& velocity){
            return 1.0/std::sqrt(1-velocity.norm_squared()); //very important --> we assume c=1, so speeds are a fraction of c
        } //in special relativity, this is what tells us how time changes at relativistic speeds

        double get_restframe_magnetic_field_squared(const Vector3D& B_field, const Vector3D& velocity){
            double lorentz_factor = get_lorentz_factor(velocity);
            double v_dot_B = velocity.dot_product(B_field);
            return B_field.norm_squared()/(lorentz_factor*lorentz_factor) + v_dot_B*v_dot_B;
        } //squared magnetic field experienced by fluid at rest

        double get_alfven_speed_squared(double Pressure, double rho_density, const Vector3D& B_field, const Vector3D& velocity){
            double epsilon = get_epsilon(Pressure, rho_density);
            double h = get_h(epsilon, Pressure, rho_density);

            double b2 = get_restframe_magnetic_field_squared(B_field, velocity);

            return b2/(rho_density*h + b2);
        } //kind of how transverse waves travel on a plucked string, alfven waves travel on B-field lines in plasma
        //this function simply finds the squared speed of those alfven waves

        double get_magnetosonic_speed_squared(double Pressure, double rho_density, const Vector3D& B_field, const Vector3D& velocity){
            double cs2 = get_relativistic_sound_squared(Pressure, rho_density);
            double va2 = get_alfven_speed_squared(Pressure, rho_density, B_field, velocity);
            return (cs2+va2-cs2*va2)/(1-cs2*va2);
        } //there are 7 waves propagating in actual ultrarelativistic environments;
          //the fast magnetosonic waves are 2 of them, and are for the most necessary details
          //the other 5 waves will be used in more advanced simulations, as you'll (hopefully) see

        

        
};

//first, let's approximate the speeds of these 2 (left&right) waves


// possible research paper directions:
// -- shock unstable regimes, unphysical oscillations, what magnetic-thermal ratio (beta) fails the solver
// -- affect of magnetic fields on dynamics: reduce B field to near 0 to reduce to hydrodynamics, to see affects of B field
// -- existence of convergence of compound wave: does it sharpen/lessen with higher resolution?