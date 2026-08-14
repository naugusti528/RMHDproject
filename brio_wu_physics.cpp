//file for physics of brio wu shock tube

//what we're going to do here is use an HLLE solver for this:
//instead of calculating every single wave speed,
//we're just going to use 1 wave going left, and another going right

#include <cmath>
#include <iostream>


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

struct Conserved{
    double D, tau;
    Vector3D S;
    Vector3D B;
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

        void get_fast_wavespeeds(double Pressure, double rho_density,
                          const Vector3D& B_field, const Vector3D& velocity,
                          double& lambda_plus, double& lambda_minus){
            double cs2 = get_relativistic_sound_squared(Pressure, rho_density);
            double va2 = get_alfven_speed_squared(Pressure, rho_density, B_field, velocity);
            double a2  = get_magnetosonic_speed_squared(Pressure, rho_density, B_field, velocity); // reuse your existing function directly

            double v2 = velocity.norm_squared();
            double vx = velocity.vector[0];

            double denom     = 1.0 - v2*a2;
            double bracket   = (1.0 - v2*a2) - (1.0 - a2)*vx*vx;
            double sqrt_term = std::sqrt(a2*(1.0 - v2)*bracket);

            lambda_plus  = ((1.0 - a2)*vx + sqrt_term) / denom;
            lambda_minus = ((1.0 - a2)*vx - sqrt_term) / denom;
        }       //lambda_plus = right-going fast wave in lab frame, lambda_minus = left-going
        
        bool cons2prim(const Conserved& U, GridCell& out, int& iterations, double tol = 1e-10, int max_iter = 50){
            double S2 = U.S.norm_squared();
            double B2 = U.B.norm_squared();
            double SdotB = U.S.dot_product(U.B);

            double W = U.tau + U.D; // crude initial guess

            auto residual = [&](double Wt, double& P_out, double& v2_out)->double{
                double v2 = (S2*Wt*Wt + SdotB*SdotB*(2.0*Wt + B2)) / (Wt*Wt*(Wt+B2)*(Wt+B2));
                if(v2 >= 1.0) v2 = 1.0 - 1e-12;
                double Lorentz = 1.0/std::sqrt(1.0-v2);
                double rho = U.D/Lorentz;
                double h = Wt/(rho*Lorentz*Lorentz);
                double P = rho*(h-1.0)*(gamma-1.0)/gamma;
                double b2 = B2*(1.0-v2) + SdotB*SdotB/(Wt*Wt);
                P_out = P; v2_out = v2;
                return Wt + B2 - P - 0.5*b2 - U.D - U.tau;
            };

            double P_dummy, v2_dummy, delta = 1e-6;
            int iter = 0;
            for(; iter < max_iter; iter++){
                double f = residual(W, P_dummy, v2_dummy);
                double f_hi = residual(W+delta, P_dummy, v2_dummy);
                double f_lo = residual(W-delta, P_dummy, v2_dummy);
                double dfdW = (f_hi - f_lo)/(2.0*delta);
                if(std::abs(dfdW) < 1e-14) break;
                double W_new = W - f/dfdW;
                if(W_new <= B2) W_new = B2 + 1e-10; // physical guard: W must exceed B^2
                if(std::abs(W_new - W) < tol){ W = W_new; iter++; break; }
                W = W_new;
            }
            iterations = iter;

            double P_final, v2_final;
            double f_final = residual(W, P_final, v2_final);
            if(std::abs(f_final) > 1e-6 || P_final <= 0.0 || v2_final >= 1.0) return false;

            double Lorentz = 1.0/std::sqrt(1.0-v2_final);
            out.rho = U.D/Lorentz;
            out.P = P_final;
            out.B = U.B;
            for(int i=0;i<3;i++)
                out.v.vector[i] = (U.S.vector[i] + (SdotB/W)*U.B.vector[i]) / (W+B2);
            out.epsilon = get_epsilon(out.P, out.rho);
            out.h = get_h(out.epsilon, out.P, out.rho);
            return true;
        }

        
};

int main(){
    Brio_Wu_Physics physics;
    Vector3D zero_v{{0,0,0}};
    Vector3D zero_B{{0,0,0}};
    double lp, lm;

    // test 1: no velocity and no magnetic field -- result should be just sound speed
    physics.get_fast_wavespeeds(0.05, 1.0, zero_B, zero_v, lp, lm);
    double cs_expected = std::sqrt(physics.get_relativistic_sound_squared(0.05, 1.0));
    std::cout << "Test 1 (no velocity, no field):\n";
    std::cout << "lambda_plus = " << lp << " (expect " << cs_expected << ")\n";
    std::cout << "lambda_minus = " << lm << " (expect " << -cs_expected << ")\n";

    // test 2: no velocity, some field -- result should be symmetric and magnitude should be greater than cs
    Vector3D left_B{{0.75, 1.0, 0.0}};
    physics.get_fast_wavespeeds(0.05, 1.0, left_B, zero_v, lp, lm);
    std::cout << "Test 2 (with field):\n";
    std::cout << "lambda_plus = " << lp << "\n";
    std::cout << "lambda_minus = " << lm << "\n";
    return 0;
}


// possible research paper directions:
// -- shock unstable regimes, unphysical oscillations, what magnetic-thermal ratio (beta) fails the solver
// -- affect of magnetic fields on dynamics: reduce B field to near 0 to reduce to hydrodynamics, to see affects of B field
// -- existence of convergence of compound wave: does it sharpen/lessen with higher resolution?