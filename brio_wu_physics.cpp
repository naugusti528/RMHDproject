//file for physics of brio wu shock tube
//what we're going to do here is use an HLLE solver for this

#include <cmath>
#include <iostream>
#include <fstream>
#include <limits>

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
    double D;    // relativistic mass-density from stationary pov  -> rest mass conservation
    double tau;  // total energy density without rest mass         -> energy density conservation
    Vector3D S;  // relativistic momentum density (spatial vector) -> momentum conservation
    Vector3D B;  // magnetic field (spatial vector)                -> magnetic flux conservation
};

struct SweepResult{
    double P_left, P_right, By_left, By_right;
    double vx_left, vx_right;  // velocity values of left&right
    double beta_left;          // plasma beta of the left state, for reference
    int floor_count;
    bool completed;            // false if the run produced NaN/Inf and had to bail early
    int max_iterations_seen;   // worst-case NR iteration count across all cells, all timesteps
    double min_pressure_seen;  // lowest pressure encountered anywhere, ever, in the run
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

        void get_fast_wavespeeds(double Pressure, double rho_density, const Vector3D& B_field, const Vector3D& velocity, double& lambda_plus, double& lambda_minus){
            double cs2 = get_relativistic_sound_squared(Pressure, rho_density);
            double va2 = get_alfven_speed_squared(Pressure, rho_density, B_field, velocity);
            double a2  = get_magnetosonic_speed_squared(Pressure, rho_density, B_field, velocity);

            double v2 = velocity.norm_squared();
            double vx = velocity.vector[0];

            double denom     = 1.0 - v2*a2;
            double bracket   = (1.0 - v2*a2) - (1.0 - a2)*vx*vx;
            double sqrt_term = std::sqrt(a2*(1.0 - v2)*bracket);

            lambda_plus  = ((1.0 - a2)*vx + sqrt_term) / denom;  // right-going fast wave in stationary pov
            lambda_minus = ((1.0 - a2)*vx - sqrt_term) / denom;  // left-".... ... ..."
        }

        Conserved prim2cons(const GridCell& prim){ // translating physics quantities into those that the evolution equations operate on
            double v2 = prim.v.norm_squared();
            double Lorentz = 1.0/std::sqrt(1.0-v2);
            double eps = get_epsilon(prim.P, prim.rho);
            double h = get_h(eps, prim.P, prim.rho);
            double W = prim.rho*h*Lorentz*Lorentz;
            double B2 = prim.B.norm_squared();
            double vdotB = prim.v.dot_product(prim.B);
            double b2 = get_restframe_magnetic_field_squared(prim.B, prim.v);

            Conserved U;
            U.D = prim.rho*Lorentz;
            for(int i=0;i<3;i++)
                U.S.vector[i] = (W+B2)*prim.v.vector[i] - vdotB*prim.B.vector[i];
            U.tau = W + B2 - prim.P - 0.5*b2 - U.D;
            U.B = prim.B;
            return U;
        }
        
        bool cons2prim(const Conserved& U, GridCell& out, int& iterations, double tol = 1e-10, int max_iter = 50){ // can the evolution equation quantities reflect measurable physics quantities?
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
                if(W_new <= 0) W_new = 0.5*W; // physical guard: W must be positive
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

        // this tells us the fluxes of everything, from mass and energy to momentum and magnetic field
        Conserved physical_flux(const GridCell& state){
            double v2 = state.v.norm_squared();
            double Lorentz = 1.0/std::sqrt(1.0-v2);
            double vdotB = state.v.dot_product(state.B);
            double b0 = Lorentz*vdotB;

            Vector3D b; // comoving field, lab-frame spatial components
            for(int i=0;i<3;i++)
                b.vector[i] = state.B.vector[i]/Lorentz + b0*state.v.vector[i];

            double b2 = get_restframe_magnetic_field_squared(state.B, state.v);
            double p_tot = state.P + 0.5*b2;

            Conserved U = prim2cons(state); // reuse — gives D, S, tau, B
            double vx = state.v.vector[0];
            double Bx = state.B.vector[0];

            Conserved F;
            F.D = U.D * vx;
            for(int i=0;i<3;i++)
                F.S.vector[i] = U.S.vector[i]*vx + (i==0 ? p_tot : 0.0) - b.vector[i]*Bx/Lorentz;
            F.tau = U.S.vector[0] - U.D*vx;
            F.B.vector[0] = 0.0;
            F.B.vector[1] = vx*state.B.vector[1] - state.v.vector[1]*Bx;
            F.B.vector[2] = vx*state.B.vector[2] - state.v.vector[2]*Bx;
            return F;
        }

        // function for computing numerical flux at grid cell level
        Conserved hlle_flux(const GridCell& left, const GridCell& right){
            double lp_L, lm_L, lp_R, lm_R;
            get_fast_wavespeeds(left.P, left.rho, left.B, left.v, lp_L, lm_L);
            get_fast_wavespeeds(right.P, right.rho, right.B, right.v, lp_R, lm_R);

            double S_L = std::min({lm_L, lm_R, 0.0});
            double S_R = std::max({lp_L, lp_R, 0.0});

            Conserved U_L = prim2cons(left);
            Conserved U_R = prim2cons(right);
            Conserved F_L = physical_flux(left);   // <- need this next, see below
            Conserved F_R = physical_flux(right);

            if(S_L >= 0.0) return F_L;
            if(S_R <= 0.0) return F_R;

            Conserved F_hlle;
            F_hlle.D = (S_R*F_L.D - S_L*F_R.D + S_L*S_R*(U_R.D - U_L.D)) / (S_R - S_L);
            F_hlle.tau = (S_R*F_L.tau - S_L*F_R.tau + S_L*S_R*(U_R.tau - U_L.tau)) / (S_R - S_L);
            for(int i=0;i<3;i++)
                F_hlle.S.vector[i] = (S_R*F_L.S.vector[i] - S_L*F_R.S.vector[i] + S_L*S_R*(U_R.S.vector[i]-U_L.S.vector[i])) / (S_R - S_L);
            for(int i=0;i<3;i++)
                F_hlle.B.vector[i] = (S_R*F_L.B.vector[i] - S_L*F_R.B.vector[i] + S_L*S_R*(U_R.B.vector[i]-U_L.B.vector[i])) / (S_R - S_L);
            return F_hlle;
        }
};

void run_unit_tests(Brio_Wu_Physics& physics, GridCell& left_state){
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
    
    // test 3: round trip (prim -> consv -> prim) -- result should show that prim var's are preserved
    left_state.rho = 1.0;
    left_state.P = 0.05;
    left_state.v = Vector3D{{0,0,0}};
    left_state.B = Vector3D{{0.75,1.0,0.0}};
    Conserved U = physics.prim2cons(left_state);
    GridCell recovered;
    int iters;
    bool ok = physics.cons2prim(U, recovered, iters);
    std::cout << "\nRound-trip test:\n";
    std::cout << "success = " << ok << ", iterations = " << iters << "\n";
    std::cout << "rho: " << recovered.rho << " (expect " << left_state.rho << ")\n";
    std::cout << "P: " << recovered.P << " (expect " << left_state.P << ")\n";
    std::cout << "vx: " << recovered.v.vector[0] << " (expect " << left_state.v.vector[0] << ")\n";

    // test 4a: flux test with physical flux values -- result should show accurate flux values
    Conserved F = physics.physical_flux(left_state);
    std::cout << "\nphysical_flux test:\n";
    std::cout << "F_D = " << F.D << " (expect 0)\n";
    std::cout << "F_tau = " << F.tau << " (expect 0)\n";
    std::cout << "F_Sx = " << F.S.vector[0] << " (expect 0.26875)\n";
    std::cout << "F_Sy = " << F.S.vector[1] << " (expect -0.75)\n";
    // test 4b: flux test with hlle flux -- left-left pair should show identical results
    Conserved F_same = physics.hlle_flux(left_state, left_state);
    std::cout << "\nhlle_flux (identical states) test:\n";
    std::cout << "F_D = " << F_same.D << " (expect 0)\n";
    std::cout << "F_Sx = " << F_same.S.vector[0] << " (expect 0.26875)\n";
    std::cout << "F_Sy = " << F_same.S.vector[1] << " (expect -0.75)\n";
    // test 4c: flux test with hlle flux -- actual left-right test
    GridCell right_state;
    right_state.rho = 0.125;
    right_state.P = 0.005;
    right_state.v = Vector3D{{0,0,0}};
    right_state.B = Vector3D{{0.75,-1.0,0.0}};
    Conserved F_riemann = physics.hlle_flux(left_state, right_state);
    std::cout << "\nhlle_flux (Brio-Wu L/R) test:\n";
    std::cout << "F_D = " << F_riemann.D << "\n";
    std::cout << "F_Sx = " << F_riemann.S.vector[0] << "\n";
}

SweepResult run_simulation(double P_left, double P_right, double By_left, double By_right, double vx_left, double vx_right, int N = 400, double CFL = 0.4, double t_end = 0.2){
    Brio_Wu_Physics physics;
    double x_min = -0.5, x_max = 0.5;
    double dx = (x_max - x_min) / N;
    double t = 0.0;

    GridCell left_state, right_state;
    left_state.rho = 1.0;   left_state.P = P_left;
    left_state.v = Vector3D{{vx_left,0,0}};  left_state.B = Vector3D{{0.75, By_left, 0.0}};
    right_state.rho = 0.125; right_state.P = P_right;
    right_state.v = Vector3D{{vx_right,0,0}}; right_state.B = Vector3D{{0.75, By_right, 0.0}};

    std::vector<GridCell> W(N);
    std::vector<Conserved> U(N);
    for(int i=0;i<N;i++){
        double x = x_min + (i+0.5)*dx;
        W[i] = (x < 0.0) ? left_state : right_state;
        U[i] = physics.prim2cons(W[i]);
    }

    int floor_count = 0;
    bool completed = true;
    int max_iterations_seen = 0;
    double min_pressure_seen = std::numeric_limits<double>::max();

    // safety check to ensure the program stops when it stalls
    int step_count = 0; int max_steps = 100000;

    while(t < t_end){
        double max_speed = 0.0;
        for(int i=0;i<N;i++){
            double lp, lm;
            physics.get_fast_wavespeeds(W[i].P, W[i].rho, W[i].B, W[i].v, lp, lm);
            max_speed = std::max({max_speed, std::abs(lp), std::abs(lm)});
        }
        if(!(max_speed > 0.0) || std::isnan(max_speed)){ completed = false; break; }
        double dt = CFL * dx / max_speed;
        if(t + dt > t_end) dt = t_end - t;

        std::vector<Conserved> F(N+1);
        for(int i=1;i<N;i++) F[i] = physics.hlle_flux(W[i-1], W[i]);
        F[0] = physics.physical_flux(W[0]);
        F[N] = physics.physical_flux(W[N-1]);

        std::vector<Conserved> U_new(N);
        for(int i=0;i<N;i++){
            U_new[i].D   = U[i].D   - (dt/dx)*(F[i+1].D   - F[i].D);
            U_new[i].tau = U[i].tau - (dt/dx)*(F[i+1].tau - F[i].tau);
            for(int k=0;k<3;k++){
                U_new[i].S.vector[k] = U[i].S.vector[k] - (dt/dx)*(F[i+1].S.vector[k] - F[i].S.vector[k]);
                U_new[i].B.vector[k] = U[i].B.vector[k] - (dt/dx)*(F[i+1].B.vector[k] - F[i].B.vector[k]);
            }
        }

        for(int i=0;i<N;i++){
            int iters;
            bool ok = physics.cons2prim(U_new[i], W[i], iters);
            if(!ok) floor_count++;
            max_iterations_seen = std::max(max_iterations_seen, iters);
            min_pressure_seen = std::min(min_pressure_seen, W[i].P);
            U[i] = U_new[i];
        }
        t += dt;
        step_count++;
        if(step_count > max_steps) {completed=false; break;}
    }

    double beta_left = (2.0*P_left) / (left_state.B.norm_squared());

    return SweepResult{P_left, P_right, By_left, By_right, vx_left, vx_right, beta_left, floor_count, completed, max_iterations_seen, min_pressure_seen};
}

int main(){
    
    Brio_Wu_Physics physics;
    GridCell left_state;
    left_state.rho = 1.0;
    left_state.P = 0.05;
    left_state.v = Vector3D{{0,0,0}};
    left_state.B = Vector3D{{0.75,1.0,0.0}};

    run_unit_tests(physics, left_state);

    std::vector<double> pressure_scales = {1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001};
    // each scale multiplies the baseline P_left=0.05, P_right=0.005 -- lower scale = lower beta
    std::vector<double> velocities = {-0.9, -0.5, 0.0, 0.5, 0.9, 0.99};

    std::ofstream outfile("beta_velocity_asymmetric_sweep_results.csv");
    outfile << "P_left,P_right,beta_left,vx_left,vx_right,floor_count,completed,max_iterations_seen,min_pressure_seen\n";

    for(double scale : pressure_scales){
        for(double vL : velocities){
            for(double vR : velocities){
                SweepResult r = run_simulation(0.05*scale, 0.005*scale, 1.0, -1.0, vL, vR);
                outfile << r.P_left <<","<< r.P_right <<","<< r.beta_left <<","<< r.vx_left <<","<< r.vx_right <<","<< r.floor_count <<","<< r.completed <<","<< r.max_iterations_seen <<","<< r.min_pressure_seen<<"\n";
            }
        }
    }
    outfile.close();

    return 0;
}


// possible research paper directions:
// -- shock unstable regimes, unphysical oscillations, what magnetic-thermal ratio (beta) fails the solver
// -- affect of magnetic fields on dynamics: reduce B field to near 0 to reduce to hydrodynamics, to see affects of B field
// -- existence of convergence of compound wave: does it sharpen/lessen with higher resolution?