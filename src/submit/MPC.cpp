#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration (from Proejcct Q&A) dlo - tunable parameters
size_t N = 10; 			// Number of calculated trajectory path.  value of 20 takes longer to compute and tend to be off track
double dt = 0.1;		// one second into the future.  Too small takes too long to calc

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//

// This is the length from front to CoG (Center of Gravity) that has a similar radius.
//  production self driving would have more tunable parameters related to CoG
const double Lf = 2.67; 	// dlo - tunable param, car of diff length - front to center of gravity of the car

// Both the reference cross track and orientation errors are 0.
// The reference velocity is set to 100 mph
double ref_cte = 0;
double ref_epsi = 0;
double ref_v = 100;  // dlo - top speed, QA crash on 200MPH (make coef higher to slow down more at high speed limit
//  dlo - cte 0 means on the line and epsi 0 means align with 0 error

// dlo - fllowing 8 lines from QA and quiz
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

// *****************************************************************
// dlo - function from quiz mpc_to_line/MPC.cpp, plus code from SDC QA
class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, 
		// `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;  // dlo - initialize cost function to 0

    // Reference State Cost =- dlo - from quiz solution
    // TODO: Define the cost related the reference state and
    // any anything you think may be beneficial.
		// dlo - from class room quiz
    for (int i = 0; i < N; i++) {
			// dlo - SDC  QA 2000 instead of 3000, higher value, more attentuation (error correction)
			// dlo - cost function - higher values (weight)  lower the error, 2000 cte gives max speed of 94MPH
      fg[0] += 2000*CppAD::pow(vars[cte_start + i]  - ref_cte, 2);  // dlo - cte cost constraints
      fg[0] += 2000*CppAD::pow(vars[epsi_start + i] - ref_epsi, 2); // dlo - steering cost constrainsts
      fg[0] += CppAD::pow(vars[v_start + i] - ref_v, 2); // dlo - throttle cost contraints
    }
		// Minimize the use of actuators.
    for (int i = 0; i < N - 1; i++) {
      fg[0] += 5*CppAD::pow(vars[delta_start + i], 2);
      fg[0] += 5*CppAD::pow(vars[a_start + i], 2);
      // dlo - try adding penalty for speed + steer - important consideration for speed
      fg[0] += 600*CppAD::pow(vars[delta_start + i] * vars[v_start+i], 2);
    }
		// calc steeting angle
		// Minimize the value gap between sequential actuations.
    for (int i = 0; i < N - 2; i++) {
      fg[0] += 200*CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
      fg[0] += 10*CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
    }

    //
    // Setup Constraints
    // NOTE: In this section you'll setup the model constraints.
    // Initial constraints
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
		// dlo - from class room quiz
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints
		// dlo - values at t and t + 1, code from the classroom quiz
    for (int t = 0; t < N -1 ; t++) {
			// time t 
      AD<double> x0 = vars[x_start + t ];
      AD<double> y0 = vars[y_start + t ];
      AD<double> psi0 = vars[psi_start + t ];
      AD<double> v0 = vars[v_start + t ];
      AD<double> cte0 = vars[cte_start + t ];
      AD<double> epsi0 = vars[epsi_start + t ];

			// Only consider the actuation at time t .
      AD<double> a = vars[a_start + t ];
      AD<double> delta = vars[delta_start + t ];

			// time t + 1
      AD<double> x1 = vars[x_start + t + 1];
      AD<double> y1 = vars[y_start + t + 1];
      AD<double> psi1 = vars[psi_start + t + 1];
      AD<double> v1 = vars[v_start + t + 1];
      AD<double> cte1 = vars[cte_start + t + 1];
      AD<double> epsi1 = vars[epsi_start + t + 1];

      if (t > 0) {   // use previous actuations (to account for latency)
        a = vars[a_start + t - 1];
        delta = vars[delta_start + t - 1];
      }
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2));

			// Here's `x` to get you started.
			// The idea here is to constraint this value to be 0.
			//
			// NOTE: The use of `AD<double>` and use of `CppAD`!
			// This is also CppAD can compute derivatives and pass
			// these to the solver.
			//
			// TODO: Setup the rest of the model constraints
			//  dlo - from class quiz
			// We add 1 to each of the starting indices due to cost being located at
			// index 0 of `fg`.
			// This bumps up the position of all the other values.
			// dlo dt - rate of change
      fg[2 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);  // dlo - CppAD lib to do the polynomial algebra
      fg[2 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);		// dlo - 
			// dlo -  steering t + 1 = psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt  // steering angle at t+1
			// dlo -  simply       psi1 = psi0 + v0 / Lf * delta * dt  // steering angle at t+1
      fg[2 + psi_start + t] = psi1 - (psi0 - v0/Lf * delta * dt); // dlo - steering angle constraint
      fg[2 + v_start + t] = v1 - (v0 + a * dt); // dlo - velocity constraint
      fg[2 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));  // dlo - Cross track error contstraint
      fg[2 + epsi_start + t] = epsi1 - ((psi0 - psides0) - v0/Lf * delta * dt);	 // dlo - steering angle error contstraint
    }
  }
	// end TODO ******************* for FG_eval()
};

// *****************************************************************
//  dlo - function is in mpc_to_line class quiz, additional code from SDC QA
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}
vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
	// dlo - parse the state vector
	/// dlo -  CPPAD simplify gradient and linear algebra  math (so does tensorflow)
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];
  // 4 * 10 + 2 * 9
  size_t n_vars = N * 6 + (N - 1) * 2; // dlo - 6 state values + cte and epsi (delta and acceleration error)
	//    N - 1 - ignore the last end point
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6; // dlo - 6 equations and 10 predictions = 60 constraints

	// dlo get codes from SDC QA video ***********
  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);  // double vector
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.

  // Set the initial variable values - dlo - additional code  from quiz mpc_to_line but not in SDC QA
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

	// dlo - next 3 set of values from SDC QA and quiz
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
	// dlo - 0.4... in SDC QA  but not in quiz, but ok due to QA typo
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;  // dlo - max breaking 
    vars_upperbound[i] = 1.0;  // dlo - max acceleration, about 5 Meter/sec in simulator in max accleration
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
	// initial state values when car starts
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";  // dlo - ok for now, increase if calc takes more time (N ...) for the simulator
	// dlo - simulator is sending about 20 frames/sec, don't wait too long to compute the solution, related to N value (10)

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  //std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;

  result.push_back(solution.x[delta_start]);  // dlo - steering
  result.push_back(solution.x[a_start]);  // dlo - throttle

  for (int i = 0; i < N-1; i++) {
    result.push_back(solution.x[x_start + i + 1]);
    result.push_back(solution.x[y_start + i + 1]);
  }

  return result;
}
