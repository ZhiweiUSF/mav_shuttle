#include "gurobi_c++.h"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <cmath>
#include <algorithm>

using namespace std;

/*
 * This script is calls Gurobi to solve the discrete model in 
 * Chen, Z., Li, X., & Zhou, X. (2019). Operational design for shuttle systems with modular vehicles under oversaturated traffic: Discrete modeling method. Transportation Research Part B: Methodological, 122, 1-19.
 * Author: Zhiwei Chen. 
 * email: zhiweic@usf.edu
 */

struct input_parameters {
	/* A structure encapsulating all input parameters in the model */
	double exchange_rate; // coefficient to convert RMB to US dollars
	double fixed_operational_cost; // c_f
	double variable_operational_cost;  // c_v
	double power_index_in_operational_cost;  // alpha
	int capacity_per_unit; // c
	double length_per_time_interval; // delta
	double unit_time_waiting_cost; // w
	int min_headway; // h_underbar
	int min_headway_interval; // j_underbar
	vector<int> vehicle_formation;  // set I
	vector<double> vehicle_capacity; // c*i, for all i
	vector<double> vehicle_operational_cost; // f_i, 
	vector<double> passenger_demand; // a_j
	int num_time_indexes; // J
	int num_vehicle_formation; // I
	vector<double> cumulative_demand; // A_j
	vector<double> oversaturated_queue; // deltaq_j 
};

template <typename T>
void print_one_dimensional_vector(vector<T> one_d_vector){
	/* A function to print one dimensional vector on the screen, used for testing the code */
	for(int i = 0; (unsigned) i < one_d_vector.size(); ++i){
		cout << one_d_vector[i] << "  ";
	}
	cout << endl;
}

bool read_file(const string file_path, vector<double>& vector_to_read) {
	/* A function read the passenger demand file.
	 * The file should be with an extension .csv. The data should be just a row. 
	 */
	ifstream in_file(file_path);

	if(!in_file) {
		cout << "No file has been open!" << endl;
		return false;
	} else {
		cout << "Start reading file from  " << file_path << endl;
		string line;
		while (getline(in_file, line, ',')) {
			vector_to_read.push_back(atof(line.c_str()));
		}

		in_file.close();
		return true;
	}
}


template <typename T>
void print_results (string file_path, vector<T> vector_to_print) {
	/* A function to write a one dimensional vector into a csv file */
	ofstream outstream;
	outstream.open(file_path);
	for (int j = 0; j < vector_to_print.size(); j++) {
		outstream << vector_to_print[j] << endl;
	}
	outstream << endl;
	outstream.close();
}

void compute_cumulative_passenger_demand(vector<double> passenger_demand, vector<double>& cumulative_demand, int num_time_indexes) {
	/* A function to compute the cumulative passenger demand */
	cumulative_demand[0] = passenger_demand[0];
	for (int j = 1; j < num_time_indexes; ++j) {
		cumulative_demand[j] = cumulative_demand[j-1] + passenger_demand[j];
	}
}

void create_input_parameters (input_parameters& input_params) {
	/* A function to generate all parameters. These are from the Beijing case study in the paper.
	 * You will need to do it here if you want to change any parameters.
	 */

	input_params.exchange_rate = 6.6; // unitless
	double unit_electricity_cost = 0.8; // RMB/kWh
	double length_of_line = 18.964; // km

	input_params.fixed_operational_cost = 2.049 * unit_electricity_cost * length_of_line; // RMB
	input_params.variable_operational_cost = 0.155 * unit_electricity_cost * length_of_line; // RMB
	input_params.power_index_in_operational_cost = 0.5; // unitless
	input_params.capacity_per_unit = 226; // passengers/unit
	input_params.unit_time_waiting_cost = 40.0 / 60.0; // RMB/minute
	input_params.min_headway = 3; // minutes
	
	input_params.vehicle_formation = vector<int>{1, 2, 3, 4, 5, 6};  // 6 types of vehicles, 0 means not dispatching anything
	input_params.vehicle_capacity.resize(6); // passengers
	input_params.vehicle_operational_cost.resize(6); // RMB
	for (int i = 0; i < input_params.vehicle_formation.size(); i++) {
		input_params.vehicle_capacity[i] = input_params.vehicle_formation[i] * input_params.capacity_per_unit;
		input_params.vehicle_operational_cost[i] = input_params.fixed_operational_cost + input_params.variable_operational_cost * pow(input_params.vehicle_capacity[i],
			 input_params.power_index_in_operational_cost);
	}
	//print_one_dimensional_vector(input_params.vehicle_operational_cost);

	input_params.length_per_time_interval = 3.0; // minute
	input_params.min_headway_interval = ceil(input_params.min_headway/input_params.length_per_time_interval); // time intervals
	read_file("../data/1.csv", input_params.passenger_demand);
	//print_one_dimensional_vector(input_params.passenger_demand);

	input_params.num_time_indexes = input_params.passenger_demand.size();
	input_params.num_vehicle_formation = input_params.vehicle_formation.size();
	input_params.cumulative_demand.resize(input_params.num_time_indexes);
	compute_cumulative_passenger_demand(input_params.passenger_demand, input_params.cumulative_demand, input_params.num_time_indexes);
	//print_one_dimensional_vector(input_params.cumulative_demand);
}

void build_obj(GRBModel& model, GRBVar* vars_x, GRBVar* vars_q, vector<double> x_coefficient, vector<double> q_coefficient) {
	/* A function to build the objective function */
	GRBLinExpr obj = 0;
	for (int i = 0; i < x_coefficient.size(); ++i) {
		obj += x_coefficient[i] * vars_x[i];
	}
	for (int i = 0; i < q_coefficient.size(); ++i) {
		obj += q_coefficient[i] * vars_q[i];
	}
	model.setObjective(obj);
	//return 0;
}

void build_constrs(GRBModel& model, GRBVar* vars_x, GRBVar* vars_q, input_parameters input_params) {
	/* A function to build constraints */
	// min headway requirement, Eq.(2)
	for (int j = input_params.min_headway_interval - 1; j < input_params.num_time_indexes; ++j) { //input_params.min_headway_interval - 1
		GRBLinExpr lhs = 0;
		for (int j_prime = j - input_params.min_headway_interval + 1; j_prime <= j; ++j_prime) { //
			for (int i = 0; i < input_params.num_vehicle_formation; ++i) {
				lhs += 1 * vars_x[ i+j_prime*input_params.num_vehicle_formation ];
			}
		}
		
		model.addConstr(lhs, '<', 1);
	}
	// queueing conservation 
	model.addConstr(vars_q[0], '=', 0); // Eq.(5)
	model.addConstr(vars_q[input_params.num_time_indexes-1], '=', 0); //Eq.(6)
	for (int j = 1; j < input_params.num_time_indexes; ++j) { //Eq.(3), q_j >=0 by default so Eq.(4) needs not to be written
		GRBLinExpr lhs = vars_q[j] - vars_q[j-1];
		for (int i = 0; i < input_params.num_vehicle_formation; ++i) {
			lhs += input_params.vehicle_capacity[i] * vars_x[ i+j*input_params.num_vehicle_formation ];
		}
		model.addConstr(lhs, '>', input_params.passenger_demand[j]);
	}
}


void get_results(GRBModel& model, GRBVar* vars_x, GRBVar* vars_q, input_parameters input_params, double& opt_obj, string file_path) {
	/* a function to extract the optimal solution */

	// find optimal objective value
	opt_obj = model.get(GRB_DoubleAttr_ObjVal);
	// add constant term into objective function
	opt_obj = (opt_obj + input_params.cumulative_demand[input_params.num_time_indexes - 1]*input_params.unit_time_waiting_cost*input_params.length_per_time_interval/2)/input_params.exchange_rate;

	// find optimal dispatch time and vehicle formation
	vector<int> optimal_dispatch_time;
	vector<int> optimal_vehicle_formation;
	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {

		for (int j = 0; j < input_params.num_time_indexes; ++j) {
			for (int i = 0; i < input_params.num_vehicle_formation; ++i) {
				if ( vars_x[ i+j*input_params.num_vehicle_formation ].get(GRB_DoubleAttr_X) > 0 ) {
					optimal_dispatch_time.push_back(j*input_params.length_per_time_interval);
					optimal_vehicle_formation.push_back(i+1);
				}
			}
		}
	}
	
	// write resutls to csv
	print_results(file_path+"opt_time_gu.csv", optimal_dispatch_time);
	print_results(file_path+"opt_formation_gu.csv", optimal_vehicle_formation);
}


int main(int argc, char* argv[]) {
	/*** input ***/
	input_parameters input_params;
	create_input_parameters(input_params);

	/*** create variables  ***/
	vector<double> x_coefficient;
	vector<double> q_coefficient;
	int num_x;
	for (int j = 0; j < input_params.num_time_indexes; ++j) {
		for (int i = 0; i < input_params.num_vehicle_formation; ++i) {
			x_coefficient.push_back(input_params.vehicle_operational_cost[i]);
		}
	}
	for (int j = 0; j < input_params.num_time_indexes; ++j) {
		q_coefficient.push_back(input_params.unit_time_waiting_cost*input_params.length_per_time_interval);
	}

	/*** model  ***/
	GRBEnv* env = 0;
	try {
		// start timing
		clock_t start_time = clock(); 
		// set model parameters
		env = new GRBEnv();
		GRBModel model = GRBModel(*env);
		model.set(GRB_DoubleParam_TimeLimit, 3600);
		model.set(GRB_IntParam_OutputFlag, 0);
		model.set(GRB_DoubleParam_MIPGap, 0.000001);
		// add variables
		GRBVar* vars_x = new GRBVar[x_coefficient.size()];
		GRBVar* vars_q = new GRBVar[q_coefficient.size()];
		vars_x = model.addVars(x_coefficient.size(), GRB_BINARY);
		vars_q = model.addVars(q_coefficient.size(), GRB_CONTINUOUS);
		// add obj
		build_obj(model, vars_x, vars_q, x_coefficient, q_coefficient);
		// add constrs
		build_constrs(model, vars_x, vars_q, input_params);
		// optimize
		model.optimize();
		// end timing
		clock_t end_time = clock(); 
		cout << "The computation time is " << (end_time - start_time) / double(CLOCKS_PER_SEC) << " seconds" << endl;
		// get results
		double opt_obj;
		get_results(model, vars_x, vars_q, input_params, opt_obj, "../data/");
		
		cout << "The optimal objective value is " << opt_obj << " dollars" << endl;

		} catch(GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		} catch(...) {
			cout << "Exception during optimization" << endl;
		}
		delete env; 

	return 0;
}
