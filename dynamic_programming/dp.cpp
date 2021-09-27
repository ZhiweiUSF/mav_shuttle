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
 * This script implements the Dynamic Programmingm (DP) algorithm in 
 * Chen, Z., Li, X., & Zhou, X. (2019). Operational design for shuttle systems with modular vehicles under oversaturated traffic: Discrete modeling method. Transportation Research Part B: Methodological, 122, 1-19.
 * Author: Zhiwei Chen, 08/31/2021
 * email: zhiweic@usf.edu
 */


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

int main() {

	//********************************************* input **************************************************//
	double epsilon = 0.1;
	double exchange_rate = 6.6; // coefficient to convert RMB to US dollars
	double unit_electricity_cost = 0.8; // RMB/kWh
	double length_of_line = 18.964; // km
	double fixed_operational_cost = 2.049 * unit_electricity_cost * length_of_line; // c_f, RMB
	double variable_operational_cost= 0.155 * unit_electricity_cost * length_of_line; // c_v, RMB
	double power_index_in_operational_cost = 0.5;  // alpha
	int capacity_per_unit = 226; // c, passengers/unit
	double unit_time_waiting_cost = 40.0 / 60.0; // w, RMB/minute
	int min_headway = 3; // h_underbar, minutes

	vector<int> vehicle_formation{0, 1, 2, 3, 4, 5, 6}; // set I, 6 types of vehicles, 0 means not dispatching anything
	vector<double> vehicle_capacity(7); // c*i, for all i 
	vector<double> vehicle_operational_cost(7); // f_i
	for (int i = 1; i < vehicle_formation.size(); ++i) {
		vehicle_capacity[i] = vehicle_formation[i] * capacity_per_unit;
		vehicle_operational_cost[i] = fixed_operational_cost + variable_operational_cost * pow(vehicle_capacity[i],
			 power_index_in_operational_cost);
	}
	//print_one_dimensional_vector(vehicle_operational_cost);

	
	double length_per_time_interval = 1.5; // delta, minute
	int min_headway_interval = ceil(min_headway/length_per_time_interval); // j_underbar
	vector<double> passenger_demand; // a_j
	read_file("../data/2.csv", passenger_demand);
	//print_one_dimensional_vector(passenger_demand);
	
	int num_vehicle_formation = vehicle_formation.size(); // I
	int num_time_indexes = passenger_demand.size(); // J

	vector<double> cumulative_demand(num_time_indexes); // A_j
	cumulative_demand[0] = passenger_demand[0];
	for (int j = 1; j < num_time_indexes; ++j) {
		cumulative_demand[j] = cumulative_demand[j-1] + passenger_demand[j];
	}
	//print_one_dimensional_vector(cumulative_demand);


	//*********************************** detect oversaturated traffic **********************************//
	vector<double> demand_in_min_headway_interval(num_time_indexes); // for each j, this vector stores the demand from j-min_headway_interval+1 to j
	vector<int> start_point_set; // set of starting points of oversaturated periods
	vector<int> end_point_set; // set of ending points of oversaturated periods
	vector<double> oversaturated_queue(num_time_indexes); // oversaturated queue, delta q_j

	// impletation of algorithm 1 in the paper to find the starting and ending points of each oversaturated period
	for (int j = min_headway_interval; j < num_time_indexes; ++j) {
		demand_in_min_headway_interval[j] += passenger_demand[j];
		for (int k = j-min_headway_interval+1; k < j; ++k) {
			demand_in_min_headway_interval[j] += passenger_demand[k];
		}
	}

	for (int j = min_headway_interval; j < num_time_indexes; ++j) {
		if (demand_in_min_headway_interval[j] > vehicle_capacity[num_vehicle_formation-1]) {
			start_point_set.push_back(j);
			int tem;
			for (int m = j+1; m < num_time_indexes; ++m) {
				if ( (cumulative_demand[m]-cumulative_demand[j]) / (m-j) <= 
						vehicle_capacity[num_vehicle_formation-1]/min_headway_interval ) {
					end_point_set.push_back(m);
					tem = m;
					break;
				} // end if
			} // end for m
			j = tem;
		} // end if
	} // end for j

	//print_one_dimensional_vector(start_point_set);
	//print_one_dimensional_vector(end_point_set);

	// compute oversaturated queue
	for (int i = 0; i < start_point_set.size(); ++i) { 
		for (int j = start_point_set[i]; j < end_point_set[i]; ++j) {
			oversaturated_queue[j] = cumulative_demand[j] - cumulative_demand[start_point_set[i]] -
						 (j-start_point_set[i])*length_per_time_interval*vehicle_capacity[num_vehicle_formation-1]/min_headway;
			if (oversaturated_queue[j] < 0){
				oversaturated_queue[j] = 0;
			}
		}
	}



	//print_one_dimensional_vector(oversaturated_queue);
	clock_t start_time = clock(); // start timing

 	// ******************************************* dp parameters ************************************************//
	vector<vector<int> > formation_record(num_time_indexes); // record of vehicle formation dispatched
	vector<vector<int> > back_tracking_info(num_time_indexes); // record information for backtracking to find the optimal solution
	vector<double> queue; // queue: number of passenger left at station, state variable q_j in the paper
	vector<double> cost;  // cost: sum of operational cost and passenger waiting cost
	vector<int> headway;  // number of time intervals from the previosu dispatch to time j+1 and min dispatch interval, state variable h_hat_j in the paper
	vector<int> bounded_headway;  //  minimum between the # time intervals from the previosu dispatch to time j+1 and min dispatch interval, state variable h_j in the paper

	vector<double> queue_prev; // queue: number of passenger left at station, state variable q_j in the paper
	vector<double> cost_prev;  // cost: sum of operational cost and passenger waiting cost
	vector<int> bounded_headway_prev;  // minimum between the # time intervals from the previosu dispatch to time j+1 and min dispatch interval, state variable h_j in the paper
	vector<int> headway_prev;  // # time intervals from the previosu dispatch to time j+1 and min dispatch interval, state variable h_hat_j in the paper	


	// ****************************************** dp initialization *******************************************//
	// create variables to temporarily store intermedate information
	vector<int> current_formation_record; // vehcle formation for all solutions in the current iteration
	vector<int> current_back_tracking_info; // info in the current iteration that will be used for backtracking
	double queue_state, cost_state;
        int headway_state, bounded_headway_state;
	double upper_bound_0, upper_bound_1, upper_bound_2, upper_bound, var1, var2; // these are used to compute queue bound

	// at time 0, enumerate all possible vehicle formations to dispatch
	for (int i = 0; i < num_vehicle_formation; ++i) {
		// update states and costs using state transition functions
		queue_state = (passenger_demand[0] - vehicle_capacity[i] < 0) ? 0 : passenger_demand[0] - vehicle_capacity[i];
		bounded_headway_state = (i==0) ? min_headway_interval : 1;
		headway_state = (i==0) ? -2000000 : 0;
		cost_state = queue_state*unit_time_waiting_cost*length_per_time_interval+vehicle_operational_cost[i];

		// add solutions
		if (queue_prev.size() == 0) { // if there have no solutions been found, just add this new solution
			queue_prev.push_back(queue_state);
			cost_prev.push_back(cost_state);
			bounded_headway_prev.push_back(bounded_headway_state);
			headway_prev.push_back(headway_state);
			current_formation_record.push_back(i);
			current_back_tracking_info.push_back(0);
		} else {		     // if there are already some solutions found
			int state_already_exist = 0; // a variable to denote whether the state (q_j, h_j) of the new solution has already been found before
			for (int i1 = 0; i1 < queue_prev.size(); ++i1) { // enumerate all existing solutions to see if we have already found a solution with the same state of the new solutions
				if (abs(queue_prev[i1] - queue_state) < epsilon && bounded_headway_prev[i1] == bounded_headway_state) {  // if yes, we see if the new solution produces lower cost
					state_already_exist = 1;
					if (cost_prev[i1] > cost_state) { // we add the new solution only if the new solution produces lower cost, i.e., it is better
						cost_prev[i1] = cost_state;
						current_formation_record[i1] = i;
						current_back_tracking_info[i1] = 0;
					}
				}
			} // end for
			if (state_already_exist == 0) {	// if we do not find the state of the new solution among existing solutions, add this new state
				queue_prev.push_back(queue_state);
				cost_prev.push_back(cost_state);
				bounded_headway_prev.push_back(bounded_headway_state);
				headway_prev.push_back(headway_state);
				current_formation_record.push_back(i);
				current_back_tracking_info.push_back(0);
			}
			
		} // end  if
	} // end for i
	formation_record[0] = current_formation_record;  // store current_formation_record from the current iteration
	back_tracking_info[0] = current_back_tracking_info;  // store current_back_tracking_info from the current iteration


	// ********************************************* dp body ************************************************//
	for (int j = 1; j < num_time_indexes; ++j) {
		// clear vectors at the beginning of each iteration
		queue.clear(); 
		cost.clear();  
		headway.clear();  
		bounded_headway.clear();  
		current_formation_record.clear();
		current_back_tracking_info.clear();

		// compute queue bound, see Eq.(10)
		var2 = (1 + floor( (num_time_indexes-j) / min_headway_interval ) ) * vehicle_capacity[num_vehicle_formation-1];
		var1 = capacity_per_unit+oversaturated_queue[j];
		upper_bound_1 = (var1 < var2) ? var1 : var2; // second branch in Eq.(10)

		var1 = (num_vehicle_formation+1)*capacity_per_unit + oversaturated_queue[j];
		upper_bound_0 = (var1 < var2) ? var1 : var2; // third branch in Eq.(10)

		var1 = oversaturated_queue[j];
		upper_bound_2 = (var1 < var2) ? var1 : var2; // first branch in Eq.(10)

		for (int k = 0; k < queue_prev.size(); ++k) {   // enumerate all previous solutions
			if (bounded_headway_prev[k] < min_headway_interval) { // no vehicles can be dispatched
				// update states and costs using state transition functions
				queue_state = (queue_prev[k] + passenger_demand[j] > 0) ? queue_prev[k] + passenger_demand[j]: 0;
				bounded_headway_state = bounded_headway_prev[k] + 1;
				headway_state = headway_prev[k] + 1;
				cost_state = cost_prev[k] + queue_prev[k]*unit_time_waiting_cost*length_per_time_interval;
				// if the resultant queue is less than or equal to the upper bound, check if we need to add it, the process is the same as that in the initialization process 
				if (queue_state <= upper_bound_0) { 
					if (queue.size() == 0) {
						queue.push_back(queue_state);
						cost.push_back(cost_state);
						bounded_headway.push_back(bounded_headway_state);
						headway.push_back(headway_state);
						current_formation_record.push_back(0);
						current_back_tracking_info.push_back(k);
					} else {
						int state_already_exist = 0;
						for (int i1 = 0; i1 < queue.size(); ++i1) {
							if (abs(queue[i1] - queue_state) < epsilon && bounded_headway[i1] == bounded_headway_state) { 
								state_already_exist = 1;
								if (cost[i1] > cost_state) {
									cost[i1] = cost_state;
									current_formation_record[i1] = 0;
									current_back_tracking_info[i1] = k;
								} // end if cost[i1]	
							} // end if abs
						} // end for
						if (state_already_exist == 0) {	
							queue.push_back(queue_state);
							cost.push_back(cost_state);
							bounded_headway.push_back(bounded_headway_state);
							headway.push_back(headway_state);
							current_formation_record.push_back(0);
							current_back_tracking_info.push_back(k);
						} // end if state_already_exist == 0
					} // end if queue.size() == 0
				} // end if queue_sate
			} else { // vehicles can be dispatched
				for (int i = 0; i < num_vehicle_formation; ++i) {  // enumerate all possible solutions, including 0 - not dispatching anything
					// update state variables and costs according to the state transition equations
					queue_state = (queue_prev[k] + passenger_demand[j] - vehicle_capacity[i]>0) ? queue_prev[k] + passenger_demand[j] - vehicle_capacity[i] : 0;
					bounded_headway_state = (i==0) ? bounded_headway_prev[k] + 1 : 1;
					headway_state = headway_prev[k] + 1;
					cost_state = cost_prev[k] + queue_prev[k]*unit_time_waiting_cost*length_per_time_interval+vehicle_operational_cost[i];
					// compute upper bound to passenger queue
					upper_bound;
					if (i == 0) {
						upper_bound = upper_bound_0;
					} else {
						if (headway_state > min_headway_interval) {
							upper_bound = upper_bound_2;
						} else {
							upper_bound = upper_bound_1;
						}
					}
					// update headway state and bounded headway state
					if (i != 0) {headway_state = 0;}
					if (bounded_headway_state > min_headway_interval) {bounded_headway_state = min_headway_interval;}
					// if the resultant queue is less than or equal to the upper bound, check if we need to add it, the process is the same as that in the initialization process 
					if (queue_state <= upper_bound) {
						if (queue.size() == 0) {
							queue.push_back(queue_state);
							cost.push_back(cost_state);
							bounded_headway.push_back(bounded_headway_state);
							headway.push_back(headway_state);
							current_formation_record.push_back(i);
							current_back_tracking_info.push_back(k);
						} else {
							int state_already_exist = 0;
							for (int i1 = 0; i1 < queue.size(); ++i1) {
								if (abs(queue[i1] - queue_state) < epsilon && bounded_headway[i1] == bounded_headway_state) { 
									state_already_exist = 1;
									if (cost[i1] > cost_state) {
										cost[i1] = cost_state;
										current_formation_record[i1] = i;
										current_back_tracking_info[i1] = k;
									} // end if cost[i1]	
								} // end if abs
							} // end for
							if (state_already_exist == 0) {	
								queue.push_back(queue_state);
								cost.push_back(cost_state);
								bounded_headway.push_back(bounded_headway_state);
								headway.push_back(headway_state);
								current_formation_record.push_back(i);
								current_back_tracking_info.push_back(k);
							} // end if state_already_exist == 0
						} // end if queue.size() == 0
					} // end if

				} // end for i	
			} // end if
		} // end for  k

		// pass information from the current iteration to variablels with "prev"
		formation_record[j] = current_formation_record;
		back_tracking_info[j] = current_back_tracking_info;	
		queue_prev = queue;
		cost_prev = cost;
		bounded_headway_prev = bounded_headway;
		headway_prev = headway;
	} // end  for j 


	//*********************************** back track to find optimal solution **********************************//
	// find optimal solution objective value and all dispatches
	double opt_obj;  // optimal objective value
	int optimal_solution_index; // position index of the solution at each iteration
	for (int j = 0; j < cost_prev.size(); ++j) {
		if (queue_prev[j] == 0) {
			opt_obj = cost_prev[j];
			optimal_solution_index = j;
		}
	}
	// add constant term into objective function
	opt_obj = (opt_obj + cumulative_demand[num_time_indexes - 1]*unit_time_waiting_cost*length_per_time_interval/2)/exchange_rate;
	// find optimal dispatch and vehicle formations
	vector<int> optimal_dispatch_time;
	vector<int> optimal_vehicle_formation;
	for (int j = num_time_indexes - 1; j > 0; --j)  {
		if (formation_record[j][optimal_solution_index] != 0) {
			optimal_dispatch_time.push_back(j*length_per_time_interval);
			optimal_vehicle_formation.push_back(formation_record[j][optimal_solution_index]);
		}
		optimal_solution_index = back_tracking_info[j][optimal_solution_index];
		
	}
	// reverse the time and type for dispatches
	reverse(optimal_dispatch_time.begin(), optimal_dispatch_time.end());
	reverse(optimal_vehicle_formation.begin(), optimal_vehicle_formation.end()); 
	//print_one_dimensional_vector(optimal_dispatch_time_index);
	//print_one_dimensional_vector(optimal_vehicle_formation);

	// write resutls to csv
	print_results("../data/opt_time.csv", optimal_dispatch_time);
	print_results("../data/opt_formation.csv", optimal_vehicle_formation);


	clock_t end_time = clock(); // end timing
	cout << "The computation time is " << (end_time - start_time) / double(CLOCKS_PER_SEC) << " seconds" << endl;
	cout << "The optimal objective value is " << opt_obj << " dollars" << endl;	
	

}
