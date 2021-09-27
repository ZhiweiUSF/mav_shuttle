%%% This script implements the Continuum Approximation (CA) algorithm in
%%% Chen, Z., Li, X., & Zhou, X. (2020). Operational design for shuttle systems with modular vehicles under oversaturated traffic: Continuous modeling method. Transportation Research Part B: Methodological, 132, 76-100.
%%% Author: Zhiwei Chen
%%% Email: zhiweic@usf.edu

clear
clc

%% input parameters
% All parameters (unless stated) have the same name and meaning as those in the dynamic
% programming algorithm.  Values are from the Beijing case study
exchange_rate = 6.6;
unit_electricity_cost = 0.8;
length_of_line = 18.964;
fixed_operational_cost = 2.049 * unit_electricity_cost * length_of_line;
variable_operational_cost = 0.155 * unit_electricity_cost * length_of_line;
capacity_per_unit = 226;
length_per_time_interval = 3;
power_index_in_operational_cost = 0.5;
unit_time_waiting_cost = 40/60;
vehicle_formation = 1 : 6;
vehicle_capacity = vehicle_formation * capacity_per_unit;
min_headway = 3;
vehicle_operational_cost = fixed_operational_cost + variable_operational_cost * vehicle_capacity .^ power_index_in_operational_cost; %dispatch cost of each vehicle (Yuan)
min_headway_interval = ceil(min_headway/length_per_time_interval);
epsilon = 0.03; % acceptable error for discretization, see below for its usage

passenger_demand = csvread('../data/1.csv');
time_point_set = [0:1:length(passenger_demand)-1]; % set of discretized time points
num_time_indexes = length(time_point_set);
num_vehicle_formation = length(vehicle_formation);


%% preprocessing
% compute cumulative passenger demand
cumulative_demand = zeros(1,num_time_indexes);
cumulative_demand(1) = passenger_demand(1);
for j = 2 : num_time_indexes
    cumulative_demand(j) = cumulative_demand(j-1) + passenger_demand(j);
end
% compute oversaturated queue (see Algorithm 1 in paper)
oversaturated_queue = zeros(1,num_time_indexes);
for j = 1 : num_time_indexes
    tem = cumulative_demand(:,1:j) - vehicle_capacity(end)/min_headway * ([1:1:j]-1) * length_per_time_interval;
    oversaturated_queue(j) = cumulative_demand(j) - vehicle_capacity(end)/min_headway*(j-1)*length_per_time_interval - min(tem);
end
% compute preferred virtual arrival demand, which is called revised demand
% for short
revised_cumulative_demand = cumulative_demand - oversaturated_queue; % preferred virtual arrival demand
revised_demand = diff(revised_cumulative_demand); %dequeued arrival demand rate
revised_demand = [0 revised_demand];

tic % start timing


%% unit-time problem
%solve unit-time problem for each local neighborhood (i.e., each time
%point)

% define the following varaibles to store optimal solution
local_opt_vehicle_formation = []; %optimal vehicle formation
local_opt_dispatch_headway = []; %optimal dispatch headway
local_opt_cost = []; % optimal cost

for j = 1 : num_time_indexes % solve the unit-time problem for each time point
    
    current_best_cost = inf;
    current_best_headway = 0;
    current_best_formation = 0;
    
    for i = 1 : num_vehicle_formation % enumerate all possible vehicle formations to be dispatched to find the one with the min cost
        
        cost = inf;
        % use Eq.(17) and (18) to analytically solve dispatch headway
        if vehicle_capacity(i) >= roundn(revised_demand(j) /length_per_time_interval * min_headway, 0)
            if (2 * vehicle_operational_cost(i)) / (unit_time_waiting_cost * (vehicle_capacity(i))) > min_headway
                if revised_demand(j)/length_per_time_interval  < unit_time_waiting_cost * (vehicle_capacity(i)) ^ 2 / (2 * vehicle_operational_cost(i))
                    headway = sqrt(2 * vehicle_operational_cost(i) / unit_time_waiting_cost / (revised_demand(j)/length_per_time_interval));
                else
                    headway = (vehicle_capacity(i)) / (revised_demand(j) / length_per_time_interval);
                end
            else
                if revised_demand(j) / length_per_time_interval <= (2 * vehicle_operational_cost(i)) / (unit_time_waiting_cost * min_headway ^ 2)
                    headway = sqrt(2 * vehicle_operational_cost(i) / unit_time_waiting_cost / (revised_demand(j)/length_per_time_interval));
                else
                    headway = min_headway;
                end
            end
            % compute cost using Eq.(13)
            cost = roundn(vehicle_operational_cost(i)/headway + unit_time_waiting_cost*(revised_demand(j)/length_per_time_interval)*headway/2, -2);
        end
        
        %update the optimal solution for this time point
        if cost < current_best_cost
            current_best_cost = cost;
            current_best_headway = headway;
            current_best_formation = i;
        end
        
    end
    
    % store the optimal solution for each time point
    local_opt_cost = [local_opt_cost, current_best_cost];
    local_opt_dispatch_headway = [local_opt_dispatch_headway, current_best_headway];
    local_opt_vehicle_formation = [local_opt_vehicle_formation, current_best_formation];
end


%% plotting
%plot the results for each neighborhood, just for testing purposes
figure(1);clf;hold all;
plot(time_point_set, local_opt_cost);
figure(2);clf;hold all;
plot(time_point_set(2:end), local_opt_dispatch_headway(2:end));
plot([0,15],[0,15],':')
figure(3);clf;hold all;
plot(time_point_set, local_opt_vehicle_formation);


%% discretization (Algorithm 2 in the paper)
% discrete the headway: find the time point for each dispatch
current_time_index = length(time_point_set); % initilize current time point index
current_time = (current_time_index-1)*length_per_time_interval; % initialize current time value
opt_num_dispatch = 1; % number of dispatches, initialize as 1
opt_time_index = [length(time_point_set)]; % time index for each dispatch in the optimal solution 
first_time_index_without_demand = find(local_opt_dispatch_headway~=0); % the first point where demand is not 0
first_time_index_without_demand = first_time_index_without_demand(1);

while current_time > min_headway % keep looping until current time is less than the min dispatch headway (i.e., no dispatches can be found)
    error = inf;
    current_time_index = opt_time_index(opt_num_dispatch);
    
    % The following for loop is essentially to draw the 45 degree line as stated in the paper,
    % i.e., to find a point j between (first_time_index_without_demand, current_time_index-1)
    % such that h1 = current_time-t1. Since we search over a discrete time
    % horizon, a point perfectly satisfying this condition may not exist.
    % Thus, we define a parameter epsilon to find a point that 'almost'
    % satifies this condition
    for j = current_time_index-1 : -1 : first_time_index_without_demand
        t1 = (j - 1) * length_per_time_interval; % time value for this time index
        h1 = local_opt_dispatch_headway(j); % headway value for this time index
        
        if (current_time_index-j) >= min_headway_interval && abs(h1/(current_time-t1)-1) < error
            error = abs(h1/(current_time-t1)-(1)); %update error
            temtp = j; %temporary storage of time index
        end
        if error < epsilon
            break %jump out of the loop if we've already found a time index whose error is less than epsilon
        end
    end
    
    opt_num_dispatch = opt_num_dispatch + 1; % We have found one dispatch, increment dispatch index
    current_time = (current_time_index-1) * length_per_time_interval; % update current_time
    opt_time_index(opt_num_dispatch) = temtp; % record the newly found dispatch time
    if opt_time_index(opt_num_dispatch) == opt_time_index(opt_num_dispatch-1)
        break %if the time point can not be updated anymore, exit
    end
end
opt_time_index = opt_time_index(1:end-1); % delete the last one since it repeats the second last one
opt_time_index = fliplr(opt_time_index); % reverse the order
opt_time = (opt_time_index - 1) * length_per_time_interval; % convert time index to time values
opt_headway = diff(opt_time_index) * length_per_time_interval; % compute the headway
opt_headway = [(opt_time_index(1)-1)*length_per_time_interval opt_headway]; % add the headway of the first dispatch

% discretize the vehicle formation
starting_time_index = first_time_index_without_demand;
opt_formation = [local_opt_vehicle_formation(1)];
for j = 1 : length(opt_time_index)
    ending_time_index = opt_time_index(j);
    opt_formation(j) = sum(local_opt_vehicle_formation(:,starting_time_index:ending_time_index))/(ending_time_index-starting_time_index+1); %weighted sum
    starting_time_index = ending_time_index + 1;
end
opt_formation = round(opt_formation);
opt_formation(opt_formation > num_vehicle_formation) = num_vehicle_formation;

computation_time = toc; % end_timing

%%
%compute the objective value by discretization
opt_queue = zeros(1,1); 
opt_vehicle_operational_cost = zeros(1,1);
opt_passenger_waiting_cost = zeros(1,1);
opt_vehicle_occupancy = zeros(1,1);
z = 0;

for j = 2 : length(time_point_set)
    index = 0;
    for k = 1 : length(opt_time)
        if roundn((j - 1) * length_per_time_interval, -2) == roundn(opt_time(k),-2)
            z = z + 1;
            index = index + 1;
            opt_queue(j) = max(0, opt_queue(j-1) + passenger_demand(j) - vehicle_capacity(opt_formation(k)));
            opt_vehicle_operational_cost(j) =  vehicle_operational_cost(opt_formation(k));
            opt_passenger_waiting_cost(j) = (opt_queue(j-1) + passenger_demand(j)/2) * unit_time_waiting_cost * length_per_time_interval;
            opt_vehicle_occupancy(j) = min(opt_queue(j-1) + passenger_demand(j), vehicle_capacity(opt_formation(k)))/vehicle_capacity(opt_formation(k)) * 100;
        end
    end
    
    if index == 0
        opt_queue(j) = opt_queue(j-1) + passenger_demand(j);
        opt_vehicle_operational_cost(j) = 0;
        opt_passenger_waiting_cost(j) =  (opt_queue(j-1) + passenger_demand(j)/2) * unit_time_waiting_cost * length_per_time_interval;
        opt_vehicle_occupancy(j) = 0;
    end
end

% compute optimal objectiva value from the continuous solution
opt_obj_continuous = sum(local_opt_cost(local_opt_cost~=Inf)*length_per_time_interval)/exchange_rate ...
    + sum(oversaturated_queue)*unit_time_waiting_cost*length_per_time_interval/exchange_rate; %optimal objective value by integral
% compute optimal objectiva value from the discretized solution
opt_obj_discrete = (sum(opt_vehicle_operational_cost) + sum(opt_passenger_waiting_cost))/exchange_rate;

