%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field Dimensions - Maximum x and y dimensions of the network area (in meters)
xm = 100; % x dimension
ym = 100; % y dimension
% Coordinates of the Sink - Position of the base station
sink.x = 0.5 * xm;
sink.y = 0.5 * ym;
% Number of Nodes in the field
n = 100;
% Optimal Election Probability of a node to become cluster head
p = 0.1;
% Energy Model (all values in Joules)
% Initial Energy of nodes
Eo = 0.5;
% Energy consumption parameters for various operations
ETX = 50 * 0.000000001; % Energy consumed during transmission
ERX = 50 * 0.000000001; % Energy consumed during reception
Efs = 10 * 0.000000000001; % Energy for free space model
Emp = 0.0013 * 0.000000000001; % Energy for multi-path fading model
EDA = 5 * 0.000000001; % Energy consumed during data aggregation
% Maximum number of rounds for the simulation
rmax = 3000;
% Threshold for number of alive nodes (10% of initial nodes)
alive_threshold = 0.1 * n;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
% Computation of do - Calculating the optimal distance for energy calculations
do = sqrt(Efs / Emp);

% Creation of the random Sensor Network
S = struct('xd', num2cell(rand(1, n) * xm), ...
           'yd', num2cell(rand(1, n) * ym), ...
           'G', num2cell(zeros(1, n)), ...
           'type', num2cell(repmat('N', 1, n)), ...
           'E', num2cell(repmat(Eo, 1, n)), ...
           'remaining_energy', num2cell(repmat(Eo, 1, n)));  % Added field for remaining energy

figure(1);
for i = 1:n
    plot(S(i).xd, S(i).yd, 'o');
    hold on;
end
S(n + 1).xd = sink.x;
S(n + 1).yd = sink.y;
plot(S(n + 1).xd, S(n + 1).yd, 'x');
% Initialize statistics arrays
STATISTICS.DEAD = zeros(rmax, 1);
STATISTICS.DEAD_N = zeros(rmax, 1);
STATISTICS.CLUSTERHEADS = zeros(rmax, 1);
STATISTICS.PACKETS_TO_BS = zeros(rmax, 1);
STATISTICS.PACKETS_TO_CH = zeros(rmax, 1);
STATISTICS.ENERGY = zeros(rmax, 1);
% First Iteration - Initializing the first iteration
figure(1);
% Counter for CHs - Counts the number of cluster heads
countCHs = 0;
% Counter for CHs per round - Counts the number of cluster heads per round
rcountCHs = 0;
% Initializing cluster counter
cluster = 1;
% Initialize counters
rcountCHs = rcountCHs + countCHs;
% Flag for the first dead node
flag_first_dead = 0;
% Total initial energy
total_initial_energy = n * Eo;
% Run simulation until 10% of nodes are alive
r = 0;
while n - STATISTICS.DEAD > 10 
    % Reset cluster head selection every epoch
    if (mod(r, round(1/p)) == 0)
        for i = 1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end
    
    hold off;
    % Number of dead nodes
    dead = 0;
    % Number of dead Normal Nodes
    dead_n = 0;
    % Counter for bits transmitted to Base Station and to Cluster Heads
    packets_TO_BS = 0;
    packets_TO_CH = 0;
        % Energy consumption in the network
    total_energy_consumed = 0;
    figure(1);
    for i = 1:n
        % Checking if there is a dead node
        if (S(i).E <= 0)
            plot(S(i).xd, S(i).yd, 'red .');
            dead = dead + 1;
            dead_n = dead_n + 1;
            hold on;    
        else
            S(i).type = 'N';
            plot(S(i).xd, S(i).yd, 'o');
            hold on;
            
           distance_to_sink = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
           energy_spent_during_CH_selection = ETX * 4000 * distance_to_sink^2;  % Adjust based on your model

        end
    end
    plot(S(n + 1).xd, S(n + 1).yd, 'x');
    % Recording statistics about dead nodes
    STATISTICS.DEAD(r + 1) = dead;
    STATISTICS.DEAD_N(r + 1) = dead_n;
    % When the first node dies
    if (dead == 1 && flag_first_dead == 0)
        first_dead = r;
        flag_first_dead = 1;
    end

    countCHs = 0;
    cluster = 1;
    for i = 1:n
        if (S(i).E > 0)
            temp_rand = rand;
            if ((S(i).G) <= 0)
                % Election of Cluster Heads
                if (temp_rand <= (p / (1 - p * mod(r, round(1 / p)))))
                    countCHs = countCHs + 1;
                    packets_TO_BS = packets_TO_BS + 1;
                    STATISTICS.PACKETS_TO_BS(r + 1) = packets_TO_BS;
                    
                    S(i).type = 'C';
                    S(i).G = round(1 / p) - 1;
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    plot(S(i).xd, S(i).yd, 'k*');

                    distance = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    X(cluster) = S(i).xd;
                    Y(cluster) = S(i).yd;
                    cluster = cluster + 1;
                    
                    % Calculation of Energy dissipated (update remaining energy)
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * 4000 + Emp * 4000 * (distance^4));
                        S(i).remaining_energy = S(i).remaining_energy - ((ETX + EDA) * 4000 + Emp * 4000 * (distance^4)); 
                    else
                        S(i).E = S(i).E - ((ETX + EDA) * 4000 + Efs * 4000 * (distance^2));
                        S(i).remaining_energy = S(i).remaining_energy - ((ETX + EDA) * 4000 + Efs * 4000 * (distance^2)); 
                    end
                end     
            end
        end
    end
    STATISTICS.CLUSTERHEADS(r + 1) = cluster - 1;
    % Election of Associated Cluster Head for Normal Nodes
    for i = 1:n
        if (S(i).type == 'N' && S(i).E > 0)
            if (cluster - 1 >= 1)
                min_dis = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
                min_dis_cluster = 1;
                for c = 1:cluster - 1
                    temp = sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2);
                    if (temp < min_dis)
                        min_dis = temp;
                        min_dis_cluster = c;
                    end
                end
                % Energy dissipated by associated Cluster Head
                if (min_dis > do)
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - (ETX * 4000 + Emp * 4000 * (min_dis^4)); 
                    S(C(min_dis_cluster).id).remaining_energy = S(C(min_dis_cluster).id).remaining_energy - (ETX * 4000 + Emp * 4000 * (min_dis^4)); 
                else
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - (ETX * 4000 + Efs * 4000 * (min_dis^2)); 
                    S(C(min_dis_cluster).id).remaining_energy = S(C(min_dis_cluster).id).remaining_energy - (ETX * 4000 + Efs * 4000 * (min_dis^2)); 
                end
                % Energy dissipated (update remaining energy)
                if (min_dis > 0)
                    S(i).E = S(i).E - (ETX * 4000 + EDA * 4000); 
                    S(i).remaining_energy = S(i).remaining_energy - (ETX * 4000 + EDA * 4000); 
                    packets_TO_CH = packets_TO_CH + 1;
                end
                S(i).min_dis = min_dis;
                S(i).min_dis_cluster = min_dis_cluster;
            end
        end
    end
    STATISTICS.PACKETS_TO_CH(r + 1) = packets_TO_CH;
    
    % Calculate total energy consumption in the network
    for  i = 1:n
        total_energy_consumed = total_energy_consumed + (Eo - S(i).E);
    end
    STATISTICS.ENERGY(r + 1) = total_energy_consumed / total_initial_energy;
    % Display round information
    fprintf('Round %d: Alive Nodes = %d, Dead Nodes = %d\n', r, n - dead, dead);
    hold on;
    rcountCHs = rcountCHs + countCHs;
    r = r + 1;
end

% Plot Network Lifetime
figure;
plot(0:r-1, n - STATISTICS.DEAD(1:r));
xlabel('Time (rounds)');
ylabel('Number of alive nodes');
title('Network Lifetime');
% Plot Energy Consumption
figure;
plot(0:r-1, STATISTICS.ENERGY(1:r));
xlabel('Time (rounds)');
ylabel('Total energy consumption (normalized)');
title('Energy Consumption');
% Plot Cluster Head Distribution
figure;
plot(0:r-1, STATISTICS.CLUSTERHEADS(1:r));
xlabel('Time (rounds)');
ylabel('Number of cluster heads');
title('Cluster Head Distribution');

% Extract node IDs and remaining energy
node_ids = 1:n;  % Assuming n represents the number of nodes
remaining_energy_values = [S(1:n).remaining_energy];

% Create the histogram
figure;
bar(node_ids, remaining_energy_values);
xlabel('Node ID');
ylabel('Remaining Energy (Joules)');
title('Remaining Energy Distribution per Node');
