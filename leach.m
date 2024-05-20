% LEACH Protocol 
% Clearing workspace and command window
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Field Dimensions - Maximum x and y dimensions of the network area (in meters)
xm = 100; % x dimension
ym = 100; % y dimension
% Coordinates of the Sink - Position of the base station
sink.x = 0.5 * xm;
sink.y = 1.75 * ym;
% Number of Nodes in the field
n = 100;
% Optimal Election Probability of a node to become cluster head
p = 0.05;
% Energy Model (all values in Joules)
Eo = 0.5; % Initial Energy of nodes
ETX = 50 * 0.000000001; % Energy consumed during transmission
ERX = 50 * 0.000000001; % Energy consumed during reception
Efs = 10 * 0.000000000001; % Energy for free space model
Emp = 0.0013 * 0.000000000001; % Energy for multi-path fading model
EDA = 5 * 0.000000001; % Energy consumed during data aggregation
% Values for Heterogeneity 
m = 0.05; % Percentage of nodes that are advanced
% Parameter alpha for energy model
a = 0.1;
% Maximum number of rounds for the simulation
rmax = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
% Computation of do - Calculating the optimal distance for energy calculations
do = sqrt(Efs / Emp);

% Creation of the random Sensor Network
figure(1);
hold on;
for i = 1:n
    S(i).xd = rand * xm;
    S(i).yd = rand * ym;
    S(i).G = 0;
    S(i).type = 'N'; % initially there are no cluster heads only nodes
    if i > m * n
        S(i).E = Eo;
        S(i).ENERGY = 0;
        plot(S(i).xd, S(i).yd, 'o');
    else
        S(i).E = Eo * (1 + a);
        S(i).ENERGY = 1;
        plot(S(i).xd, S(i).yd, '+');
    end
end
S(n + 1).xd = sink.x;
S(n + 1).yd = sink.y;
plot(S(n + 1).xd, S(n + 1).yd, 'x');
hold off;

% Network Lifetime Graph Initialization
alive_nodes = zeros(rmax + 1, 1); % Store number of alive nodes per round
average_energy = zeros(rmax + 1, 1); % Store average energy per round
packets_TO_BS = zeros(rmax + 1, 1); % Store number of packets to base station per round
packets_TO_CH = zeros(rmax + 1, 1); % Store number of packets to cluster head per round
pdr = zeros(rmax + 1, 1); % Store packet delivery ratio per round
total_energy_consumed = zeros(rmax + 1, 1); % Store total energy consumed per round
cluster_heads = zeros(rmax + 1, 1); % Store number of cluster heads per round

figure(2);
hold on;
xlabel('Rounds');
ylabel('Number of Alive Nodes');
title('Network Lifetime');
h_alive_nodes = plot(alive_nodes, 'g-'); % Initialize empty plot line (green) for alive nodes

figure(3);
hold on;
xlabel('Rounds');
ylabel('Average Node Energy (Normalized)');
title('Average Node Energy');
h_avg_energy = plot(average_energy, 'b-'); % Initialize empty plot line (blue) for avg energy

figure(4);
hold on;
xlabel('Rounds');
ylabel('Packet Delivery Ratio (%)');
title('Packet Delivery Ratio');
h_pdr = plot(pdr, 'm-'); % Initialize empty plot line (magenta) for PDR

figure(5);
hold on;
xlabel('Rounds');
ylabel('Total Energy Consumption (Normalized)');
title('Energy Consumption');
h_total_energy = plot(total_energy_consumed, 'r-'); % Initialize empty plot line (red) for energy consumption

figure(6);
hold on;
xlabel('Rounds');
ylabel('Number of Cluster Heads');
title('Cluster Head Distribution');
h_cluster_heads = plot(cluster_heads, 'k-'); % Initialize empty plot line (black) for cluster heads

figure(7);
hold on;
xlabel('Rounds');
ylabel('Throughput (Packets Delivered)');
title('Throughput');
h_throughput = plot(packets_TO_BS, 'c-'); % Initialize empty plot line (cyan) for throughput

hold off;

% First Iteration - Initializing the first iteration
flag_first_dead = 0; % Initialize flag for the first dead node
DEAD = zeros(rmax + 1, 1); % Store number of dead nodes per round
DEAD_A = zeros(rmax + 1, 1); % Store number of dead advanced nodes per round
DEAD_N = zeros(rmax + 1, 1); % Store number of dead normal nodes per round
alive_nodes_percentage = 100; % Initialize to 100%
r = 0; % Initialize round counter
first_dead = inf; % Initialize round of first dead node (set to infinity)
while alive_nodes_percentage > 10

    % Increment round counter
    r = r + 1;
    disp(['Current round: ', num2str(r)]);
    
    % LEACH Protocol with Average Node Energy

    % ... (rest of the code remains the same until here) ...

    if mod(r, round(1/p)) == 0
        for i = 1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end

    % Number of dead nodes
    dead = 0;
    dead_a = 0;
    dead_n = 0;
    % Counter for bits transmitted to Base Station and to Cluster Heads
    packets_TO_BS(r + 1) = 0;
    packets_TO_CH(r + 1) = 0;
    % Clear the current figure and hold the plot for animation
    figure(1);
    clf;
    hold on;
    % Plot the nodes
    for i = 1:n
        if S(i).E <= 0
            plot(S(i).xd, S(i).yd, 'r.', 'MarkerSize', 10); % Dead node
            dead = dead + 1;
            if S(i).ENERGY == 1
                dead_a = dead_a + 1;
            else
                dead_n = dead_n + 1;
            end
            % Update first_dead if this is the first dead node
            if r < first_dead
                first_dead = r;
            end
        else
            if S(i).ENERGY == 0
                plot(S(i).xd, S(i).yd, 'bo', 'MarkerSize', 8); % Normal node
            else
                plot(S(i).xd, S(i).yd, 'g+', 'MarkerSize', 8); % Advanced node
            end
        end
    end
    % Plot the sink
    plot(S(n + 1).xd, S(n + 1).yd, 'kx', 'MarkerSize', 12, 'LineWidth', 2); % Sink
    % Elect Cluster Heads
    countCHs = 0;
    cluster = 1;
    for i = 1:n
        if S(i).E > 0 && S(i).G <= 0
            if rand <= (p / (1 - p * mod(r, round(1/p))))
                countCHs = countCHs + 1;
                packets_TO_BS(r + 1) = packets_TO_BS(r + 1) + 1;
                S(i).type = 'C';
                S(i).G = round(1/p) - 1;
                plot(S(i).xd, S(i).yd, 'k*', 'MarkerSize', 10); % Cluster head
                distance = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
                if distance > do
                    S(i).E = S(i).E - ((ETX + EDA) * 4000 + Emp * 4000 * distance^4);
                else
                    S(i).E = S(i).E - ((ETX + EDA) * 4000 + Efs * 4000 * distance^2);
                end
                cluster = cluster + 1;
            end
        end
    end

    % Associate normal nodes with the nearest cluster-head
    for i = 1:n
        if S(i).type == 'N' && S(i).E > 0
            min_dis = inf;
            min_dis_cluster = 0;
            for c = 1:n
                if S(c).type == 'C'
                    distance = sqrt((S(i).xd - S(c).xd)^2 + (S(i).yd - S(c).yd)^2);
                    if distance < min_dis
                        min_dis = distance;
                        min_dis_cluster = c;
                    end
                end
            end
            if min_dis_cluster > 0
                S(i).cl = min_dis_cluster;
                packets_TO_CH(r + 1) = packets_TO_CH(r + 1) + 1;
                if min_dis > do
                    S(i).E = S(i).E - (ETX * 4000 + Emp * 4000 * min_dis^4);
                else
                    S(i).E = S(i).E - (ETX * 4000 + Efs * 4000 * min_dis^2);
                end
            end
        end
    end

    % Update the energy of the cluster-heads for receiving data from member nodes
    for c = 1:n
        if S(c).type == 'C'
            for i = 1:n
                if S(i).cl == c
                    distance = sqrt((S(i).xd - S(c).xd)^2 + (S(i).yd - S(c).yd)^2);
                    if distance > do
                        S(c).E = S(c).E - (ERX * 4000 + Emp * 4000 * distance^4 + EDA * 4000); % Include energy for data aggregation (EDA)
                    else
                        S(c).E = S(c).E - (ERX * 4000 + Efs * 4000 * distance^2 + EDA * 4000); % Include energy for data aggregation (EDA)
                    end
                end
            end
        end
    end
    hold off;

    % Calculate Average Node Energy
    total_energy = 0;
    for i = 1:n
        if S(i).E > 0
            total_energy = total_energy + S(i).E;
        end
    end
    average_energy(r + 1) = total_energy / (n - dead);

    % Update plot data
    alive_nodes(r + 1) = n - dead;
    cluster_heads(r + 1) = countCHs;
    pdr(r + 1) = packets_TO_BS(r + 1) / (packets_TO_CH(r + 1) + packets_TO_BS(r + 1));
    total_energy_consumed(r + 1) = total_energy;

    % Display some statistics on the current round
    disp(['Number of dead nodes: ', num2str(dead)]);
    disp(['Number of alive nodes: ', num2str(n - dead)]);
    disp(['Number of Cluster Heads: ', num2str(countCHs)]);
    
    % Update Network Lifetime Graph Data
    alive_nodes_percentage = ((n - dead) / n) * 100;
    DEAD(r + 1) = dead;
    DEAD_A(r + 1) = dead_a;
    DEAD_N(r + 1) = dead_n;

    % Update plots
    if ishghandle(h_alive_nodes)
        set(h_alive_nodes, 'YData', alive_nodes); % Update alive nodes line data
    end
    if ishghandle(h_avg_energy)
        set(h_avg_energy, 'YData', average_energy); % Update avg energy line data
    end
    if ishghandle(h_pdr)
        set(h_pdr, 'YData', pdr); % Update PDR line data
    end
    if ishghandle(h_total_energy)
        set(h_total_energy, 'YData', total_energy_consumed); % Update total energy consumption line data
    end
    if ishghandle(h_cluster_heads)
        set(h_cluster_heads, 'YData', cluster_heads); % Update cluster heads line data
    end
    if ishghandle(h_throughput)
        set(h_throughput, 'YData', packets_TO_BS); % Update throughput line data
    end

    % Check if the network is dead (all nodes depleted)
    if alive_nodes_percentage <= 10
        break;
    end
end

% Network Lifetime Analysis
disp(['Total Rounds: ', num2str(r)]);
disp(['First Node Dead at Round: ', num2str(first_dead)]);
disp(['Number of Dead Advanced Nodes: ', num2str(dead_a)]);
disp(['Number of Dead Normal Nodes: ', num2str(dead_n)]);

% Plot Dead Nodes vs Rounds
figure(8);
plot(DEAD);
xlabel('Rounds');
ylabel('Number of Dead Nodes');
title('Dead Nodes vs Rounds');

% Plot Packets Transmitted to BS vs Rounds
figure(9);
plot(1:r, packets_TO_BS(1:r));
xlabel('Rounds');
ylabel('Packets Transmitted to BS');
title('Packets Transmitted to BS vs Rounds');

% Plot Histogram of Remaining Energy
remaining_energy_data = zeros(1, n);
for i = 1:n
    remaining_energy_data(i) = S(i).E;
end
figure(10);
bar(1:n, remaining_energy_data); % Bar graph with x-axis as nodes and y-axis as remaining energy
xlabel('Node Index');
ylabel('Remaining Energy');
title('Remaining Energy of Nodes');
