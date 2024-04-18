function generateClusterEnergyHistogram(clusterModel)
% Generate a histogram of energy of clusters in the LEACH algorithm
%
%   Input:
%       clusterModel: The cluster model generated by the LEACH algorithm
%
% Mohammad Hossein Homaei, Homaei@wsnlab.org & Homaei@wsnlab.ir
% Ver 1. 10/2014

    % Extract necessary information from the cluster model
    clusterNode = clusterModel.clusterNode;
    nodeArch = clusterModel.nodeArch;

    % Initialize an array to store the energy levels of cluster heads
    clusterEnergy = zeros(1, clusterNode.countCHs);

    % Extract energy levels of cluster heads
    for i = 1:clusterNode.countCHs
        clusterEnergy(i) = nodeArch.node(clusterNode.no(i)).energy;
    end

    % Plot histogram of energy levels of cluster heads
    figure;
    histogram(clusterEnergy);
    xlabel('Energy');
    ylabel('Frequency');
    title('Histogram of Energy Levels of Cluster Heads');
end