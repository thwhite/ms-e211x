%% MS&E Project Code
% Thomas White 11/26/21
clear all;
close all;
clc;

%% Problem Setup 

rng(343, 'twister');

% Problem constraints
n_customers = 10000;
n_items = 10;

% Initialize arrays
pi_bids = zeros(n_items, n_customers);
A_requests = zeros(n_items, n_customers);
b_inventory = [];
true_prices = [7,1,3,4,4,7,4,7,3,4];

for i = 1:n_customers
    A_requests(:, i) = rand(1, 10) < 0.5; % assign each a 1 or 0 randomly
    pi_bids(:, i) = true_prices'.*A_requests(:, i) + sqrt(0.2)*randn(10, 1).*A_requests(:, i);
    % note about this - it's not the listed formulation, but it seems like
    % we need that or people will bid les than 0
end 

%% Solve the Problem Offline



%% PROBLEM 1




%% Problem 2




%% Problem 3



%% Problem 4


