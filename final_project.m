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
n_inventory = 1000;

% Initialize arrays
pi_bids = zeros(1, n_customers);
A_requests = zeros(n_items, n_customers);
b_inventory = n_inventory*ones(1, n_items);
true_prices = [7,1,3,4,4,7,4,7,3,4];

for i = 1:n_customers
    A_requests(:, i) = rand(1, 10) < 0.5; % assign each a 1 or 0 randomly
    pi_bids(i) = true_prices*A_requests(:, i) + sqrt(0.2)*randn;

end 

%% Solve the Problem Offline
[x_offline, offline_revenue, ~, ~, lambda_offline] = linprog(-pi_bids', A_requests, b_inventory, [], [], zeros(size(pi_bids)), ones(size(pi_bids)));


%% PROBLEM 1

ks = [50, 100, 200, 400, 800, 1600];
revenues = [];
for j = 1:length(ks)

    
    revenue = 0;

    % get the shadow prices
    A_requests_SLPM = A_requests(:, 1:ks(j));
    pi_bids_SLPM = pi_bids(1:ks(j));
    b_inventory_SLPM = b_inventory*(ks(j)/n_customers);
    [x,~,~,~,lambda] = linprog(-pi_bids_SLPM', A_requests_SLPM, b_inventory_SLPM, [], [], zeros(size(pi_bids)), ones(size(pi_bids)));
    
    % use the shadow prices to make future purchasing decisions
    decision = [];
    for l = (ks(j)+1):n_customers % for all the remaining customers
        if (pi_bids(l) >= lambda.ineqlin'*A_requests(:, l)) && all(b_inventory' > A_requests(:, l))
            b_inventory = b_inventory - A_requests(:, l)';
            revenue = revenue + pi_bids(l);
        end 
    end 
    revenues(j) = revenue;
    ratios(j) = revenue/offline_revenue;

end 

disp("Problem 1 Results")



%% Problem 2



ks = [1, 100, 200, 400, 800, 1600, 3200, 6400];
test_run_length = 50;

revenue_dyn = 0;
b_inventory_dyn = n_inventory*ones(1, n_items);

for m = 1:n_customers


    % set shadow prices

    if ismember(m, ks)
        % update the shadow prices
        A_requests_dyn_upd = A_requests(:, m:m+test_run_length);
        pi_bids_dyn_upd = pi_bids(m:m+test_run_length);
        b_inventory_dyn_upd = b_inventory_dyn/(m:m+test_run_length/n_customers);
        [x,~,~,~,lambda] = linprog(-pi_bids_dyn_upd', A_requests_dyn_upd, b_inventory_dyn_upd, [], [], zeros(size(pi_bids)), ones(size(pi_bids)));
    end 
    
    % use the shadow prices to make future purchasing decisions

    if (pi_bids(m) >= lambda.ineqlin'*A_requests(:, m)) && all(b_inventory_dyn' >= A_requests(:, m))
        b_inventory_dyn = b_inventory_dyn - A_requests(:, m)';
        revenue_dyn = revenue_dyn + pi_bids(m);
    end
end


ratio_dyn = revenue_dyn/offline_revenue;


%% Problem 3



%% Problem 4


