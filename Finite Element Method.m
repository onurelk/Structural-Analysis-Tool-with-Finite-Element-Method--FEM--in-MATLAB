clc;
clear all;

% Elasticity Modulus
EI = 0.667 * 10E7;

% Input the number of bars (CS)
CS = input('Enter the Number of Bars=');

% Calculate the number of nodes (NS)
NS = CS + 1;

% Initialize an array to store the x-coordinates of nodes
kord = zeros(NS, 1);

% Input the x-coordinate values for each node
for i = 1:NS
    fprintf('Enter x-coordinate value for node %d=', i)
    kord(i, 1) = input('');
end

% Extract the x-coordinates
x = kord(:, 1);

% Initialize an array for y-coordinates
y = zeros(NS, 1);

% Initialize an array to store element node connectivity
EN = zeros(CS, 2);

% Define element node connectivity
for i = 1:CS
    EN(i, 1) = i;
    EN(i, 2) = i + 1;
end

% Calculate the total number of degrees of freedom (SD)
SD = NS * 2;

% Initialize matrices for stiffness, forces, and displacements
K = zeros(SD);
xf = zeros(CS, 1);
yf = zeros(CS, 1);
L = zeros(CS, 1);
KK = zeros(SD);
F = zeros(SD, 1);

% Loop through each bar element
for e = 1:CS
    i = EN(e, :);
    
    % Define element degrees of freedom
    ESD = [(2 * i(1) - 1), (2 * i(1)), (2 * i(2) - 1), (i(2) * 2)];
    
    % Calculate the x and y components of the bar
    xf(e, 1) = x(i(2)) - x(i(1));
    yf(e, 1) = y(i(2)) - y(i(1));
    
    % Calculate the length of the bar
    L(e, 1) = sqrt(xf(e, 1) .* xf(e, 1) + yf(e, 1) .* yf(e, 1));
    
    % Calculate the element stiffness matrix
    k = (EI / (L(e, 1))^3) * [12, 6 * L(e, 1), -12, 6 * L(e, 1);
                              6 * L(e, 1), 4 * (L(e, 1))^2, -6 * L(e, 1), 2 * (L(e, 1))^2;
                              -12, -6 * L(e, 1), 12, -6 * L(e, 1);
                              6 * L(e, 1), 2 * (L(e, 1))^2, 12, -6 * L(e, 1)];
    
    % Assemble the global stiffness matrix
    K(ESD, ESD) = K(ESD, ESD) + k;
    KK(ESD, ESD) = KK(ESD, ESD) + k;
    
    % Input the distributed load on the element
    fprintf('Enter the distributed load on element %d=', e)
    q = input('');
    
    % Calculate the element force vector
    f = [-q * L(e, 1) / 2; -q * (L(e, 1))^2 / 12; -q * L(e, 1) / 2; q * (L(e, 1))^2 / 12];
    
    % Assemble the global force vector
    F(ESD, 1) = F(ESD, 1) + f;
end

% Remove rows and columns corresponding to constrained degrees of freedom
z = [1 5];
K(:, z) = [];
K(z, :) = [];
F(z, :) = [];

% Solve for displacements using inverse of stiffness matrix
U = inv(K) * F;

% Calculate the total number of degrees of freedom (again)
ss = length(U);

% Initialize arrays for global displacements and node coordinates
Ug = zeros(1, SD);
NS = zeros(1, SD);

% Assign node numbers to degrees of freedom
for i = 1:SD
    % Node Degrees of Freedom
    NS(1, i) = i;
end

% Remove rows corresponding to constrained degrees of freedom
NS(:, z) = [];

% Apply the calculated displacements to node coordinates
Ug(NS) = U;

% Calculate the deformed coordinates of nodes
for i = 1:NS
    xL(i, 1) = x(i, 1) + Ug(1, 2 * i - 1);
    yL(i, 1) = y(i, 1) + Ug(1, 2 * i);
end

% Calculate the reaction forces
R = KK(z, :) * Ug';

% Display the stiffness matrix, forces, and displacements
display(K)
display(F)
display(Ug)
This code appears to perform structural analysis on a system of bars, calculating displacements, forces, and reactions based on input parameters and boundary conditions. Please note that the code assumes some user input for various parameters.





