% This is the master file that sets parameters and calls the other files.

% Lutscher, Popovic, Shaw, Theoretical Ecology
% K types
% total population  = 2 * female population
% No need to distinguish between male and female because of simple mating structure
% track dispersal types only
% With dispersal-induced mortality and initial condition adjustment
% Gives individual speeds, speed without mutation and speed with mutation
%clear
clear

% parameters
b = 4;  % offspring per female
n = 5; % maximum local density
tau = 1; % mate search rate
dmin = 0.1; % minimal dispersal parameter
dmax = 0.5; %maximal dispersal parameter
mu = 0.03;  % mutation probability for simulation run
mui = mu;  % mutation probability for  initial conditions
m = 0.3; % mortality factor: survival is exp(-d(i)*m) 
K = 2; % number of types
P = [b,n,tau,dmin,dmax,mu,mui,m,K];

% use the following to calculate the speed of a single type
% There are K equally spaced types between dmin and dmax

speedsingle = SingleSpeciesAlleeSpreadKFUNCTION(P)

% use the following to calculate the joint speed of
% K equally spaced types between dmin and dmax with mutation parameter mu

speedwith = AllSpreadKTypesMortalityFUNCTION(P)

% use the following to calculate the joint speed of
% K equally spaced types between dmin and dmax with mutation set to zero

P(6) = 0; % set mutation to zero in the run
speedwithout = AllSpreadKTypesMortalityFUNCTION(P)