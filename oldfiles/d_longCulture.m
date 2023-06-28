%% Title: Model of free/biofilm cells with two types of predator
% Author: S williams

% This version of the code is intended to allow for logistic-growth of the
% biofilm and the free swimming bacteria. 
% I have also included cell mortality for each of the different species.

% The equations for this section are outlined in the overleaf.
% See section: MATHEMATICAL MODELLING > Steady state model

%% Constants
% System
t0 = 0; % Start time.
t1 = 1000000; % End time.
%N = 10; % Number of samples.
timeSteps = linspace(t0,t1)'; % Sample points of ODE solver.
%timeSteps = linspace(t0,t1,N)';

%% Initial conditions
N0 = zeros(4,1); % Store for intial conditions.
N0(1) = 0.1; % Free cells.
N0(2) = 10; % Biofilm.

%predators = "N"; % No predators.
%predators = "TS"; % Both predators.
% predators = "T"; % Only B predators.
 predators = "S"; % Only F predators.

if (predators == "TS") % Both predators
    N0(3) = 0.01; % 0.5 if free pred % Free predators.
    N0(4) = 0.01; % 1.9 if bf pred % Biofilm predators.
elseif (predators == "T") % Just biofilm predator
    N0(3) = 0; % 0.5 if free pred % Free predators.
    N0(4) = 0.01; % 1.9 if bf pred % Biofilm predators.
elseif (predators == "S") % Just free predator
    N0(3) = 0.01; % 0.5 if free pred % Free predators.
    N0(4) = 0; % 1.9 if bf pred % Biofilm predators.
else % No predators
    N0(3) = 0; % 0.5 if free pred % Free predators.
    N0(4) = 0; % 1.9 if bf pred % Biofilm predators.
end

%% Main

[t,N] = ode89(@(t,y)... % Variables.
              odefun(t,y), ... % Function (declared at bottom).
              timeSteps, ... % Time steps to be sampled by solver.
              N0); % Initial conditions.

%% Plotting tools

tst = 1;
tend = 100;

plot(timeSteps(tst:tend),N(tst:tend,1),'lineWidth',5); % Carbon conc.
hold on;
plot(timeSteps(tst:tend),N(tst:tend,2),'lineWidth',5); % free prey.
plot(timeSteps(tst:tend),N(tst:tend,3),'lineWidth',5); % biofilm prey.
plot(timeSteps(tst:tend),N(tst:tend,4),'lineWidth',5); % free pre.

%% ODE functions

% ODE terms.
% y(1) - Free cells.
% y(2) - Biofil cells.
% y(3) - Predator 1, free cell type.
% y(4) - Predator 2, biofilm type.

% Population dynamics ODE.
function dNdt = odefun(~,y)

    rF = 0.21; rB = 0.007; % Growth rates.
    Kf = 10; Kb = 10; % Carrying capacities.
    chiMax=0.05; chiMin=0.0005; % Range of BF growth.
    a=0.01; % Biofilm dependence strength.
    chiFB = (a*chiMax + y(2)*chiMin)/(a+y(2)); % F to B rate.
    chiBF = 0.005; % B to F rate.
    HF = 1;HB = 0.1; % Predation half saturations
    gS = 0.12;gT = 0.09; % Grazing rates.
    epsS = 0.5;epsT = 0.33; % Predator Growth efficiency.
    muF = 0.1; muB = 0.01; muS = 0.05; muT = 0.1;
    
    Fgrow = rF * y(1) * (1 - y(1)/Kf); % F growth.
    Bgrow = rB * y(2) * (1 - y(2)/Kb); % B growth.
    FtoB = chiFB * y(1); % Migration F to B.
    BtoF = chiBF * y(2); % Migration B to F.
    Spred = gS * y(1)/(y(1)+HF) * y(3); % Predation of F by S.
    Tpred = gT * y(2)/(y(2)+HB) * y(4); % Predation of B by T.
    Fmort = muF * y(1); % F intrinsic death rate.
    Bmort = muB * y(2); % B intrinsic death rate.
    Smort = muS * y(3); % S intrinsic death rate.
    Tmort = muT * y(4); % T intrinsic death rate.

    dNdt = zeros(4,1);

    dNdt(1) = Fgrow - FtoB + BtoF - Spred - Fmort; % Changes in F.
    dNdt(2) = Bgrow + FtoB - BtoF - Tpred - Bmort; % Changes in B.
    dNdt(3) = epsS*Spred - Smort; % Changes in S.
    dNdt(4) = epsT*Tpred - Tmort; % Changes in T.

end