%% Title: Model of free/biofilm cells with two types of predator
% Author: S williams

% This code doesn't really work or do anything interesting...

%% Conditions (Modify this part for different run-types)

runtype = "early"; % or "late" % type of run.
predators = "N"; % Enter S for free, T for bf, TS for both, otherwise none will be present.

%% Constants
% System
t0 = 0; % Start time
if (runtype == "early")
    t1 = 72; 
else
    t1 = 96;
end
N = 1000; % Number of samples
timeSteps = linspace(t0,t1,N)'; % Sample points of ODE solver.

%% Initial conditions
N0 = zeros(6,1); % Store for intial conditions.
N0(1) = 200; % Carbon concentration in media.
N0(2) = 1.5; % Free cells.
if (runtype == "early")
    N0(3) = 0; % No biofilm present.
else
    N0(3) = 5; % Established biofilm.
end
if (predators == "TS") % Both predators
    N0(4) = 0.5; % 0.5 if free pred % Free predators.
    N0(5) = 1.9; % 1.9 if bf pred % Biofilm predators.
elseif (predators == "T") % Just biofilm predator
    N0(4) = 0; % 0.5 if free pred % Free predators.
    N0(5) = 1.9; % 1.9 if bf pred % Biofilm predators.
elseif (predators == "S") % Just free predator
    N0(4) = 0.5; % 0.5 if free pred % Free predators.
    N0(5) = 0; % 1.9 if bf pred % Biofilm predators.
else % No predators
    N0(4) = 0; % 0.5 if free pred % Free predators.
    N0(5) = 0; % 1.9 if bf pred % Biofilm predators.
end

%% Main

opts = odeset('NonNegative',1:6); % Ensure all solutions are non-negative.

[t,N] = ode89(@(t,y)... % Variables.
              odefun(t,y), ... % Function (declared at bottom).
              timeSteps, ... % Time steps to be sampled by solver.
              N0, ...
              opts); % Initial conditions.

%% Plotting tools

plot(timeSteps,N(:,1),'lineWidth',5); % Carbon conc.
hold on;
plot(timeSteps,N(:,2),'lineWidth',5); % free prey.
plot(timeSteps,N(:,3),'lineWidth',5); % biofilm prey.
plot(timeSteps,N(:,4),'lineWidth',5); % free pre.
plot(timeSteps,N(:,5),'lineWidth',5); % biofilm pred.

%% saver

% str = string([N(1,2),N(1,3),N(1,4),N(1,5),N(end,2),N(end,3),N(end,4),N(end,5)]);
% writelines(' ','output.txt',WriteMode="append")
% writelines(str(:),'output.txt',WriteMode="append"); % Display the final populations.

%% ODE functions

% ODE terms.
% y(1) - carbon conc.
% y(2) - Free cells.
% y(3) - Biofil cells.
% y(4) - Predator 1, free cell type.
% y(5) - Predator 2, biofilm type.

% Population dynamics ODE.
function dNdt = odefun(~,y)
    % Default parameters, from paper.
    eb=0.2; % Carbon useage parameters.
    rpl=0.21; rbf=0.007; % Cell growth rates.
    Hpl=1; Hbf=1; Hpa=1; Ham=0.1; % Cell carbon useage hill function.
    SoV=8/2.5; % Well growth surface area:volume ratio. % /3 in model §, /2.5 in exp't §.
    gpa=0.12;gam=0.09; % predator grazing rates.
    epa=0.5;eam=0.33; % predator growth rates.
    chibf=0.005; chimax=0.05; chimin=0.0005; 
    a=0.01; chipl = (a*chimax+chimin*y(3))/(a+y(3)); 
    % Store for ODE values
    dNdt = zeros(5,1);
    % Terms
    %t1 = rpl*(y(1)/(y(1)+Hpl))*y(2); % Growth of free cells.
    t1 = rpl*y(2)*sign(y(1)); % Linearised Growth of free cells.
    %t2 = rbf*SoV*(y(1)/(y(1)+Hbf))*y(3); % Growth of biofilm.
    t2 = rbf*SoV*y(3)*sign(y(1)); % Linearised Growth of biofilm.
    t3 = chipl*SoV*y(2); % free -> biofilm migration.
    t4 = chibf*SoV*y(3); % biofilm -> free migration.
    %t5 = gpa*y(4)*y(2)/(y(2)+Hpa); % grazing of free.
    t5 = gpa*y(4); % grazing of free.
    %t6 = gam*y(5)*y(3)/(y(3)+Ham); % grazing of biofilm.
    t6 = gam*y(5); % grazing of biofilm.
    % Calculate ODE term values
    dNdt(1) = -(1/eb)*( t1 + t2 );
    dNdt(2) = t1 - t5 - t3 + t4;
    dNdt(3) = t2/SoV - t6 + t3/SoV - t4/SoV;
    dNdt(4) = epa*t5;
    dNdt(5) = eam*t6;
    dNdt(6) = (1/eb - 1)*(t1+t2) + (1-epa)*t5 + (1-eam)*t6; % Waste.
end