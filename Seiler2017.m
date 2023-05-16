%% Title: Basic model of equilibirum reached by free and biofilm vibrio 
% Author: S williams

%% Constants
% System
t0 = 0; % Start time
t1 = 10; % End time
N = 100; % Number of samples
timeSteps = linspace(t0,t1,N)'; % Sample points of ODE solver.
% Inidivuals

% pair-wise

%% Initial conditions
N0 = zeros(5,1);
N0(1) = 1; % Carbon access.
N0(2) = 1; % Free cells.
N0(3) = 0; % Biofilm.
N0(4) = 0.1; % Free predators.
N0(5) = 0.1; % Biofilm predators.

%% Main

[t,N] = ode45(@(t,y)... % Variables.
              odefun(t,y), ... % Function (declared at bottom).
              timeSteps, ... % Time steps to be sampled by solver.
              N0); % Initial conditions.

%% Plotting tools

plot(timeSteps,N(:,1),'lineWidth',5);
hold on;
plot(timeSteps,N(:,2),'lineWidth',5);
plot(timeSteps,N(:,3),'lineWidth',5);
plot(timeSteps,N(:,4),'lineWidth',5);
plot(timeSteps,N(:,5),'lineWidth',5);

%% ODE functions

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
    SoV=8/3; % Well growth surface area:volume ratio.
    gpa=0.12;gam=0.09; % predator grazing rates.
    epa=0.5;eam=0.33; % predator growth rates.
    chibf=0.005; chimax=0.05; chimin=0.0005; a=0.01; chipl = (a*chimax+chimin*y(3))/(a+y(3)); % 

    dNdt = zeros(5,1);

    t1 = rpl*(y(1)/(y(1)+Hpl))*y(2); % Growth of free cells.
    t2 = rbf*SoV*(y(1)/(y(1)+Hbf))*y(3); % Growth of biofilm.
    t3 = chipl*SoV*y(2); % free -> biofilm migration.
    t4 = chibf*SoV*y(3); % biofilm -> free migration.
    t5 = gpa*y(4)*y(2)/(y(2)+Hpa); % grazing of free.
    t6 = gam*y(5)*y(3)/(y(3)+Ham); % grazing of biofilm.

    dNdt(1) = -(1/eb)*( t1 + t2 );
    dNdt(2) = t1 - t5 - t3 + t4;
    dNdt(3) = t2 - t6 + t3 - t4;
    dNdt(4) = epa*t5;
    dNdt(5) = eam*t6;

end