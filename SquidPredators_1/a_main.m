%% Title: Model of free/biofilm cells with various predator types
% Author: S williams

%% Preamble

clear all
close all

%% Constants

model = struct; % Store for all of the parameters in the model.

% Default parameters.
model.eb=0.2; % Carbon use efficiency.
model.rpl=0.21; model.rbf=0.007; % Cell growth rates.
model.Hpl=1; model.Hbf=1; model.Hpa=1; model.Ham=1; % Cell carbon useage hill function.
model.SoV=8/3; % Well growth surface area:volume ratio. % /3 in model §, /2.5 in exp't §.
model.gpa=0.12; model.gam=0.09; % (DEFAULT: gpa=0.12;gam=0.09;)predator grazing rates.
model.epa=0.5; model.eam=0.33; % predator growth efficiency.
model.chibf=0.005; model.chimax=0.05; model.chimin=0.0005; model.a=0.01; % Migration rates.
model.rcy1 = 0; model.rcy2 = 0; % Efficiency of recycling of predator waste.

%% Initial conditions

t0 = 0; t1 = 72; % Start and endpoints for the solver.

T0 = 0.1; % See e_tuningT0 for origin of value.
S0 = 0.1; % See e_tuningS0 for origin of value.

N0 = zeros(9,1); % Store for intial conditions.
N0(1) = 1; % Carbon concentration in media.
N0(2) = 0.1; % Free cells.
N0(3) = 0; % No intial biofilm.
N0(4) = 0; % Initial predator S.
N0(5) = 0; % Initial predator T.

%% Main loop

[N] = solvePop(N0,t0,t1,model);

%% Plotting tools

%% Functions

% ODE terms.
% y(1) - carbon conc.
% y(2) - Free cells.
% y(3) - Biofil cells.
% y(4) - Predator 1, free cell type.
% y(5) - Predator 2, biofilm type.
% y(6) - waste F.
% y(7) - waste B.
% y(8) - waste S.
% y(9) - waste T.

% Model.
% This struct contains all the parameters.

% Population dynamics ODE.
function dNdt = odefun(~,y,model)

% Default parameters, from paper.
eb=model.eb; % Carbon useage parameters.
rpl=model.rpl; rbf=model.rbf; % Cell growth rates.
Hpl=model.Hpl; Hbf=model.Hbf; Hpa=model.Hpa; Ham=model.Ham; % Cell carbon useage hill function.
SoV=model.SoV; % Well growth surface area:volume ratio. % /3 in model §, /2.5 in exp't §.
gpa=model.gpa;gam=model.gam; % (DEFAULT: gpa=0.12;gam=0.09;)predator grazing rates.
epa=model.epa;eam=model.eam; % predator growth efficiency.
chibf=model.chibf; chimax=model.chimax; chimin=model.chimin; 
a=model.a; chipl = (a*chimax+chimin*y(3))/(a+y(3)); % Migration rates.
muF=model.muF; muB=model.muB; muS=model.muS; muT=model.muT; % Species mortality.
rcy1 = model.rcy1; rcy2 = model.rcy2; % Predator recycling;

% Store for ODE values.
dNdt = zeros(9,1);

% Standard terms
t1 = rpl*(y(1)/(y(1)+Hpl))*y(2); % Growth of free cells.
t2 = rbf*SoV*(y(1)/(y(1)+Hbf))*y(3); % Growth of biofilm.
t3 = (chipl*SoV)*y(2); % free -> biofilm migration.
t4 = chibf*(SoV*y(3)); % biofilm -> free migration.
t5 = gpa*y(4)*y(2)/(y(2)+Hpa); % grazing of free.
t6 = gam*y(5)*y(3)/(y(3)+Ham); % grazing of biofilm.

% Waste stores
w1 = (1/eb - 1)*(t1); % Waste prey free.
w2 = (1/eb - 1)*(t2); % Waste prey biofilm.
w3 = (1-epa)*t5; % Waste predator S.
w4 = (1-eam)*t6*SoV; % Waste predator T.

% Mortality stores
m1 = muF*y(2);
m2 = muB*y(3);
m3 = muS*y(4);
m4 = muT*y(5);

% Calculate ODE term values
dNdt(1) = -(1/eb)*( t1 + t2 ) + rcy1*w3 + rcy2*w4;
dNdt(2) = t1 - t5 - t3 + t4 - m1;
dNdt(3) = t2/SoV - t6 + t3/SoV - t4/SoV - m2;
dNdt(4) = epa*t5 - m3;
dNdt(5) = eam*t6 - m4;

% Calculate waste ODE
dNdt(6) = w1 + m1; 
dNdt(7) = w2 + m2; 
dNdt(8) = w3 + m3;
dNdt(9) = w4 + m4;

end

function [N] = solvePop(N0,t0,t1,model)

    opts = odeset('NonNegative',1:9); % Ensure all solutions are non-negative.

    [~,N] = ode89(@(t,y)... % Variables.
        odefun(t,y,model), ... % Function (declared at bottom of the script).
        [t0 t1], ... % Time range for which to solve.
        N0, ...  % Initial conditions.
        opts); % Extra conditions for the solver.

end