%% Title: Model of free/biofilm cells with two types of predator
% Author: S williams

% This script is a modification of c_Seiler_fig2.m
% We focus entirely on modifying the S population here.
% In Seiler_fig2 given the 'email' initial conditions the S population
% quickly finish off the F (overcompensate). In the 'non-email' initial
% conditions they underperform and the F are largely unaffected by their
% presence.
% In this code the aim is to interpolate between these two cases to find initial conditions
% which gives a match to the values which aregiven for the early case of
% the paper.
% Following which we will then modify to the 'late' biofilm case and see
% how these values perform (when scaled by the experimental case).

%% Conditions (Modify this part for different run-types)

runs = 100; % Number of initial conditions to sample

residuals = zeros(runs,2); % Store for the metrics of the outcomes.
popfin = zeros(runs,2);
IVPs = linspace(0,0.1,runs); % Initial conditions to be run.

% System
t0 = 0; % Start time
t1 = 72; % End time

N0 = zeros(7,1); % Store for intial conditions.
N0(1) = 200; % Carbon concentration in media.
N0(2) = 1.5; % Free cells.
N0(3) = 0; % No intial biofilm.
N0(6) = 0; % Initially no waste.

for run = 1:runs % Loop on the different initial conditions.

    %% IVP input

    N0(4) = IVPs(run); % Free predator corresponding to the S0 of the run.
    N0(5) = 0; % No biofilm predators.

    %% Main ODE solver step

    [F,B] = solvePop(N0,t0,t1);

    %% Output processing
    popfin(run,1) = F; % Get -> (freecells - expected(freecells))
    popfin(run,2) = B; % Get -> (biofilmcells - expected(biofilmcells))
    residuals(run,1) = F - 1.1162; % Get -> (freecells - expected(freecells))
    residuals(run,2) = B - 0.257; % Get -> (biofilmcells - expected(biofilmcells))

end % End of loop on runs.

%% Process the residuals found in all the runs.

quality = sqrt(residuals(:,1).^2 + residuals(:,2).^2); % Quantify the fit quality.
%quality = abs(residuals(:,1)); % Quantify the fit quality.
%quality = abs(residuals(:,2)); % Quantify the fit quality.
[~,I] = min(quality); % Get the "best" fit by minimising the error.

% Rerun the ODE system with the control.
N0(4) = 0;
[F,B] = solvePop(N0,t0,t1);
a = [F B]; % top figs, total prey composition.

% Rerun the ODE system with the best conditions.
N0(4) = IVPs(I);
[F,B] = solvePop(N0,t0,t1);
b = [F B]; % top figs, total prey composition.

bar([a;b],'stacked')

%% Plot varrying the initial predator population size

bar(IVPs,[popfin(:,1) popfin(:,2)], 'stacked')

%% FUNCTIONS 
%% ODE functions

% ODE terms.
% y(1) - carbon conc.
% y(2) - Free cells.
% y(3) - Biofil cells.
% y(4) - Predator 1, free cell type.
% y(5) - Predator 2, biofilm type.
% y(6) - waste.

% Population dynamics ODE.
function dNdt = odefun(~,y)

% Default parameters, from paper.
eb=0.2; % Carbon useage parameters.
rpl=0.21; rbf=0.007; % Cell growth rates.
Hpl=1; Hbf=1; Hpa=1; Ham=0.1; % Cell carbon useage hill function.
SoV=8/3; % Well growth surface area:volume ratio. % /3 in model §, /2.5 in exp't §.
gpa=0.12;gam=0.09; % (DEFAULT: gpa=0.12;gam=0.09;)predator grazing rates.
epa=0.5;eam=0.33; % predator growth efficiency.
chibf=0.005; chimax=0.05; chimin=0.0005;
a=0.01; chipl = (a*chimax+chimin*y(3))/(a+y(3)); % Migration rates.
% Store for ODE values
dNdt = zeros(5,1);
% Terms
t1 = rpl*(y(1)/(y(1)+Hpl))*y(2); % Growth of free cells.
t2 = rbf*SoV*(y(1)/(y(1)+Hbf))*y(3); % Growth of biofilm.
t3 = (chipl*SoV)*y(2); % free -> biofilm migration.
t4 = chibf*(SoV*y(3)); % biofilm -> free migration.
t5 = gpa*y(4)*y(2)/(y(2)+Hpa); % grazing of free.
t6 = gam*y(5)*y(3)/(y(3)+Ham); % grazing of biofilm.
% Calculate ODE term values
dNdt(1) = -(1/eb)*( t1 + t2 );
dNdt(2) = t1 - t5 - t3 + t4;
dNdt(3) = t2/SoV - t6 + t3/SoV - t4/SoV;
dNdt(4) = epa*t5;
dNdt(5) = eam*t6;
dNdt(6) = (1/eb - 1)*(t1+t2) ; % Waste prey.
dNdt(7) = (1-epa)*t5 + (1-eam)*t6*SoV; % Waste predator.

end

%% method 1 via eqn 1 romanova 2009.

% Function to caclulate the biovolume from the carbon conc.

function V = biovolume(in,type) % 'in' is ugC either (v) per ml or (s) per cm^2.

if type == 'v'
    vol = 3; % volume.
    in = in * vol; % total ug carbon.
end
if type == 's'
    surf = 8; % surface area.
    in = in * surf; % total ug carbon.
end

fgpug = 10^9; % fg per ug.
fgC = in * fgpug; % Total grams carbon in fg.
fgCpcell = 10^9*1.5/(6*10^6);
v0 = (fgCpcell/133.754)^(1/0.436); % um^3 per cell.
N = fgC/fgCpcell; % Number of cells.
V = N*v0 * 10^(-9); % get final volume in mm^3.

end

%% Solve ODE function (refactor)

function [F, B] = solvePop(N0,t0,t1)

    opts = odeset('NonNegative',1:6); % Ensure all solutions are non-negative.

    [~,N] = ode89(@(t,y)... % Variables.
        odefun(t,y), ... % Function (declared at bottom of the script).
        [t0 t1], ... % Time range for which to solve.
        N0, ...  % Initial conditions.
        opts); % Extra conditions for the solver.
    F = biovolume(N(end,2),'v');
    B = biovolume(N(end,3),'s');

end
