%% Title: Model of free/biofilm cells with two types of predator
% Author: S williams

%% Early time fig2 (bar 2/3 optimised, bar 4 predictive).

t0 = 0; t1 = 72; % Start and endpoints for the solver.

pops = zeros(4,2); % Store for the outputs.

N0 = zeros(7,1); % Store for intial conditions.
N0(1) = 200; % Carbon concentration in media.
N0(2) = 1.5; % Free cells.
N0(3) = 0; % No intial biofilm.
N0(6) = 0; % Initially no waste.

T0 = 0.1919; % See e_tuningT0 for origin of value.
S0 = 0.1263; % See e_tuningS0 for origin of value.

% Rerun the ODE system with the control.
N0(4) = 0; N0(5) = 0;
[pops(1,1),pops(1,2)] = solvePop(N0,t0,t1);

% Rerun the ODE system with bf predators.

N0(4) = 0; N0(5) = T0;
[pops(2,1),pops(2,2)] = solvePop(N0,t0,t1);

% Rerun the ODE system with free predators.
N0(4) = S0; N0(5) = 0;
[pops(3,1),pops(3,2)] = solvePop(N0,t0,t1);

% Rerun the ODE system with both predators.
N0(4) = 24*S0/27; N0(5) = 2127*T0/3160; % Values scaled by experiments.
[pops(4,1),pops(4,2)] = solvePop(N0,t0,t1);

bar(pops,'stacked')

%% Late time fig2 (bar 2/3/4 predictive).

t0 = 0; t1 = 96; % Start and endpoints for the solver.

pops = zeros(4,2); % Store for the outputs.

N0 = zeros(7,1); % Store for intial conditions.
N0(1) = 200; % Carbon concentration in media.
N0(2) = 1.5; % Free cells.
N0(3) = 5; % No intial biofilm.
N0(6) = 0; % Initially no waste.

T0 = 0.1919; % See e_tuningT0 for origin of value. Scaled by experiments.
S0 = 0.1263; % See e_tuningS0 for origin of value. Scaled by experiments.

% Rerun the ODE system with the control.
N0(4) = 0; N0(5) = 0;
[pops(1,1),pops(1,2)] = solvePop(N0,t0,t1);

% Rerun the ODE system with bf predators.

N0(4) = 0; N0(5) = 1847*T0/3160;
[pops(2,1),pops(2,2)] = solvePop(N0,t0,t1);

% Rerun the ODE system with free predators.
N0(4) = 40*S0/27; N0(5) = 0;
[pops(3,1),pops(3,2)] = solvePop(N0,t0,t1);

% Rerun the ODE system with both predators.
N0(4) = 41*0.1263/27; N0(5) = 1853*T0/3160;
[pops(4,1),pops(4,2)] = solvePop(N0,t0,t1);

bar(pops,'stacked')

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

    [~,N] = ode45(@(t,y)... % Variables.
        odefun(t,y), ... % Function (declared at bottom of the script).
        [t0 t1], ... % Time range for which to solve.
        N0, ...  % Initial conditions.
        opts); % Extra conditions for the solver.
    F = biovolume(N(end,2),'v');
    B = biovolume(N(end,3),'s');

end