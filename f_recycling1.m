%% Title: Model of free/biofilm cells with two types of predator
% Author: S williams

%% Initialize parameters

model = struct; % Store for all of the parameters in the model.

% Default parameters, from paper.
model.eb=0.2; % Carbon useage parameters.
model.rpl=0.21; model.rbf=0.007; % Cell growth rates.
model.Hpl=1; model.Hbf=1; model.Hpa=1; model.Ham=0.1; % Cell carbon useage hill function.
model.SoV=8/3; % Well growth surface area:volume ratio. % /3 in model §, /2.5 in exp't §.
model.gpa=0.12; model.gam=0.09; % (DEFAULT: gpa=0.12;gam=0.09;)predator grazing rates.
model.epa=0.5; model.eam=0.33; % predator growth efficiency.
model.chibf=0.005; model.chimax=0.05; model.chimin=0.0005; model.a=0.01; % Migration rates.

%% Main loop 

Nruns = 100; % Number of loops to run
cond = linspace(0,1,Nruns); % Conditions for the main loop
outputs = zeros(Nruns,1);

for run = 1:Nruns

    % New parameters, normalised values.
    model.rcy1=0; model.rcy2=cond(run);

    %% Early time fig2 (bar 2/3 optimised, bar 4 predictive).

    t0 = 0; t1 = 72; % Start and endpoints for the solver.

    T0 = 0.1919; % See e_tuningT0 for origin of value.
    S0 = 0.1263; % See e_tuningS0 for origin of value.

    N0 = zeros(9,1); % Store for intial conditions.
    N0(1) = 200; % Carbon concentration in media.
    N0(2) = 1.5; % Free cells.
    N0(3) = 0; % No intial biofilm.
    N0(4) = 0; % Initial predator S.
    N0(5) = 0; % Initial predator T.

    %% Solve the ODE.

    [~,~,N] = solvePop(N0,t0,t1,model);

    %% Solution processing
    outputs(run) = N(end,2)+N(end,3)*model.SoV; % Total prey biomass

end % End of main loop

%% Plotting tools 1

plot(outputs)
hold on

%% Plotting tools 2

% plot(N(:,1)) % Carbon conc.

% plot(N(:,2)) % Free prey.
% plot(N(:,3)) % Biofilm prey.
% plot(N(:,4)) % Free predator.
% plot(N(:,5)) % biofilm predator.

% plot(N(:,6)) % Free prey waste.
% plot(N(:,7)) % Biofilm prey waste.
% plot(N(:,8)) % Free predator waste.
% plot(N(:,9)) % Biofilm predator waste.

%% ODE functions

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

% Recycling efficiencies.
rcy1 = model.rcy1; rcy2 = model.rcy2;

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

% Calculate ODE term values
dNdt(1) = -(1/eb)*( t1 + t2 ) + rcy1*w3 + rcy2*w4;
dNdt(2) = t1 - t5 - t3 + t4;
dNdt(3) = t2/SoV - t6 + t3/SoV - t4/SoV;
dNdt(4) = epa*t5;
dNdt(5) = eam*t6;

% Calculate waste ODE
dNdt(6) = w1 ; 
dNdt(7) = w2; 
dNdt(8) = w3;
dNdt(9) = w4;

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

function [F, B, N] = solvePop(N0,t0,t1,model)

    opts = odeset('NonNegative',1:9); % Ensure all solutions are non-negative.

    [~,N] = ode89(@(t,y)... % Variables.
        odefun(t,y,model), ... % Function (declared at bottom of the script).
        [t0 t1], ... % Time range for which to solve.
        N0, ...  % Initial conditions.
        opts); % Extra conditions for the solver.

    F = biovolume(N(end,2),'v');
    B = biovolume(N(end,3),'s');

end