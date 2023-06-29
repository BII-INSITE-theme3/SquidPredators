%% Title: Model of free/biofilm cells with two types of predator
% Author: S williams

% This code is a rudimentary sensistivity analysis.

% This can be run in two ways: 
% (Option 1), examines the impact on the total prey biomass.
% (Option 2), examples the ratio of free to biofilm prey.
% Each has a corresponding set of plotting tools to visualise.
% Array "fields" contains the parameter types plotted on the x-axis.
% By modifying N0(4) and N0(5) you can change predator populations.

%% Initialize model (base) parameters

model = struct; % Store for all of the parameters in the model.

% Default parameters, from paper.
model.eb=0.2; % Carbon useage parameters.
model.rpl=0.21; model.rbf=0.007; % Cell growth rates.
model.Hpl=1; model.Hbf=1; model.Hpa=1; model.Ham=0.1; % Cell carbon useage hill function.
model.gpa=0.12; model.gam=0.09; % (DEFAULT: gpa=0.12;gam=0.09;)predator grazing rates.
model.epa=0.5; model.eam=0.33; % predator growth efficiency.
model.chibf=0.005; model.chimax=0.05; model.chimin=0.0005; model.a=0.01; % Migration rates.
%model.rcy1=0; model.rcy2=0; % Recycling of predator waste.
model.SoV=8/3; % Well growth surface area:volume ratio. % /3 in model §, /2.5 in exp't §.

fields = fieldnames(model);

%% Constants

% Runtime constants
t0 = 0; t1 = 96; % Start and endpoints for the solver.

% Predator concentrations
T0 = 0.1919; % See e_tuningT0 for origin of value.
S0 = 0.1263; % See e_tuningS0 for origin of value.

% Initial conditions
N0 = zeros(9,1); % Store for intial conditions.
N0(1) = 200; % Carbon concentration in media.
N0(2) = 1.5; % Free cells.
N0(3) = 0; % No intial biofilm.
N0(4) = S0; % Initial predator S.
N0(5) = T0; % Initial predator T.

%% Main loop 

Nruns1 = 15; % Number of loops to run
Nruns2 = 3; % Number of loops to run
fmod = 0.1; % Percentage modifier for parameter
modifier = linspace(1-fmod,1+fmod,Nruns2); % Get modified parameters as normalised values
outputs1 = zeros(Nruns1,Nruns2); % Store for the outputs
outputs2 = zeros(Nruns1,Nruns2); % Store for the outputs

for run1 = 1:Nruns1
    for run2 = 1:Nruns2

        % Set the modified parameters
        modeltemp = model; % Store the "base" model
        modeltemp.(fields{run1}) = modeltemp.(fields{run1})*modifier(run2); % Modify the parameter of interest

        % Solve the ODE system
        [~,~,N] = solvePop(N0,t0,t1,modeltemp); % Solve for the modified system

        % Solution processing
        outputs1(run1,run2) = N(end,2)+N(end,3); % (option 1) Total prey biomass
        outputs2(run1,run2) = N(end,2)/N(end,3);  % (option 2) Relative prey composition

    end
end

%% Plotting tools 1 - total prey biomass

plot(outputs1(:,1)./outputs1(:,2),'*','lineStyle','none')
hold on
plot(outputs1(:,3)./outputs1(:,2),'*','lineStyle','none')

xlabel('parameter');
ylabel('outcome relative to contol');
lab1 = ['-' num2str(100*fmod) '%'];
lab2 = ['-' num2str(100*fmod) '%'];
legend(lab1,lab2);

%% Plotting tools 2 - relative prey composition

hold on
plot(outputs2(:,1),'*','lineStyle','none')
plot(outputs2(:,2),'*','lineStyle','none')
plot(outputs2(:,3),'*','lineStyle','none')

xlabel('parameter');
ylabel('F/B at 72hr');
lab1 = ['-' num2str(100*fmod) '%'];
lab2 = 'control';
lab3 = ['-' num2str(100*fmod) '%'];
legend(lab1,lab2,lab3);

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
%rcy1 = model.rcy1; rcy2 = model.rcy2; % Recycling efficiencies.

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
dNdt(1) = -(1/eb)*( t1 + t2 ); %+ rcy1*w3 + rcy2*w4;
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