%% Title: Model of free/biofilm cells with two types of predator
% Author: S williams

%% Conditions (Modify this part for different run-types)

runtype = ["early","late"]; % type of run.
nrun = 2;
predators = ["N", "T","S","TS"]; % Enter S for free, T for bf, TS for both, otherwise none will be present.
npred = 4;
email = "y";

%% loop on conditions

for run = 1%nrun

    F = zeros(4,1);
    B = zeros(4,1);

    for preds = 1:npred

        runtype1 = runtype(run);
        predators1 = predators(preds);

        % System
        t0 = 0; % Start time
        if (runtype1 == "early")
            t1 = 72;
        else
            t1 = 96;
        end
        %N = 1000; % Number of samples
        timeSteps = linspace(t0,t1)'; % Sample points of ODE solver.

        N0 = zeros(7,1); % Store for intial conditions.
        N0(1) = 200; % Carbon concentration in media.
        N0(2) = 1.5; % Free cells.

        %% IVP based on exp measurements.
        if email == "n"
            if (runtype1 == "early")
                N0(3) = 0; % No biofilm present.
                if (predators1 == "TS") % Both predators
                    N0(4) = predMassF(24);   N0(5) = predMassBF(2127);
                elseif (predators1 == "T") % Just biofilm predator
                    N0(4) = 0;               N0(5) = predMassBF(3160);
                elseif (predators1 == "S") % Just free predator
                    N0(4) = predMassF(27);   N0(5) = 0;
                else % No predators
                    N0(4) = 0;               N0(5) = 0;
                end
            else % "Late" biofilm
                N0(3) = 5; % Established biofilm.
                if (predators1 == "TS") % Both predators
                    N0(4) = predMassF(41);   N0(5) = predMassBF(1853);
                elseif (predators1 == "T") % Just biofilm predator
                    N0(4) = 0;               N0(5) = predMassBF(1847);
                elseif (predators1 == "S") % Just free predator
                    N0(4) = predMassF(40);   N0(5) = 0;
                else % No predators
                    N0(4) = 0;               N0(5) = 0;
                end
            end

        %% IVP based on email.
        elseif email == "y"
            if (runtype1 == "early")
                N0(3) = 0; % No biofilm present.
                if (predators1 == "TS") % Both predators
                    N0(4) = 0.5;   N0(5) = 1.9;
                elseif (predators1 == "T") % Just biofilm predator
                    N0(4) = 0;               N0(5) = 1.9;
                elseif (predators1 == "S") % Just free predator
                    N0(4) = 0.5;   N0(5) = 0;
                else % No predators
                    N0(4) = 0;               N0(5) = 0;
                end
            else % "Late" biofilm
                N0(3) = 5; % Established biofilm.
                if (predators1 == "TS") % Both predators
                    N0(4) = 0.5;             N0(5) = 1.9;
                elseif (predators1 == "T") % Just biofilm predator
                    N0(4) = 0;               N0(5) = 1.9;
                elseif (predators1 == "S") % Just free predator
                    N0(4) = 0.5;             N0(5) = 0;
                else % No predators
                    N0(4) = 0;               N0(5) = 0;
                end
            end
        end

        %%

        N0(6) = 0; % Initially no waste.

        opts = odeset('NonNegative',1:6); % Ensure all solutions are non-negative.

        [t,N] = ode89(@(t,y)... % Variables.
            odefun(t,y), ... % Function (declared at bottom).
            timeSteps, ... % Time steps to be sampled by solver.
            N0, ...  % Initial conditions.
            opts);

        % save vals for this run type
        F(preds) = N(end,2);
        B(preds) = N(end,3);

        %N(end,4)/N(1,4) % Free predator growth
        %N(end,5)/N(1,5) % Biofilm predator growth
        %hold on;

    end

    F = biovolume(F,'v');
    B = biovolume(B,'s');

    pops = [F B]; % top figs, total prey composition.

    bar(pops,'stacked')

end

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

%% Calculate the initial predator concentrations.

function M = predMassF(CV) % Free predator mass.

% Here, we assume that the cilliate is a spheroid.

rmaj = 320/2; % Major axis.
rmin = (50/120)*rmaj; % Minor axis.
Vpred = (4/3)*pi*rmaj*rmin^2; % Volume of cell in um^3.
M = 133.754*Vpred^(0.438); % Romanova 2009, get carbon mass per cell.
M = M * CV; % Convert from Mass/cell to mass/ml by multiplying cells/ml.
fudge1 = 1; % 15 is good
M = M*fudge1;
M = M * 10^-9; % Convert to ugC from fgC.

end

function M = predMassBF(CV) % Biofilm predator mass.

% Here, we assume that the ameoba are spherical...

r = 23; % Cell radius (assume spherical).
Vpred = (4/3)*pi*r^3; % Volume of cell in um^3.
M = 133.754*Vpred^(0.438); % Romanova 2009, get carbon mass per cell.
M = M * CV; % Convert from Mass/cell to mass/ml by multiplying cells/ml.
M = M * 10^-9; % Convert to ugC from fgC.
fudge2 = 1; % 20 is good
M = M*fudge2;
M = M * (3/8); % Multiply by vol and divide by area to convert from /ml to /cm^2.

end
