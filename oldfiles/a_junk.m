%% Junk file.

% Note: I like to always make a script like this. Each section is just some
% random spitballing of ideas. They usually arent important enough to be
% comprehensively commented (or, if they were they made their way into main files, 
% where they get full documentation). 

%%

x = 8;

cpc = x*10^9/(6*10^6);

a = (40+x)/cpc;
b = (cpc/133)^0.453;
c = 3;
disp(a*b*c)

%%

% Here, we assume that the PT cells are spheroids.

predMassF(30)
predMassBF(3000)

function M = predMassF(CV)

rmaj = 120; % Major axis.
rmin = 50; % Minor axis.
Vpred = (4/3)* pi*rmaj * rmin^2; % Volume of cell in um^3.
M = 133.754*Vpred^(0.436); % Romanova 2009, get carbon mass per cell.
M = M * CV; % Convert from Mass/cell to mass/ml by multiplying cells/ml.
M = M * 10^-9; % Convert to ugC from fgC.

end

function M = predMassBF(CV)

r = 18;
Vpred = (4/3)*pi*r^3; % Volume of cell in um^3.
M = 133.754*Vpred^(0.436); % Romanova 2009, get carbon mass per cell.
M = M * CV; % Convert from Mass/cell to mass/ml by multiplying cells/ml.
M = M * 10^-9; % Convert to ugC from fgC.

end

%%


