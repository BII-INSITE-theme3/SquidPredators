%% Title: Model of free/biofilm cells with two types of predator
% Author: S williams

% This code required another code to be run beforehand so you have access
% to the array N(:,3), which describes the biofilm population's dynamics.

%% binding rate free -> biofilm plotter

a = 0.01; % Strength of dependence
chimax = 0.05; % Maximum binding rate
chimin = 0.0005; % Minimum binding rate
chipl = (a*chimax+chimin.*N(:,3))./(a+N(:,3)); % Effective binding rate
t = linspace(t0,t1,length(N(:,3))); % Time array
plot(t,chipl) % plot result

% Note: this should be run after b_Seilier2017.m with "early" and "late"
% conditions as a comparison. 
% The results show that the form of the detachment rate is largely unused.
