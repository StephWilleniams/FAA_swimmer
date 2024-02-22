
clear all 
close all

%%

pass = load("hists_extrema/1_1_10.mat");
pass = pass.h1p;

%% Calculate the histograms

n = 320; % Corresponds to (100/(2*12*320)) swimmer lengths per bin.
edges = linspace(-8.5,8.5,n+1);
centers = edges(1:end-1) + (edges(2)-edges(1))/2;

h1a = histcounts(act(:,4),edges,'Normalization','pdf');
h1p = histcounts(pass(:,4),edges,'Normalization','pdf');

%% Get the positional distributions

%% Main body

smWindow = 2; % Small window of smoothing, multiply value in n comment for effective bin width.

% Symmetrise the distributions
sym_distance = (100/12) + centers(1:length(centers)/2);
sym_active = h1a(1:length(h1a)/2) + h1a(end:-1:1+length(h1a)/2);
symp_passive = h1p(1:length(h1p)/2) + h1p(end:-1:1+length(h1a)/2);

plot(sym_distance,smooth(sym_active,smWindow),'LineWidth',4,'Color',[0.2 0.8 0.2],'LineStyle','-'); hold on
plot(sym_distance,smooth(symp_passive,floor(smWindow*(12/10))),'LineWidth',4,'Color',[0.2 0.2 1],'LineStyle','-');

% Additional images and interaction line manually added.


%% Inset (unsmoothed)

plot(centers,h1a,'LineWidth',4,'Color',[0.2 0.8 0.2],'LineStyle','-'); hold on
plot(centers,h1p,'LineWidth',4,'Color',[0.2 0.2 1],'LineStyle','-');

%xlim([-0.2,3])
%ylim([0,4])