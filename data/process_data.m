
clear all 
close all

n = 320;
xx = 8.333-linspace(0,8,n);

%%

act = load("1_1_output_active.txt");
pass = load("1_1_output_passive.txt");

% h1a = histcounts(abs(act(abs(act(:,4)),4)),2000);
% h1p = histcounts(abs(pass(abs(pass(:,4)),4)),2000);
h1a = histcounts(act(:,4),n,'Normalization','pdf');
h1p = histcounts(pass(:,4),n,'Normalization','pdf');

%%

clear act pass

%%

act = load("40_12_output_active.txt");
pass = load("40_12_output_passive.txt");

h2a = histcounts(act(:,4),n,'Normalization','pdf');
h2p = histcounts(pass(:,4),n,'Normalization','pdf');

%% Actives

smWindow = 2;

c = jet(4);

%plot(xx,smooth(h2p,floor(smWindow*(12/10))),'LineWidth',4,'Color',c(1,:),'LineStyle','-');hold on
plot(xx,smooth(h1a,smWindow),'LineWidth',4,'Color',[0.2 0.8 0.2],'LineStyle','-'); hold on
%plot(xx,smooth(h2a,smWindow),'LineWidth',4,'Color',c(3,:),'LineStyle','-'); 
plot(xx,smooth(h1p,floor(smWindow*(12/10))),'LineWidth',4,'Color',[0.2 0.2 1],'LineStyle','-');

%xlim([-0.2,3])
%ylim([0,4])
axis equal