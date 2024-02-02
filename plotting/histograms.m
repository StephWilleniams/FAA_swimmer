%%

clear all
close all

%%

S1 = readmatrix(['../goodData/puller_noKicks_dt_over_2/output_active.txt']);
S2 = readmatrix(['../goodData/puller_noKicks/output_active.txt']);
S3 = readmatrix(['../goodData/puller_noKicks_dt_times_2/output_active.txt']);
%S = readmatrix(['../goodData/dt_over_2_pusher_noKicks/output_active.txt']);

%%

%histogram(S(:,3),100)
numBins = 100;
H1 = histcounts(S1(1000:end,4),numBins,'BinLimits',[-8.3,8.3],'Normalization','pdf');
H2 = histcounts(S2(1000:end,4),numBins,'BinLimits',[-8.3,8.3],'Normalization','pdf');
H3 = histcounts(S3(1000:end,4),numBins,'BinLimits',[-8.3,8.3],'Normalization','pdf');

H1_sym = H1(1:length(H1)/2) + H1(end:-1:1+length(H1)/2);
H2_sym = H2(1:length(H2)/2) + H2(end:-1:1+length(H2)/2);
H3_sym = H3(1:length(H3)/2) + H3(end:-1:1+length(H3)/2);
plot(H1_sym,'LineWidth',5);hold on;
plot(H2_sym,'LineWidth',5);
plot(H3_sym,'LineWidth',5);

%%

% a = find(abs(S2(:,4))<9);
% histogram(S2(a,4),numBins,'BinLimits',[-64,64],'Normalization','pdf')
