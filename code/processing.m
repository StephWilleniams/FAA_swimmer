%%

clear all
close all
 
%%

num = 4;

S = readmatrix([num2str(num) '_outputs/output_active.txt']);
C = readmatrix([num2str(num) '_outputs/output_passive.txt']);

%%

%histogram(S(:,3))
histogram(S(:,4),100)