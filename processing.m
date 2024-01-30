%%

clear all
close all
 
%%

S = readmatrix('outputs/output_active.txt');
C = readmatrix('outputs/output_passive.txt');

%%

%histogram(S(:,3))
histogram(S(:,4))