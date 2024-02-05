%% Script to plot the outputs from output code.

%%

clear all
close all

%% pol 1

%rSeg = zeros(3,1);
rA = 1;
rC = 0.5;
r0 = (4*sqrt(2)/(5*sqrt(2)+2));
rSeg = [r0 * rA, ...
       r0 * rA/sqrt(2), ...
       r0 * rA/2];
%sDV = zeros(3,1);
sDV = [rSeg(1) - rA, ...
    rSeg(1)*2 - rA, ...
    rSeg(1)*2 + rSeg(2) - rA];

%% pol -1

%rSeg = zeros(3,1);
rA = 1;
rC = 0.5;
r0 = (4*sqrt(2)/(5*sqrt(2)+2));
rSeg = [r0 * rA, ...
       r0 * rA/sqrt(2), ...
       r0 * rA/2];
%sDV = zeros(3,1);
sDV = -[rSeg(1) - rA, ...
    rSeg(1)*2 - rA, ...
    rSeg(1)*2 + rSeg(2) - rA];

%%

S = readmatrix('outputs/1_1_outputs/output_active.txt');
C = readmatrix('outputs/1_1_outputs/output_passive.txt');

%%

unTime = unique(C(:,1));
theta = linspace(0,2*pi,100);

for i = 1:1:16000%1600:1:2001%:10:1000%length(unTime)

    acts = find(S(:,1) == unTime(i));
    pass = find(C(:,1) == unTime(i));

    nS = length(acts);
    nC = length(pass);
    for n = 1:nC
        plot(C(pass(n),3) + rC*cos(theta),C(pass(n),4) + rC*sin(theta),'k','LineWidth',5)
        hold on
    end

    for n = 1:nS
        for seg = 1:3
            plot( S(acts(n),3)+sDV(3-seg+1)*cos(S(acts(n),5))+rSeg(3-seg+1)*cos(theta) , S(acts(n),4)+sDV(3-seg+1)*sin(S(acts(n),5))+rSeg(3-seg+1)*sin(theta),'g','LineWidth',5)
        end
    end
    xlim([-30,30])
    ylim([-30,30])
    pause(0.5)
    hold off;

end

%%

n = 1;

a = find( S(:,2) == n );

plot(S(a,6))
