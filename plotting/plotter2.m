%% Script to plot the outputs from output code.

%%

clear all
close all
 
%%

S = readmatrix('../code/outputs/1_1_outputs/output_active.txt');
C = readmatrix('../code/outputs/1_1_outputs/output_passive.txt');

%%

num = 1;

S = readmatrix([num2str(num) '_outputs/output_active.txt']);
C = readmatrix([num2str(num) '_outputs/output_passive.txt']);

%% pol 1

%rSeg = zeros(3,1);
rA = 1;
rC = 0.5;
r0 = (4*sqrt(2)/(5*sqrt(2)+2))/2;
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

unTime = unique(S(:,1));
theta = linspace(0,2*pi,100);

for i = 3000:5:4000%1600:1:2001%:10:1000%length(unTime)

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
            if S(acts(n),6) == 0
                plot( S(acts(n),3)+sDV(3-seg+1)*cos(S(acts(n),5))+rSeg(3-seg+1)*cos(theta) , S(acts(n),4)+sDV(3-seg+1)*sin(S(acts(n),5))+rSeg(3-seg+1)*sin(theta),'g','LineWidth',5)
            else
                plot( S(acts(n),3)+sDV(3-seg+1)*cos(S(acts(n),5))+rSeg(3-seg+1)*cos(theta) , S(acts(n),4)+sDV(3-seg+1)*sin(S(acts(n),5))+rSeg(3-seg+1)*sin(theta),'r','LineWidth',5)
            end
        end
    end
    xlim([5,20])
    ylim([-9,-4])
    pause(0.1)
    hold off;

end
