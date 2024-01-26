A = readmatrix('output_active.txt');

un = unique(A(:,1));

close all

rA = 1;

rL = (4*sqrt(2)/(5*sqrt(2)+2))*rA;
rM = rL/sqrt(2);
rS = rL/2;
p = 1;

xp1 = x + p*(rL - rA)*[cos(x(3)),sin(x(3)),0];
xp2 = x + p*(rL*2 - rA)*[cos(x(3)),sin(x(3)),0];
xp3 = x + p*(rL*2 + rM - rA)*[cos(x(3)),sin(x(3)),0];

for iq = 1:1:length(A(:,1))
    a = find(A(:,1) == un(iq));
    theta = linspace(0,2*pi,100);
    R = 1;
    for ii = 1:length(a)
        plot(A(a(ii),3)+R*cos(theta),A(a(ii),4)+R*sin(theta))
        hold on
        plot(A(a(ii),3) + p*(rL - R)*cos(A(a(ii),5))        + rL*cos(theta), A(a(ii),4) +  p*(rL - R)*sin(A(a(ii),5))        + rL*sin(theta));
        plot(A(a(ii),3) + p*(rL*2 - R)*cos(A(a(ii),5))      + rM*cos(theta), A(a(ii),4) +  p*(rL*2 - R)*sin(A(a(ii),5))      + rM*sin(theta));
        plot(A(a(ii),3) + p*(rL*2 + rM - R)*cos(A(a(ii),5)) + rS*cos(theta), A(a(ii),4) +  p*(rL*2 + rM - R)*sin(A(a(ii),5)) + rS*sin(theta));
        xlim([-10,10])
        ylim([-10,10])
        %axis equal
    end
    pause(0.1)
    hold off
end

%%

st = 4;
et = 1000;

plot(A(1:2:end,3),A(1:2:end,4))