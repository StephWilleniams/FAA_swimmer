
%% Get the lengths, visual confirmation this is correct.

theta = linspace(0,2*pi);

x1 = 1;
x2 = 0;
x3 = 0;

x = [x1,x2,x3];

c = 0.5;

m = (5+sqrt(2))/2;

s1 = x + c*(2-m)*[cos(x3),sin(x3),0];
s2 = s1 + 2*c*[cos(x3),sin(x3),0];
s3 = s2 + c*sqrt(2)*[cos(x3),sin(x3),0];

scatter(x(1),x(2),1000)
hold on;

scatter(s1(1),s1(2));
plot(s1(1) + 2*c*cos(theta),s1(2) + 2*c*sin(theta));

scatter(s2(1),s2(2));
plot(s2(1) + sqrt(2)*c*cos(theta),s2(2) + sqrt(2)*c*sin(theta));

scatter(s3(1),s3(2));
plot(s3(1) + c*cos(theta),s3(2) + c*sin(theta));
axis equal

%% 

rA = 1;
p = -1;

rL = (4*sqrt(2)/(5*sqrt(2)+2))*rA;
rM = rL/sqrt(2);
rS = rL/2;

x = [0,0,0];

xp1 = x + p*(rL - rA)*[cos(x(3)),sin(x(3)),0];
xp2 = x + p*(rL*2 - rA)*[cos(x(3)),sin(x(3)),0];
xp3 = x + p*(rL*2 + rM - rA)*[cos(x(3)),sin(x(3)),0];

theta = linspace(0,2*pi);
hold on
plot(xp1(1) + rL*cos(theta'),xp1(2) + rL*sin(theta'));
plot(xp2(1) + rM*cos(theta'),xp2(2) + rM*sin(theta'));
plot(xp3(1) + rS*cos(theta'),xp3(2) + rS*sin(theta'));
