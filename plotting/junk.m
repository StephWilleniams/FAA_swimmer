
theta = pi;
c = cos(theta);
s = sin(theta);
fpar = 0.5;
fperp = 2;
ft = fpar*[c*c, c*s; s*c, s*s] + fperp*([1,0;0,1]-[c*c, c*s; s*c, s*s]);
F = inv(ft)*[1;0];
disp(F)

%%



a = (1+1+0.5+1/sqrt(2))/2;

fricPar = 2*pi/ (log(a) - 0.207  + (0.980/a) - (0.133/a^2));
fricPerp = 4*pi/(log(a) + 0.839  + (0.185/a) + (0.233/a^2)); 
fR = pi*a^2/(3* (log(a) - 0.662  + (0.917/a) - (0.050/a^2))); 

%%

dt = 0.0005;
k=0:0.1:20000;
out = 1-(exp(-dt*k));
plot(k,out);

%%


in = -21;
yT = 10;
yB = -10;

((in+10)-floor((in+10)/(yT-yB))*(yT-yB))-10
