
clear all;

i = 1;
j = 1;
smwin = 5;

histo = load([ 'hists_extrema/' num2str(i) '_'  num2str(j) '_10.mat']);
comp = histo.h1p;
sym_comp = comp(1:length(comp)/2)+comp(end:-1:length(comp)/2+1);

noise = 0.5;
r = linspace(0,100);

%t = 2 - 1.5*exp(-r/20);
t = (4.7*exp(-r/2.5)+1.5) - 5*exp(-r/2) + (0.5 - 0.5*exp(-r/10));
plot(r,t); hold on
t = 1./t;

m1 = -0.1-25*exp(-r/5) + 25*exp(-r/3.5) + 1*exp(-((r-25)/2).^2) -1*exp(-((r-30)/2).^2)  ;

%m1 = m1 + noise*randn(length(m1),1)';
%m1 = -0.1-25*exp(-r/5) + 25*exp(-r/3.5);
plot(r,m1)

m2 = (7*exp(-r/20)+6).^2 - (7+6)^2*exp(-r/4);
%m2 = m2 + noise^2*randn(length(m2),1)';
plot(r,sqrt(m2));

D0 = 2;
Deff = (D0 + 0.5*t.*m2);
%Veff = (t.*m1 - 0.5*gradient(t.*m2)); 

%%

alpha = 12;

plot(alpha*t.*m1+gradient(Deff))

%%

expo = zeros(length(r),1);
P = zeros(length(r),1);

alpha = 1;

expo(1) = alpha*t(1).*(m1(1)./Deff(1));

for ii = 2:length(r)
    expo(ii) = expo(ii-1) + alpha*t(ii)*(m1(ii)/Deff(ii));
end

for ii = 1:length(P)
    P(ii) = (1/Deff(ii))*exp(expo(ii));
end

sx = linspace(0,1,160);
r = linspace(0,1,length(P));

%plot(expo); hold on; plot(exp(expo));
plot(r,P/sum(P)); hold on
ylim([0,0.02])
%plot(sx,sym_comp(1:end)/sum(sym_comp(1:end)))
