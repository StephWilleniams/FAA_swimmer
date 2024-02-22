
clear all
close all

i = 40;
j = 10;
smwin = 1;

histo = load([ 'hists_extrema/' num2str(i) '_'  num2str(j) '_40.mat']);
comp = histo.h1p;
sym_comp = comp(1:length(comp)/2)+comp(end:-1:length(comp)/2+1);

path = [ num2str(i) '_'  num2str(j) '_wait.mat'];
in = load(path);
t = in.tmean;
t = smooth(t,5);

path = [ num2str(i) '_'  num2str(j) '_m1.mat'];
in = load(path);
m1 = in.y1mean;
m1 = smooth(m1,5);

path = [ num2str(i) '_'  num2str(j) '_m2.mat'];
in = load(path);
m2= in.y2mean;
m2 = smooth(m2,5);

a = find(t == 0);
b = find(m1 == 0);
c = find(m2 == 0);

d = min([a(1),b(1),c(1)]);

t = 1./flip(t(1:d));
m1 = -flip(m1(1:d));
m2 = flip(m2(1:d));

Deff = 0.05 + 0.5*t.*m2;

p = zeros(length(t),1);
expo = zeros(length(t),1);

alpha = 12;

expo(1) = alpha*t(1)*m1(1)/Deff(1);
p(1) = (1/Deff(1))*exp(t(1)*m1(1)/Deff(1));

for i = 2:length(t)

    expo(i) = expo(i-1) + alpha*t(i)*m1(i)/Deff(i);
    p(i) = (1/Deff(i))*exp(t(i)*m1(i)/Deff(i));

end

px = linspace(0,1,length(p));
sx = linspace(0,1,length(sym_comp(20:end)));

plot(px,p,'LineWidth',5); hold on
plot(sx,sym_comp(20:end),'LineWidth',5)
ylim([0,0.5])

% plot(expo)

%%
