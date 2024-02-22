
% Non-local distribution calculator.

Nruns = 10;
n = 320; % Corresponds to (100/(2*12*320)) swimmer lengths per bin.
edges = linspace(-8.5,8.5,n+1);

counter = 0;

for i = [1,40]
    for j = [1,12]
        hp = zeros(length(edges)-1,Nruns);
        for run = 1:Nruns
            counter = counter+1;
            pass = load(['outputs_extrema/' num2str(j) '_' num2str(i) '_' num2str(counter) '_outputs/output_passive.txt']);
            hp(:,run) = histcounts(pass(:,4),edges,'Normalization','pdf');
        end
        h1p = mean(hp(:,:),2);
        save(['outputs_extrema/' num2str(i) '_' num2str(j) '.mat'],'h1p');
    end
end

%% This assumes that the previous code has been correctly run.

n = 320; % Corresponds to (100/(2*12*320)) swimmer lengths per bin.
edges = linspace(-10,10,n+1);
centers = edges(1:end-1) + (edges(2)-edges(1))/2;
sym_distance = 100+12*centers(1:length(centers)/2);

smWindow = 1; % Small window of smoothing, multiply value in n comment for effective bin width.

h1 = load(['hists_extrema/1_1_10.mat']);
h1p = h1.h1p;
symp_passive = h1p(1:length(h1p)/2) + h1p(end:-1:1+length(h1p)/2);
plot(sym_distance,smooth(symp_passive,floor(smWindow*(12/10))),'LineWidth',4,'LineStyle','-'); hold on

h2 = load(['hists_extrema/1_10_20.mat']);
h2p = h2.h1p;
symp_passive = h2p(1:length(h2p)/2) + h2p(end:-1:1+length(h2p)/2);
plot(sym_distance,smooth(symp_passive,floor(smWindow*(12/10))),'LineWidth',4,'LineStyle','-');

h3 = load(['hists_extrema/40_1_30.mat']);
h3p = h3.h1p;
symp_passive = h3p(1:length(h3p)/2) + h3p(end:-1:1+length(h3p)/2);
plot(sym_distance,smooth(symp_passive,floor(smWindow*(12/10))),'LineWidth',4,'LineStyle','-');

h4 = load(['hists_extrema/40_10_40.mat']);
h4p = h4.h1p;
symp_passive = h4p(1:length(h4p)/2) + h4p(end:-1:1+length(h4p)/2);
plot(sym_distance,smooth(symp_passive,floor(smWindow*(12/10))),'LineWidth',4,'LineStyle','-');

xlim([0,100])

savefig('f2a');
saveas(gcf,'f2a.png');

%% Not symmetric

n = 320; % Corresponds to (100/(2*12*320)) swimmer lengths per bin.
edges = linspace(-8.5,8.5,n+1);
centers = edges(1:end-1) + (edges(2)-edges(1))/2;
sym_distance = (100/12) + centers(1:length(centers)/2);

smWindow = 1; % Small window of smoothing, multiply value in n comment for effective bin width.

h1 = load(['hists_extrema/1_1_10.mat']);
h1p = h1.h1p;
plot(centers,h1p,'LineWidth',4,'LineStyle','-'); hold on

h2 = load(['hists_extrema/1_10_20.mat']);
h2p = h2.h1p;
plot(centers,h2p,'LineWidth',4,'LineStyle','-');

h3 = load(['hists_extrema/40_1_30.mat']);
h3p = h3.h1p;
plot(centers,h3p,'LineWidth',4,'LineStyle','-');

h4 = load(['hists_extrema/40_10_40.mat']);
h4p = h4.h1p;
plot(centers,h4p,'LineWidth',4,'LineStyle','-');

savefig('f2a');
saveas(gcf,'f2a.png');

%% Accum test

h1 = h1.h1p;
h2 = h2.h1p;
h3 = h3.h1p;
h4 = h4.h1p;

%%

ind = 1:160;
valh1 = zeros(160,1);
valh2 = zeros(160,1);

h1(h1<0.005)=0;
h4(h4<0.005)=0;

for ii = ind

valh1(ii) = (sum(h1(1:ii)) + sum(h1(end-ii+1:end)))/sum(h1);
valh2(ii) = (sum(h4(1:ii)) + sum(h4(end-ii+1:end)))/sum(h4);

end

symp_passive1 = h1(1:length(h1)/2) + h1(end:-1:1+length(h1)/2);
symp_passive2 = h4(1:length(h4)/2) + h4(end:-1:1+length(h4)/2);

plot(10*symp_passive1,'LineWidth',5); hold on
plot(10*symp_passive2,'LineWidth',5)

yy = valh2./valh1;
yy(isnan(yy)) = 0;

plot(yy,'LineWidth',5)

