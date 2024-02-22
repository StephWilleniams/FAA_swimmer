% Non-local distribution calculator.

Nruns = 10;
n = 320; % Corresponds to (100/(2*12*320)) swimmer lengths per bin.
edges = linspace(-8.5,8.5,n+1);

for i = 1:40
    for j = 1:12

        pass = load(['outputs/' num2str(j) '_' num2str(i) '_outputs/output_passive.txt']);
        hp = zeros(length(edges)-1,1);
        hp = histcounts(pass(:,4),edges,'Normalization','pdf');
        save(['outputs/' num2str(i) '_' num2str(j) '_hist.mat'],'hp');

    end
end

%% The following assumes that the previous code was run correctly.

k = (1:10)*7.3/10; s = 0.1 + (15-0.1)*(1:40)/39;

n = 320; % Corresponds to (100/(2*12*320)) swimmer lengths per bin.
edges = linspace(-8.5,8.5,n+1);
centers = edges(1:end-1) + (edges(2)-edges(1))/2;

accum = zeros(10,40);

h1 = load(['hists_all/' num2str(1) '_' num2str(1) '_1.mat']);
h = h1.hp;
sym_distance = (100/12) + centers(1:length(centers)/2);
sym_passive = h(1:length(h)/2) + h(end:-1:1+length(h)/2);
comparison = sum(sym_passive(1:30))/sum(sym_passive);

for i = 1:10
    for j = 1:40

        h1 = load(['hists_all/' num2str(j) '_' num2str(i) '_1.mat']);
        h = h1.hp;
        sym_distance = (100/12) + centers(1:length(centers)/2);
        sym_passive = h(1:length(h)/2) + h(end:-1:1+length(h)/2);
        accum(i,j) = (sum(sym_passive(1:30))/sum(sym_passive))/comparison;

    end
end

imagesc(s,k,accum)
savefig('f2b');
saveas(gcf,'f2b.png');
