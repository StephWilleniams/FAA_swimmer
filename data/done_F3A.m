
clear all

nbins = 100; % Number of bins
ntimes = 100; % Number of times in the bins
bins = linspace(0,100,nbins+1);
%bins = linspace(-8.75,8.75,nbins+1);
edges = linspace(0,30,ntimes+1);
cents = 100-(bins(1:end-1) + (bins(2)-bins(1))/2);
%cents = (bins(1:end-1) + (bins(2)-bins(1))/2);

% Script to calculate the wait times of the particles
counter = 0; % Counter, needed to order the files.

for i = [1,40] % Loop on values of kick freq
    for j = [1,10] % Loop on values of kick str

        all_times = [];

        for k = 1:10 % Loop over sim repeats 

            counter = counter+1; % Increment thecounter for loading the correct file

            load(['jumps_extrema/' num2str(j) '_' num2str(i) '_' num2str(counter) '_jumps.mat']); % Load the file
            
            times = zeros(length(jumps(:,1))-1,3);
            times(:,1) = jumps(2:end,2)-jumps(1:end-1,1);
            times(:,2) = jumps(2:end,6)-jumps(1:end-1,5);
            times(:,3) = 12*(jumps(1:end-1,6));
            all_times = [all_times;times];

        end
        
        %all_times = all_times(all_times(:,2)>0.2,:);
        all_times = all_times(all_times(:,1)>0,:); % Remove the times at which the particle index changes.

        all_times1 = zeros(size(all_times));
        all_times1(:,1) = 0.0005*all_times(:,1);
        all_times1(:,3) = abs(all_times(:,3));
        %all_times1(:,2) = (all_times(:,2));
        all_times1 = all_times1( all_times1(:,1) > min(all_times1(:,1))*10 ,:);
        
        h1 = zeros(nbins,ntimes);
        tmean = zeros(nbins,1);
        
        for bin = 1:nbins
            a = find( all_times1(:,3) > bins(bin) );
            b = find( all_times1(a,3) < bins(bin+1) ); b = a(b);
            h1(bin,:) = histcounts(all_times1(b,1),edges,'Normalization','pdf');
            if length(b)>20 % Arbitrary data requirement to stop weird stuff
                tmean(bin) = mean(all_times1(b,1));
            else
                tmean(bin) = 0;
            end
            %plot(h1(bin,3:end));hold on;
        end
        
        tmean(isnan(tmean)) = 0;

        %imagesc(h1(:,3:end))
        plot(cents,smooth(tmean,3),'LineWidth',4); hold on;

        save([num2str(i) '_' num2str(j) '_wait'],'tmean');

    end
end

savefig('f3a');
saveas(gcf,'f3a.png');
