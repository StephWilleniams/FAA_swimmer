%% Get the mixed moment of the distribution 

clear all

nbins = 80; % Number of bins
ndisps = 1000; % Number of points in the bins
bins = linspace(0,100,nbins+1);
disp_max = 100;
edges = linspace(-disp_max,disp_max,ndisps+1);
dispcents = edges(1:end-1) + (edges(2)-edges(1))/2;
bincents = bins(1:end-1) + (bins(2)-bins(1))/2;

% Script to calculate the wait times of the particles
counter = 0; % Counter, needed to order the files.
for i = [1,40] % Loop on values of kick freq
    for j = [1,10] % Loop on values of kick str

        all_y1 = [];

        for k = 1:10 % Loop over sim repeats 
            counter = counter+1; % Increment thecounter for loading the correct file
            load(['jumps_extrema/' num2str(j) '_' num2str(i) '_' num2str(counter) '_jumps.mat']); % Load the file
            y1 = zeros(length(jumps(:,1)),3);
            y1(:,1) = 12*(jumps(:,4)-jumps(:,3)); % x-displacement
            y1(:,2) = 12*(jumps(:,6)-jumps(:,5)); % y-displacement
            y1(:,3) = 12*jumps(:,5); 
            all_y1 = [all_y1;y1];

        end

        a = find(all_y1(:,3) < 0);
        all_y1(a,1:2) = -all_y1(a,1:2); % Flip events all into top half of channel.
        
        all_y1 = all_y1(abs(all_y1(:,1))<50,:);
        

        h1 = zeros(nbins,ndisps);
        mixmomean = zeros(nbins,1);
        
        for bin = 1:nbins
            a = find( all_y1(:,3) > bins(bin) );
            b = find( all_y1(a,3) < bins(bin+1) ); b = a(b);

            mixmos = zeros(length(b),1);

            mixmos = (all_y1(b,1) - mean(all_y1(b,1))).*(all_y1(b,2) - mean(all_y1(b,2)))/(std(all_y1(b,1) - mean(all_y1(b,1)))*std(all_y1(b,2) - mean(all_y1(b,2))));

            h1(bin,:) = histcounts(mixmos,edges,'Normalization','pdf');
            if length(b)>20
                mixmomean(bin) = mean(mixmos);
            else
                mixmomean(bin) = 0;
            end
            %plot(h1(bin,3:end));hold on;
        end

        mixmomean(isnan(mixmomean))=0;

        plot(bincents,mixmomean,'LineWidth',4); hold on;

    end

end

savefig('f4a');
saveas(gcf,'f4a.png');