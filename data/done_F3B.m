%% Get the first moment of the distribution 

clear all

nbins = 100; % Number of bins
ndisps = 1000; % Number of times in the bins
bins = linspace(0,100,nbins+1);

disp_max = 1;
edges = linspace(-disp_max,disp_max,ndisps+1);
dispcents = edges(1:end-1) + (edges(2)-edges(1))/2;

% Script to calculate the wait times of the particles
counter = 0; % Counter, needed to order the files.
for i = [1,40] % Loop on values of kick freq
    for j = [1,10] % Loop on values of kick str

        all_y1 = [];

        for k = 1:10 % Loop over sim repeats 
            counter = counter+1; % Increment thecounter for loading the correct file
            load(['jumps_extrema/' num2str(j) '_' num2str(i) '_' num2str(counter) '_jumps.mat']); % Load the file
            y1 = zeros(length(jumps(:,1)),2);
            y1(:,1) = 12*(jumps(:,6)-jumps(:,5));
            y1(:,2) = 12*(jumps(:,5)); 
            all_y1 = [all_y1;y1];

        end

        all_y1 = all_y1(abs(all_y1(:,1))>0.5,:);
        all_y1 = all_y1(abs(all_y1(:,1))<50,:);
        
        all_y1(all_y1(:,2)<0,:) = -all_y1(all_y1(:,2)<0,:);
        
        h1 = zeros(nbins,ndisps);
        y1mean = zeros(nbins,1);
        
        for bin = 1:nbins
            a = find( all_y1(:,2) > bins(bin) );
            b = find( all_y1(a,2) < bins(bin+1) ); b = a(b);
            h1(bin,:) = histcounts(all_y1(b,1),edges,'Normalization','pdf');

            if length(b)>20
                y1mean(bin) = mean(all_y1(b,1));
            else
                y1mean(bin) = 0;
            end

            %plot(h1(bin,3:end));hold on;
        end
        
        bincents = bins(1:end-1) + (bins(2)-bins(1))/2;
        plot(100-bincents,smooth(y1mean,3),'LineWidth',4); hold on;
        %imagesc(dispcents,bincents,h1)

        save([num2str(i) '_' num2str(j) '_m1'],'y1mean');

    end
end

savefig('f3b');
saveas(gcf,'f3b.png');
