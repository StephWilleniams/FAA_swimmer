%% Get the second parallel moment of the distribution 

clear all

% Script to calculate the wait times of the particles
counter = 0; % Counter, needed to order the files.
for i = [1,40] % Loop on values of kick freq
    for j = [1,10] % Loop on values of kick str

        all_x1 = [];

        for k = 1:10 % Loop over sim repeats 
            counter = counter+1; % Increment thecounter for loading the correct file
            load(['jumps_extrema/' num2str(j) '_' num2str(i) '_' num2str(counter) '_jumps.mat']); % Load the file
            y1 = zeros(length(jumps(:,1)),2);
            y1(:,1) = 12*(jumps(:,4)-jumps(:,3));
            y1(:,2) = 12*(jumps(:,5)); 
            all_x1 = [all_x1;y1];

        end
        
        nbins = 80; % Number of bins
        ndisps = 10000; % Number of times in the bins
        bins = linspace(0,100,nbins+1);

        disp_max = 100;
        edges = linspace(-disp_max,disp_max,ndisps+1);
        dispcents = edges(1:end-1) + (edges(2)-edges(1))/2;
        
        a = find(all_x1(:,2)<0);
        all_x1(a,1:2) = -all_x1(a,1:2);
        all_x1 = all_x1(abs(all_x1(:,1))<100,:);
        %all_x1(:,1) = all_x1(:,1).^2;
        
        h1 = zeros(nbins,ndisps);
        x2mean = zeros(nbins,1);
        
        for bin = 1:nbins
            a = find( all_x1(:,2) > bins(bin) );
            b = find( all_x1(a,2) < bins(bin+1) ); b = a(b);
            h1(bin,:) = histcounts(all_x1(b,1),edges,'Normalization','pdf');
            if length(b)>10
                x2mean(bin) = mean( (all_x1(b,1) - mean(all_x1(b,1))).^2 );
            else
                x2mean(bin) = 0;
            end
            %plot(h1(bin,3:end));hold on;
        end
        
        bincents = 100 - (bins(1:end-1) + (bins(2)-bins(1))/2);
        plot(bincents,smooth(sqrt(x2mean),3),'LineWidth',4); hold on;
        %imagesc(dispcents,bincents,h1)

    end
end

savefig('f3d');
saveas(gcf,'f3d.png');