
clear all
close all

%%
pass = load("outputs_10K/12_40/output_passive.txt");

%%

c = jet(10);

for part_ind = 18

    npts = 200; % Number of points along the trajectory

    pos = zeros(npts,3);
    t0 = 500; tend = t0+(50*npts)-50;
    times = t0:50:tend;

    for n = 1:npts
        t = times(n);
        a = find(pass(:,1) == t);
        b = find(pass(a,2) == part_ind);
        b = a(b);
        pos(n,1:2) = pass(b,3:4);
    end

    pos(1:end-1,3) = sqrt((pos(2:end,1) - pos(1:end-1,1)).^2 + (pos(2:end,2) - pos(1:end-1,2)).^2);

    %plot(pos(:,1),pos(:,2),'Color',c(1,:),linewidth=3); hold on
    %scatter(pos(1,1),pos(1,2))

    x = pos(:,1)';
    y = pos(:,2)';
    z = zeros(length(pos(:,1)),1)';
    lineColor = (smooth(12*pos(:,3)/(0.0005*50),3))';

    surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);

    ylim([6.5,8.6])

    hold on;

end

%%

x = linspace(0, 2*pi, 1920); % HDTV resolution.
y = sin(x);
z = zeros(size(x));
lineColor = x;  % This is the color, it varies with x in this case.
% Plot the line with width 8 so we can see the colors well.
surface([x;x], [y;y], [z;z], [lineColor;lineColor],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 8);
grid on;
