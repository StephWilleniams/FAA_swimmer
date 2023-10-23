A = readmatrix('output.txt');

un = unique(A(:,1));

%%

close all

for iq = 1:1:length(un)
    a = find(A(:,1) == un(iq));
    theta = linspace(0,2*pi,100);
    R = 1;
    for ii = 1:length(a)
        plot(A(a(ii),3)+R*cos(theta),A(a(ii),4)+R*sin(theta))
        xlim([-30,30])
        ylim([-30,30])
        hold on
    end
    pause(0.5)
    hold off
end

%%

st = 4;
et = 1000;

plot(A(1:2:end,3),A(1:2:end,4))