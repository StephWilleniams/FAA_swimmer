%%

S = readmatrix(['1_0_outputs/output_active.txt']);

%%

i = 1;

S = readmatrix(['outputs/1_' num2str(i) '_outputs/output_active.txt']);
C = readmatrix(['outputs/1_' num2str(i) '_outputs/output_passive.txt']);

%%

un = unique(S(:,1));
theta = linspace(0,2*pi,100);

%rSeg = zeros(3,1);
rA = 1;
rC = 0.5;
r0 = (4*sqrt(2)/(5*sqrt(2)+2));
rSeg = [r0 * rA, ...
       r0 * rA/sqrt(2), ...
       r0 * rA/2];
%sDV = zeros(3,1);
sDV = -[rSeg(1) - rA, ...
    rSeg(1)*2 - rA, ...
    rSeg(1)*2 + rSeg(2) - rA];

n = 30;
for i = 1:n:length(un)

    acts = find(S(:,1) == un(i));

    for np = 1
        for seg = 1:3
            if S(acts(np),6) == 0
                plot( S(acts(np),3)+sDV(3-seg+1)*cos(S(acts(np),5))+rSeg(3-seg+1)*cos(theta) , S(acts(np),4)+sDV(3-seg+1)*sin(S(acts(np),5))+rSeg(3-seg+1)*sin(theta),'g','LineWidth',5);
            else
              plot( S(acts(np),3)+sDV(3-seg+1)*cos(S(acts(np),5))+rSeg(3-seg+1)*cos(theta) , S(acts(np),4)+sDV(3-seg+1)*sin(S(acts(np),5))+rSeg(3-seg+1)*sin(theta),'r','LineWidth',5);
            end
            xlim([-60,60])
            ylim([-10,10])
            hold on;
        end
    end
    pause(0.1)
    hold off;

end


%%

nindmax = 10;
nparts=5;

store = [];
dx =[];

for i = 1:nindmax
    
    S = readmatrix(['outputs/1_' num2str(i) '_outputs/output_active.txt']);

    untime = unique(S(:,1));

    for na = 1:nparts

        a = find( S(:,2) == na);
        b = find(S(a,6)==1);
        
        if ~isempty(b)
            store = [store;b(1)*0.0005];
        else 
            disp('empty')
        end

        if (S(a(b(1)),3)-S(a(1),3) < 0)
            dx = [dx;S(a(b(1)),3)-S(a(1),3)];
        else
            dx = [dx;0];
        end

    end

end

plot(dx*10);
hold on
plot(-store*50);

%%

i = 1;

S = readmatrix(['outputs/1_' num2str(i) '_outputs/output_active.txt']);

plot(S(:,3),S(:,4))
% plot(S(:,5))
% plot(S(:,6))


%%

kickFreq = 0.0005*1;

for i = 1:1/0.0005
    if ceil( rand-(exp(-kickFreq)) ) > 0
        ceil( rand-(exp(-kickFreq)) ) * cos(theta);
    end
end

%%

theta = pi/2+0.00001;
ceil( 1-(exp(-kickFreq)) ) * cos(theta)/abs(cos(theta));

%%

C


