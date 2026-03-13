% plot output from vmmc StickySquares simulation

clear;

ptime = 0.0;  % how many seconds to pause for after each frame
kpl = 1;   % how many frames to skip when plotting. 1 = plot each frame. 


filehead = 'hier'; 
%filehead = 'crosstalk';


statsfile = [filehead,'_stats.txt'];
datafile = [filehead,'_traj.txt'];




% open the file and extract data
fid = fopen(datafile,'r');
tline = fgetl(fid);
disp(tline);
params = fscanf(fid,'%d %f %d\n\n',3);  % get the number of particles
n = params(1);  % number of particles
L = params(2);  % axis for plotting
n0 = params(3);  % base # of particles
x0 = fscanf(fid,'0 %f %f %f\n',[3,n]);  % get the first set of coordinates


% set up figure
figure(1)
clf
cmap = cmocean('thermal');  % linear map (light on top)  
axis([0,L,0,L]);
axis square
set(gca,'nextplot','replacechildren');
squares(x0(1:2,:)',cmap, n0);
drawnow
title('press a key to continue');
pause;

% loop through steps
nsteps = 0;
x = x0;
while(~isempty(x))
    
    % plot
    if(nsteps >= 0 && (mod(nsteps,kpl) == 0 && ~isinf(kpl)) || ((kpl == inf && nsteps == nt-1)))
        x = x(1:2,:)';
        figure(1)
        squares(x,cmap, n0);
        title(['t = ',num2str(nsteps)]);
        drawnow
        pause(ptime);
    end

    % get next point
    nsteps = nsteps + 1;
    fscanf(fid,'%d\n\n',1);
    x = fscanf(fid,'0 %f %f %f\n',[3,n]);
end

disp(['nsteps = ',num2str(nsteps)]);

% close file
fclose(fid);




% Read in statistics
if(0)
    % f0 = which fragment size to plot, 
    % maxf0 = max # of copies of this fragment size (limits for y axis)
    stats = readmatrix(statsfile,'NumHeaderLines',1);
    time = stats(:,1);
    energy = stats(:,2);
    eavg = cumsum(energy)./[1:length(energy)]';

    % plot energy versus time
    figure(2)
    clf
    plot(time,energy,time,eavg);
    xlabel('nsteps');
    ylabel('Energy');
    legend({'energy';'time avg'},'Location','Best');

    % find maximum fragment size
    fragmenthist = stats(:,3:end);
    nt = size(fragmenthist,1);
    s = zeros(nt,1);
    for it=1:nt
        s(it) = find(fragmenthist(it,:) > 0, 1, 'last');
    end

    % plot fragment size versus time
    figure(3)
    clf
    subplot(1,2,1)
    plot(time,stats(:,f0+2));
    set(gca,'ylim',[0 maxf0+0.1]);
    xlabel('nsteps');
    ylabel('# of completed fragments');
    subplot(1,2,2)
    plot(time,s);
    set(gca,'ylim',[0 f0+0.1]);
    xlabel('nsteps');
    ylabel('max fragment size');
end





% draw squares
% south-east corners is n x 2 vector
function squares(corners,cmap,n0)
iftext = 0;
xc = [0,1,1,0];
yc = [0,0,1,1];

cla
hold on
n = size(corners,1);
for i=1:n
    x0 = corners(i,1);
    y0 = corners(i,2);
    xpl=x0+xc;
    ypl=y0+yc;
    idx = mod(i,n0);
    if(idx == 0) idx = n0; end
    fill(xpl,ypl,cmap(ceil(idx*256/n0),:),'Linewidth',1);
    if(iftext) 
        text(x0+0.5,y0+0.5,num2str(idx-1),'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end
hold off
end






