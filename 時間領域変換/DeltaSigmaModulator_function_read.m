% Number of points per signal
N = 800;

% Figure
close('all');
Color = 0.7*[1 1 1];
Figure = figure('Color',Color);
hold('on');
grid('on');
box('on');

% Second Y axes
[Axes,~,~] = plotyy(0,0,0,0);
linkaxes(Axes,'xy');
YTicks = linspace(0,1,11);
YTickLabels = arrayfun(@(y)sprintf('%G',y),YTicks-0.5,'UniformOutput',false);
set(Axes(1),...
    'YColor',  'k',...
    'Fontsize', 8);
set(Axes(2),...   
    'YColor',     'y',...
    'YLim',       [-1 +1],...
    'Ytick',      YTicks,...
    'Yticklabel', YTickLabels,...
    'Fontsize',   8);
grid(Axes(2),'on');

% Default plots
plot([1 N],0.5*[1 1],'r:');
h(1) = plot(0,'y');
h(2) = plot(0,'y');
h(3) = plot(0,'b');
h(4) = stairs(0,'g');

% Axes
set(gca,...
    'Color',    Color,...
    'Xtick',    0:10:N,...
    'Xlim',     [1,N],...
    'Ytick',    -1:0.5:1,...
    'Ylim',     [-1.1 +1.1],...
    'Fontsize', 8);
title('\Delta\Sigma modulation','Fontweight','Light','Fontsize',9);
xlabel(Axes(1),'Sample number','Fontsize',8);
ylabel(Axes(1),'DC','Fontsize',8);
ylabel(Axes(2),'AC','Fontsize',8);
legend([h(3) h(4)],{'Input signal (DC+AC)','Output signal'},'Color',Color,'Fontsize',8);

% Full screen
drawnow;
warning('off','all');
jFrame = get(Figure,'JavaFrame');
jFrame.setMaximized(true);
warning('on','all');

% AC test
%fileID = fopen('Impulse response for dispersion compensation.txt','r');
fileID = fopen('Impulse response for dispersion compensation.txt');
C = textscan(fileID,'%6f %6f');
fclose(fileID);
whos
C
cosG = C{1}{linspace(1,800,N)}

% Delta Sigma modulator
DSM = DeltaSigmaModulator('Oversampling',1);


    
    % Delta sigma modulator reset
    set(DSM,...
        'Sigma',          0,...
        'PreviousOutput', 0);        
           
    
    % Delta sigma modulation
    [Signal,SignalDS] = DSM.update(cosG);
    
    % Plot updates
    set(h(1),'Xdata',[1 N],...
             'Ydata',[1 1]*(cosG));    
    set(h(2),'Xdata',[1 N],...
             'Ydata',[1 1]*(SinG));
    set(h(3),'Ydata',Signal);
    set(h(4),'Ydata',SignalDS);
    
    fileID = fopen('DSM_data.txt','w');
    fprintf(fileID,'%6s %12s\n','cosG','output');
    fprintf(fileID,'%4f %4f\n',Signal,SignalDS);
    %   SignalとSignalDSに横軸の定義がないため１００個のデータがセットで定義されている。
    fclose(fileID);

