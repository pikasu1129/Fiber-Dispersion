% Number of points per signal
N = 100;

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

% Default plots :color select
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
    'Ytick',    0:0.5:1,...
    'Ylim',     [-0.1 +1.1],...
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
DC = 1/2;
AC = 1/2;
SINWT = sin(linspace(0,2*pi,N));

% Delta Sigma modulator
DSM = DeltaSigmaModulator('Oversampling',1);


    
    % Delta sigma modulator reset
    set(DSM,...
        'Sigma',          0,...
        'PreviousOutput', 0);        
           
    % Delta sigma modulation
    [Signal,SignalDS] = DSM.update(DC + AC * SINWT);
    
    % Plot updates
    set(h(1),'Xdata',[1 N],...
             'Ydata',[1 1]*(DC+AC));    
    set(h(2),'Xdata',[1 N],...
             'Ydata',[1 1]*(DC-AC));
    set(h(3),'Ydata',Signal);
    set(h(4),'Ydata',SignalDS);
    
    fileID = fopen('DSM_AC_data.txt','w');
    fprintf(fileID,'%6s %12s\n','DC+AC','output');
    fprintf(fileID,'%4f %4f\n',Signal,SignalDS);
    %   SignalとSignalDSに横軸の定義がないため１００個のデータがセットで定義されている。
    fclose(fileID);

