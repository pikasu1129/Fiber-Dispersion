% Number of points per signal
N = 30;

% Figure
close('all');
Color = 0.7*[1 1 1];
Figure = figure('Color',Color);
hold('on');
grid('on');
box('on');

% Default plots
h(1) = plot(0,'b');
h(2) = stairs(0,'g');

% Axes
set(gca,...
    'Color',    Color,...
    'Xtick',    0:10:N,...
    'Xlim',     [1,N],...
    'Ytick',    0:0.1:1,...
    'Ylim',     [-0.1 1.1],...
    'Fontsize', 8);
title('\Delta\Sigma modulation','Fontweight','Light','Fontsize',9);
xlabel('Sample number','Fontsize',8);
ylabel('Value','Fontsize',8);
legend({'Input signal (DC)','Output signal'},'Color',Color,'Fontsize',8);

% Full screen
drawnow;
warning('off','all');
jFrame = get(Figure,'JavaFrame');
jFrame.setMaximized(true);
warning('on','all');

% DC test
DC = 0:1e-3:1;
DSM = DeltaSigmaModulator('Oversampling',N);

for dc = 1:numel(DC)
    
    % Delta sigma modulator reset
    set(DSM,...
        'Sigma',          0,...
        'PreviousOutput', 0);        
       
    % Delta sigma modulation
    [Signal,SignalDS] = DSM.update(DC(dc));
    
    % Plot update
    set(h(1),'Ydata',Signal);
    set(h(2),'Ydata',SignalDS);
    
    % Pause
    pause(0.01);
    
end



