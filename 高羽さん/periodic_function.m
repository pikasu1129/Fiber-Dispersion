close('all');
%Color = 0.7*[1 1 1];
%set(gca,...
%    'Color',    Color,...
%    'Xtick',    0:10:N,...
%    'Xlim',     [1,N],...
%    'Ytick',    -1:0.5:1,...
%    'Ylim',     [-1.1 +1.1],...
%    'Fontsize', 8);


x = linspace(-400*10^(-12),400*10^(-12),800);
L = 80000;
ramda = 1.55*10^(-6);
c = 3*10^(8);
D = -0.000016;
p = 1/2;
%b = -2*10^(-26);
b = -D*(ramda^2)/(2*pi*c);
cosH_cmp = cos((x.^2)/(2*b*L)-pi/4);
sinH_cmp = sin((x.^2)/(2*b*L)-pi/4);
figure(1);
plot(x,cosH_cmp,'LineWidth',5);

hold on
grid on
figure(2);
plot(x,sinH_cmp,'LineWidth',5);
hold on
disp (b);
DSM = DeltaSigmaModulator('Oversampling',1);


    
    % Delta sigma modulator reset
    set(DSM,...
        'Sigma',          0,...
        'PreviousOutput', 0);        
           
    
    % Delta sigma modulation
    [Signal,SignalDS] = DSM.update(p*cosH_cmp+p);
    figure(1);
    plot(x,SignalDS);
    hold on
    
DSM = DeltaSigmaModulator('Oversampling',1);


    
    % Delta sigma modulator reset
    set(DSM,...
        'Sigma',          0,...
        'PreviousOutput', 0);        
           
    
    % Delta sigma modulation
    [Signal,SignalDS] = DSM.update(p*sinH_cmp+p);
    figure(2);
    plot(x,SignalDS);
    
    fileID = fopen('DSM_Re_data.txt','w');
    fprintf(fileID,'%12s\n','signal');
    fprintf(fileID,'%4f\n',SignalDS);
    fclose(fileID);