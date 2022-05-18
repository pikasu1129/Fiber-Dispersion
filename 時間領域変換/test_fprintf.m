x = 0:0.1:1;
A = [x; exp(x)];

fileID = fopen('exp.txt','w');
fprintf(fileID,'%6s %12s\r\n','x','exp(x)');
fprintf(fileID,'%6.2f %12.8f\r\n',A);
fclose(fileID);