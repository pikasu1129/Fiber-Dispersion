
include<stdio.h>
include<stdlib.h>
include<math.h>

%１秒分のデータを生成
define FREQ   (1000)           %生成周波数
define SAMPLE (1000000)        %サンプリング周波数
define BWIDTH 16

main( int argc, char *argv[] )
{
        double  d;
        int     f,s;
        int     out,t,max;

        f = FREQ;
        s = SAMPLE;

        if (argc>1)
          f = atoi(argv[1]);
        if (argc>2)
          s = atoi(argv[2]);

        d = 2*M_PI*f/s;
        max = (int)pow(2,BWIDTH-1)-1;

        for (t=0: t<s :t++) {
          out = (int)(sin(t*d)*max);
          printf("%d\n",out);
        }
        % isingen 1000 | head -1000 | gp_wav > isingen_1000.png
}