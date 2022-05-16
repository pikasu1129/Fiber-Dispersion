
%include<stdio.h>
%include<stdlib.h>
%include<math.h>

%‚P•b•ª‚Ìƒf[ƒ^‚ð¶¬
define FREQ   (1000)           %¶¬Žü”g”
define SAMPLE (1000000)        %ƒTƒ“ƒvƒŠƒ“ƒOŽü”g”
define BWIDTH 16

int main( int argc, char *argv[] )

        double  d;
        int     f;
        int     s;
        int     out;
        int     t;
        int     max;

        f = FREQ;
        s = SAMPLE;

        if argc>1
          f = atoi(argv[1]);
        elseif (argc>2)
          s = atoi(argv[2]);
        end
          
        d = 2*M_PI*f/s;
        max = pow(2,BWIDTH-1)-1;

        for (t=0: t<s :t++) 
          out = sin(t*d)*max;
          printf("%d\n",out);
        
        isingen 1000 | head -1000 | gp_wav > isingen_1000.png
        return 0;
