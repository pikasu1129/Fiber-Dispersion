function Ellipke



%
%�@�R�����g��łꍇ�́A���Ŏn�߂�B
%�@
%
%
M = 0:0.01:1;
[K,E] = ellipke(M);
plot(M,K,M,E)
grid on
xlabel('M')
title('Complete Elliptic Integrals of First and Second Kind')
legend('First kind','Second kind')
