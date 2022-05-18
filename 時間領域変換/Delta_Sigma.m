function Ellipke



%
%　コメントを打つ場合は、％で始める。
%　
%
%
M = 0:0.01:1;
[K,E] = ellipke(M);
plot(M,K,M,E)
grid on
xlabel('M')
title('Complete Elliptic Integrals of First and Second Kind')
legend('First kind','Second kind')
