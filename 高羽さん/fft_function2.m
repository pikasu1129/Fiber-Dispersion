x = linspace(-1 * 10^(-10),1 * 10^(-10),100);

y = heaviside(x) - 2 * heaviside(x - 5.51 * 10^(-12)) + 2 * heaviside(x - 1.49 * 10^(-11)) -2 * heaviside(x - 2.35 * 10^(-11)) +  2 * heaviside(x - 7.42 * 10^(-11)) -  2 * heaviside(x - 8.02 * 10^(-11)) +  2 * heaviside(x - 8.74 * 10^(-11)) -  2 * heaviside(x - 9.09 * 10^(-11)) +  heaviside(x - 9.41 * 10^(-11));
figure(1);
plot(x, y)
Y = fft(y,10000);
plot(Y)
disp(Y);

real_Y = real(Y);
plot(real_Y);

figure(2);
imag_Y = imag(Y);
plot(imag_Y);

y =( heaviside(x) - 2 * heaviside(x - 5.51 * 10^(-12)) + 2 * heaviside(x - 1.49 * 10^(-11)) -2 * heaviside(x - 2.35 * 10^(-11)) +  2 * heaviside(x - 7.42 * 10^(-11)) -  2 * heaviside(x - 8.02 * 10^(-11)) +  2 * heaviside(x - 8.74 * 10^(-11)) -  2 * heaviside(x - 9.09 * 10^(-11)) +  heaviside(x - 9.41 * 10^(-11)));
figure(3);
plot(x, y)
Y = fft(y,10000);
plot(Y)
disp(Y);

real_Y = real(Y);
plot(real_Y);


hold on
imag_Y = imag(Y);
plot(imag_Y);
hold off


y =-( heaviside(x) - 2 * heaviside(x - 5.51 * 10^(-12)) + 2 * heaviside(x - 1.49 * 10^(-11)) -2 * heaviside(x - 2.35 * 10^(-11)) +  2 * heaviside(x - 7.42 * 10^(-11)) -  2 * heaviside(x - 8.02 * 10^(-11)) +  2 * heaviside(x - 8.74 * 10^(-11)) -  2 * heaviside(x - 9.09 * 10^(-11)) +  heaviside(x - 9.41 * 10^(-11)));
figure(4);
plot(x, y)
Y = fft(y,10000);
plot(Y)
disp(Y);

real_Y = real(Y);
plot(real_Y);

hold on
imag_Y = imag(Y);
plot(imag_Y);
hold off