x = linspace(-1 * 10^(-10),1 * 10^(-10),100);
a = 0.7723 * 10^(-8)        % deltaT / deltaL = 2.18 / 28.226 
y = heaviside(x) - 2 * heaviside(x - a * 1.654 * 10^(-3)) + 2 * heaviside(x - a * 4.466 * 10^(-3)) -2 * heaviside(x - a * 7.064 * 10^(-3)) +  2 * heaviside(x - a * 22.266 * 10^(-3)) -  2 * heaviside(x - a * 24.057 * 10^(-3)) +  2 * heaviside(x - a * 26.206 * 10^(-3)) -  2 * heaviside(x - a * 27.258 * 10^(-3)) +  heaviside(x - a * 28.226 * 10^(-3));
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

y = heaviside(x) - 2 * heaviside(x - a * 1.654 * 10^(-3)) + 2 * heaviside(x - a * 4.466 * 10^(-3)) -2 * heaviside(x - a * 7.064 * 10^(-3)) +  2 * heaviside(x - a * 22.266 * 10^(-3)) -  2 * heaviside(x - a * 24.057 * 10^(-3)) +  2 * heaviside(x - a * 26.206 * 10^(-3)) -  2 * heaviside(x - a * 27.258 * 10^(-3)) +  heaviside(x - a * 28.226 * 10^(-3));
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

