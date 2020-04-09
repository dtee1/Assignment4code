clear;
close all;
% David Talson
% Student number: 101022690

global C G b
G = zeros(5,5);
C = zeros(5,5);
b = zeros(5,1);

res(1,2,1);
cap(1,2,0.25);
res(2,0,2);
res(3,0,44.23);
res(4,5,0.1);
res(5,0,1000);
ind(2,3,0.2);
vol(1,0,1);
vcvs(4,0,3,0,(100/44.23)+(200/44.23^2)+(300/44.23^3));

Vin = linspace(-10, 10, 100);

for i = 1:100
    b(7) = Vin(i);
    nodeVoltages = G\b;
    vnode3(i) = nodeVoltages(3);
    vOut(i) = nodeVoltages(5);
    gain(i) = vOut(i)/Vin(i);
end

figure(1)
plot(Vin, vnode3,'LineWidth',2);
hold on
plot(Vin, vOut,'LineWidth',2);
xlabel('Input Voltage')
ylabel('Voltage')
legend('Voltage at node 3', 'Output voltage (Vo)');
grid on
hold off

F = logspace(0,8,5000);

for i=1:length(F)
    w = 2*pi*F(i);
    s = 1i*F(i);
    A = G + (s*C);   

    nodeVoltages = A\b;
    Vout(i) = nodeVoltages(5);
    gain(i) = 20*log(abs(nodeVoltages(5)));
end

figure(2);
semilogx(F, Vout,'LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('Output voltage (V)');
title('Frequency Response')
grid on

figure(3);
semilogx(F, gain,'LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
title('Gain');
grid on

cnew = 0.25 + 0.05*randn(1,1000);
for i = 1:1000
    s = 1i*pi;
    
    C(1,1) = cnew(i);
    C(2,2) = cnew(i);
    C(1,2) = cnew(i)*-1;
    C(2,1) = cnew(i)*-1;
    c(i)=cnew(i);
    A = G +(s*C);
    nodeVoltages = A\b;
    gainnew(i) = 20*log10(abs(nodeVoltages(5)));
end

figure(4)
histogram(gainnew);
title('Gain Histogram')
grid on

Vin1 = 0;
flag = 1;
time = 0;
for i = 10:10:1000
    if i >= 30 & flag
       Vin1 = 1;
       flag = 0;
    end  
    index = round(i/10);
   
    Vinnew(index)=Vin1;

    time(index) = i/1000;
    b(7) = Vin1;
    w = 2*pi*(1/(i/1000));
    s = 1i*w;
    A = G + (s*C);   

    nodeVoltages = A\b;
    Voutnew(index) = nodeVoltages(5);  
    
end

figure(5);
plot(time,Voutnew,'LineWidth',2);
hold on
plot(time,Vinnew,'LineWidth',2)
xlabel('Time (ms)');
ylabel('Voltage(V)');
title('AC voltage response for step function');
legend('Output voltage','Input voltage')
grid on
hold off

f = 0.1*(-(100/2):(100/2-1))./100;
VinFft = abs(fftshift(fft(Vinnew)));
VoutFft = abs(fftshift(fft(Voutnew)));

figure(6)
plot(f,VoutFft,'LineWidth',2)
hold on
plot(f,VinFft,'LineWidth',2)
title('Fourier transform of step function')
xlim([-0.04 0.04])
grid on 
hold off

Vin2 = 0;

flag = 1;
count = 1;
time = 0;
for i=10:10:1000
    
    Vin2 = sin(2*pi*(1/0.1)*(i/1000));
    index = round(i/10);
    Vinnew(index)=Vin2;
    if i<10
        index = 1;
    end
    time(index) = i/1000;
    b(7) = Vin2;
    w = 2*pi*(1/0.1);
    s = 1i*w;
    A = G + (s*C);   

    nodeVoltages = A\b;
    Voutnew(index) = nodeVoltages(5); 
    
end
figure(7);
plot(time,Voutnew,'LineWidth',2);
hold on
plot(time,Vinnew,'LineWidth',2)
xlabel('Time (ms)');
ylabel('Voltage(V)');
title('AC voltage response for sine signal')
legend('Output voltage','Input voltage')
grid on
hold off

f = 0.1*(-(100/2):(100/2-1))./100;
VinFft = abs(fftshift(fft(Vinnew)));
VoutFft = abs(fftshift(fft(Voutnew)));

figure(8)
plot(f,VoutFft,'LineWidth',2)
hold on
plot(f,VinFft,'LineWidth',2)
title('Fourier transform of sine signal')
legend('Output','Input')
xlim([-0.04 0.04])
grid on
hold off

Vin3 = 0;
time = 0;
flag = 1;
for i=1:1000
    
   
    Vin3 =  exp((-((i/1000) - 0.06)^2)/(2*(0.03^2)));
    
     
    index=round(i);
    
    Vinnew(index)=Vin3;
    time(index) = i/1000;
    b(7) = Vin3;
    b(8) = (100/44.23)+(200*Vin3^2/44.23^2)+(300*Vin3^3/44.23^3);
    w = 2*pi*(1/(i/1000));
    s = 1i*w;
    A = G + (s*C);   

    nodeVoltages = A\b;
    Voutnew(index) = abs(nodeVoltages(5));    
end

figure(9);
plot(Voutnew,'LineWidth',2);
hold on
plot(Vinnew,'LineWidth',2)
xlabel('Time (ms)');
ylabel('Voltage(V)');
title('AC voltage response for gaussian input pulse')
legend('Output','Input')
grid on
hold off

f = (-(1000/2):(1000/2-1))./1000;
VinFft = abs(fftshift(fft(Vinnew)));
VoutFft = abs(fftshift(fft(Voutnew)));

figure(10)
plot(f,VoutFft,'LineWidth',2)
hold on
plot(f,VinFft,'LineWidth',2)
grid on 
title('Fourier transform of gaussian input pulse and output')
xlim([-0.1 0.1])
legend('Output','Input')
hold off

cap(3,0,0.00001);
cur(3,0,0.01);

flag = 1;
for i=1:1000
    
    Vin3 =  exp((-((i/1000) - 0.06)^2)/(2*(0.03^2)));    
    
    index=round(i);
    Vinnew(index)=Vin3;
    time(index) = i/1000;
    
    currNoise = 0.01*normrnd(0,1);
    b(3) = currNoise;
    b(7) = Vin3;
    w = 2*pi*(1/(i/1000));
    s = 1i*w;
    A = G + (s*C);   

    nodeVoltages = A\b;
    Voutnew(index) = abs(nodeVoltages(5));    
end
figure(11)
plot(Voutnew,'LineWidth',2);
hold on
plot(Vinnew,'LineWidth',2)
xlabel('Time (ms)');
ylabel('Voltage(V)');
title('AC voltage response for gaussian input pulse')
legend('Output','Input')
grid on
hold off

f = (-(1000/2):(1000/2-1))./1000;
VinFft = abs(fftshift(fft(Vinnew)));
VoutFft = abs(fftshift(fft(Voutnew)));

figure(12)
plot(f,VoutFft,'LineWidth',2)
hold on
plot(f,VinFft,'LineWidth',2)
grid on 
title('Fourier transform of gaussian input pulse and output')
xlim([-0.1 0.1])
legend('Output','Input')
hold off