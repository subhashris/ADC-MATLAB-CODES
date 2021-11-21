%Frequency Modulation
clc; close all; clear all;

Am = input("Enter the amplitude of the message signal");
fm = 5e3;


kf = 10e3;
delf = kf*Am;

beta = delf/fm;
%beta = 0.1;
disp(beta);
Ac = sqrt(2);
fc = 200e3;

fs = 16*fc;
t = 0:1/fs:1000e-3;

m = Am*cos(2*pi*fm*t);
c = Ac*cos(2*pi*fc*t);
s = Ac*cos(2*pi*fc*t + beta*sin(2*pi*fm*t)); % from single tone modulation

y = s;

[ps, f] = pspectrum(y, fs, 'FrequencyResolution', 100);
figure;
plot(f/1e3, 10*log10(ps)); grid on
title('Power Spectrum')
xlabel('Frequency(Hz)'); ylabel('Power(dB)')
axis([10 400 -60 10])

mout = fmdemod(y, fc, fs, kf*Am);

figure;
subplot(5,1,1); plot(t(1:1000), m(1:1000)); xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Modulating signal')
subplot(5,1,2); plot(t(1:1000), c(1:1000)); xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Carrier signal')
subplot(5,1,3); plot(t(1:1000), s(1:1000)); xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Modulated signal')
subplot(5,1,4); plot(t(1:1000), y(1:1000)); xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Received signal')
subplot(5,1,5); plot(t(1:1000), mout(1:1000)); xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Demodulated signal')
