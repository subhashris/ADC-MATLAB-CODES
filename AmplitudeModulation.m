%Amplitude Modulation
clc; close all; clear all;

fm = 5e3; fc = 100e3;
ka = 0.250;
Ac = 1;
Am = input("Enter the amplitude of the message signal");

mu = ka * Am;

fs = 16 * fc;
t = 0:1/fs:1000e-3;

m = Am * cos(2*pi*fm*t);
c = Ac * cos(2*pi*fc*t);
s = (1 + ka*m).*c;

s = awgn(s, 30, 'measured');



[ps, f] = pspectrum(s, fs, 'FrequencyResolution', 100);
figure;
plot(f, 10*log10(ps))
title('Power Spectrum')
xlabel('Frequency(Hz)'); ylabel('Power(dB)')
axis([0 200e3 -70 10])

%y = 1+ka*m;
[yupper,ylower] = envelope(s);
figure;
subplot(4,1,1); plot(t(1:1000), m(1:1000)); 
xlabel('Time(s)');
ylabel('Amplitude(V)');
title('Message signal')
subplot(4,1,2); plot(t(1:1000), c(1:1000));
xlabel('Time(s)'); 
ylabel('Amplitude(V)'); 
title('Carrier signal')
subplot(4,1,3); plot(t(1:1000), s(1:1000)); 
xlabel('Time(s)'); 
ylabel('Amplitude(V)'); 
title('Modulated signal')
subplot(4,1,4); plot(t(1:1000), yupper(1:1000));
xlabel('Time(s)'); 
ylabel('Amplitude(V)'); 
title('Demodulated signal')
