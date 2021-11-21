
%Sampling
% clc; clear all; close all;
fo = 2000; fc = 2*fo; 
Fs = input("Enter the sampling frequency");
Nord = input("Enter the order of the filter");

delf = fo/20;
To = 1/fo;

Ts = 1/Fs;
Fsdas = 8*Fs;
Dutycycle = 20/100;
T1 = (1/Fsdas)*Dutycycle;
No = Fsdas/Fs;
N1 = round(No*Dutycycle);

t = 0:1/Fsdas:5000*To;
N = length(t);

Vm = 2; x = Vm*cos(2*pi*fo*t)+4;
pulse = [ones(1,N1) zeros(1,No-N1)];

%Natural Sampling
deltrain = zeros(1,N-No+1); deltrain(1:No:end) = 1;
p = conv(deltrain,pulse);
v = p.*x;

[b,a] = butter(Nord,fc/(Fsdas/2));
[H,w] = freqz(b,a, Fsdas, Fsdas);
y = filter(b,a,v);

figure;
subplot(4,1,1); plot(t,x); axis([min(t) max(t), min(x)-0.1 max(x)+0.1]);
xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Analog Message Signal');
subplot(4,1,2); plot(p); %axis([min(t) max(t), min(x)-0.1 max(p)+0.1]);
xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Sampling Pulse Train');
subplot(4,1,3); plot(v); %axis([min(t) max(t), min(x)-0.1 max(v)+0.1]);
xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Sampled Signal');
subplot(4,1,4); plot(y); %axis([min(t) max(t), min(x)-0.1 max(y)+0.1]);
xlabel('Time(s)'); ylabel('Amplitude(V)'); title('Reconstructed Signal');


[px,fx] = pspectrum(x, Fsdas, 'FrequencyResolution', delf);
[pv,fv] = pspectrum(v, Fsdas, 'FrequencyResolution', delf);
[py,fy] = pspectrum(y, Fsdas, 'FrequencyResolution', delf);
figure;
subplot(4,1,1); plot(fx,10*log10(px)); axis([-1 4*Fs, -50 max(10*log10(px)+1)]); grid on
xlabel('Frequency(Hz)'); ylabel('Power(dB)'); title('Analog Message Signal');
subplot(4,1,2); semilogx(w,20*log10(abs(H))); axis([-1 4*Fs, -50 max(abs(H)+1)]); grid on
xlabel('Frequency(Hz)'); ylabel('Power(dB)'); title('Sampling Pulse Train');
subplot(4,1,3); plot(fv,10*log10(pv)); axis([-1 3*Fs, -50 max(10*log10(pv)+1)]); grid on
xlabel('Frequency(Hz)'); ylabel('Power(dB)'); title('Sampled Signal');
subplot(4,1,4); plot(fy,10*log10(py)); axis([-1 3*Fs, -50 max(10*log10(py)+1)]); grid on
xlabel('Frequency(Hz)'); ylabel('Power(dB)'); title('Reconstructed Signal');
