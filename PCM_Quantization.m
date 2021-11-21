%Quantization in PCM
%Sinusoidal signal as input
clc;
clear all;
close all;
fm = 1e3;
fs = 1024*fm;
Am = 2;
n = input('Enter the no. of bits');

t = 0:1/fs:2/fm;
x = Am*cos(2*pi*fm*t);
N = length(t);
x_power = x*x'/ N;

for i = 1:length(n)
    delta = 2*Am/(2^n(i)-1);
    partition = -Am+delta/2:delta:Am-delta/2;
    codebook = -Am:delta:Am;

    [index,xq] = quantiz(x, partition, codebook);
    qe = x - xq;
    qe_avg = mean(qe);
    qe_power = qe*qe'/N;
    SQNR = x_power/qe_power;
    SQNRdB(i) = 10*log10(SQNR);
    SQNRdB_th(i) = 6*n(i)+1.72;
end

plot(t,x,'b-',t,xq,'r--'); xlabel('Time(s)'); ylabel('Amplitude(V)');legend('input signal x(t)', 'quantised signal x_q(t)'); title("Time Domain Representation");
figure; plot(n,SQNRdB,'b-',n, SQNRdB_th,'r-*'); xlabel('Number of bits(n)'); ylabel('SQNR(dB)'); title("No. of bits vs SQNR(dB)");

%Random signal as input(Gaussian Distribution)
% t = 0:1/fs:2/fm;
% x = zeros(1, length(t));
% sigx2 = 4;
% LF = 4;
% sigx = sqrt(sigx2);
% Am = LF*sigx;
% x = awgn(x, sigx2);
% x_power = x*x';
% 
% for j = 1:length(n)
%     delta = 2*Am/(2^n(j)-1);
%     partition = -Am+delta/2:delta:Am-delta/2;
%     codebook = -Am:delta:Am;
% 
%     [index,xq] = quantiz(x, partition, codebook);
%     qe = x - xq;
%     qe_avg = mean(qe);
%     qe_power = qe*qe';
%     SQNR = x_power/qe_power;
%     SQNRdB(j) = 10*log10(SQNR);
%     SQNRdB_th(j) = 6*n(j)-7.2;
% end
% 
% plot(t,x,'b-',t,xq,'r--'); xlabel('Time(s)'); ylabel('Amplitude(V)'); legend('input signal x(t)', 'quantised signal x_q(t)');  title("Time Domain Representation");
% figure; plot(n,SQNRdB,'b-',n, SQNRdB_th,'r-*'); xlabel('Number of bits(n)'); ylabel('SQNR(dB)');  title("No. of bits vs SQNR(dB)");

