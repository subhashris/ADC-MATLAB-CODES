%Line coding techniques
clc; clear all; close all;
Rb = 1e3;
Nb = 100;
Tb = 1/Rb;
Bo = Rb/2;
Ntow = 3;
Ts = Tb/100; fs = 1/Ts';
N = Tb/Ts;
figure(1);

n = 0:Nb-1; b = randi([1 2], 1, Nb)-1;


%Unipolar NRZ
A= 2
t = (0:N-1)*Ts
v = ones(1,N); a = A*b;
subplot(3,1,1); plot(t,v); axis([-Tb,(Nb+1)*Tb, -0.3, max(v)+0.3]);
xlabel("Time(s)");ylabel("p(t)");title("Basic pulse");

%Polar NRZ
% A = 2
% t=(0:N-1)*Ts
% v=ones(1,N); a=A*b; a(b==0)=-A;
% subplot(3,1,1);plot(t,v);axis([-Tb,(Nb+1)*Tb,-0.3, 1.3])
% xlabel("Time(s)");ylabel("p(t)");title("Basic pulse");

%Manchester NRZ
% t=(0:N-1)*Ts
% v = [ones(1,N/2) -ones(1,N/2)];a = b; a (b==0)=-1;
% subplot(3,1,1); plot(t,v); axis([min(t)-Tb,max(t)+Tb,min(v)-0.3,max(v)+0.3]);
% xlabel("Time(s)");ylabel("p(t)");title("Basic pulse");

%Nyquist pulse
% a=b; a(a==0)=-1;
% t=Ntow*Tb:Ts:Ntow*Tb; v=sinc(2*Bo*t);
% subplot(3,1,1);plot(t,v);axis([min(t)-Tb,max(t)+Tb,min(v)-0.3,max(v)+0.3]);
% xlabel("Time(s)");ylabel("p(t)");title("Basic pulse");

%Raised Cosine Pulse
% alp=input("Enter the value of roll off factor");
% a=b;a(a==0)=-1;
% k=16*Bo^2*alp^2;
% delta=Tb*1e-10;
% t=-Ntow*Tb+delta:Ts:Ntow*Tb;
% v=sinc(2*Bo*t).*cos(2*pi*alp*Bo*t)./(1-16*Bo^2*alp^2*t.^2);
% subplot(3,1,1);plot(t,v);
% axis([min(t)-Tb,max(t)+Tb, min(v)-0.3,max(v)+0.3]);
% xlabel("Time(s)");ylabel("p(t)");title("Basic pulse");

subplot(3,1,2);stem(n*Tb,b);axis([-Tb,(Nb+1)*Tb, -1.3,1.3])
imp_train = zeros(1,(Nb-1)*N+1); imp_train((1:N:end))=a;
xlabel("Time(s)");ylabel("del_T(t)");title("Impulse train");

x=conv(imp_train,v);
t=0:Ts:(length(x)-1)*Ts;
subplot(3,1,3);plot(t,x);axis([min(t)-Tb,max(t)+Tb, min(x)-0.3,max(x)+0.3]);
xlabel("Time(s)");ylabel("x(t)"); title("Transmitted Signal");

[ps,f]=pspectrum(x,fs,'FrequencyResolution',100);
figure;
plot(f,10*log10(ps));axis([0 10*Rb -70 10])
% plot(f,ps);axis([0 10*Rb -0.1 10])
xlabel("Frequency(Hz)"); ylabel("Power(dB)"); title("Power Spectrum (Frequency domain)");

% %Eye diagram
% x=awgn(x,100,'measured')
% subplot(3,1,3);plot(t,x);
% axis([min(t)-Tb,max(t)+Tb,min(x)-0.3,max(x)+0.3])
% figure;hold on;
% x1=x((Ntow+2)N:end-3(Ntow+2));
% t=(0:2*N-1)*Ts
% for ii=Ntow+2:length(x0/N-Ntow-2
%     plot(t,x((ii-1)N+N/2+1:(ii+1)(N)+N/2))
% end
