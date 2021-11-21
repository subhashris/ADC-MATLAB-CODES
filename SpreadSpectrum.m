%Spread Spectrum
clc
close all
clear all
Rb=1e3;     %bit rate
Nb=input("Enter no. of bits: ");    %10,1000
b=2*randi([0 1],Nb,1)-1; %Random binary bits
PNSeq=[0 0 1 1 1 0 1];
LPN=length(PNSeq);
Rc=LPN*Rb;  %chip rate
OSRc=10;    %4,10
OSRb=LPN*OSRc;

fs=OSRb*Rb;
pulsec=ones(1,OSRc);
pulseb=ones(1,OSRb);
PNSeq(PNSeq==0)=-1;
chipseq=zeros(1,(LPN-1)*OSRc+1); chipseq(1:OSRc:end)=PNSeq; chipseq=conv(chipseq,pulsec);
tc=(0:LPN*OSRc-1)/(Rc*OSRc);


bitseq1=zeros(1,(Nb-1)*OSRb+1); bitseq1(1:OSRb:end)=b;
bitseq_raw=conv(bitseq1,pulseb);
tb=(0:Nb*OSRb-1)/(Rb*OSRb);
c=[];
for ii=1:Nb
    c=[c chipseq];
end

xseq=bitseq_raw.*c;
y = xseq.*c;
%xseq=conv(bitseq1,chipseq);
figure;subplot(4,1,1);plot(tb,bitseq_raw); axis([0 max(tb) -1.1 1.1]); 
xlabel("Time(s)"); ylabel("Amplitude b(t)"); title("Message signal");
subplot(4,1,2);plot(tb,c); axis([0 max(tb) -1.1 1.1]);
xlabel("Time(s)"); ylabel("Amplitude c(t)"); title("Pseudo-Random Noise(PRN) signal");
subplot(4,1,3);plot(tb,xseq); axis([0 max(tb) -1.1 1.1]);
xlabel("Time(s)"); ylabel("Amplitude x(t)"); title("Spreaded signal(Transmitted)");
subplot(4,1,4);plot(tb,y); axis([0 max(tb) -1.1 1.1]);
xlabel("Time(s)"); ylabel("Amplitude y(t)"); title("Received signal");

[pb,f]=pwelch(bitseq_raw,[],[],[],fs);
[px,f]=pwelch(xseq,[],[],[],fs);
figure; subplot(2,1,1);plot(f,10*log10(pb));
xlabel("Frequency(Hz)"); ylabel("Power(dB)"); title("Frequency spectrum of Message signal");
subplot(2,1,2);plot(f,10*log10(px));
xlabel("Frequency(Hz)"); ylabel("Power(dB)"); title("Frequency spectrum of Spreaded signal");

