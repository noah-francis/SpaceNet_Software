%%
clc;clear;close all;
%%
[y,Fs]=audioread('SDRSharp_20201111_092402Z_106300000Hz_IQ.wav');
L=length(y);
t=linspace(0,L/Fs,L);

subplot(2,1,1)
plot(t,y(:,1),'b')
grid on
ylabel('I')
subplot(2,1,2)
plot(t,y(:,2),'r')
grid on
ylabel('Q')
xlabel('Time [s]')

%% FFT BABYYYYY
I = fft(y(:,1));
P2 = abs(I/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L+106300000;
figure;
semilogy(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%%
I = fft(y(:,2));
P2 = abs(I/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L+106300000;
figure;
semilogy(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')