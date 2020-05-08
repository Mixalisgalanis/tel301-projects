clear all;
close all;
 
T=0.01;
over=10;
Ts=T/over;
Fs=1/Ts;
A=4;
a=0.5;
Nf=2048;
 
%??????? ?1
 
%?????????? ?????? SRRC ?(t)
[phi, t] = srrc_pulse(T, Ts, A, a); 
 
%??????????? ??????????????? Fourier
y=abs(fftshift(fft(phi,Nf)*Ts));
 
%???????? ??????????, ?????????? ????????? ????????
F = [-Fs/2:Fs/Nf:Fs/2-Fs/Nf];
 
%???????? ?????????? ?????????? ?????????
figure;
semilogy(F,y.^2,'b')
xlabel('Frequency (Hz)');
title('????????? ????????? ????????? (semilogy)')
legend('|?(F)|^2')
grid on;
 
%??????? ?2
N=50;
 
%?????????? N ??????????? ??? ?????????? bits
b = (sign(randn(N, 1)) + 1)/2;
 
%????????? ??? bits ?? 2-PAM ???????
X=bits_to_2PAM(b);
 
%???????? over-1 ????????? ?????? ??? ????????
X_delta = 1/Ts*upsample(X, over);
 
%??????? ????? ?????? ??? ????????? 
tx = [-A*T:Ts:(A*T+N*T-Ts)];
 
%???????? ??? X_delta ?? ?? ?(t)
x=conv(X_delta,phi)*Ts;
 
%??????? ?3
Px=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
 
figure;
subplot(2,1,1);
plot(F,Px)
grid on;
xlabel('Frequency (Hz)');
title('????????????? ??? X(t)(plot)');
 
subplot(2,1,2);
semilogy(F,Px)
grid on;
xlabel('Frequency (Hz)');
title('????????????? ??? X(t)(semilogy)');
 
K=100;
Px_new=zeros(K,Nf);
%??????????????? ??? ???????? ?????????? ??? ? ???????????
for i=1:K 
    
    %?????????? N ??????????? ??? ?????????? bits
    b = (sign(randn(N, 1)) + 1)/2;
    
    %????????? ??? bits ?? 2-PAM ???????
    X=bits_to_2PAM(b);
 
    %???????? over-1 ????????? ?????? ??? ????????
    X_delta = 1/Ts*upsample(X, over);
    
    %??????? ????? ?????? ??? ????????? 
    tx = [-A*T:Ts:(A*T+N*T-Ts)];
 
    %???????? ??? X_delta ?? ?? ?(t)
    x=conv(X_delta,phi)*Ts;
 
    Px_new(i,:)=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
end;
 
%??????????? ???????????? ?????????? ?????????? ??????
Px_exp=sum(Px_new)./K;
 
%??????????? ?????????? ?????????? ?????????? ?????? ??? 
%??? ?????? ???? 
Px_th=(var(X)^2/T).*(y.^2);
 
figure;
semilogy(F,Px_exp,'b')
hold on;
semilogy(F,Px_th,'r')
xlabel('Frequency (Hz)');
title('????????? ????????? ??????, ??? ?=100 ??? ?=50');
legend('???????????','?????????');
 
N=100;
K=1000;
Px_new=zeros(K,Nf);
%??????????????? ??? ???????? ?????????? ??? ? ???????????
for i=1:K 
    
    %?????????? N ??????????? ??? ?????????? bits
    b = (sign(randn(N, 1)) + 1)/2;
    
    %????????? ??? bits ?? 2-PAM ???????
    X=bits_to_2PAM(b);
 
    %???????? over-1 ????????? ?????? ??? ????????
    X_delta = 1/Ts*upsample(X, over);
    
    %??????? ????? ?????? ??? ????????? 
    tx = [-A*T:Ts:(A*T+N*T-Ts)];
 
    %???????? ??? X_delta ?? ?? ?(t)
    x=conv(X_delta,phi)*Ts;
 
    Px_new(i,:)=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
end;
 
%??????????? ???????????? ?????????? ?????????? ??????
Px_exp=sum(Px_new)./K;
 
figure;
semilogy(F,Px_exp,'b')
hold on;
semilogy(F,Px_th,'r')
xlabel('Frequency (Hz)');
title('????????? ????????? ??????, ??? ?=1000 ??? ?=100');
legend('???????????','?????????');
 
%??????? ?4
 
N=50;
 
%?????????? N ??????????? ??? ?????????? bits
b = (sign(randn(N, 1)) + 1)/2;
 
%????????? ??? bits ?? 4-PAM ???????
X=bits_to_4PAM(b);
 
%???????? over-1 ????????? ?????? ??? ????????
X_delta = 1/Ts*upsample(X, over);
 
%??????? ????? ?????? ??? ????????? 
tx = [-A*T:Ts:(A*T+N*T-Ts)];
 
%???????? ??? X_delta2 ?? ?? ?(t)
x=conv(X_delta,phi)*Ts;
 
Px=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
 
K=1000;
Px_new=zeros(K,Nf);
%??????????????? ??? ???????? ?????????? ??? ? ???????????
for i=1:K 
    
    %?????????? N ??????????? ??? ?????????? bits
    b = (sign(randn(N, 1)) + 1)/2;
    
    %????????? ??? bits ?? 2-PAM ???????
    X=bits_to_4PAM(b);
 
    %???????? over-1 ????????? ?????? ??? ????????
    X_delta = 1/Ts*upsample(X, over);
    
    %??????? ????? ?????? ??? ????????? 
    tx = [-A*T:Ts:(A*T+N*T-Ts)];
 
    %???????? ??? X_delta ?? ?? ?(t)
    x=conv(X_delta,phi)*Ts;
 
    Px_new(i,:)=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
end;
 
%??????????? ?????????? ?????????? ?????????? ?????? ??? 
%??? ?????? ???? 
Px_th=(5/T).*(y.^2);
 
%??????????? ???????????? ?????????? ?????????? ??????
Px_exp=sum(Px_new)./K;
 
figure;
semilogy(F,Px_exp)
hold on;
semilogy(F,Px_th)
xlabel('F (Hz)');
title(['A.4 - 4-PAM Periodogram (Logarithmic) with K =  ' num2str(K)]);
legend('Theoretical P_X','Estimated P_X');
 
%??????? ?5
T=0.02;
over=20;
Ts=T/over;
Fs=1/Ts;
 
%?????????? ?????? SRRC ?(t)
[phi, t] = srrc_pulse(T, Ts, A, a); 
 
%??????????? ??????????????? Fourier
y=abs(fftshift(fft(phi,Nf)*Ts));
 
%???????? ??????????, ?????????? ????????? ????????
F = [-Fs/2:Fs/Nf:Fs/2-Fs/Nf];
 
N=50;
 
%?????????? N ??????????? ??? ?????????? bits
b = (sign(randn(N, 1)) + 1)/2;
 
%????????? ??? bits ?? 2-PAM ???????
X=bits_to_2PAM(b);
 
%???????? over-1 ????????? ?????? ??? ????????
X_delta = 1/Ts*upsample(X, over);
 
%??????? ????? ?????? ??? ????????? 
tx = [-A*T:Ts:(A*T+N*T-Ts)];
 
%???????? ??? X_delta ?? ?? ?(t)
x=conv(X_delta,phi)*Ts;
 
Px=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
 
figure;
subplot(2,1,1);
plot(F,Px)
grid on;
xlabel('Frequency (Hz)');
title('????????????? ??? X(t)(plot)');
 
subplot(2,1,2);
semilogy(F,Px)
grid on;
xlabel('Frequency (Hz)');
title('????????????? ??? X(t)(semilogy)');
 
K=100;
Px_new=zeros(K,Nf);
%??????????????? ??? ???????? ?????????? ??? ? ???????????
for i=1:K 
    
    %?????????? N ??????????? ??? ?????????? bits
    b = (sign(randn(N, 1)) + 1)/2;
    
    %????????? ??? bits ?? 2-PAM ???????
    X=bits_to_2PAM(b);
 
    %???????? over-1 ????????? ?????? ??? ????????
    X_delta = 1/Ts*upsample(X, over);
    
    %??????? ????? ?????? ??? ????????? 
    tx = [-A*T:Ts:(A*T+N*T-Ts)];
 
    %???????? ??? X_delta ?? ?? ?(t)
    x=conv(X_delta,phi)*Ts;
 
    Px_new(i,:)=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
end;
 
%??????????? ???????????? ?????????? ?????????? ??????
Px_exp=sum(Px_new)./K;
 
%??????????? ?????????? ?????????? ?????????? ?????? ??? 
%??? ?????? ???? 
Px_th=(var(X)^2/T).*(y.^2);
 
figure;
semilogy(F,Px_exp,'b')
hold on;
semilogy(F,Px_th,'r')
xlabel('Frequency (Hz)');
title('????????? ????????? ??????, ??? ?=100 ??? ?=50');
legend('???????????','?????????');
 
K=1000;
Px_new=zeros(K,Nf);
%??????????????? ??? ???????? ?????????? ??? ? ???????????
for i=1:K 
    
    %?????????? N ??????????? ??? ?????????? bits
    b = (sign(randn(N, 1)) + 1)/2;
    
    %????????? ??? bits ?? 2-PAM ???????
    X=bits_to_2PAM(b);
 
    %???????? over-1 ????????? ?????? ??? ????????
    X_delta = 1/Ts*upsample(X, over);
    
    %??????? ????? ?????? ??? ????????? 
    tx = [-A*T:Ts:(A*T+N*T-Ts)];
 
    %???????? ??? X_delta ?? ?? ?(t)
    x=conv(X_delta,phi)*Ts;
 
    Px_new(i,:)=(abs(fftshift(fft(x,Nf)*Ts)).^2)/(length(tx)*Ts);
end;
 
%??????????? ???????????? ?????????? ?????????? ??????
Px_exp=sum(Px_new)./K;
 
figure;
semilogy(F,Px_exp,'b')
hold on;
semilogy(F,Px_th,'r')
xlabel('Frequency (Hz)');
title('????????? ????????? ??????, ??? ?=1000 ??? ?=100');
legend('???????????','?????????');
