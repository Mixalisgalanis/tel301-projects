%C1

a=0.5;
A=5;
T=0.1;
over=10;
Ts=T/over;
N=50;

%dimiourgia N tixaiwn bits
b = (sign(randn(N, 1)) + 1)/2; 

%================================================

%C2

%metatropi twn bits se 2-PAM simvola
X=bits_to_2PAM(b);

%prosthiki over-1 midenikwn metaksi twn simvolwn
X_delta = 1/Ts*upsample(X, over);

%orsimos tou aksona tou xronou
t=0:Ts:((N+N*(over-1))-1)*Ts;

figure(1)
plot(t,X_delta)
xlabel('t(s)')
ylabel('X-delta')
title('C - X-delta')

%dimiourgia palmou Ö(t)
[phi,t1]=srrc_pulse(T, Ts, A, a);
%orismos aksona xronou tis sineliksis
tx = [t(1)+t1(1):Ts:t(end)+t1(end)];
%sineliksi tou X_delta me to ö(t)
x=conv(X_delta,phi)*Ts;

figure(2)
plot(tx,x)
xlabel('t(s)')
title('C - X(t): Convolution of X-delta & SRRC Pulse')
ylabel('Conv')
%dimiourgia tou ö(-t) palmou
phi1=fliplr(phi);
t2=-fliplr(t1);
%sineliksi tou x me to ö(-t)
z=conv(x,phi1)*Ts;
%orismos aksona xronou tis sineliksis
tz= [tx(1)+t2(1):Ts:tx(end)+t2(end)];

figure(3)
plot(tz,z)
xlabel('t(s)');
ylabel('Conv');
title('C - Z(t): Convolution of X & Phi(-t)');

figure(4)
plot(tz,z)
hold on;
%emfanisi se koino diagramma to Z(t) me to X
stem((0:N-1)*T,X,'filled');
xlabel('t(s)')
legend('Z(kT)', 'X(kT)')

%================================================