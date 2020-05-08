%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC INFORMATION                                   %
% Course: Telecommunication Systems I - Excersize 3   %
% Deadline: 20-Dec-18                                 %
% FullName: Mixalis Galanis                           %
% Academic ID: 2016030036                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing things up
close all
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful Functions                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting Preperation Functions
set_diag_labels = @(tit, xlab, ylab){title(tit); xlabel(xlab); ylabel(ylab);};
set_diag_layout = @(fig, rows, cols, selected){figure(fig); subplot(rows,cols,selected)};
%Some Useful Plotting Functions (continuous & discrete)
c_plot = @(t, f){ plot(t,f);};
c_semilogy = @(t, f){ semilogy(t,f);};
d_plot = @(t, f){ stem(t,f);};
%Some Useful Signal Processing Functions
t_conv = @(t1, t2, dt) (t1(1) + t2(1)):dt:(t1(end) + t2(end));  %Convolution time
gen_N_bits = @(N) ((sign(randn(N,1)) + 1)/2);                   %Generates N Random Bits
F_T = @(X,Nf,Ts) (abs(fftshift(fft(X,Nf)*Ts)));                 %Fourier Transform of X
P_X = @(t_X, X,Nf,Ts) ((F_T(X,Nf,Ts).^2)/(length(t_X)*Ts));     %Periodogram of X


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A1 - Creating N Bits Sequence                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 50;
bit_seq = gen_N_bits(3*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A2 - Generating 8-PSK bits                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = 1; %TODO Figure out A's value
bit_8_psk = bits_to_PSK_8(bit_seq, A); %N length



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A3 - Filtering Xi & Xq with SRRC Pulses                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Settings for plotting the diagrams
r = 2; c = 2; fig = 1; sub_fig = 1;

Nf = 2^11;
DT = 1;

Xi = squeeze(bit_8_psk(:,1));
t_Xi = 0:length(Xi);
Xq = squeeze(bit_8_psk(:,2));
t_Xq = 0:length(Xq);

%Building up parameters for srrc_pulse function
T = 0.01;       %Nyquist parameter (>0)
over = 10;      %Oversampling factor (>0)
Ts = T/over;    %Sampling period (>0)
Fs = 1/Ts;      %Sampling frequency (>0)
Hd = 4;         %Half duration of the pulse (>0)
a = 0.5;        %Roll-off factor(0<a<1)

%Calling srrc_pulse function to store phi and t variables
[phi, t_phi] = srrc_pulse(T, Ts, Hd, a);

%Convoluting Xi, Xq with phi
Xi_delta = Fs * upsample(Xi,over);            %Inserts zeros in between bits
Xq_delta = Fs * upsample(Xq,over);            %Inserts zeros in between bits

Xi_conv = conv(Xi_delta, phi)*Ts;
Xq_conv = conv(Xq_delta, phi)*Ts;

t_Xi_conv = t_conv(t_phi, t_Xi, Ts);
t_Xq_conv = t_conv(t_phi, t_Xq, Ts);

%Calculating Periodgrams
F_ES = (-Fs/2) :(Fs/Nf): (Fs/2 - Fs/Nf);
PXi = P_X(t_Xi_conv, Xi_conv,Nf,Ts);
PXq = P_X(t_Xq_conv, Xq_conv,Nf,Ts);

%Plotting output diagrams
%Xi_conv
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_Xi_conv,Xi_conv);
set_diag_labels("A.3 - Convolution of Xi,n with SRRC pulse","t(s)","");
grid on;
%Xq_conv
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_Xq_conv,Xq_conv);
set_diag_labels("A.3 - Convolution of Xq,n with SRRC pulse","t(s)","");
grid on;
%PXi
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PXi);
set_diag_labels("A.3 - Periodgram of Xi,n (Logarithmic)","F [Hz]","");
grid on;
%PXq
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PXq);
set_diag_labels("A.3 - Periodgram of Xq,n (Logarithmic)","F [Hz]","");
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A4 - Multiplication of Xi, Xq by carrier                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = fig + 1; sub_fig = 1; %New Figure
F0 = 200;

%Calculating Multiplication by Carrier
Xi_t = Xi_conv .* (2*cos(2*pi*F0*t_Xi_conv));
Xq_t = Xq_conv .* ((-2)*sin(2*pi*F0*t_Xq_conv));

%Calculating Periodgrams
PXi_t = P_X(t_Xi_conv, Xi_t,Nf,Ts);
PXq_t = P_X(t_Xq_conv, Xq_t,Nf,Ts);

%Plotting output diagrams
%Xi_t
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_Xi_conv,Xi_t);
set_diag_labels("A.4 - Xi_t (multiplied by carrier)","t(s)","");
grid on;
%Xq_t
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_Xq_conv,Xq_t);
set_diag_labels("A.4 - Xq_t (multiplied by carrier)","t(s)","");
grid on;
%PXi
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PXi_t);
set_diag_labels("A.4 - Periodgram of Xi_t,n (Logarithmic)","F [Hz]","");
grid on;
%PXq
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PXq_t);
set_diag_labels("A.4 - Periodgram of Xq_t,n (Logarithmic)","F [Hz]","");
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A5 - Multiplication of Xi, Xq by carrier                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = fig + 1; sub_fig = 1; r = 2; c = 2; %New Figure

X = Xi_t + Xq_t;

PX = P_X(t_Xi_conv, X,Nf,Ts);

%Plotting output diagrams
%X
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_Xi_conv,X);
set_diag_labels("A.5 - Xi_t (multiplied by carrier)","t(s)","");
grid on;
%PX
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PX);
set_diag_labels("A.5 - Periodgram of Xi_t,n (Logarithmic)","F [Hz]","");
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A7 - White Gassuan Noise                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR = 10; %Setting up signal to noise ratio
sw = Fs * (10^(-(SNR/10)));
wn= randn(1,length(X)) .* sqrt(sw);
X_wn = X + wn;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A8 - Multiplication of X by carrier                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = fig + 1; sub_fig = 1; r = 2; c = 2; %New Figure

XI = X .* cos(2*pi*F0*t_Xi_conv);
XQ = X .* cos(2*pi*F0*t_Xq_conv);

%Calculating Periodgrams
PXI = P_X(t_Xi_conv, XI,Nf,Ts);
PXQ = P_X(t_Xq_conv, XQ,Nf,Ts);

%Plotting output diagrams
%XI
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_Xi_conv,XI);
set_diag_labels("A.8 - XI (multiplied by carrier)","t(s)","");
grid on;
%XQ
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_Xq_conv,XQ);
set_diag_labels("A.8 - XQ (multiplied by carrier)","t(s)","");
grid on;
%PXI
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PXI);
set_diag_labels("A.8 - Periodgram of XI (Logarithmic)","F [Hz]","");
grid on;
%PXQ
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PXQ);
set_diag_labels("A.8 - Periodgram of XQ (Logarithmic)","F [Hz]","");
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A9 - Filtering XI & XQ with SRRC Pulses                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = fig + 1; sub_fig = 1; r = 2; c = 2; %New Figure

XI_conv = conv(XI, phi)*Ts;
XQ_conv = conv(XQ, phi)*Ts;

t_XI_conv = t_conv(t_phi, t_Xi, Ts);
t_XQ_conv = t_conv(t_phi, t_Xq, Ts);

%Calculating Periodgrams
F_ES = (-Fs/2) :(Fs/Nf): (Fs/2 - Fs/Nf);
PXI_2 = P_X(t_XI_conv, XI_conv,Nf,Ts);
PXQ_2 = P_X(t_XQ_conv, XQ_conv,Nf,Ts);

%Plotting output diagrams
%XI_conv
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_XI_conv,XI_conv);
set_diag_labels("A.9 - Convolution of XI with SRRC pulse","t(s)","");
grid on;
%XQ_conv
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_XQ_conv,XQ_conv);
set_diag_labels("A.9 - Convolution of XQ with SRRC pulse","t(s)","");
grid on;
%PXI_2
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PXI_2);
set_diag_labels("A.9 - Periodgram of XI Filtered (Logarithmic)","F [Hz]","");
grid on;
%PXQ_2
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,PXQ_2);
set_diag_labels("A.9 - Periodgram of XQ Filtered (Logarithmic)","F [Hz]","");
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A10 - Using Scatterplot to sample plotted signals        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cout = 1;

for i = 1:over:length(t_XI_conv)
    X_wn(count,1) = XI_conv;
    X_wn(count,2) = XQ_conv;
    count = count + 1;
end
Y = [XI_conv XQ_conv];
scatterplot(Y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A11 - Detecting PSK Sequence                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[est_X, est_bit_seq] = detect_PSK_8(Y,A);


