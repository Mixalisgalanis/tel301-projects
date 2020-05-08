%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC INFORMATION                                   %
% Course: Telecommunication Systems I - Excersize 2   %
% Deadline: 22-Nov-18                                 %
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
% A.1 - Creating SRRC Pulses and Displaying Energy Spectrum%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Settings for plotting the diagrams
r = 6; c = 2; fig = 1; sub_fig = 1;

%Building up parameters for srrc_pulse function
T = 0.01;       %Nyquist parameter (>0)
over = 10;      %Oversampling factor (>0)
Ts = T/over;    %Sampling period (>0)
Fs = 1/Ts;      %Sampling frequency (>0)
A = 4;          %Half duration of the pulse (>0)
a = 0.5;        %Roll-off factor(0<a<1)

%Calling srrc_pulse function to store phi and t variables
[phi, t] = srrc_pulse(T, Ts, A, a);

%Calculating Fourier Transformation (fft & fftshift)
Nf = 2^11;                                  %Number of Samples
X_ES = (F_T(phi,Nf,Ts)).^2;                 %Energy Spectrum
F_ES = (-Fs/2) :(Fs/Nf): (Fs/2 - Fs/Nf);    %Frequency Vector

%Customizing & Plotting Energy spectrum Logarithmically
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES,X_ES);
set_diag_labels("A.1 - Energy Spectrum of SRRC pulse (Logarithmic)","F [Hz]","log|X(F)|^2");
grid on;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A.2 - Energy Spectrum of 2-PAM Sequence Symbols          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 50;         %Number of bits

%Generating 2PAM Sequence
b = gen_N_bits(N);   %Creating N bits
X_seq = bits_to_2PAM(b);        %Using bits_to_2PAM to generate 2-PAM symbols

t_Delta = 0 : Ts : ((N + N*(over - 1)) - 1)*Ts; %Defining time for X_Delta
X_delta = Fs * upsample(X_seq,over);            %Inserts zeros in between bits

%Convoluting X(t) = X_delta with phi
X_conv_time = t_conv(t_Delta,t,Ts); %Setting up convolution length
X = conv(X_delta,phi)*Ts;           %Calculating actual convolution

%Customizing & Plotting X(t)
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_plot(X_conv_time, X);
set_diag_labels('A.2 - X(t): Convolution of X-delta & SRRC Pulse',"t (s)","Conv");
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A.3 - Comparing Theoretical & Estimated Energy Spectrums %
%       of 2-PAM Symbols                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating Periodogram
PX = P_X(X_conv_time, X,Nf,Ts);

%Customizing & Plotting X(t) Linearly
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_plot(F_ES, PX);
set_diag_labels('A.3 - 2-PAM Periodogram (Linear)',"F (Hz)","Conv");
grid on;

%Customizing & Plotting X(t) Logarithmically
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES, PX);
set_diag_labels('A.3 - 2-PAM Periodogram (Logarithmic)',"F (Hz)","Conv");
grid on;

K = [100,1000];
%For every element of K (repeats)
for i=1:length(K)
    Sum_PX = 0;
    
    %For every repeat
    for j=1:K(i)
        b = gen_N_bits(N);   %Creating N bits
        X_seq = bits_to_2PAM(b);        %Using bits_to_2PAM to generate 2-PAM symbols

        t_Delta = 0 : Ts : ((N + N*(over - 1)) - 1)*Ts; %Defining time for X_Delta
        X_delta = Fs * upsample(X_seq,over);            %Inserts zeros in between bits
         
        %Convoluting X(t) = X_delta with phi
        X_conv_time = t_conv(t_Delta,t,Ts); %Setting up convolution length
        X = conv(X_delta,phi)*Ts;           %Calculating actual convolution

        %Summing up in every repeat to create average
        Sum_PX = Sum_PX + P_X(X_conv_time, X,Nf,Ts);   
    end
    
    Average_PX = Sum_PX./K(i);    %Estimated Energy Spectrum
    S_X = (var(X_seq)^2 / T) .* ((F_T(phi,Nf,Ts)).^2); %Theoretical Energy Spectrum
    
    %Customizing & Plotting Theoretical & Estimated P_X Logarithmically
    sub_fig = sub_fig + 1;
    set_diag_layout(fig, r, c, sub_fig);
    c_semilogy(F_ES, S_X);          %Theoretical P_X
    set_diag_labels(['A.3 - 2-PAM Periodogram (Logarithmic) with K =  ' num2str(K(i))],"F (Hz)","Conv");
    hold on;
    c_semilogy(F_ES, Average_PX);  %Estimated P_X
    hold off;
    legend('Theoretical P_X','Estimated P_X')
    grid on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A.4 - Energy Spectrum of 4-PAM Symbols                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating 4PAM Sequence
N = 100;
b = gen_N_bits(N);      %Generate bits
Xseq = bits_to_4PAM(b); %Using bits_to_4PAM to generate 4-PAM symbols

t_Delta = 0 : Ts : ((N + N*(over - 1)) - 1)*Ts; %Defining time for X_Delta
X_delta = Fs * upsample(X_seq,over);            %Inserts zeros in between bits

t_x = t_conv(t,t_Delta,Ts); %Setting up convolution length
X = conv(X_delta,phi)*Ts;   %Calculating actual convolution

%Constructing Periodogram
PX = P_X(t_x, X, Nf, Ts);

%Calculating Theoretical Energy Spectrum
S_X = (var(X_seq)^2 / T) .* ((F_T(phi,Nf,Ts)).^2); %Theoretical Energy Spectrum

%Calculating Estimation of Energy Spectrum
K = 1000;
Sum_PX = 0;
for i=1:K
    b = gen_N_bits(N);      %Generate bits (N remains the same)
    Xseq = bits_to_4PAM(b); %Using bits_to_4PAM to generate 4-PAM symbols

    t_Delta = 0 : Ts : ((N + N*(over - 1)) - 1)*Ts; %Defining time for X_Delta
    X_delta = Fs * upsample(X_seq,over);            %Inserts zeros in between bits
    
    t_x = t_conv(t,t_Delta,Ts); %Setting up convolution length
    X = conv(X_delta,phi)*Ts;   %Calculating actual convolution

    Sum_PX = Sum_PX + P_X(X_conv_time, X,Nf,Ts);
end

Average_PX = Sum_PX ./ K;   %Estimated Energy Spectrum with K = 1000

%Customizing & Plotting Theoretical & Estimated P_X Logarithmically
sub_fig = sub_fig + 2;
set_diag_layout(fig, r, c/2, sub_fig/2);
c_semilogy(F_ES, S_X);          %Theoretical P_X
set_diag_labels(['A.4 - 4-PAM Periodogram (Logarithmic) with K =  ' num2str(K)],"F (Hz)","P_X");
hold on;
c_semilogy(F_ES, Average_PX);  %Estimated P_X
hold off;
legend('Theoretical P_X','Estimated P_X')
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A.5 - Comparing Theoretical & Estimated Energy Spectrums %
%       of 2-PAM Symbols (T' = 2T)                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Changing parameters
T = 2*T;
over = 2*over;
N = 50;         %Number of bits

[phi,t] = srrc_pulse(T, Ts, A, a); %Creating new pulse with updated parameters

%Generating 2PAM Sequence
b = gen_N_bits(N);   %Creating N bits
X_seq = bits_to_2PAM(b);        %Using bits_to_2PAM to generate 2-PAM symbols

t_Delta = 0 : Ts : ((N + N*(over - 1)) - 1)*Ts; %Defining time for X_Delta
X_delta = Fs * upsample(X_seq,over);            %Inserts zeros in between bits

%Convoluting X(t) = X_delta with phi
X_conv_time = t_conv(t_Delta,t,Ts); %Setting up convolution length
X = conv(X_delta,phi)*Ts;           %Calculating actual convolution

%Calculating Periodogram
PX = P_X(X_conv_time, X,Nf,Ts);

%Customizing & Plotting X(t) Linearly
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_plot(F_ES, PX);
set_diag_labels('A.5 - 2-PAM Periodogram (Linear) with T* = 2T',"F (Hz)","Conv");
grid on;

%Customizing & Plotting X(t) Logarithmically
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_semilogy(F_ES, PX);
set_diag_labels('A.5 - 2-PAM Periodogram (Logarithmic) with T* = 2T',"F (Hz)","Conv");
grid on;

K = [100,1000];
%For every element of K (repeats)
for i=1:length(K)
    Sum_PX = 0;
    
    %For every repeat
    for j=1:K(i)
        b = gen_N_bits(N);          %Creating N bits
        X_seq = bits_to_2PAM(b);    %Using bits_to_2PAM to generate 2-PAM symbols

        t_Delta = 0 : Ts : ((N + N*(over - 1)) - 1)*Ts; %Defining time for X_Delta
        X_delta = Fs * upsample(X_seq,over);            %Inserts zeros in between bits
         
        %Convoluting X(t) = X_delta with phi
        X_conv_time = t_conv(t_Delta,t,Ts); %Setting up convolution length
        X = conv(X_delta,phi)*Ts;           %Calculating actual convolution

        Sum_PX = Sum_PX + P_X(X_conv_time, X,Nf,Ts);   %Summing up in every repeat to create average
    end
    
    Average_PX = Sum_PX./K(i);    %Estimated Energy Spectrum
    S_X = (var(X_seq)^2 / T) .* ((F_T(phi,Nf,Ts)).^2); %Theoretical Energy Spectrum
    
    %Customizing & Plotting Theoretical & Estimated P_X Logarithmically
    sub_fig = sub_fig + 1;
    set_diag_layout(fig, r, c, sub_fig);
    c_semilogy(F_ES, S_X);          %Theoretical P_X
    set_diag_labels(['A.5 - P_x(F) 2-PAM Periodogram (Logarithmic) with K =  ' num2str(K(i)) 'and T* = 2T'],"F (Hz)","Conv");
    hold on;
    c_semilogy(F_ES, Average_PX);  %Estimated P_X
    hold off;
    legend('Theoretical P_X','Estimated P_X')
    grid on;
end