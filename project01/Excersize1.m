%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC INFORMATION                                   %
% Course: Telecommunication Systems I - Excersize 1   %
% Deadline: 1-Nov-18                                  %
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
t_conv = @(t1, t2, dt) (t1(1) + t2(1)):dt:(t1(end) + t2(end)); %convolution time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A.1 - Creating pulses of SRRC (Square Root Raised Cosine)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Settings for plotting the diagrams
r = 2; c = 2; fig = 1; sub_fig = 1;

%Building up parameters for srrc_pulse function
T = 0.1;        %Nyquist parameter (>0)
over = 10;      %Oversampling factor (>0)
Ts = T/over;    %Sampling period (>0)
Fs = 1/Ts;      %Sampling frequency (>0)
A = 5;          %Half duration of the pulse (>0)
a = [0,0.5,1];  %3 Roll-off factors in a vector (0<a<1)

%Plotting 3 srrc_pulses with each roll-off factor in the same diagram
for i=1:length(a)
    %Calling srrc_pulse function to store phi and t variables
    [phi, t] = srrc_pulse(T, Ts, A, a(i));
    
    %Plotting phi in t time in a figure with a title and x,y labels
    set_diag_layout(fig, r, c, sub_fig);
    set_diag_labels("A.1 - Srrc pulses with different roll-off factors","t","phi(t)");
    c_plot(t, phi);
    
    grid on;    %Grid Enabled
    hold on;    %Any Next plots are applied in this diagram
end
legend('a = 0','a = 0.5','a = 1')
hold off;   %Stop plotting on the same diagram



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A.2 - Fourier Transformation of previous srrc_pulses     %
%       and Display of their Energy Spectrum               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plotting (Linear) Amplitude Spectrums for F.T of the above 3 srrc_pulses
sub_fig = sub_fig + 1;
for i=1:length(a)
    %Calling srrc_pulse function to get phi and t variables
    [phi, t] = srrc_pulse(T, Ts, A, a(i));
    
    NFFT = 2^nextpow2(length(phi));     %Length of the signal for the F.T
    X_seq = abs(fftshift(fft(phi,NFFT)*Ts));     %F.T (Aligned to center with fftshift)

    t_Delta = (-Fs/2) :(Fs/NFFT): (Fs/2 - Fs/NFFT);       %Frequency Vector
    
    %Plotting Energy spectrum Linearly
    set_diag_layout(fig, r, c, sub_fig);
    set_diag_labels("A.2 - Amplitude Spectrum of srrc pulses (Linear)","F [Hz]","|X(F)|^2");
    c_plot(t_Delta,X_seq);
    
    grid on;    %Grid Enabled
    hold on;    %Any Next plots are applied in this diagram
end
legend('a = 0','a = 0.5','a = 1')
hold off;   %Stop plotting on the same diagram

%Plotting (Logarithmic) Energy Spectrums for F.T of the above 3 srrc_pulses
sub_fig = sub_fig + 1;
for i=1:length(a)
    %Calling srrc_pulse function to get phi and t variables
    [phi, t] = srrc_pulse(T, Ts, A, a(i));
    
    NFFT = 2^nextpow2(length(phi));     %Length of the signal for the F.T
    X_seq = abs(fftshift(fft(phi,NFFT)*Ts));     %F.T (Aligned to center with fftshift)

    t_Delta = (-Fs/2) :(Fs/NFFT): (Fs/2 - Fs/NFFT);       %Frequency Vector
    
    %Plotting Energy spectrum Logarithmically
    set_diag_layout(fig, r, c, sub_fig);
    set_diag_labels("A.2 - Amplitude Spectrum of srrc pulses (Logarithmic)","F [Hz]","log|X(F)|^2");
    c_semilogy(t_Delta,X_seq);  
    
    grid on;    %Grid Enabled
    hold on;    %Any Next plots are applied in this diagram
end
legend('a = 0','a = 0.5','a = 1')
hold off;   %Stop plotting on the same diagram



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A.3 - Creating pulses of SRRC (Square Root Raised Cosine)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_fig = sub_fig + 1;
C1 = T / 10^1;
C2 = T / 10^3;

for i=1:length(a)
    BW = (1 + a(i))/(2 * T);	%Calculating Bandwidth
    
    %Calling srrc_pulse function to store phi and t variables
    [phi, t] = srrc_pulse(T, Ts, A, a(i));
    
    NFFT = 2^nextpow2(length(phi));     %Length of the signal for the F.T
    X_seq = abs(fftshift(fft(phi,NFFT)*Ts));     %F.T (Aligned to center with fftshift)

    t_Delta = (-Fs/2) :(Fs/NFFT): (Fs/2 - Fs/NFFT);       %Frequency Vector
    
    %Plotting Energy spectrum Logarithmically
    set_diag_layout(fig, r, c, sub_fig);
    set_diag_labels("A.2 - Amplitude Spectrum of srrc pulses (Logarithmic)","F [Hz]","log|X(F)|^2");
    c_semilogy(t_Delta,X_seq);  
    
    grid on;    %Grid Enabled
    hold on;    %Any Next plots are applied in this diagram
end
c_plot(t_Delta,C1);
c_plot(t_Delta,C2);
legend('a = 0','a = 0.5','a = 1', 'c_1 = 10^-^2', 'c_2 = 10^-^4');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B.1 - SRRC & ORTHONORMALITY                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%New excersize, new figure
r = 4; c = 3; fig = fig+1; sub_fig = 1;

%Setting up k
k = [0, 1];

%For every a value
for i=1:length(a)
    %Calling srrc_pulse function to store phi and t variables
    [phi, t] = srrc_pulse(T, Ts, A, a(i));  %get phi and t from srrc_pulse
    
    %For every k value
    for j=1:length(k)
        %Calculating Phi_moved
        phi_moved = [zeros(1, (1/Ts)*k(j)*T) phi(1:end - (1/Ts)*k(j)*T)];
        %Preparing Graphs
        set_diag_layout(fig, r, c, sub_fig);
        set_diag_labels(['B.1 - SRRC with a = ' num2str(a(i)) ' and k = ' num2str(k(j))],'t (s)','Phi(t)');
        grid on;
        %Plotting Phi(t) & Phi_moved(t) In Same Diagram
        hold on;
        c_plot(t, phi);
        c_plot(t, phi_moved);
        legend('phi(t)','phi(t-kT)');
        hold off;
        sub_fig = sub_fig + c;
        
        
        %Calculating multiplication of Phi & Phi_moved        
        mult_phis = phi .* phi_moved; 
        %Preparing Graphs
        set_diag_layout(fig, r, c, sub_fig);
        set_diag_labels(['B.1 - Mult of Phi(t), Phi(t-kT) with a = ' num2str(a(i)) ' and k = ' num2str(k(j))],"t (s)","Phi(t)");
        %Plotting Multiplication of Phi(t) & Phi_moved(t)                
        grid on;
        c_plot(t, mult_phis);
        title(['B.1 - Mult of Phi(t), Phi(t-kT) with a = ' num2str(a(i)) ' and k = ' num2str(k(j))]);
        legend('mult(phi(t),phi(t-kT))');
        sub_fig = sub_fig + c;
        
        %Calculating Integral (area)
        area = sum(mult_phis)*Ts   
    end
    sub_fig = sub_fig - 4*c + 1 ;
    hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C - 2-PAM BASIC SYSTEM ZONE (USING 2-PAM MODULATION)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%New excersize, new figure
r = 2; c = 2; fig = fig+1; sub_fig = 1;

a = a(2);       %a becomes 0.5, other variables remain the same
N = 50;         %Number of bits

b = (sign(randn(N,1)) + 1)/2;   %Creating N bits
X_seq = bits_to_2PAM(b);        %Using bits_to_2PAM to generate PAM symbols

t_Delta = 0 : Ts : ((N + N*(over - 1)) - 1)*Ts; %Defining time for X_Delta
X_delta = Fs * upsample(X_seq,over);            %Inserts zeros in between bits

%Customizing and Plotting X_delta
set_diag_layout(fig, r, c, sub_fig);
c_plot(t_Delta, X_delta);
set_diag_labels('C - X-delta',"t (s)","X-Delta");
grid on;


[phi, t] = srrc_pulse(T,Ts,A,a);        %Gathering SRRC Pulse
X_conv_time = t_conv(t_Delta, t, Ts);   %Calculating convolution time
X = conv(X_delta, phi)*Ts;              %Calculating convolution X

%Customizing and Plotting X
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, sub_fig);
c_plot(X_conv_time, X);
set_diag_labels('C - X(t): Convolution of X-delta & SRRC Pulse',"t (s)","Conv");
grid on;

%Calculating Phi(-t)
phi_flipped = fliplr(phi);  %Flipping the phi signal
t_flipped = -fliplr(t);      %Flipping time axis

Z = conv(X,phi_flipped)*Ts;                        %Calculating convolution Z
Z_conv_time = t_conv(X_conv_time,t_flipped,Ts);    %Calculating convolution time

%Customizing and Plotting Z
sub_fig = sub_fig + 1;
set_diag_layout(fig, r, c, [sub_fig, sub_fig+1]);
c_plot(Z_conv_time, Z);
set_diag_labels('C - Z(t): Convolution of X & Phi(-t)',"t (s)","Conv");
grid on;
hold on;

d_plot(0:(N-1)*T,X);    %Show Z and X at the same time
legend('Z(kT)', 'X(kT)');
hold off;



