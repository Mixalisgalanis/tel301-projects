function [psk_bits] = bits_to_PSK_8(b, A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psk_bits = bits_to_PSK_8(b)                                                   %
% OUTPUT                                                                        %    
%      psk_bits: a vector of 8-PSK symbols                                      %
% INPUT                                                                         %      
%      b: sequence (vector) of binary bits                                      %
%      A: Amplitude                                                             %
%                                                                               %
%    M. Galanis, Dec. 2018                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = length(b)/3;
psk_bits = zeros(L,2);
Xn = [0 0];
n = 1;

for i = 1:3:length(b)
    if ((b(i) == 0) && (b(i+1) == 0) && (b(i+2) == 0))
        m = 0;
    elseif ((b(i) == 0) && (b(i+1) == 0) && (b(i+2) == 1))
        m = 1;
    elseif ((b(i) == 0) && (b(i+1) == 1) && (b(i+2) == 1))
        m = 2;
    elseif ((b(i) == 0) && (b(i+1) == 1) && (b(i+2) == 0))
        m = 3;
    elseif ((b(i) == 1) && (b(i+1) == 1) && (b(i+2) == 0))
        m = 4;
    elseif ((b(i) == 1) && (b(i+1) == 1) && (b(i+2) == 1))
        m = 5;
    elseif ((b(i) == 1) && (b(i+1) == 0) && (b(i+2) == 1))
        m = 6;
    elseif ((b(i) == 1) && (b(i+1) == 0) && (b(i+2) == 0))
        m = 7;
    end
    Xn(1) = A*cos((pi*m)/4);
    Xn(2) = A*sin((pi*m)/4);
    psk_bits(n, 1) = Xn(1);
    psk_bits(n, 2) = Xn(2);
    n = n + 1;
end
return

