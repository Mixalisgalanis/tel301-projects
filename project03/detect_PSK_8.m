function [est_X,est_bit_seq] = detect_PSK_8(Y,A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [est_X,est_bit_seq] = detect_PSK_8(Y,A)                                       %
% OUTPUT                                                                        %
%      est_X: estimated symbols                                                 %
%      est_bit_seq: estimated bit sequence                                      %
% INPUT                                                                         %      
%      Y: PSK Sequence                                                          %
%      A: Amplitude                                                             %
%                                                                               %
%    M. Galanis, Dec. 2018                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
est_X = zeros(1, length(Y));
m_length = 2^3;

possible_values_i = zeros(1,m_length);
for m = 0: m_length - 1
   possible_values_i(m + 1) = A*cos((m*pi)/4); 
end

possible_values_q = zeros(1,m_length);
for m = 0 : m_length - 1
   possible_values_q(m + 1) = A*sin((m*pi)/4); 
end

for i = 1 : length(Y)
distances = zeros(1,m_length);
    for m = 1 : m_length
        if (i < lenth(Y) /2)
            distances(m) = sqrt((Y(i,1) - possible_values_i(m))^2 + (Y(i,1) - possible_values_i(m))^2);
        else
            distances(m) = sqrt((Y(i,1) - possible_values_q(m))^2 + (Y(i,1) - possible_values_q(m))^2);
        end
    end
    
    index = 1;
    for m = 1 : m_length
        if (distances(m) == min(distances))
            est_X(i) = distances(m);
            index = m;
        end
    end
    
    b = de2bi(est_X(index));
    est_bit_seq(3 * i) = b(1);
    est_bit_seq(3 * i + 1) = b(2);
    est_bit_seq(3 * i + 2) = b(3);
end
return
