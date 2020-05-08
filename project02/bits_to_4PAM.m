function [pam_bits] = bits_to_4PAM(b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pam_bits = bits_to_4PAM(b)                                                    %
% OUTPUT                                                                        %    
%      pam_bits: a vector of 4-PAM symbols X                                    %
% INPUT                                                                         %      
%      b: sequence (vector) of binary bits                                      %
%                                                                               %
%    M. Galanis, Nov. 2018                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pam_bits = zeros(1,length(b)/2);

%Equivalent bits are converted like so:
%00 -----> +3
%01 -----> +1
%11 -----> -1
%10 -----> -3
for i = 1:2:length(b)
    if (b(i) == 0 && b(i+1) == 0);     pam_bits((i+1)/2) = +3;
    elseif (b(i) == 0 && b(i+1) == 1); pam_bits((i+1)/2) = +1;
    elseif (b(i) == 1 && b(i+1) == 1); pam_bits((i+1)/2) = -1;
    elseif (b(i) == 1 && b(i+1) == 0); pam_bits((i+1)/2) = -3;
    end
end
return
