function [pam_bits] = bits_to_2PAM(b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pam_bits = bits_to_PAM(b)                                                     %
% OUTPUT                                                                        %    
%      pam_bits: a vector of 2-PAM symbols X                                    %
% INPUT                                                                         %      
%      b: sequence (vector) of binary bits                                      %
%                                                                               %
%    M. Galanis, Oct. 2018                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pam_bits = zeros(1,length(b));

%Equivalent bits are converted from 0 to +1 and from 1 to -1
for i = 1:length(b)
   if (b(i) == 0); pam_bits(i) = 1; else; pam_bits(i) = -1; end
end
return
