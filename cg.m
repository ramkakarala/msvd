function x=cg(v)
% usage
%        x=cg(v)
% Computes the coding gain for vector v
% which is the ratio of arithmetic to gemetric
% means
% rk 1 June 2000
%                                 kakarala@labs.agilent.com
% 
% Agilent GIVES NO EXPRESS OR IMPLIED WARRANTY OF ANY KIND AND 
% ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR PURPOSE ARE DISCLAIMED.
% Agilent SHALL NOT BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
% OR CONSEQUENTIAL DAMAGES ARISING OUT OF ANY USE OF THIS SOFTWARE.

am=mean(v);
gm=prod(v)^(1/length(v));
if gm>0
   x=am/gm;
else
   x=NaN;
end;
