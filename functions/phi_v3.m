% function data = phi(sp,ph0,ph1,ph2)
%	it applies ph0, ph1 and ph2 to columns of sp
%
%   ph1 defined such that ph1 = {-180, 0, +180} @ {-f,0,+f}
%   ph2 defined such that ph2 = {+180, 0, +180} @ {-f,0,+f}
%   for a spectral width of 2f going from -f to +f
%
function [data, phase] = phi_v3(sp,ph0,ph1,ph2,offset)

%if nargin < 5 ph3 = 0; end
if nargin < 4 ph2 = 0; end
if nargin < 3 ph1 = 0; end
if nargin < 2 ph0 = 0; end

[m n]=size(sp);
a=(-m/2:1:m/2-1)'/m;
%b=((-m/2:1:m/2-1)'.^2)/(m^2/2);
b=(((-m/2-offset):1:(m/2-offset)-1)'.^2)/((m)^2/2);
phase=ph0+a*ph1+b*ph2;
data=sp.*(exp(i*(phase*pi/180))*ones(1,n));
end