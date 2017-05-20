function yInt = simpInt(func,minx,maxx,he,i)
%SIMPINT integrates using Simpson's rule
% Input: 
% func - Function handle to function to be integrated
yInt = ((maxx-minx)/6) * (func(minx,he,i) + 4*func((maxx+minx)/2,he,i) + func(maxx,he,i));
end