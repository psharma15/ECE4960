%% Pragya Sharma, ps847, 23th Feb 2017.
% This code scales i-th row by a and adds it to j-th row of square full
% rank matrices. 
function [Anew] = rowScaleF(A,i,j,a)
Anew = A;
Anew(j,:) = a*A(i,:)+A(j,:);
end