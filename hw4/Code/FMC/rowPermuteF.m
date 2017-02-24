%% Pragya Sharma, ps847, 23th Feb 2017
% This code interchanges rows i and j, where i<j for full matrix
function [Anew] = rowPermuteF(A,i,j)
Anew = A;
Anew(i,:) = A(j,:);
Anew(j,:) = A(i,:);
end