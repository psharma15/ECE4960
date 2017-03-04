%% Calculating pi with variable digits precision
% Pragya Sharma, ps847
% 11-02-2017

% Chudnovsky algorithm
% Reference: https://en.wikipedia.org/wiki/Approximations_of_pi
syms k pi_inv
k1 = 13591409;
k2 = 545140134;
k3 = 640320;
pi_inv = 0;
get_prec1 = 30;
get_prec2 = 20000;
Niter = 20; 
for i = 1:Niter 
    k = (i - 1);
    pi_inv = pi_inv + ((((-1)^k) * factorial(6*k) * (k1 + (k2*k)))/...
        (factorial(3*k) * (factorial(k))^3 * (k3)^((3*k)+1.5)));
end
pi_inv = 12*pi_inv;
pi_approx = 1/pi_inv;

tic;
pi1 = vpa(pi_approx,get_prec1);
fprintf('Claculating pi with %d digits precision. ',get_prec1)
display(pi1);
t1 = toc;
disp(['Time to calculate: ',num2str(t1),'secs']);
fprintf('-----------------------------------------------------\n');
tic;
pi2 = vpa(pi_approx,get_prec2);
fprintf('Claculating pi with %d digits precision. ',get_prec2)
display(pi2);
t2 = toc;
disp(['Time to calculate: ',num2str(t2),'secs']);
