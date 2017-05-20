function y4 = func4(x,he,i)
ds = [-1/he,  1/he];          
s = [(1 - x/he), x/he];       % Shape function (Hat)
y4 = exp(2*(x+(i-1)*he))*s';
end