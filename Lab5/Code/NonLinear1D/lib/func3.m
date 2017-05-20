function y3 = func3(x,he,i)
ds = [-1/he,  1/he];          
s = [(1 - x/he), x/he];       % Shape function (Hat)
y3 = s'*ds*s(1);
end