function y2 = func2(x,he,i)
ds = [-1/he,  1/he];          
s = [(1 - x/he), x/he];       % Shape function (Hat)
y2 = s'*ds*s(2);
end