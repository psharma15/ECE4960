function y1 = func1(x,he,i)
s = [(1 - x/he), x/he];       % Shape function (Hat)
y1 = s'*s;
end