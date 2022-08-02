function [output] = movavgvar(signal,windowSize)
% windowSize = 5;
output = signal;
for pp = 2:windowSize
    b = (1/pp)*ones(1,pp);
    a = 1;
    tmp = filter(b,a,signal);
    output(pp:(end-pp)) = tmp(pp:(end-pp));
end
end