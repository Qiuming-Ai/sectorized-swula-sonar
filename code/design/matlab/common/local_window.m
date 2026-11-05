function w = local_window(name, n)
%LOCAL_WINDOW Toolbox-free symmetric window generator.
if n < 1 || fix(n) ~= n
    error('n must be a positive integer.');
end
if n == 1
    w = 1;
    return;
end
k = (0:n-1).';
switch lower(string(name))
    case "rect"
        w = ones(n, 1);
    case "hann"
        w = 0.5 - 0.5*cos(2*pi*k/(n-1));
    case "hamming"
        w = 0.54 - 0.46*cos(2*pi*k/(n-1));
    otherwise
        error('Unknown window: %s', name);
end
end
