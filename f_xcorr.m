function [c] = f_xcorr(sig1,sig2,maxlag)

m = numel(sig1(:,1));
maxlagDefault = m-1;
mxl = min(maxlag,maxlagDefault);
m2 = findTransformLength(m);

X = fft(sig1,m2,1);
Y = fft(sig2,m2,1);

c1 = ifft(X.*conj(Y),[],1,'symmetric');

c = [c1(m2 - mxl + (1:mxl),:); c1(1:mxl+1,:)];

cxx0 = sum(sig1.^2,1);
cyy0 = sum(sig2.^2,1);
scaleCoeffCross = sqrt(cxx0.*cyy0);

c = c./scaleCoeffCross;

%-------------------------------------------------------------------------
function m = findTransformLength(m)
m = 2*m;
while true
    r = m;
    for p = [2 3 5 7]
        while (r > 1) && (mod(r, p) == 0)
            r = r / p;
        end
    end
    if r == 1
        break;
    end
    m = m + 1;
end
