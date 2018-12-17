function y = threshf(x,sorh,t,alpha)
%%% THRESHf Perform soft or hard or trimmed thresholding. 
%%% Y = threshf(X,SORH,T) returns soft (if SORH = 's')
%%% or hard (if SORH = 'h') T-thresholding  of the input 
%%% vector or matrix X. T is the threshold value.
%%%
%%% Y = threshf(X,'s',T) returns Y = SIGN(X).(|X|-T) for |X|>T else 0. soft 
%%% thresholding is shrinkage.
%%%
%%% Y = threshf(X,'h',T) returns Y = X for (|X|>T) else 0. hard
%%% thresholding is cruder.
%%%
%%% Y = threshf(X,'t',T,alpha) returns Y = X*[|X|^alpha-T^alpha]/|X|^alpha for (|X|>T), 
%%% else 0. trimmed thresholding is something between soft and hard thresholding.
%%%
%%% Author : B.K. SHREYAMSHA KUMAR 
%%% Created on 01-12-2009.
%%% Updated on 01-12-2009.

switch sorh
    case 's'
        tmp = (abs(x)-t);
        tmp = (tmp+abs(tmp))/2;
        y   = sign(x).*tmp;

    case 'h'
        y   = x.*(abs(x)>t);

    case 't'
        x=x+eps;
        tmp = abs(x).^alpha-t^alpha;
        tmp = (tmp+abs(tmp))/2;
        y = x.*tmp./abs(x).^alpha;

    otherwise
        error('Invalid argument value.')
end
