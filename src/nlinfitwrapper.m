% "nlinfit" wrapper function
%
% mikael.mieskolainen@cern.ch, 2019

function yhat = nlinfitwrapper(param, data)

% input data is global variable in costfunc1

[~,yhat,~] = costfunc1(param, data);

end