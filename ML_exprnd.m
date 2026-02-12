% y=ML_exprnd(mu,sigma,T,N)
%    mu = mean
% Sigma = Variance
% 
function y=ML_exprnd(mu,sigma,T,N)

lambda=1/sqrt(sigma); % rate

y=mu+exprnd(1/lambda,T,N)-1/lambda;