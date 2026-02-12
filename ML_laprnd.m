
function yy=ML_laprnd(mu,sigma,T,N,sk,findk)

yy=NaN(T,N);
if nargin == 4
    for nn=1:N
        y1=ML_exprnd(0,0.5*sigma(nn),T,1);
        y2=ML_exprnd(0,0.5*sigma(nn),T,1);
        yy(:,nn)=mu(nn)+y1-y2;
    end
else        % asymmetric laplace distribution    
    if nargin ==5; findk=0; end
    for nn=1:N
        if findk==1
            syms k
            kappa=double(abs(vpasolve(( ( 2*(1-k^6) ) / ( (k^4+1)^1.5 ) ) ==sk(nn),k )));    
        else kappa=sk(nn);
        end
        num=1+kappa^4;
        den=sigma(nn)*kappa^2;
        lambda=sqrt(num/den);
    
        y1=ML_exprnd(0,(kappa*lambda)^(-2),T,1);
        y2=ML_exprnd(0,(kappa/lambda)^2,T,1);    
        yy(:,nn)=mu(nn)+y1-y2;    
    end
end

  