% BL_SimulationsI0_2026_aux - 
% 
%     T - number of obs
%     N - number of variables
%     r - number of common factors
%     q - number of common shocks
%   tau - correlation idiosyncratic component
%    mu - max eigenvalue Companion matrix VAR factors
% delta - Autocorrelation idiosyncratic component
% theta - signal to noise ratio
%   stT - shocks distribution (0,1,2,3,4) - Normal, Student-t, Laplace, Asymmetric Laplace, Hansen (1994) Skew-t
%  iota - Mess-up EM initialization
%    pi - Perturb generated factors (0,1)
% 
function [EM, PCA, DGP]=BL_SimulationsI0_2026_aux(T,N,r,q,tau,mu,delta,theta, stT,iota,perturb)


try isempty(stT); if isempty(stT); stT=0; end; catch; stT=0; end
try isempty(stT); if isempty(stT); stT=0; end; catch; stT=0; end
Tb = T + 10;                                                             % Precompute burn-in length

%%%%% ============================================== %%%%%
%%%%% ==== ------------------------------------ ==== %%%%%
%%%%% ==== ---                              --- ==== %%%%%
%%%%% ==== ---  PART 1: Simulate the model  --- ==== %%%%%
%%%%% ==== ---                              --- ==== %%%%%
%%%%% ==== ------------------------------------ ==== %%%%%
%%%%% ============================================== %%%%%


%%% ========================= %%%
%%%  Simulate common factors  %%%
%%% ========================= %%%

if stT==0   % ------------------------------------------------------------- % common shocks 
    u=randn(Tb,q);                                                          % Normally distributed
elseif stT==1; u=trnd(4,Tb,q);                                              % Student-t shocks
elseif stT==2; u=ML_laprnd(zeros(q,1),ones(q,1),Tb,q);                      % Laplace shocks    
elseif stT==3; u=ML_laprnd(zeros(q,1),ones(q,1),Tb,q,0.9+0.2*rand(q,1));    % Asymmetric Laplace shocks   
elseif stT==4                                                               % Hansen (1994) Skew-t
    nu = 4 + 8*rand(1,q);                                                   % The nu of the skew T is ~ U(3,13)
    gamma = -.1 + .2*rand(1,q);                                             % The Gamma of the skew T is ~ U(-.1,.1)    
    u=zeros(Tb,q); for ii=1:q; u(:,ii)= skewtrnd(nu(ii),gamma(ii),Tb,1); end % -------------------- 
else
    u=randn(Tb,q); 
end         % ------------------------------------------------------------- %                          
uu=ML_Standardize(u);                                                       % standardize common shocks 

            % ------------------------------------------------------------- % VAR Common Factors
A = diag(ML_uniform(q, 0.5, 0.8));
off_diag_vals = ML_uniform(q*(q-1), 0.1, 0.4);
A(~eye(q)) = off_diag_vals;
mut = max(abs(eig(A)));                                                     % Max Root of the polynomial
A = (mu / mut) * A;                                                         % Fix max eigenvalue
            % ------------------------------------------------------------- % 

            % ------------------------------------------------------------- % Generate common factors
ff0 = zeros(Tb, q);
ff0(1,:) = uu(1,:);
for tt = 2:Tb
    ff0(tt,:) = (A * ff0(tt-1,:)')' + uu(tt,:);
end
            % ------------------------------------------------------------- %
if r > q
    ff1=ML_lag(ff0,(r/q)-1,0);                                                  % take lags of the common factors
    ff1=ff1(end-T+1:end,:);                                                     % keeps only the last T observations 
else
    ff1 = ff0(end-T+1:end, :);
end

%%% =================================== %%%
%%%  Simulate idiosyncratic components  %%%
%%% =================================== %%%
if delta==0     % --------------------------------------------------------- % Autocorrelation idiosyncratic component
    d=zeros(N,1);                                                           % d = 0
elseif delta==1; d=[zeros(round(N/2),1); .3*ones(N-round(N/2),1)];          % [d=0 ; d = {0, 0.3}]
elseif delta>2; d=ML_uniform(N,0,delta-2);                                  % d ~ U(0,delta)
else; d=ML_uniform(N,delta,1-2*delta);                                      % d ~ U(delta,1-2*delta)
end             % --------------------------------------------------------- %

TT1=diag(rand(N,1)+0.5);                                                    % variance idiosyncratic shocks
if tau==0; TT2=eye(N); else; TT2=toeplitz(tau.^(0:N-1)); end                % covariance matrix idiosyncratic shocks
TT2=TT2-triu(TT2,11)-tril(TT2,-11);                                         % keep only 10 correlations

sTT1_cTT2 = sqrt(TT1) * cholcov(TT2);                                       % Precompute Cholesky decomposition
if stT==0   % ------------------------------------------------------------- % idiosyncratic shocks 
    v=randn(Tb,N)*sTT1_cTT2;                                                % Normally distributed
elseif stT==1; v=trnd(4,Tb,N)*sTT1_cTT2;                                    % Student-t shocks
elseif stT==2; v=ML_laprnd(zeros(N,1),ones(N,1),Tb,N)*sTT1_cTT2;            % Laplace shocks
elseif stT==3                                                               %
    v=ML_laprnd(zeros(N,1),ones(N,1),Tb,N,0.9+0.2*rand(N,1))*sTT1_cTT2;     % Asymmetric Laplace shocks                                                     
elseif stT==4                                                               % Hansen (1994) Skew-t    
    nu = 3 + 10*rand(1,N);                                                  % The nu of the skew T is ~ U(3,13)
    gamma = -.15 + .3*rand(1,N);                                            % The Gamma of the skew T is ~ U(-.15,.15)
    vv=zeros(Tb,N); for ii=1:N; vv(:,ii)= skewtrnd(nu(ii),gamma(ii),Tb,1); end % -------------------- 
    v=vv* sTT1_cTT2;                                                        % --------------------                                   
end         % ------------------------------------------------------------- %


e=v; for ii=1:N; e(:,ii)=filter(1,[1 -d(ii)],v(:,ii)); end                  % Idiosyncratic component
e = e(11:end, :);                                                           % Remove burn-in
e=ML_center(e);                                                             % center idiosyncratic component

%%% ================================= %%%
%%%  Common components and variables  %%%
%%% ================================= %%%
LL1=1+randn(N,r);                                                           % factor loadings
chi=ML_center( ff1*LL1' );                                                  % common components

% ------------------------------------------------------------------------- % fix signal to noise ratio
temp=var(chi)./var(e);                                                      % ratio of variances
theta2=theta+ML_center(ML_uniform(N,theta-.25,theta+.25));                  % draw signal to noise ratio
scaling_factors = (1 ./ sqrt(theta2 .* temp'))';
chi = chi .* repmat(scaling_factors, T, 1);                                 % impose signal to noise ratio
% ------------------------------------------------------------------------- %


%%% =============================== %%%
%%%  Identify factors and loadings  %%%
%%% =============================== %%%
[W,M] = eigs(cov(chi), r,'LM'); W=W*diag(sign(W(1,:))');                    % eigenvalue-eigenvectors of chi      
lambda=W*sqrt(M);                                                           % identified factor loadings
f=chi*W/sqrt(M);                                                            % identified common factors
if perturb==1; f=f*chol(eye(r)+T^-0.5)'; end                                     % perturb common factors
chi=f*lambda';
x=chi+e;                                                                    % variables
H1=NaN(r+1,r); for rr=1:r;  H1(:,rr)=ML_ols(f(:,rr),ff1,1); end             % Rotation matrix
const=H1(1,:); H1=H1(2:end,:); H2=eye(r)/H1;                                % ---------------
AA=H1'*[A zeros(q,r-q);eye(r-q) zeros(r-q,r-q)]*H2';
G=H1'*[eye(q); zeros(r-q,r-q)];    
const2=((eye(r)-AA)*const')';



%%% =================== %%%
%%%  Save DGP results   %%%
%%% =================== %%%
DGP.x=x;
DGP.lambda=lambda;
DGP.f=f;
DGP.chi=chi;
DGP.xi=e;
DGP.G=G;
DGP.A=AA;
DGP.const=const2;




%%%%% ============================================== %%%%%
%%%%% ==== ------------------------------------ ==== %%%%%
%%%%% ==== ---                              --- ==== %%%%%
%%%%% ==== ---  PART 2: Estimate the model  --- ==== %%%%%
%%%%% ==== ---                              --- ==== %%%%%
%%%%% ==== ------------------------------------ ==== %%%%%
%%%%% ============================================== %%%%%

p=1;
pr=p*r;
iter=50; 
tresh = 1e-4;

standardizzo=1;     % ----------------------------------------------------- % Determine whether to standardize...
if standardizzo==1                                                          % ... the data prior to estimation
    [xx,mx,sx]=ML_Standardize(x);                                           %
else                                                                        %
    xx=x; mx=zeros(1,N); sx=ones(1,N);                                      %
end                                                                         %
MX=repmat(mx,T,1); SX=repmat(sx,T,1);                                       % Useful objects
DGP.xx=xx; DGP.mx=mx; DGP.sx=sx;                                            % store variances and mean of the data

%%% ======================================== %%%
%%%  Estimate factors and loadings with PCA  %%%
%%% ======================================== %%%
Gamma=cov(xx);                                                              % covariance matrix of the standardize data
[W,M] = eigs(Gamma, r,'LM'); W=W*diag(sign(W(1,:))');                       % eigenvalue-eigenvectors
L0=W*sqrt(M);                                                               % estimate factor loadings
f0=xx*W/sqrt(M);                                                            % estimate common factors
for rr=1:r          % ----------------------------------------------------- % Fix the sign of the estimated ...
    rho=corr(f(:,rr),f0(:,rr));                                             % ... factors. These are just some ...
    if rho<0; f0(:,rr)=-f0(:,rr); L0(:,rr)=-L0(:,rr); end                   % ... weird and perverse cases
end                 % ----------------------------------------------------- %
e0=xx-f0*L0';                                                               % "residual" PCA
[~, uu0,AL0]=ML_VAR(f0,p,0);                                                % VAR on the Static Factors
[~, G0]=ML_edynfactors2(uu0,q);                                             % Common Shocks
chi0=MX+(f0*L0').*SX;                                                       % common component PCA
xi0=x-chi0;                                                                 % idiosyncratic component PCA

PCA.F=f0; PCA.Lambda=L0; PCA.chi=chi0; PCA.xi=xi0; % ---------------------- % Saves PCA estimates

if iota~=0      % ========================================================= % Mess-up EM initialization
    W=W+randn(N,r)*cholcov(toeplitz(iota.^(0:r-1)));                        % contaminates the eigenvector
    W=W/sqrt(diag(diag(W'*W)));                                             % ----------------------------
    L0=W*sqrt(M);                                                           % estimate factor loadings
    f0=xx*W/sqrt(M);                                                        % estimate common factors
    for rr=1:r          % ------------------------------------------------- % Fix the sign of the estimated ...
        rho=corr(f(:,rr),f0(:,rr));                                         % ... factors. These are just some ...
        if rho<0; f0(:,rr)=-f0(:,rr); L0(:,rr)=-L0(:,rr); end               % ... weird and perverse cases
    end                 % ------------------------------------------------- %
    e0=xx-f0*L0';                                                           % "residual" PCA
    [~, uu0,AL0]=ML_VAR(f0,p,0);                                            % VAR on the Static Factors
    [~, G0]=ML_edynfactors2(uu0,q);                                         % Common Shocks
    chi0=MX+(f0*L0').*SX;                                                   % common component PCA
    xi0=x-chi0;                                                             % idiosyncratic component PCA
end             % ========================================================= %

%%% ==================================== %%%
%%%  Expectation Maximization Algorithm  %%%
%%% ==================================== %%%
AL0b=cat(3,AL0,zeros(r,r));                                                 % add one lag to avoid computing PtTm
[F0,lambda0,A0,P0,Q0,R0]=ML_SS_DFM_I0(f0,L0,AL0b,G0,e0,N,p+1,r);            % State-Space representation
[xitT,PtT1,~,~,~,~,A1,L1,R1,Q1,G1]=...                                      % EM algorithm
    ML_efactorML2026(xx,F0,P0,A0,lambda0,R0,Q0,r,q,p,iter,tresh,0);         % ------------
PtT1=PtT1(1:pr,1:pr,:);A1=A1(1:pr,1:pr); Q1=Q1(1:pr,1:pr);                  % Eliminate the extra lag
f1=xitT(:,1:r); L1=L1(:,1:pr);                                              % factors
jr=diag(sign(L1(1,1:r))'); L1(:,1:r)=L1(:,1:r)*jr; f1=f1*jr;                % Loadings of first variable are positive (normalization)
for rr=1:r          % ----------------------------------------------------- % Fix the sign of the estimated ...
    rho=corr(f(:,rr),f1(:,rr));                                             % ... factors. These are just some ...
    if rho<0; f1(:,rr)=-f1(:,rr); L1(:,rr)=-L1(:,rr); end                   % ... weird and perverse cases
end                 % ----------------------------------------------------- %
chi1=MX+(f1*L1(:,1:r)').*SX;                                                % Common component
xi1=x-chi1;                                                                 % idiosyncratic component


% -------------------------------------------------------------------- %
% Verify identification by computing MSE of distances from diag matrix %
% -------------------------------------------------------------------- %
eig_L1L1 = eig(L1' * L1 / N); diag_L1L1 = diag(L1' * L1 / N);
I_L1 = mean((sort(eig_L1L1) - sort(diag_L1L1)).^2);             % factor loadings
I_f1=mean( ( eig(f1'*f1/T) - 1 ).^2 );                                      % factors


%%% ==================== %%%
%%%  Save EM results    %%%
%%% ==================== %%%
EM.F=f1;
EM.Lambda=L1;
EM.chi=chi1;
EM.xi=xi1;
EM.R=R1;
EM.G=G1;
EM.A=A1;
EM.Q=Q1;
EM.PtT=PtT1;
EM.Gdag=(G1'*G1)\G1'; 
EM.xitT=xitT;
EM.I_f1=I_f1;
EM.I_L1=I_L1;





% ===========================================================================
% ===========================================================================
% ---------------------------------------------------------------------------
% 
% FUNCTIONS
% 
% ---------------------------------------------------------------------------
% ===========================================================================
% ===========================================================================






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_Standardize - Standardize Variables
%
% [y M s] = ML_Standardize(x)
%
% Written by Matteo Luciani (matteoluciani@yahoo.it)
% =========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, M, s] = ML_Standardize(x)

T=size(x,1);
s = nanstd(x);
M = nanmean(x);
ss = ones(T,1)*s;
MM = ones(T,1)*M;
y = (x-MM)./ss;





% ===========================================================================
% ML_uniform - Generate U(a,b) of size T
% x=ML_uniform(T,a,b)
% 
% Written by Matteo Luciani (matteo.luciani@ulb.ac.be)
% ===========================================================================

function x=ML_uniform(T,a,b)
 x=a+(b-a).*rand(T,1);
 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% =========================================================================
% ML_center - Demean variables
% CENTER XC = center(X)
%   Centers each column of X.
%   J. Rodrigues 26/IV/97, jrodrig@ulb.ac.be
% =========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XC = ML_center(X)
T = size(X,1);
XC = X - ones(T,1)*(sum(X)/T); % Much faster than MEAN with a FOR loop





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_lag - Given x(t) produces x(t-j0) .... x(t-k)
%
% xx=ML_lag(x,k,j0)
% by default j0 is 1
%
% eg x=[x1 x2]
% xx=ML_lag(x,2)=[x1_{t-1} x1_{t-2} x2_{t-1} x2_{t-2}]
%
% Matteo Luciani (matteoluciani@yahoo.it)
% =========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xx=ML_lag(x,k,j0)
[T, N] = size(x);

if nargin<3; j0=1; end;    
n=1;

if N*(k+1-j0)==0
    xx=[];
elseif T==1
    xx=x;
else
    xx=zeros(T-k,N*(k+1-j0));
    for i=1:N
        for j=j0:k
            xx(:,n)=x(k+1-j:T-j,i);
            n=n+1;
        end;
    end;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_ols - OLS Estimation
% 
% [beta,u,v,esu,r2,espqr]=ML_ols(y,x,det);
%   Inputs:
%       y   = Endogenous Variable (vector)
%       x   = Exogenous Variables (matrix)
%       det = 0 no constant
%       det = 1 constant
%       det = 2 time trend
%       det = 3 constant + time trend
%   Outputs:
%       beta  = estimated coefficient
%       u     = residuals
%       v     = varcov matrix of the estimates
%       esu   = Residual Variance
%       r2    = R-Squared
%       espqr = estimates standard errors
%
% Written by Matteo Luciani (matteoluciani@yahoo.it)
% This is a modified version of the codes available on Fabio Canova webpage
% =========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [beta,u,v,esu,r2,espar,yhat]=ML_ols(y,x,det)
T = size(x,1);

    % Add deterministic components
    if det == 1
        x = [ones(T,1), x];
    elseif det == 2
        x = [(1:T)', x];
    elseif det == 3
        x = [ones(T,1), (1:T)', x];
    end
k=size(x,2);                                                                % number of parameters
xx = (x' * x) \ eye(k);                                                     % inv(x'x)
beta=xx*(x'*y);                                                             % ols coeff
yhat=x*beta;                                                                % fitted values
u=y-yhat;                                                                   % residuals
uu=u'*u;                                                                    % SSR
esu=uu/(T-k);                                                               % Residual Variance
yc=y-mean(y);                                                               % Centered variables
r2=1-uu/(yc'*yc);                                                           % R2
v=esu*xx;                                                                   % varcov matrix of the estimates
espar=sqrt(diag(v));                                                        % Standard errors





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
%  ML_SS_DFM_I0 - Build State-Space representation of stationary DFM
% 

function [F0,lambda0,A0,P0,Q0,R0]=ML_SS_DFM_I0(f0,L0,AL0,G0,xi0,N,p,r,IDEN)

if nargin <9; IDEN=1; end 
pr=p*r;
Q0 = zeros(pr,pr); Q0(1:r,1:r)=G0*G0';                                      % Variance of the VAR residuals
lambda0 = [L0 zeros(N,r*(p-1))];                                            % Loadings
R0 = diag(diag(cov(xi0)));                                                  % covariance idio
F0=[]; for pp=1:p; F0=cat(1,F0,f0(p+1-pp,:)'); end                          % initial conditions factors
A0=ML_VAR_companion_matrix(AL0);                                            % VAR companion form aka VAR(1)

if IDEN==1; P0=eye(pr);                                                     % Impose initial variance = identity matrix
else; P0 = reshape(inv(eye(pr^2)-kron(A0,A0))*Q0(:),pr,pr);                 % initial variance
end
% =========================================================================






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_efactorML2026 - ML estimation I0 DFM, straight EM no initialization, P00=I (for simulations)
% 
% [xitT,PtT,PtTm,xitt,xittm,Ptt,Pttm,A,Lambda,R,Q]=...
%     ML_efactorML4(y,F0,P0,A,Lambda,R,Q,q,p,maxiter,tresh,cc)% 
% 
%       Y - data with deterministic component
%      F0 - initial values for the states
%      P0 - initial variance for the states
%       A - VAR(1) for the states
%  Lambda - Loadings
%       R - covariance matrix errors
%       Q - covariance matrix shocks
%       r - number of factors
%       q - number of shocks
%       p - number of lags VAR factors
% maxiter - max number iterations
%   tresh - threshold for stopping algorithm
%      cc - # of obs eliminated before estimating parameters
% 

function [xitT,PtT,xitt,xittm,Ptt,Pttm,A,Lambda,R,Q,G]=...
    ML_efactorML2026(y,F0,P0,A,Lambda,R,Q,r,q,p,maxiter,tresh,cc)


OPTS.disp = 0;
[T, N]=size(y);                                                             % size of the panel
pr=p*r; pr1=(p+1)*r;

%%% =================== %%%
%%%  Start with E-Step  %%%
%%% =================== %%%
[xitt,xittm,Ptt,Pttm,loglik]=ML_KalmanFilter2026(F0,P0,y,A,Lambda,R,Q);        % Kalman Filter
[xitT,PtT]=ML_KalmanSmoother2026(y,A,xitt',xittm',Ptt,Pttm,Lambda,R);          % Kalman Smoother
F00=xitT(1,:)'; P00=eye(pr1);                                               % initial conditions  
F=xitT(:,1:r);                                                              % factors
loglik1=loglik;                                                             % likelihood

yyc=ML_center(y)';                                                      % endogenous variables
yy=y';

for jj=1:maxiter         
    
    %%%% ======================== %%%%
    %%%% ====                ==== %%%%
    %%%% ====    M - Step    ==== %%%%
    %%%% ====                ==== %%%%
    %%%% ======================== %%%%
    
    %%% =========================== %%%
    %%% 1: Estimate factor Loadings %%%
    %%% =========================== %%%
    
    xx=ML_center(F)';                                                       % exogenous variables, aka the factors        
    num=zeros(N,r); den=zeros(r,r);                                         % preallocates
    for tt=cc+1:T  % ------------------------------------------------------ % \sum_{t=1}^T
        num=num+yyc(:,tt)*xx(:,tt)';                                         % y_t F_{t|T}'
        den=den+xx(:,tt)*xx(:,tt)'+PtT(1:r,1:r,tt);                         % F_{t|T}F_{t|T}'+P_{t|T}
    end    
    L=num/den;                                                              % factor loadings           
    Lambda(1:N,1:r) = L;                                                    % store factor loadings

    %%% ======================================================= %%%
    %%% 2: Estimate parameters for law of motion of the Factors %%%
    %%% ======================================================= %%%  
    ff=cat(2,zeros(r,1),F(2:T,:)');                                         % F_t
    xx=cat(2,ML_lag(xitT(:,1:pr),1)',zeros(pr,1));                          % F_{t-1}
    EF1=zeros(r,pr); EF1F1=zeros(pr,pr); EF=zeros(r,r);                     % initialize
    for tt=cc+2:T  % ------------------------------------------------------ % \sum_{t=1}^T
        EF1=EF1+ff(:,tt)*xx(:,tt-1)'+PtT(1:r,r+1:pr1,tt);                   % E(F_t F_{t-1})
        EF1F1=EF1F1+xx(:,tt-1)*xx(:,tt-1)'+PtT(1:pr,1:pr,tt-1);             % E(F_{t-1} F_{t-1})
        EF=EF+ff(:,tt)*ff(:,tt)'+PtT(1:r,1:r,tt);                           % E(F_t F_t)
    end   
    A(1:r,1:pr)=EF1/EF1F1;                                                  % parameter VAR factors        
    Q(1:r,1:r) = (EF - A(1:r,1:pr)*EF1') / (T-cc);                          % covariance
    if r>q; [P,M] = eigs(Q(1:r,1:r),q,'lm',OPTS); Q(1:r,1:r) = P*M*P'; end  % in case the factors are singular
    
    %%% ============================== %%%
    %%% 3: Covariance Prediction Error %%%
    %%% ============================== %%%    
    xx=xitT';
    R=zeros(N,N);                                                           % preallocates
    for tt=cc+1:T  % ------------------------------------------------------ % \sum_{t=1}^T
        eta=yy(:,tt)-Lambda*xx(:,tt);                                       % prediction error
        R = R+ eta*eta'+ Lambda*PtT(:,:,tt)*Lambda';                        % E(eta_t eta_t')
    end
    R = diag(diag(R/(T-cc)));
    
    %%%% ======================== %%%%
    %%%% ====                ==== %%%%
    %%%% ====    E - Step    ==== %%%%
    %%%% ====                ==== %%%%
    %%%% ======================== %%%%            
    [xitt,xittm,Ptt,Pttm,loglik]=ML_KalmanFilter2026(F00,P00,y,A,Lambda,R,Q);  % Kalman Filter
    [xitT,PtT]=ML_KalmanSmoother2026(y,A,xitt',xittm',Ptt,Pttm,Lambda,R);      % Kalman SMoother
    F=xitT(:,1:r);                                                          % jj-th step factors            
    F00=xitT(1,:)'; P00=eye(pr1);                                           % initial conditions    

    %%%% ================================= %%%%
    %%%% ====                         ==== %%%%
    %%%% ====    Check convergence    ==== %%%%
    %%%% ====                         ==== %%%%
    %%%% ================================= %%%%
    delta_loglik = abs(loglik - loglik1);                                   % |logL(t) - logL(t-1)|
    avg_loglik = (abs(loglik) + abs(loglik1) + 10^(-3))/2;                  % avg = (|logL(t)| + |logL(t-1)|)/2
    if (delta_loglik / avg_loglik) < tresh; break; end                      % convergence if |f(t) - f(t-1)| / avg < threshold
    if jj>1; if loglik - loglik1<0; break; end; end                         % stop algorithm in case likelihood decrease 
    loglik1=loglik;                                                         % store log-likelihood    
end

if r>q; G = P*sqrt(M); else; G=eye(r); end 
% =========================================================================





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_KalmanFilter2026 - Kalman Filter Estimation of a Factor Model
% 
%  THE MODEL
% 
%   X_t = beta1 + beta2*t + C S_t + e_t,        e_t ~ N(0,R), R is diagonal
%   S_t = mu + A S_{t-1} + u_t,                 u_t ~ N(0,Q)
% 
%       S_t|X^{t-1} ~ N( S_{t|t-1} , P_{t|t-1} )
%       X_t|X^{t-1} ~ N( X_{t|t-1} , H_{t|t-1} )
% 
% 
%  THE PROCEDURE
% 
% [xitt,xittm,Ptt,Pttm,loglik]=ML_KalmanFilter(initx,initV,x,A,C,R,Q,mu,beta)
% 
% INPUTS:
%   x - the data
%   C - the observation matrix 
%   A - the system matrix
%   Q - the system covariance 
%   R - the observation covariance
%   initx - the initial state (column) vector 
%   initV - the initial state covariance 
%   mu   - constant in transition equation (optional)
%   beta - constant and linear trend in observation equation (optional)
% OUTPUTS:
%   xittm = S_{t|t-1}
%    Pttm = P_{t|t-1}
%    xitt = S_{t|t} 
%     Ptt = P_{t|t}
%  loglik = value of the log-likelihood
% 
% Matteo Luciani (matteoluciani@yahoo.it)

function [xitt,xittm,Ptt,Pttm,loglik]=ML_KalmanFilter2026(initx,initV,x,A,C,R,Q,mu,beta)

[T,N]=size(x);                                                              % Number of Observations
r=size(A,1);                                                                % Number of states
xittm=[initx zeros(r,T)];                                                   % S_{t|t-1}
xitt=zeros(r,T);                                                            % S_{t|t} 
Pttm=cat(3,initV,zeros(r,r,T));                                             % P_{t|t-1}
Ptt=zeros(r,r,T);                                                           % P_{t|t}
loglik=0;                                                                   % Initialize the log-likelihood
y=x';                                                                       % transpose for convenience

if nargin<8; mu=zeros(r,1); end
if nargin<9; beta=zeros(1,2); end

for j=1:T
    
    %%% ============= %%%
    %%% Updating Step %%%
    %%% ============= %%%   
    X=C*xittm(:,j) + beta*[1;j];                                            % X_{t|t-1} - Prediction
    H=C*Pttm(:,:,j)*C'+R;                                                   % H_{t|t-1} - Conditional Variance of the Observable   
    try 
        Hinv = H\eye(N);    
    catch
        Hinv = inv(H+ eye(N)*1e-10);    
    end
    
    K = Pttm(:,:,j)*C'*Hinv;                                                % Kalman gain
    e = y(:,j) - X;                                                         % error (innovation)
    xitt(:,j)=xittm(:,j)+K*e;                                               % S_{t|t}
    Ptt(:,:,j) = (eye(r) - K*C)*Pttm(:,:,j);                                % P_{t|t}   

    %%% =============== %%%
    %%% Prediction Step %%%
    %%% =============== %%%
    xittm(:,j+1)=mu+A*xitt(:,j);                                            % S_{t|t-1} - States
    Pttm(:,:,j+1)=A*Ptt(:,:,j)*A'+Q;                                        % P_{t|t-1} - Conditional Variance of the State 
    loglik = loglik + 0.5*(log(det(Hinv))  - e'*Hinv*e);                    % Log Likelihood   
    
end

xitt=xitt';
xittm=xittm';
% =========================================================================





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_KalmanSmoother2026 - Kalman Smoother for a Factor Model
% 
%  THE MODEL
%  
%   X_t = C S_t + e_t,              e_t ~ N(0,R), R is diagonal
%   S_t = A S_{t-1} + u_t,          u_t ~ N(0,Q)
% 
%   S_t|X^{t-1} ~ N( S_{t|t-1} , P_{t|t-1} )
%   X_t|X^{t-1} ~ N( X_{t|t-1} , H_{t|t-1} )
% 
%  
%  THE PROCEDURE
%  
% [xitT,PtT,PtTm]=ML_KalmanSmoother2026(x,A,xitt,xittm,Ptt,Pttm,C,R)
% 
% INPUTS:
%       A - the system matrix
%    xitt - S_{t|t}
%   xittm - S_{t|t-1}
%     Ptt - P_{t|t}
%    Pttm - P_{t,t-1|t-1}
%       C - the observation matrix 
%       R - the observation covariance
% OUTPUTS:
%    xitT = S_{t|T}
%     PtT = P_{t|T} 
%    PtTm = P_{t,t-1|T} 
%     
% This use the formulas in Chapter 4 of Durbin and Koopman "Time Series Analysis by
% State Space Methods" (second edition)
% Matteo Luciani (matteoluciani@yahoo.it)

function [xitT,PtT]=ML_KalmanSmoother2026(x,A,xitt,xittm,Ptt,Pttm,C,R)

[T]=size(xitt,2);                                                           % Number of Observations
[n, r]=size(C);                                                             % Number of Variables and Number of States
Pttm=cat(3,zeros(r,r,1),Pttm(:,:,1:end-1));                                 % P_{t|t-1}, remove the last observation because it has dimension T+1
xittm=[zeros(r,1) xittm(:,1:end-1)];                                        % S_{t|t-1}, remove the last observation because it has dimension T+1
xitT=[zeros(r,T)  xitt(:,T)];                                               % S_{t|T} 
PtT=cat(3,zeros(r,r,T),Ptt(:,:,T));                                         % P_{t|T} 
y=[zeros(n,1) x'];                                                          % transpose for convenience

rr=zeros(r,T+1); N=zeros(r,r,T+1); L=zeros(r,r,T+1);                        % Preallocates

CPCR = C*Pttm(:,:,T+1)*C'+R;                                                % useful matrices
iCPCR = eye(n)/CPCR;                                                        % ---------------
CiCPCRC = C'*iCPCR *C;                                                      % ---------------
L(:,:,T+1) = A - A*Pttm(:,:,T+1)*CiCPCRC;                                   % ---------------

for tt = T:-1:2    
    CPCR = C*Pttm(:,:,tt)*C'+R;                                             % useful matrices
    iCPCR = eye(n)/CPCR;                                                    % ---------------
    CiCPCRC = C'*iCPCR *C;                                                  % ---------------    
    rr(:,tt-1) = C'*iCPCR * ( y(:,tt) - C*xittm(:,tt) ) + L(:,:,tt)'*rr(:,tt);    
    N(:,:,tt-1) = CiCPCRC + L(:,:,tt)'*N(:,:,tt)*L(:,:,tt);    
    L(:,:,tt) = A - A*Pttm(:,:,tt)*CiCPCRC;             
    xitT(:,tt) = xittm(:,tt) + Pttm(:,:,tt)*rr(:,tt-1);                     % state        
    PtT(:,:,tt) = Pttm(:,:,tt) - Pttm(:,:,tt)*N(:,:,tt-1)*Pttm(:,:,tt);     % covariance of the state 
end

PtT(:,:,1)=[];
xitT=xitT(:,2:end)';
% =========================================================================





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_VAR - Estimates a VAR(k) for a vector of variables y
% if y is a single variable it estimates an AR model with OLS
%
% [A,u,AL,CL]=ML_VAR(y,det,jj);
%   y  = vector of endogenous variables
%   k  = number of lags
%   det = 0 no constant
%   det = 1 constant
%   det = 2 time trend
%   det = 3 constant + time trend
%   A  = Matrix of coefficients for the reduced form 
%   u  = Vector of Residuals
%
%  y_t=A(L)y_{t-1}+u_t
%  y_t=C(L)u_t
% 
% written by Matteo Luciani (matteoluciani@yahoo.it)

function [A,u,AL,CL,C1]=ML_VAR(y,k,det,s)
[T, N] = size(y);

%%% Building Up the vector for OLS regression %%%
yy=y(k+1:T,:);
xx=NaN(T-k,N*k);
for ii=1:N    
    for jj=1:k
        xx(:,k*(ii-1)+jj)=y(k+1-jj:T-jj,ii);
    end
end

%%% OLS Equation-By-Equation %%%
if det==0; ll=0; elseif  det==3; ll=2; else; ll=1; end
A=NaN(N*k+ll,N); u=NaN*yy;
for ii=1:N
    [A(:,ii),u(:,ii)]=ML_ols(yy(:,ii),xx,det);
end

At=A; if det==3; At(1:2,:)=[]; elseif det==1||det==2; At(1,:)=[]; end
AL=NaN(N,N,k); for kk=1:k; AL(:,:,kk)=At(1+kk-1:k:end,:)'; end  


%%% Impulse Responses %%%
if nargin<4;s=20; end
CL(:,:,1) = eye(N);
for ss=2:s
    CL(:,:,ss) = 0;
    for ii = 1:min(ss-1,k)        
        temp3=AL(:,:,ii)*CL(:,:,ss-ii);        
        CL(:,:,ss)=CL(:,:,ss)+temp3;
    end
end

C1=eye(N); for ii=1:k; C1=C1-AL(:,:,ii); end; C1=inv(C1);                     % Long Run Multipliers
% =========================================================================





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_edynfactors2 - Estimation of Dynamic Factors Innovation 
% 
% Principal Components of the Residuals of a VAR(p) estimated on r Static Factors
% Refined version of ML_edynfactors in which svd rather than eigs is used
% [etahat, G]=ML_edynfactors2(y,q)
% The Model:
%       X(t) = lambda*F(t) + xsi(t)
%       F(t) = A(L)*F(t-1) + epsilon
%       epsilon = G*eta
%       where G is r by q
% Outputs:
%       etahat = estimates of dynamic Factors Innovations
%       G = estimates of G
% Inputs:
%       y = epsilon
%       q = number of dynamic factors
%

% Written by Matteo Luciani (matteoluciani@yahoo.it)

function [eta, G]=ML_edynfactors2(y,q)
opt.disp=0;
N=size(y,2);
sigma=cov(y);                                                               % Variance Covariance Matrix of VAR Residuals
if q<N; [P,M]=eigs(sigma,q,'LM',opt); else [P,M]=eig(sigma); end            % Eigenvalue eigenvectors decomposition
eta=y*P*(M^-.5);                                                            % Dynamic Factor Innovations            
G=P(:,1:q)*(M^.5);
% =========================================================================





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================
% ML_VAR_companion_matrix - build the companion matrix for a VAR
% PHI=ML_VAR_companion_matrix(AL,det);
% Y(t)=AL(:,:,1)*Y(t-1)+...+AL(:,:,p)*Y(t-p)+u(t)
% For a VAR(3) with k variables:
% | AL(:,:,1)   AL(:,:,2)   AL(:,:,3)   |
% | eye(k,k)    zeros(k,k)  zeros(k,k) |
% | zeros(k,k)  eye(k,k)    zeros(k,k) |
%

% Written by Matteo Luciani (matteoluciani@yahoo.it)

function PHI=ML_VAR_companion_matrix(A)

s=size(A);
if length(s)==1; PHI=A;                             % AR(1)
elseif length(s)==2;                        
    if s(2)==1; p=s(1,1);                           % Autoregressive model
        if  p==2; PHI=[A'; 1 0];                    % AR(2)   
        else PHI=diag(ones(1,p-1),-1); PHI(1,:)=A'; % AR(p)
        end
    else PHI=A; end                                 % VAR(1)
else    
    k=s(1,1); p=s(1,3);
    PHI=[];
    for i=1:p; PHI=[PHI A(:,:,i)]; end;
    PHI=[PHI;eye(k*(p-1),k*(p-1)) zeros(k*(p-1),k)];
end
% =========================================================================