function MM=BL_SimulationsI0_2026_Eff_aux(T,N,r,q,tau,mu,delta,theta, stT,perturb)


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

if     floor(tau)==1; tau=tau-1; cut=1; 
elseif floor(tau)==2; tau=tau-2; cut=2;
elseif floor(tau)==3; tau=tau-3; cut=3;    
else;                            cut=0; 
end
TT1=diag(rand(N,1)+0.5);                                                    % variance idiosyncratic shocks
                % --------------------------------------------------------- % covariance matrix idiosyncratic shocks
if tau==0; TT2=eye(N);                                                      % Idio not correlated 
elseif tau==100; TT1=eye(N); TT2=BaiLiao(T,N);                                  % this is Bai and Liao
else; TT2=toeplitz(tau.^(0:N-1));                                           % correlated idio
end                                                                         %
if cut >0
    if cut==1; nc=5;                                                        % keep only 5 correlations
    elseif cut==2; nc=round(N^(1/2));                                       % keep only n^(1/2) correlations
    elseif cut==3; nc=round(N^(1/5));                                       % keep only n^(1/5) correlations
    end
    TT2=TT2-triu(TT2,nc)-tril(TT2,-nc); 
end                          
                % --------------------------------------------------------- %
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
LL1=1+randn(N,r);                                                             % factor loadings
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
lambda=W*sqrt(M);                                                          % identified factor loadings
f=chi*W/sqrt(M);                                                           % identified common factors
if perturb==1; f=f*chol(eye(r)+T^-0.5)'; end                                     % perturb common factors
chi=f*lambda';
x=chi+e;                                                                    % variables



%%%%% ====================================== %%%%%
%%%%% ==== ---------------------------- ==== %%%%%
%%%%% ==== ---                      --- ==== %%%%%
%%%%% ==== ---  PART 2: Efficiency  --- ==== %%%%%
%%%%% ==== ---                      --- ==== %%%%%
%%%%% ==== ---------------------------- ==== %%%%%
%%%%% ====================================== %%%%%
GG = sqrt(TT1)*(TT2)*sqrt(TT1);                                             % Covariance Matrix idiosyncratic
DD = diag(diag(GG));                                                        % Diagonal elements
OD = GG-DD;                                                                 % extra diagonal elements

iDD=eye(N)/DD;                                                              % Useful matrices
LL=lambda'*lambda;                                                          % ---------------
iLL=eye(r)/LL;                                                              % ---------------
LiddL=lambda'*iDD*lambda;                                                   % ---------------
iLiddL=eye(r)/LiddL;                                                        % ---------------
LiDDggL=lambda'*iDD*GG*iDD*lambda;                                          % ---------------
    
WPC = N * ( iLL * lambda'*GG*lambda * iLL);                                 % asymptotic covariance PCA
WEM = N * (iLiddL * LiDDggL * iLiddL);                                      % asymptotic covariance EM

%%% Alternative version that imposes our identification %%%
iLL2=eye(r)/M;
PSIxi=eye(r)/(W'*iDD*W);                                                    % From proposition 3
PHIxi=W'*iDD*GG*iDD*W;                                                      % From proposition 3
WEM2 = sqrt(iLL2) * PSIxi * PHIxi * PSIxi *sqrt(iLL2);                      % asymptotic covariance EM
WPC2 = sqrt(iLL2) * (W' * GG * W) * sqrt(iLL2);                               % From comment 2 to Proposition 5

AB = N * iLiddL * LiDDggL * iLL;                                        % Decompose Asymptotic variances
C1 = WPC+WEM-AB-AB';
C2 = AB-WEM;
WiDD=W'*iDD;
Psi = eye(r) / (WiDD*W);
ERR = Psi * WiDD*OD*W - Psi * WiDD*OD*WiDD' * Psi;

MM = real([max(max(ERR)) min(min(ERR)) sort(real(eig(ERR)))' sort(eig(C1+C2+C2'))' sort(eig(WPC-WEM),'descend')' sort(eig(WPC2-WEM2),'descend')']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate idiosyncratic shocks as in Bai and Liao (2016) JoE
function TT2=BaiLiao(T,N)

e=randn(T,N);
ABC=[ones(N,1) .7*randn(N,3)];

u(:,1)=e(:,1);
u(:,2)=e(:,[2 1])*ABC(2,1:2)';
u(:,3)=e(:,[3 2 1])*ABC(3,1:3)';

for ii=4:N; u(:,ii)=e(:,ii:-1:ii-3)*ABC(ii,:)'; end

TT2=corr(u);                                                                % Correlation idiosyncratic shocks
alpha = 1e-4;  TT2 = (1-alpha)*TT2 + alpha*eye(N);                          % Mild shrinkage toward identity matrix for regularization

% Sigma=eye(N);
% 
% autoval=-1;
% kk=0;
% while autoval<0; kk=kk+1;
% 
%     if kk==10; disp(kk); break; end
% 
%     abc=0.7*ML_Standardize(randn(N,3));   
% 
%     for ii=1:N   
%         Sigma(ii,ii)=1+abc(ii,:)*abc(ii,:)';        
%         if ii>3
%             Sigma(ii,ii-1)=[  1 abc(ii-1,1:2)]* abc(ii,:)';
%             Sigma(ii,ii-2)=[0 1 abc(ii-2,1)]  * abc(ii,:)';
%             Sigma(ii,ii-3)=[0 0 1]            * abc(ii,:)';
%         end
%     end
% 
%     Sigma(2,1)=[  1 0 0]* abc(2,:)';
%     Sigma(3,2)=[  1 abc(2,1:2)]    * abc(3,:)';
%     Sigma(3,1)=[0 1 abc(1,1)]  * abc(3,:)';
% 
%     TT2=Sigma+Sigma'-diag(diag(Sigma));
%     autoval=min(eig(TT2));
% end


% if tau==1
%     abc=0.7*ML_Standardize(randn(N,3));
%     vv=v;
%     ii=0; v(:,ii+1)=vv(:,ii+1);
%     ii=1; v(:,ii+1)=vv(:,ii+1)+vv(:,ii:ii)*abc(ii,1)';
%     ii=2; v(:,ii+1)=vv(:,ii+1)+vv(:,ii-1:ii)*abc(ii,1:2)';
%     for ii=3:N-1
%         v(:,ii+1)=vv(:,ii+1)+vv(:,ii-2:ii)*abc(ii,:)';
%     end   
% end



