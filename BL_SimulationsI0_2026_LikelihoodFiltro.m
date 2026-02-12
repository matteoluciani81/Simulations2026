% BL_SimulationsI0_2026_LikelihoodFiltro - Study LL when initialization is not correct

clear; close all; clc
Cartella='/cmc/home/m1mxl04/Documents/';
addpath([Cartella 'ML']); addpath([Cartella 'ML/ML_Package']);
addpath('/cmc/home/m1mxl04/Documents/BL/code')
ML_graph_options

ST=0;


% ========================================================================= % Parameters for simulations
T=100;
N=100;
mu=.7;
delta=2.5;
tau=0.5;
theta=0.5; % signal to noise ratio
r=4;
q=4;
stT=3;
iota=0;
perturb=0;
p=1;
pr=p*r;
iter=50; 
tresh=10^(-4);

lgd={'$\ell(X_{nt}; \hat{{\varphi}}_n^{(k)})$','$\log\Delta \ell_k$'};
Filename=['LogLikSim_2026_' num2str(q) num2str(r)];

%%%%% ============================================== %%%%%
%%%%% ==== ------------------------------------ ==== %%%%%
%%%%% ==== ---                              --- ==== %%%%%
%%%%% ==== ---  PART 1: Simulate the model  --- ==== %%%%%
%%%%% ==== ---                              --- ==== %%%%%
%%%%% ==== ------------------------------------ ==== %%%%%
%%%%% ============================================== %%%%%

Tb = T + 10;                                                             % Precompute burn-in length
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
H1=NaN(r+1,r); for rr=1:r;  H1(:,rr)=ML_ols(f(:,rr),ff1,1); end             % Rotation matrix\
const=H1(1,:); H1=H1(2:end,:); H2=eye(r)/H1;                                % ---------------
AA=H1'*[A zeros(q,r-q);eye(r-q) zeros(r-q,r-q)]*H2';
G=H1'*[eye(q); zeros(r-q,r-q)];    
const2=((eye(r)-AA)*const')';
g=f; for tt=2:T; g(tt,:)=const2+(AA*g(tt-1,:)')'+uu(tt+10,:)*G'; end 




%%%%% ============================================== %%%%%
%%%%% ==== ------------------------------------ ==== %%%%%
%%%%% ==== ---                              --- ==== %%%%%
%%%%% ==== ---  PART 2: Estimate the model  --- ==== %%%%%
%%%%% ==== ---                              --- ==== %%%%%
%%%%% ==== ------------------------------------ ==== %%%%%
%%%%% ============================================== %%%%%

[xx,mx,sx]=ML_Standardize(x);                                               % standardize the data
MX=repmat(mx,T,1); SX=repmat(sx,T,1);                                       % Useful objects
DGP.xx=xx; DGP.mx=mx; DGP.sx=sx;                                            % store variances and mean of the data

% %%% Studia filtro %%%
% [F9,lambda9,A9,P9,Q9,R9]=ML_SS_DFM_I0(f,lambda,cat(3,AA,zeros(r,r)),G,e,N,p+1,r);              % State-Space representation
% [xitt9,xittm9,Ptt9,Pttm9]=ML_KalmanFilter2(F9,P9,xx,A9,lambda9,R9,Q9);        % Kalman Filter
% [~,PtT9]=ML_KalmanSmoother3(xx,A9,xitt9',xittm9',Ptt9,Pttm9,lambda9,R9);       % Kalman Smoother

Gamma=cov(xx);                                                              % covariance matrix of the standardize data
ijk=0;
for iota=0:.1:.9; ijk=ijk+1; % --------------------------------------------- % For different levels of miss-specification

    % ===================================================================== % Estimate factors and loadings with PCA
    [W,M] = eigs(Gamma, r,'LM'); W=W*diag(sign(W(1,:))');                   % eigenvalue-eigenvectors
    if iota~=0
        W=W+randn(N,r)*cholcov(toeplitz(iota.^(0:r-1)));                    % contaminates the eigenvector
        W=W/sqrt(diag(diag(W'*W)));                                         % ----------------------------
    end
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
    
    % ===================================================================== % Expectation Maximization Algorithm  %%%    
    AL0b=cat(3,AL0,zeros(r,r));                                             % add one lag to avoid computing PtTm
    [F0,lambda0,A0,P0,Q0,R0]=ML_SS_DFM_I0(f0,L0,AL0b,G0,e0,N,p+1,r);        % State-Space representation
    [xitT,PtT1,~,~,Ptt1,Pttm1,A1,L1,R1,Q1,G1,LogLik]=...                    % EM algorithm
        ML_efactorML4d3_MB(xx,F0,P0,A0,lambda0,R0,Q0,r,q,p,20,10^(-20),0,0);

    ll=1;
    delta_loglik = abs(LogLik(2:end,ll) - LogLik(1:end-1,ll));              % |logL(t) - logL(t-1)|
    avg_loglik = (abs(LogLik(2:end,ll)) + abs(LogLik(1:end-1,ll)) + 10^(-3))/2; % Average likelihood
    criterio=log(delta_loglik./avg_loglik);                                 % convergence criterion

    band=.1*( ML_diff(ML_minmax(LogLik(:,ll))') );
    band2=.1*ML_diff(ML_minmax(criterio)');

    MxVer(:,ijk)=LogLik(1:20,ll);                                           % Store LL
    Crit(:,ijk)=criterio;                                                   % Store information criterion
    
    % ===================================================================== % Chart likelihood and convergence criterion
    set(groot,'defaultAxesTickLabelInterpreter','latex');          
    axes('Parent',figure,'FontSize',12); ML_FigureSize,
    yyaxis right;
    pl2=plot(1:20,criterio,'-s','linewidth',1.5);
    set(gca,'ycolor','k');
    ylabel(lgd{2},'interpreter','latex','FontSize',16)
    ylim(ML_minmax(criterio)+[-band2 band2])

    yyaxis left;
    pl1=plot(0:20, LogLik(:,ll),'-o','linewidth',1.5);
    set(gca,'ycolor','k');
    ylabel(lgd{1},'interpreter','latex','FontSize',16)
    xlabel('iteration($k$)','interpreter','latex','FontSize',16)
    ylim(ML_minmax(LogLik(:,ll))+[-band band])
    xlim([-.5 20.5])    
    legend([pl1 pl2],lgd,'location','best','interpreter','latex','FontSize',16)    
    if ST==1        
        exportgraphics(gcf,[Filename 'iota' num2str(10*iota) '.pdf'],'ContentType','vector');
    end
end



% ========================================================================= % Chart likelihood and convergence criterion
ZZ=MxVer(:,2:10)./repmat(MxVer(:,1),1,9)-1;                                        % pctdev from good initialization

figure  % mock figure
plot(ZZ,'LineWidth',2); legend(num2str((0.1:.1:.9)'))

ZZ=MxVer(:,4:3:10)./repmat(MxVer(:,1),1,3)-1;

mzz=ML_round2(ML_min(ZZ),1,0,5);
tzz=ML_closest([1 2 2.5 5 10],-mzz/10,1); % tick space
ytick=ML_ytick(mzz,0,tzz);
set(groot,'defaultAxesTickLabelInterpreter','latex');
axes('Parent',figure,'FontSize',12); ML_FigureSize,
pl3=plot(ZZ,'-o','linewidth',1.5);
set(gca,'xtick',1:3:21,'XTickLabel',num2str((0:2:20)'),'ytick',ytick)
xlabel('iteration($k$)','interpreter','latex','FontSize',16)
xlim([.5 20.5])
ylim([ML_round2(ML_min(ZZ),1,0,5) 0.5 ])
gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',[.8 .8 .8],'linewidth',1);
legend(pl3,{'$\iota=0.3$','$\iota=0.6$','$\iota=0.9$'},'location','SE','interpreter','latex','FontSize',16)
if ST==1
    ax = gca;
    exportgraphics(gcf,[Filename 'CompareIota.pdf'],'ContentType','vector');
end
