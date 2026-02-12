% BL_AsymptoticVarianceChiPCA - Estimates asymptotic variances for PCA estimator
% 
% [SigmaChi,WF,VL]=BL_AsymptoticVarianceChiPCA(x,F,L) 
% 
%    x - DATA (STANDARDIZED, if the model was estimated on standardized data)
%    F - estimated factors
%    L - estimated loadings
%    m - M=2m+1 are the number of autocovariances used in the kernel
% 
% To compute the covariance of chi you have to standardize the data first
% the user then need to multiply for the variance outside the function
% 
function [SigmaChi,WF,VL]=BL_AsymptoticVarianceChiPCA_sim(x,sx,F,L,m,nv)

if nargin==5; nv=5; end                                                     % default number of variances to compute

[N,r]=size(L); T=size(F,1);                                                 % relevant sizes                 
ONEr = ones(1, r);
sx2 = sx .^ 2;  % Pre-compute squared variances
GammaXi = zeros(N, N, T);
CL=zeros(T,N,nv); CF=CL; SigmaChi=CL;  WF=cell(nv,1); VL=WF;                    % ------------
for ii=1:nv; VL{ii}=zeros(r,r,N); end                                         % ------------

xi=x-F*L';                                                                  % estimated idiosyncratic components
for tt=1:T; GammaXi(:,:,tt) = xi(tt,:)'*xi(tt,:); end                       % TV covariance idiosyncratic

GammaF = F'*F/(T - 1);  iGf= GammaF\eye(r);                                 % usefull quantities 
mGx=mean(GammaXi,3);    R = diag(diag(mGx));                                % ------------------
LL=L'*L;                iLL=(LL/N)\eye(r);                                  % ------------------


%%% ======================================= %%%
%%%  Asymptotic variance of common factors  %%%
%%% ======================================= %%%
n=floor(N^(4/5));                                                           % number of variables to use for the HCC estimator
Weights=zeros(N,N);                                                         % construct weights
ww=[1./(1:n) zeros(1,N-n)];                                                 % -----------------
for ii=1:N; Weights(ii,ii:N)=ww(1:N+1-ii); end                              % -----------------
Weights=Weights-eye(N)+Weights';                                            % -----------------
Weights(n+1:N,:)=0; Weights(:,n+1:N)=0;                                     % -----------------
CC=mGx.*Weights;                                                            % Weighted variances

CC2=L'*CC*L / N;                                                            % useful quantity
WF{1} = iLL*CC2*iLL;                                        % ------------- % WHCC

if nv>1
    Weights(1:n,1:n)=1;
    CC=mGx.*Weights;                                                        % Compute unweighted version
    CC2=L'*CC*L / N;                                                        % useful quantity
    WF{2} = iLL*CC2*iLL;                                    % ------------- % WHCC
end

if nv>2; WF{3} = mean(diag(mGx)) *iLL; end      % ------------------------- % W*

if nv>3
    CC = L'*R*L / N;
    WF{4} = iLL*CC*iLL;                         % ------------------------- % W0
end

if nv>4
    WF{5} = zeros(r, r, T);
    diagGammaXi = zeros(N, T);
    for tt=1:T
        diagGammaXi(:, tt) = diag(GammaXi(:, :, tt));
    end
    for tt = 1:T
        CC = L' * diag(diagGammaXi(:, tt)) * L / N;
        WF{5}(:,:,tt)=iLL*CC*iLL;                           % ------------- % WHC
    end
end

%%% ======================================== %%%
%%%  Asymptotic variance of factor loadings  %%%
%%% ======================================== %%%
ONEr=ones(1,r);
mmGx=mean(diag(mGx));
for ii=1:N 
    UU = F.*(xi(:,ii)*ONEr);                                                % usefull quantity
    Gamma_k=(UU'*UU)/T; 
    if nv>2; VL{3}(:,:,ii)=iGf*mmGx; end            % --------------------- % V^*    
    if nv>3; VL{4}(:,:,ii)=iGf*mGx(ii,ii); end      % --------------------- % V^0
    if nv>4; VL{5}(:,:,ii)=iGf*Gamma_k*iGf; end     % --------------------- % V^HC
    for vv=1:m
        A=UU(vv:T,:)'*UU(1:T+1-vv,:)/T;        
        weight = 1 - vv / (m + 1);
        Gamma_k = Gamma_k + weight * (A + A');
    end
    VL{1}(:,:,ii)=iGf*Gamma_k*iGf;                  % --------------------- % V^HCC
end            
if nv>1; VL{2}=VL{1}; end


%%% ========================================= %%%
%%%  Asymptotic variance of common component  %%%
%%% ========================================= %%%

for kk = 1:min(nv, 4)
    for nn = 1:N
        % Compute F * VL{kk}(:,:,nn) * F' once for all time points
        FVF = F * VL{kk}(:, :, nn) * F';  % T x T matrix
        CL(:, nn, kk) = diag(FVF) / T;     % Extract diagonal
        
        % Scalar value, same for all time points
        LWL = L(nn, :) * WF{kk} * L(nn, :)';
        CF(:, nn, kk) = LWL / N;
        
        % Vectorized over time
        SigmaChi(:, nn, kk) = sx2(nn) * (CL(:, nn, kk) + CF(:, nn, kk));
    end
end

if nv>4
    for nn=1:N
        for tt=1:T            
            CL(tt,nn,5) = ( F(tt,:)*VL{5}(:,:,nn)*F(tt,:)' ) / T;
            CF(tt,nn,5) = ( L(nn,:)*WF{5}(:,:,tt)*L(nn,:)' ) / N;
            SigmaChi(tt, nn, 5) = sx2(nn) * (CL(tt, nn, 5) + CF(tt, nn, 5));
        end
    end
end 