% BL_AsymptoticVarianceChiEM - Estimates asymptotic variances for EM estimator
% 
% [SigmaChi,WF,VL]=BL_AsymptoticVarianceChiEM(x,F,L,PtT,m) 
% 
%    x - DATA (STANDARDIZED, if the model was estimated on standardized data)
%    F - estimated factors
%    L - estimated loadings
%  PtT -  Conditional covariance of the factors
%    m - M=2m+1 are the number of autocovariances used in the kernel
% 
% To compute the covariance of chi you have to standardize the data first
% the user then need to multiply for the variance outside the function
% 

function [SigmaChi,WF,VL,SigmaChidag]=BL_AsymptoticVarianceChiEM_sim(x,sx,F,L,PtT,m,nv,G)

try isempty(nv); if isempty(nv); nv=5; end; catch; nv=5; end                % default number of variances to compute
try isempty(G); if isempty(G); rq=0; else; rq=1; end; catch; rq=0; end     

[N,r]=size(L); T=size(F,1);                                                 % relevant sizes                 
ONEr = ones(1, r);                                                          % Useful constant
GammaXi=zeros(N,N,T);                                                       % preallocates
CL=zeros(T,N,nv); CF=CL; SigmaChi=CL;  WF=cell(nv,1); VL=WF;                % ------------
for ii=1:nv; VL{ii}=zeros(r,r,N); end                                       % ------------

xi=x-F*L';                                                                  % estimated idiosyncratic components
FtL = zeros(N, N, T);
for tt = 1:T
    FtL(:, :, tt) = L * PtT(:, :, tt) * L';
end
for tt=1:T; GammaXi(:, :, tt) = xi(tt, :)' * xi(tt, :) + FtL(:, :, tt); end % TV covariance idiosyncratic             

GammaF = F'*F/(T-1);    iGf = GammaF\eye(r);                                % usefull quantities 
mGx=mean(GammaXi,3);    R = diag(diag(mGx));    iR=R\eye(N);                % ------------------
LL=L'*L;                iLL=(LL/N)\eye(r);                                  % ------------------
LiRL=L'*iR*L;           iLiRL=(LiRL/N)\eye(r);                              % ------------------

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
iRCC = iR*CC;
CC2=L'*iRCC*iR*L / N;                                                       % useful quantity
WF{1} = iLiRL*CC2*iLiRL;                    % ----------------------------- % W^HCC

if nv>1
    Weights(1:n,1:n)=1;
    CC=mGx.*Weights;                                                        % Compute unweighted version
    iRCC = iR*CC;
    CC2=L'*iRCC*iR*L / N;                                                   % useful quantity
    WF{2} = iLiRL*CC2*iLiRL;                % ----------------------------- % W^HCC
end

if nv>2; WF{3} = mean(diag(mGx)) *iLL; end  % ----------------------------- % W^*

if nv>3; WF{4} = iLiRL; end                 % ----------------------------- % W^0

if nv>4
    WF{5} = zeros(r, r, T);
    diagGammaXi = zeros(N, T);
    for tt = 1:T
        diagGammaXi(:, tt) = diag(GammaXi(:, :, tt));
    end
    
    for tt = 1:T
        CC = L' * iR * diag(diagGammaXi(:, tt)) * iR * L / N;
        WF{5}(:,:,tt)=iLiRL*CC*iLiRL;       % ----------------------------- % W^HC
    end
end

%%% ======================================== %%%
%%%  Asymptotic variance of factor loadings  %%%
%%% ======================================== %%%

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
    VL{1}(:,:,ii)=iGf*Gamma_k*iGf;  % ------------------------------------- % V^HCC
end            
if nv>1; VL{2}=VL{1}; end

%%% ========================================= %%%
%%%  Asymptotic variance of common component  %%%
%%% ========================================= %%%

sx2 = sx .^ 2;  % Pre-compute squared variances

for kk = 1:min(nv, 4)
    for nn = 1:N
        % Vectorized over time
        FVF = F * VL{kk}(:, :, nn) * F';  % T x T matrix
        CL(:, nn, kk) = diag(FVF) / T;
        
        LWL = L(nn, :) * WF{kk} * L(nn, :)';
        CF(:, nn, kk) = LWL / N;
        
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


if rq==1    % ------------------------------------------------------------- % This is the singular case
    
    SigmaChidag = zeros(T, N, nv);
    
    Gdag=(G'*G)\G';                                                             % left inverse G
    q=size(G,2);
    GdagG = Gdag' * G';
    for nn = 1:N
        for kk = 1:min(nv, 4)
            % Pre-compute common terms
            iWF = WF{kk} \ eye(r);
            WH = (G' * iWF * G) \ eye(q);
            GWHG = G * WH * G';
            LWL = L(nn, :) * GWHG * L(nn, :)';
            CF(:, nn, kk) = LWL / N;
            
            for tt = 1:T
                fdag = F(tt, :) * GdagG;
                CL(tt, nn, kk) = (fdag * VL{kk}(:, :, nn) * fdag') / T;
                SigmaChidag(tt, nn, kk) = sx2(nn) * (CL(tt, nn, kk) + CF(tt, nn, kk));
            end
        end
    end
    
    if nv>4
        for nn=1:N
            for tt=1:T
                fdag = F(tt, :) * GdagG;
                CL(tt,nn,5) = ( fdag*VL{5}(:,:,nn)*fdag ) / T;
                
                iWF = WF{5}(:, :, tt) \ eye(r);
                WH = (G' * iWF * G) \ eye(q);
                CF(tt,nn,5) = ( L(nn,:)*G*WH*G'*L(nn,:)' ) / N;
                
                SigmaChidag(tt, nn, 5) = sx2(nn) * (CL(tt, nn, 5) + CF(tt, nn, 5));
            end
        end
    end
else; SigmaChidag=[];
end         % ------------------------------------------------------------- %



% n=1;
% 
% CC=mGx;
% for ii=1:N
%     J=randperm(N,N-n);
%     if sum(J==ii)==1; J(J==ii)=[]; end
%     CC(ii,J)=0;
% end

% CC=mGx-triu(mGx,round(n/2))-tril(mGx,-round(n/2));