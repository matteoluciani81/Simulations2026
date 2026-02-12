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
function [Z_chi, Z_F, Z_L]=BL_ComputeZ_PCA2026(DFM,DGP,vfl,nv)

if nargin==2; vfl=1; end

xx=DGP.xx;
sx=DGP.sx;
chi_dgp = DGP.chi;                                                          % Simulated common component
F=DFM.F;                                                                    % estimated factors
L=DFM.Lambda;                                                               % estimated loadings
chi=DFM.chi;                                                                % estimated common components
[N,r]=size(L); T=size(F,1);                                                 % relevant sizes
m=floor(T^.25);
[SigmaChi,WF,VL]=BL_AsymptoticVarianceChiPCA_sim2026(xx,sx,F,L,m,nv);

%%% =========================================== %%%
%%%  Standardize Normal for chi, F, and Lambda  %%%
%%% =========================================== %%%
K=size(SigmaChi,3);
Z_chi=zeros(T,N,K); Z_L=[]; Z_F=[];                                           % preallocates
chi_diff = chi - chi_dgp;  % Compute difference once
for kk = 1:K; Z_chi(:, :, kk) = chi_diff ./ sqrt(SigmaChi(:, :, kk)); end

if vfl==1
    Z_L=zeros(N,r,K); Z_F=zeros(T,r,K);                                         % preallocates
    F_diff = F - DGP.f;  % Compute difference once
    sqrtN = sqrt(N);
    for kk=1:K
        for jj=1:r
            if ndims(WF{kk}) == 3
                % Time-varying case: WF is r x r x T
                sFW = sqrt(WF{kk}(jj, jj, :));  % Extract 1 x 1 x T
                sFW = squeeze(sFW);              % Convert to T x 1
                Z_F(:, jj, kk) = sqrtN * F_diff(:, jj) ./ sFW;
            else
                % Time-invariant case: WF is r x r
                sFW = sqrt(WF{kk}(jj, jj));
                Z_F(:, jj, kk) = sqrtN * F_diff(:, jj) / sFW;
            end
        end
    end
    
    L_diff = L - DGP.lambda;  % Compute difference once
    sqrtT = sqrt(T);
    for kk=1:K
        for nn=1:N
            for jj=1:r
                sVL = sqrt(VL{kk}(jj, jj, nn));
                Z_L(nn, jj, kk) = sqrtT * L_diff(nn, jj) / sVL;
            end
        end
    end
end