
function [Z_chi, Z_F, Z_L,Z_chidag]=BL_ComputeZb(DFM,DGP,vfl,nv,G)


if nargin==2; vfl=1; end
try isempty(vfl); if isempty(vfl); vfl=1; end; catch; vfl=1; end  


try isempty(nv); if isempty(nv); nv=5; end; catch; nv=5; end                % default number of variances to compute

try isempty(G); if isempty(G); rq=0; else; rq=1; end; catch; rq=0; G=[]; end     

% Extract data
xx=DGP.xx;
sx=DGP.sx;
F=DFM.F;                                                                    % estimated factors
L=DFM.Lambda;                                                               % estimated loadings
chi=DFM.chi;                                                                % estimated common components
PtT=DFM.PtT;
[N,r]=size(L); T=size(F,1);                                                 % relevant sizes
m = floor(T^0.25);
[SigmaChi,WF,VL,SigmaChidag]=...                                            % Compute asymptotic variances
    BL_AsymptoticVarianceChiEM_sim2026(xx,sx,F,L,PtT,m,nv,G); 

%%% =========================================== %%%
%%%  Standardize Normal for chi, F, and Lambda  %%%
%%% =========================================== %%%
K=size(SigmaChi,3);
Z_chi=NaN(T,N,K); Z_L=[]; Z_F=[];                                           % preallocates

chi_diff = chi - DGP.chi;                                                   % Compute difference once
for kk = 1:K; Z_chi(:,:,kk) = chi_diff ./ sqrt(SigmaChi(:,:,kk)); end

if vfl==1
    Z_L=NaN(N,r,K); Z_F=NaN(T,r,K);                                         % preallocates
    % --------------------------------------------------------------------- % Compute differences once
    F_diff = F - DGP.f;                                                     % T x r
    L_diff = L - DGP.lambda;                                                % N x r
    sqrtN = sqrt(N);
    for kk=1:K
        for jj=1:r            
            sFW = sqrt(squeeze(WF{kk}(jj, jj, :)));                         % Extract diagonal elements for all t at once
            Z_F(:, jj, kk) = sqrtN * F_diff(:, jj) ./ sFW;
        end
    end
    
    sqrtT = sqrt(T);
    for kk=1:K
        for jj = 1:r            
            if isnumeric(VL{kk}) && ndims(VL{kk}) == 3                      % Try to vectorize over N if possible                
                sVL_all = sqrt(squeeze(VL{kk}(jj, jj, :)));  % N x 1
                Z_L(:, jj, kk) = sqrtT * L_diff(:, jj) ./ sVL_all;
            else                                                            % Fallback to loop if structure doesn't allow vectorization
                for nn = 1:N
                    sVL = sqrt(VL{kk}(jj, jj, nn));
                    Z_L(nn, jj, kk) = sqrtT * L_diff(nn, jj) / sVL;
                end
            end
        end
    end
end

if rq==1
    chidag_diff = DFM.chidag - DGP.chidag;
    Z_chidag=NaN(T,N,K);                                    % preallocates

    for kk=1:K 
        Z_chidag(:,:,kk) = chidag_diff ./ sqrt(SigmaChidag(:,:,kk));
    end
else
     Z_chidag=[];
end



