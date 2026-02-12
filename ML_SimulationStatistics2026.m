function [TR01,TR02,EE,EA,R11,R12,R21,R22]=ML_SimulationStatistics2026(F0,FF,L0,LL,chi0,chat)

m=size(FF,2);                                                               % Number of models
EE=repmat(NaN*chi0,1,1,1,m); EA=EE; TR01=NaN(1,m); TR02=TR01;             	% preallocates
R11=NaN(1,size(FF{1},2),m); R12=R11; R21=R11; R22=R11;

for jj=1:m
    num=trace((F0'*FF{jj})/(FF{jj}'*FF{jj})*(FF{jj}'*F0));                  % Trace statistics for factors
    den=trace(F0'*F0);                                                      % ----------------------------
    TR01(1,jj)=num/den;                                                     % ----------------------------
    [~,~,R11(:,:,jj)]=canoncorr(F0,FF{jj});                                 % canonical correlation factors    
    R21(:,:,jj)=diag(corr(F0,FF{jj}))';                                     % pairwise (standard) correlation
    
    num=trace((L0'*LL{jj})/(LL{jj}'*LL{jj})*(LL{jj}'*L0));                  % Trace statistics for loadings
    den=trace(L0'*L0);                                                      % ----------------------------
    TR02(1,jj)=num/den;                                                     % ----------------------------
    [~,~,R12(:,:,jj)]=canoncorr(L0,LL{jj});                                 % canonical correlation loadings
    R22(:,:,jj)=diag(corr(L0,LL{jj}))';                                     % pairwise (standard) correlation
    
    EE(:,:,1,jj)=(chat{jj}-chi0).^2;                                        % squared error
    EA(:,:,1,jj)=abs((chat{jj}-chi0));                                      % absolute error
end