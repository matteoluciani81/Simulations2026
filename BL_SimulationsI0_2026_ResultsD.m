% BL_SimulationsI0_2026_ResultsD - Simulations I0 for distribution

clear all; close all; clc
Cartella='/cmc/home/m1mxl04/Documents/';
addpath([Cartella 'ML']); addpath([Cartella 'ML/ML_Package']);

ML_graph_options
NomeFile='BL_SimulationI0_2026_D_';
FolderD=[pwd '/StatisticsD/'];
stampa=0;
salva=0;
M=2;                                                                        % number of estimators
ALPHA=norminv([.75 .84 .9 .95 .975 .99]); ALPHA=sort([0 ALPHA -ALPHA]);     % Percentiles of N(0,1) for coverage
na=length(ALPHA);
Coverage=NaN(length(ALPHA),6); Coverage2=Coverage;                          % ------------



%%% ======================================== %%%
%%%  Parameters for Data Generating Process  %%%
%%% ======================================== %%%

QR=[4 4; 2 4];                               % [q r]
NT=[75 75; 100 100; 200 200; 300 300; 500 500]; % [N T]
%           mu      tau     delta   theta   stT 
%           1       2       3       4       5  
ParamDGP=[  0.7     0       0       0.5     0;...
            0.7     0.5     2.5     0.5     0;...
            0.7     0.5     0       0.5     0;...
            0.7     0       2.5     0.5     0;...            
            0.5     0.5     2.5     0.5     0;...
            0.7     0.5     2.5     1       0;...
            0.7     0.5     2.5     0.5     1;...
            0.7     0.5     2.5     0.5     2;...
            0.7     0.5     2.5     0.5     3;...
            0.7     0.5     2.5     0.5     4];
ParamDGP(:,6)=0;   % iota
ParamDGP(:,7)=0;   % pi

% ParamDGP(5:7,:)=[];

kk=0;
for qr=1
    for nt=1:size(NT,1)
        for pp=1:size(ParamDGP,1)
            kk=kk+1;
            % kk   q   r   N   T   mu   tau   delta   theta   stT   iota   pi
            %  1   2   3   4   5    6     7       8       9    10     11   12
            PARAM(kk,:)=[kk QR(qr,:) NT(nt,:) ParamDGP(pp,:)];
        end
    end
end

nsim=size(PARAM,1);
TabellaC=cell(3,2); TabellaS=TabellaC;
for ll=1:4      % ------------------------------------------------- % Different ways of computing AsyVAR
    for vv=1:2 % ------------------------------------------------- % {EM, PCA}
        TabellaC{ll,vv}=NaN(nsim,13);
        TabellaS{ll,vv}=NaN(nsim,4);
    end
end

ii=0; jj=0;
for kk=1:size(PARAM,1)  % ------------------------------------------------- % For each DGP   
    % disp(kk)
    inizio=now;
    % --------------------------------------------------------------------- % Extract parameters
    q=PARAM(kk,2); r=PARAM(kk,3); N=PARAM(kk,4); T=PARAM(kk,5);             % Retreive parameters DGP
    mu=PARAM(kk,6); tau=PARAM(kk,7); delta=PARAM(kk,8); theta=PARAM(kk,9);  % -----------------------
    stT=PARAM(kk,10); iota=PARAM(kk,11); pi=PARAM(kk,12);
    % --------------------------------------------------------------------- % Determine distribution type
    if stT==1; STUD='T';                                                    % ... and idiosyncratic shocks
    elseif stT==2; STUD='L';                                                % ----------------------------
    elseif stT==3; STUD='AL';                                               % ----------------------------
    elseif stT==4; STUD='ST';                                               % ----------------------------
    else STUD=[]; end                                                       % ----------------------------

    tipo=['q' num2str(q) 'r' num2str(r) 'T' num2str(T) 'N' num2str(N) ...   % Name for simulation process
        'mu' num2str(10*mu) 'tau' num2str(10*tau) ...                       % ----------------------------
        'delta' num2str(10*delta) 'theta' num2str(10*theta) ...             % ----------------------------
        'iota' num2str(10*iota) 'pi' num2str(pi) STUD];                     % ----------------------------
    
    if ~exist([FolderD NomeFile tipo '.mat'], 'file')
        jj=jj+1; 
        disp(['NO     - ' tipo]);
        TabellaParNO(jj,:)=PARAM(kk,:); 
        continue;
    end
    
    disp(['YES    - ' tipo]);                                           % Display DGP
    ii=ii+1;        
    load([FolderD NomeFile tipo]);
    
    for ll=1:4      % ------------------------------------------------- % Different ways of computing AsyVAR
        for vv=1:2  % ------------------------------------------------- % {EM, PCA}
            TabellaC{ll,vv}(ii,:)=Coverage(:,vv,ll)';
            TabellaS{ll,vv}(ii,:)=StatDistr(:,vv,ll)';
        end
    end
end
TabelloneRobust=[PARAM 1-TabellaC{2,1} 1-TabellaC{2,2}];
TabelloneNoNRobust=[PARAM 1-TabellaC{4,1} 1-TabellaC{4,2}];

  % q=PARAM(kk,2); r=PARAM(kk,3); N=PARAM(kk,4); T=PARAM(kk,5);             % Retreive parameters DGP
  %   mu=PARAM(kk,6); tau=PARAM(kk,7); delta=PARAM(kk,8); theta=PARAM(kk,9);  % -----------------------
  %   stT=PARAM(kk,10); iota=PARAM(kk,11); pi=PARAM(kk,12);

J(:,1)=find(PARAM(:,4)>75&PARAM(:,6)==.7&PARAM(:,7)==0&PARAM(:,8)==0&...    % Gaussian, \tau=0, \delta=0
    PARAM(:,9)==0.5&PARAM(:,10)==0);
J(:,2)=find(PARAM(:,4)>75&PARAM(:,6)==.7&PARAM(:,7)==.5&PARAM(:,8)==2.5&... % Gaussian, \tau=0.5, \delta=0.5
    PARAM(:,9)==0.5&PARAM(:,10)==0);
J(:,3)=find(PARAM(:,4)>75&PARAM(:,6)==.7&PARAM(:,7)==.5&PARAM(:,8)==2.5&... % Asymmetric Laplace, \tau=0.5, \delta=0.5
    PARAM(:,9)==0.5&PARAM(:,10)==3);
J(:,4)=find(PARAM(:,4)>75&PARAM(:,6)==.7&PARAM(:,7)==.5&PARAM(:,8)==2.5&... % Skew-t, \tau=0.5, \delta=0.5
    PARAM(:,9)==0.5&PARAM(:,10)==4);


TEMP2=[];
for ia=[.1 .05]; TEMP1=[];              % --------------------------------- % for different values of alpha
    iq=find(ML_round(normcdf(ALPHA),3)==ia/2); iq= [iq na-iq+1];            % identify quantiles of interes
    for ll=1:4; TEMP0=[];                   % ----------------------------- % for different DGP        
        if ll==1                                                            % For gaussian case only compute non-robust covariance matrix
            for vv=1:2                      % ------------------------- % For different estimators (EM/PCA)
                AA=(1-TabellaC{4,vv}(J(:,ll),iq(1)));           % Left tail 
                BB=TabellaC{4,vv}(J(:,ll),iq(2));                               % Right tail 
                TEMP0=cat(2,TEMP0,AA+BB);
            end
            TEMP1=cat(1,TEMP1,TEMP0); TEMP0=[]; 
        end
        for vv=1:2                              % ------------------------- % For different estimators (EM/PCA)
            AA=(1-TabellaC{2,vv}(J(:,ll),iq(1)));           % Left tail (robust covariance matrix)
            BB=TabellaC{2,vv}(J(:,ll),iq(2));                               % Right tail (robust covariance matrix)
            TEMP0=cat(2,TEMP0,AA+BB);
        end        
        TEMP1=cat(1,TEMP1,TEMP0);
    end
    TEMP2=cat(2,TEMP2,1-TEMP1);
    % end
end

stop


if stampa==1
    namefile='SimulationResults2021D.xls';
    Top1={'','q','r','N','T','mu','tau','delta','theta','D'};
    
    Top3{1,1}='mean'; Top3{1,2}='std'; Top3{1,3}='skew'; Top3{1,4}='kurt';
    
    AV={'HCC_BL','HCC_BN','star','0'};
    for cc=1:4
        for vv=1:2
            if vv==1
                Foglio=['Chi_EM_' AV{cc}];
            else
                Foglio=['Chi_PCA_' AV{cc}];
            end
            disp(['Creating sheet '             Foglio])
            xlswrite(namefile, Top1,            Foglio, 'A1')
            xlswrite(namefile, PARAM,           Foglio, 'A2')
            xlswrite(namefile, Top2,            Foglio, 'K1')
            xlswrite(namefile, TabellaC{cc,vv},  Foglio, 'K2')
            xlswrite(namefile, Top3,            Foglio, 'Y1')
            xlswrite(namefile, TabellaS{cc,vv},  Foglio, 'Y2')
        end
    end
end



stop



for qr=1:1
    for nt=1:length(NT)
        J(nt,:,qr)=find(PARAM(:,4)==NT(nt,1) & PARAM(:,5)==NT(nt,2));    % identifies all DGP with a given pair of N and T            
        for cc=1:4
            for jj=1:size(J,2)            
                ZZ1(nt,:,jj,qr,cc)=TabellaC{cc,1}(J(nt,jj,qr),:);
                WW1(nt,:,jj,qr,cc)=TabellaS{cc,1}(J(nt,jj,qr),:);
            end
        end        
    end
end