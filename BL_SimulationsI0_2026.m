% BL_SimulationsI0_2026 - Simulations Stationary DFM


clear all; close all; clc
Cartella='/cmc/home/m1mxl04/Documents/';
addpath([Cartella 'ML']); addpath([Cartella 'ML/ML_Package']);

ML_graph_options
filename='BL_SimulationI0_2026';
FolderH=[pwd '/Histograms/'];
FolderC=[pwd '/StatisticsC/'];
FolderD=[pwd '/StatisticsD/'];
salva=1;
stampa=0;
override=1;
nrep=2000;
nv=4;                                                                       % number of variances to consider
omega=1;                                                                    % trimmed portion
M=2;                                                                        % number of estimators
ALPHA=norminv([.75 .84 .9 .95 .975 .99]); ALPHA=sort([0 ALPHA -ALPHA]);     % Percentiles of N(0,1) for coverage
LabelM={'EM','PCA'};

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


kk=0;
for qr=1:1
    for nt=1:size(NT,1)
        for pp=1:size(ParamDGP,1)
            kk=kk+1;
            % kk   q   r   N   T   mu   tau   delta   theta   stT   iota   pi
            %  1   2   3   4   5    6     7       8       9    10     11   12
            PARAM(kk,:)=[kk QR(qr,:) NT(nt,:) ParamDGP(pp,:)];
        end
    end
end

                                                    % Select DGP for which we want to produce the histogram
PrmG0=PARAM(PARAM(:,6)==.7 & PARAM(:,7)==.5 & PARAM(:,8)==2.5 & PARAM(:,9)==.5,:);
% PrmG1=PrmG0(PrmG0(:,10)==0 | PrmG0(:,10)==3 | PrmG0(:,10)==4,:); % Select the distributions
PrmG1=PrmG0(PrmG0(:,10)==0 | PrmG0(:,10)==1,:); % Select the distributions
% PrmG1=PrmG0;
PrmG=PrmG1(PrmG1(:,4)>200 & PrmG1(:,4)<400,:);

PARAM=PrmG;
for kk=1:size(PARAM,1)
    disp(['Processing configuration ' num2str(kk) '/' num2str(size(PARAM,1))])
    inizio=now;
    
    % --------------------------------------------------------------------- % Extract parameters
    q=PARAM(kk,2); r=PARAM(kk,3); N=PARAM(kk,4); T=PARAM(kk,5);             % Retreive parameters DGP
    mu=PARAM(kk,6); tau=PARAM(kk,7); delta=PARAM(kk,8); theta=PARAM(kk,9);  % -----------------------
    stT=PARAM(kk,10); iota=PARAM(kk,11); perturb=PARAM(kk,12);
    
    % --------------------------------------------------------------------- % Determine distribution type
    if stT==1; STUD='T';
    elseif stT==2; STUD='L';
    elseif stT==3; STUD='AL';
    elseif stT==4; STUD='ST';
    else STUD=''; 
    end
    
    % --------------------------------------------------------------------- % Create identifier string
    tipo=['q' num2str(q) 'r' num2str(r) 'T' num2str(T) 'N' num2str(N) ...   % Name for simulation process
        'mu' num2str(10*mu) 'tau' num2str(10*tau) ...                       % ----------------------------
        'delta' num2str(10*delta) 'theta' num2str(10*theta) ...             % ----------------------------
        'iota' num2str(10*iota) 'pi' num2str(perturb) STUD];                     % ----------------------------

    dograph = ismember(PARAM(kk,1),PrmG(:,1));    
    
    % --------------------------------------------------------------------- % Check if already completed
    if exist([FolderC filename '_C_' tipo '.mat'], 'file')
        disp(['Already done: ' tipo]);        
        if dograph==1 && override==1
            disp('Running simulation to produce histogram')
        else
            continue;
        end
    end
    disp(tipo);                                                             % Display DGP
    disp(datestr(now));
        
    % --------------------------------------------------------------------- % Preallocate ONLY what we need
    TR01 = NaN(nrep,2); TR02 = NaN(nrep,2);
    I_L1 = NaN(nrep,1); I_f1 = NaN(nrep,1);
    R11 = NaN(nrep,r,2); R12 = R11; R21 = R11; R22 = R11;
    EE = NaN(T,N,nrep,M); EA = NaN(T,N,nrep,M);
    Z1 = NaN(T,N,nv,nrep);
    if q==r; Z0 = NaN(T,N,nv,nrep); end
    
    tic    
    parfor bb=1:nrep
        % ----------------------------------------------------------------- % Generate data and estimate models
        [EM_bb, PCA_bb, DGP_bb] = BL_SimulationsI0_2026_aux...
            (T,N,r,q,tau,mu,delta,theta,stT,iota,perturb);
                
        F0 = DGP_bb.f; L0 = DGP_bb.lambda;  chi0 = DGP_bb.chi;              % Extract true values
        FF1 = PCA_bb.F;LL1 = PCA_bb.Lambda; chat1 = PCA_bb.chi;             % Extract PCA estimates   
        FF2 = EM_bb.F; LL2 = EM_bb.Lambda;  chat2 = EM_bb.chi;              % Extract EM estimates                
        I_L1(bb,1) = EM_bb.I_L1; I_f1(bb,1) = EM_bb.I_f1;                   % Store identification criteria
        
        % ----------------------------------------------------------------- % Compute consistency statistics
        [TR01(bb,:), TR02(bb,:), EE(:,:,bb,:), EA(:,:,bb,:), ...
            R11(bb,:,:), R12(bb,:,:), R21(bb,:,:), R22(bb,:,:)] = ...
            ML_SimulationStatistics2026(F0, {FF1,FF2}, L0, {LL1,LL2}, ...
            chi0, {chat1,chat2});
        
        % ----------------------------------------------------------------- % Compute distribution statistics
        temp = BL_ComputeZ2026(EM_bb, DGP_bb, 0, nv);
        Z1(:,:,:,bb) = temp;
        
        if q==r
            temp = BL_ComputeZ_PCA2026(PCA_bb, DGP_bb, 0, nv);
            Z0(:,:,:,bb) = temp;
        end
    end
    
    elapsed = toc;
    disp(['Simulations completed in ' num2str(elapsed/60,'%.1f') ' minutes'])
    
    tic
    disp('Computing summary statistics...')
    
    % --------------------------------------------------------------------- % Trace statistics
    TR1a = mean(TR01);
    TR1b = trimmean(TR01,omega);
    TR2a = mean(TR02);
    TR2b = trimmean(TR02,omega);
    
    % --------------------------------------------------------------------- % Canonical correlations
    R11a = squeeze(mean(R11));
    R11b = squeeze(trimmean(R11,omega));
    R12a = squeeze(mean(R12));
    R12b = squeeze(trimmean(R12,omega));
    
    % --------------------------------------------------------------------- % Standard correlations
    R21a = squeeze(mean(R21));
    R21b = squeeze(trimmean(R21,omega));
    R22a = squeeze(mean(R22));
    R22b = squeeze(trimmean(R22,omega));
    
    % --------------------------------------------------------------------- % Mean Squared Errors
    temp1a=mean(EE);                                                        % mean over time
    temp1b=mean(temp1a,2);                                                  % mean over variables    
    RMSE=squeeze(mean(sqrt(temp1b),3));                                     % mean over repetitions (RMSE)    
    MSE=squeeze(mean(temp1b,3));                                            % mean over repetitions (MSE)    
    RMSEtn=RMSE*min(sqrt(T),sqrt(N));
    MSEtn=MSE*min(T,N);
    
    % --------------------------------------------------------------------- % Trimmed Mean Squared Errors
    temp1a=trimmean(EE,omega);                                              % TM over time
    temp1b=trimmean(temp1a,omega,'',2);                                     % TM over variables
    RTMSE=squeeze(trimmean(sqrt(temp1b),omega,'',3));                       % TM over repetitions (MSE)
    TMSE=squeeze(trimmean(temp1b,omega,'',3));                              % TM over repetitions (MSE)
    RTMSEtn=RTMSE*min(sqrt(T),sqrt(N));   
    TMSEtn=TMSE*min(T,N);
        
    % --------------------------------------------------------------------- % Mean Absolute Errors
    temp1a=mean(EA);                                                        % mean over time
    temp1b=mean(temp1a,2);                                                  % mean over variables    
    MAD=squeeze(mean(temp1b,3));                                            % mean over repetitions    
       
    % --------------------------------------------------------------------- % Maximum Absolute Errors
    temp1a=max(EA);                                                         % max over time
    temp1b=max(temp1a,[],2);                                                % max over variables    
    MAX=squeeze(max(temp1b,[],3));                                          % max over repetitions           
    
    % --------------------------------------------------------------------- % Identification criteria
    I_L1 = mean(I_L1);
    I_f1 = mean(I_f1);
    
    % --------------------------------------------------------------------- % Coverage statistics
    ntb=T*N*nrep;
    aa_len = length(ALPHA);
    Coverage = NaN(aa_len, 2, nv);
    StatDistr = NaN(4, 2, nv);
    for cc=1:nv
        % ----------------------------------------------------------------- % EM statistics
        ZZ1 = reshape(Z1(:,:,cc,:), ntb, 1);
        for aa=1:aa_len
            Coverage(aa,1,cc) = mean(ZZ1 > ALPHA(aa));
        end
        StatDistr(:,1,cc) = ML_FourMoments(ZZ1);
        
        if dograph==1
            if ismember(cc,[2 4])
                BL_Histogram(ZZ1, .1, 5);
                print('-dpdf','-vector','-r600', [FolderH 'Hist_C_EM_' tipo '_' num2str(cc)]);
                print('-depsc','-vector','-r600', [FolderH 'Hist_C_EM_' tipo '_' num2str(cc)]);
                close
                ZZ11=reshape(Z1(1,1,cc,:), nrep, 1);
                save([FolderH 'Hist_C_EM_' tipo '_' num2str(cc)],'ZZ11')
            end
        end
        
        % ----------------------------------------------------------------- % PCA statistics (if q==r)
        if q==r
            ZZ0 = reshape(Z0(:,:,cc,:), ntb, 1);
            for aa=1:aa_len
                Coverage(aa,2,cc) = mean(ZZ0 > ALPHA(aa));
            end
            StatDistr(:,2,cc) = ML_FourMoments(ZZ0);
            
            if dograph==1
                if ismember(cc,[2 4])
                    BL_Histogram(ZZ0, .1, 5);
                    print('-depsc','-vector','-r600',[FolderH 'Hist_C_PCA_' tipo '_' num2str(cc)]);
                    print('-dpdf','-vector','-r600',[FolderH 'Hist_C_PCA_' tipo '_' num2str(cc)]);
                    close
                    ZZ01=reshape(Z0(1,1,cc,:), nrep, 1);
                    save([FolderH 'Hist_C_PCA_' tipo '_' num2str(cc)],'ZZ01')
                end
            end         % --------------------------------------------------------- %
        end
    end
    
    toc
    
    %%% ================== %%%
    %%%  Displays Results  %%%
    %%% ================== %%%
    
    temp={'','PCA', 'MLE'};
    temp{2,1}='RMSE';  temp{3,1}='RTMSE';  temp{4,1}='MAD';  temp{5,1}='MAX';    
    for mm=1:M
        temp{2,mm+1}=num2str(RMSE(mm),'%.3f');
        temp{3,mm+1}=num2str(RTMSE(mm),'%.3f');
        temp{4,mm+1}=num2str(MAD(mm),'%.3f');
        temp{5,mm+1}=num2str(MAX(mm),'%.3f');
    end    
    disp(temp)        
            
    temp2 = cell(aa_len, nv+1);
    for aa=1:aa_len
        temp2{aa,1}=[num2str(100*(1-normcdf(ALPHA(aa))),'%.0f') '%'];        
        for cc=1:nv
            temp2{aa,cc+1}=num2str(Coverage(aa,1,cc),'%.3f');
        end           
    end    
    disp(temp2) 
    
    temp3{1,1}='mean'; temp3{2,1}='std'; temp3{3,1}='skew'; temp3{4,1}='kurt';     
    for jj1=1:4
        for cc=1:nv
            temp3{jj1,cc+1}=num2str(StatDistr(jj1,1,cc),'%.2f'); 
        end
    end  
    disp(temp3)    
        
    toc % ================================================================= %

    
    % Save results
    if salva==1
        save([FolderC filename '_C_' tipo], ...
            'TR1a', 'TR1b','TR2a','TR2b','TR2b',...
            'MSE','MSEtn','TMSE','TMSEtn','RMSE','RMSEtn','RTMSE','RTMSEtn',...        
            'MAD','MAX', 'R11a', 'R11b','R12a','R12b',...
            'R21a', 'R21b','R22a','R22b',...
            'I_L1','I_f1');
        save([FolderD filename '_D_' tipo],'Coverage','StatDistr');
    end       
       
    
    clear trPtT trPtt trPttm Z1 Z1b Z1c ZZ ZZ2 
    clear StatDistr StatDistr2
    clear temp temp2 temp3
    
    disp(['Started:  ' datestr(inizio)]);
    disp(['Finished: ' datestr(now)]);
    disp(' ')
end

disp('All simulations completed!')

stop


q=4; r=4;  mu=.5; tau=.5; delta=.2; theta=.5; stT=1; STUD=[]; kk=0;
if stT==1; STUD='T'; else STUD=[]; end
for N=[5 10 25 50 75 100 200 300]; T=100; kk=kk+1;
    tipo=['q' num2str(q) 'r' num2str(r) 'T' num2str(T) 'N' num2str(N) ...   % Name for simulation process
        'mu' num2str(10*mu) 'tau' num2str(10*tau) ...                       % ----------------------------
        'delta' num2str(10*delta) 'theta' num2str(10*theta) STUD];          % ----------------------------    
    disp(tipo);  
    tic
    EM=BL_SimulationsI0c_aux(T,N,r,q,tau,mu,delta,theta,stT);
    toc
    PtT(:,kk)=EM.trPtT(1:20,:);
    Ptt(:,kk)=EM.trPtt(1:20,:);
    Pttm(:,kk)=EM.trPttm(1:20,:);       
end

ZZ=[Pttm(:,kk) Ptt(:,kk) PtT(:,kk)];

plot(ZZ)








 % % % ALPHA2=ALPHA([3 6:8 11]);
% % % for cc=1:5    
% % %     ZZ{1}=reshape(Z1(:,:,cc,:),ntb,1);
% % %     BL_Histogram(ZZ{1},.4,5);
% % %     % print('-dpdf','-painters','-r600',['Hist3.pdf']);
% % %     M4(cc,:)=ML_FourMoments(ZZ{1});
% % %     
% % %     for aa=1:length(ALPHA2)  % ----------------------------------------- % Coverage
% % %         aaa=ALPHA2(aa);
% % %         TCov(cc,aa)=1-sum(ZZ{1}>aaa)/ntb;                           % common
% % %     end
% % % end
% % % 
% % % 
% % % for cc=1:5    
% % %     ZZ{4}=reshape(Z0(:,:,cc,:),ntb,1);
% % %     BL_Histogram(ZZ{4},.4,5);
% % %      M4(cc,:)=ML_FourMoments(ZZ{4});
% % %     
% % %     for aa=1:length(ALPHA2)  % ----------------------------------------- % Coverage
% % %         aaa=ALPHA2(aa);
% % %         TCov(cc,aa)=1-sum(ZZ{4}>aaa)/ntb;                           % common
% % %     end
% % % end
% % % 
% % % for bb=1:nrep;
% % %     ZZ{1}=reshape(Z1(:,:,4,bb),N*T,1);
% % %     M4a(bb,:)=ML_FourMoments(ZZ{1});
% % % end