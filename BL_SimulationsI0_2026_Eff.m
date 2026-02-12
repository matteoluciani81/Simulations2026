% BL_SimulationsI0_2026_Eff - Simulations Stationary DFM


clear all; close all; clc
Cartella='/cmc/home/m1mxl04/Documents/';
addpath([Cartella 'ML']); addpath([Cartella 'ML/ML_Package']);

ML_graph_options
filename='BL_SimulationI0_2026_Eff';
FolderE=[pwd '/StatisticsE/'];
salva=1;
stampa=1;
nrep=5000;


%%% ======================================== %%%
%%%  Parameters for Data Generating Process  %%%
%%% ======================================== %%%

%       q       r
%       2       3
QR=[4 4];                               % [q r]
NT=[100 100; 200 200; 300 300; 500 500; 1000 1000]; % [N T]
%           mu      tau     delta   theta   stT 
%           1       2       3       4       5  
ParamDGP=[  0.7     0       2.5   0.5     0;...
            0.7     0.1     2.5   0.5     0;...
            0.7     0.2     2.5   0.5     0;...
            0.7     0.3     2.5   0.5     0;...
            0.7     0.4     2.5   0.5     0;...
            0.7     0.5     2.5   0.5     0;...
            0.7     1.5     2.5   0.5     0;...
            0.7     2.5     2.5   0.5     0;...
            0.7     3.5     2.5   0.5     0;...
            0.7     100     2.5   0.5     0];
            % 0.7     0       2.5     0.5     0;...            
            % 0.5     0.5     2.5     0.5     0;...
            % 0.7     0.5     2.5     1       0;...
            % 0.7     0.5     2.5     0.5     1;...
            % 0.7     0.5     2.5     0.5     2;...
            % 0.7     0.5     2.5     0.5     3;...
            % 0.7     0.5     2.5     0.5     4];
ParamDGP(:,6)=0;   % iota
ParamDGP(:,7)=1;   % pi
% ParamDGP(end-1:end,:)=[];

kk=0;
for qr=1:size(QR,1)
    for nt=1:size(NT,1)
        for pp=1:size(ParamDGP,1)
            kk=kk+1;
            % kk   q   r   N   T   mu   tau   delta   theta   stT   iota   pi
            %  1   2   3   4   5    6     7       8       9    10     11   12
            PARAM(kk,:)=[kk QR(qr,:) NT(nt,:) ParamDGP(pp,:)];
        end
    end
end


MM=zeros(nrep,18);
% PARAM=PARAM(PARAM(:,4)<1000,:);
PARAM=PARAM(ismember(PARAM(:,7),[.3]),:); %PARAM(:,7)=100;
% PARAM=PARAM(PARAM(:,8)==2.5,:)

for kk=1:size(PARAM,1)  % ------------------------------------------------- % For each DGP       
disp(['Processing configuration ' num2str(kk) '/' num2str(size(PARAM,1))])
    inizio=now;
    % --------------------------------------------------------------------- % Extract parameters
    q=PARAM(kk,2); r=PARAM(kk,3); N=PARAM(kk,4); T=PARAM(kk,5);             % Retreive parameters DGP
    mu=PARAM(kk,6); tau=PARAM(kk,7); delta=PARAM(kk,8); theta=PARAM(kk,9);  % -----------------------
    stT=PARAM(kk,10); iota=PARAM(kk,11); pi=PARAM(kk,12);
    
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
        'iota' num2str(10*iota) 'pi' num2str(pi) STUD];                     % ----------------------------
    
    
    % if exist([FolderE filename '_' tipo '.mat'], 'file')
    %     disp(['This DGP already done: ' tipo]); 
    % continue;
    % end
    disp(tipo);                                                             % Display DGP
    disp(datestr(now));    
    tic % ================================================================= %
    disp('Simulating and estimating the model')
    parfor bb=1:nrep    % ------------------------------------------------- % For each repetition
        MM(bb,:)=BL_SimulationsI0_2026_Eff_aux(T,N,r,q,tau,mu,delta,theta,stT,pi);
    end                 % ------------------------------------------------- %
    toc % ================================================================= %

    if salva==1; save([FolderE filename '_' tipo], 'MM'); end
    
    disp(['Simulation started at: ' datestr(inizio)]);
    disp(['Simulation ended at:   ' datestr(now)]); disp(' ')    
    
    Tabella(kk,:)=sum(MM(:,7:10)>0);
    Tabella2(kk,:)=mean(MM(:,3:6));
    Tabella3(kk,:)=sum(MM(:,11:14)>0);
    disp(sum(MM(:,[14 18])>0))                                                 %  Displays Results
    
end




for kk=1:size(PARAM,1)  % ------------------------------------------------- % For each DGP       
    % --------------------------------------------------------------------- % Extract parameters
    q=PARAM(kk,2); r=PARAM(kk,3); N=PARAM(kk,4); T=PARAM(kk,5);             % Retreive parameters DGP
    mu=PARAM(kk,6); tau=PARAM(kk,7); delta=PARAM(kk,8); theta=PARAM(kk,9);  % -----------------------
    stT=PARAM(kk,10); iota=PARAM(kk,11); pi=PARAM(kk,12);
    
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
        'iota' num2str(10*iota) 'pi' num2str(pi) STUD];                     % ----------------------------
    
    
    load([FolderE filename '_' tipo '.mat'])
          
    Tabella(kk,:)=sum(MM(:,7:10)>0);
    Tabella2(kk,:)=mean(MM(:,3:6));
    Tabella3(kk,:)=sum(MM(:,11:14)>0);   
    Tabella4(kk,:)=sum(MM(:,[14 18])>0)/nrep;
    
end

Tabellone=[PARAM Tabella4];

[Tabellone(1:6:end,13) Tabellone(3:6:end,13) Tabellone(2:6:end,13)]