% BL_SimulationsI0_2026_ResultsC - Simulations Stationary DFM

%ADJUST THE LABELING OF THE FILES FOR IOTA

clear all; close all; clc
Cartella='/cmc/home/m1mxl04/Documents/';
addpath([Cartella 'ML']); addpath([Cartella 'ML/ML_Package']);

ML_graph_options
NomeFile='BL_SimulationI0_2026_C_';
FolderC=[pwd '/StatisticsC/'];
stampa=1;
M=2;                                                                        % number of estimators



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
ParamDGP(:,7)=1;   % pi


NameG=['Sim2026_iota' num2str(10*ParamDGP(1,6)) 'pi' num2str(ParamDGP(1,7))]; % Name for charts

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


% PARAM(1:2,:)=[];
nsim=size(PARAM,1);
TabellaChi=NaN(nsim,12); TabellaChiTN=NaN(nsim,8);
TabellaTR=NaN(nsim,8); 
TabellaR2f=NaN(nsim,8); TabellaR2l=TabellaR2f;
TabellaIdentifica=NaN(nsim,2);
ii=0; jj=0;
for kk=1:size(PARAM,1)  % ------------------------------------------------- % For each DGP   
    disp(kk)
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
    iota=PARAM(kk,11);
    % --------------------------------------------------------------------- % Create identifier string
    tipo=['q' num2str(q) 'r' num2str(r) 'T' num2str(T) 'N' num2str(N) ...   % Name for simulation process
        'mu' num2str(10*mu) 'tau' num2str(10*tau) ...                       % ----------------------------
        'delta' num2str(10*delta) 'theta' num2str(10*theta) ...             % ----------------------------
        'iota' num2str(10*iota) 'pi' num2str(pi) STUD];                     % ----------------------------
    disp(tipo);                                                             % Display DGP
    
    if ~exist([FolderC NomeFile tipo '.mat'], 'file')
        jj=jj+1; 
        TabellaParNO(jj,:)=PARAM(kk,:); 
        continue;
    end
    
    ii = ii+1;
    load([FolderC NomeFile tipo]);
    TabellaParSI(ii,:)=PARAM(kk,:);                                     % Parameters of the simulations available
    TabellaChi(ii,:)=[RMSE' RTMSE' MAD' MAX' MSE' TMSE'];
    TabellaChiTN(ii,:)=[RMSEtn' RTMSEtn' MSEtn' TMSEtn'];
    TabellaTR(ii,:)=[TR1a TR1b TR2a TR2b];                              % Trace
    TabellaR2f(ii,:)=reshape(R21a,1,2*r);                               % R2 F
    TabellaR2l(ii,:)=reshape(R22a,1,2*r);                               % R2 L
    TabellaIdentifica(ii,:)=[I_f1 I_L1];
end


ndgp=size(ParamDGP,1);

ZZ1=NaN(size(NT,1),ndgp,size(QR,1)); WW1=ZZ1;                               % preallocates
ZZ2=NaN(size(NT,1),ndgp,size(QR,1)); WW2=ZZ2; ZZ3=ZZ2; WW3=WW2;             % ------------
ZZ4=NaN(size(NT,1),ndgp,r,size(QR,1));                                      % ------------
ZZ5=NaN(size(NT,1),ndgp,size(QR,1)); WW5=ZZ5;                               % ------------ 
for qr=1
    for nt=1:length(NT)
        try
            J(nt,:,qr)=find(PARAM(:,2)==QR(qr,1) & PARAM(:,3)==QR(qr,2) ... % identifies all DGP with a given pair of q and r     
                & PARAM(:,4)==NT(nt,1) & PARAM(:,5)==NT(nt,2));             % and N and T

            AA = TabellaChiTN(J(nt,:,qr),2)';                               % MSE EM        
            ZZ0(nt,:,qr)=AA;

            AA = TabellaChi(J(nt,:,qr),2)';                                 % MSE EM
            BB = TabellaChi(J(nt,:,qr),1)';                                 % MSE PCA

%             AA = TabellaChi(J(nt,:,qr),3)';                                     % MSE EM
%             BB = TabellaChi(J(nt,:,qr),4)';                                     % MSE PCA
            ZZ1(nt,:,qr)=AA;
            WW1(nt,:,qr)=AA./BB;

            AA=TabellaTR(J(nt,:,qr),2)';                                         % Trace factors EM
            BB=TabellaTR(J(nt,:,qr),1)';                                         % Trace factors PCA        
            ZZ2(nt,:,qr)=AA;
            WW2(nt,:,qr)=AA./BB;

            AA=TabellaTR(J(nt,:,qr),6)';                                         % Trace loadings EM
            BB=TabellaTR(J(nt,:,qr),5)';                                         % Trace loadings PCA        
            ZZ3(nt,:,qr)=AA;
            WW3(nt,:,qr)=AA./BB;       

            for rr=1:QR(qr,2); ZZ4(nt,:,rr,qr) = TabellaR2f(J(nt,:,qr),r+rr)'; end

            ZZ5(nt,:,qr)= TabellaIdentifica(J(nt,:,qr),1)';  
            WW5(nt,:,qr)= TabellaIdentifica(J(nt,:,qr),2)';
        catch; disp([qr nt])
        end
    end
end


% ------------------------------------------------------------------------- % Construct the legend for the charts
cc=PARAM(J(1,:,1),:);                           % Parameter configurations with q=r and N=T
cc(cc(:,8)>0,8)=cc(cc(:,8)>0,8)-2; 
dd=repmat('    ',ndgp,1);
dd(cc(:,10)==0,:)=repmat('N   ',sum(cc(:,10)==0),1);
dd(cc(:,10)==1,:)=repmat('t(4)',sum(cc(:,10)==1),1);
dd(cc(:,10)==2,:)=repmat('L   ',sum(cc(:,10)==2),1);
dd(cc(:,10)==3,:)=repmat('AL  ',sum(cc(:,10)==3),1);
dd(cc(:,10)==4,:)=repmat('Sk-t',sum(cc(:,10)==4),1);

for ii=1:ndgp
    Legenda{ii}=['$\mu=$' num2str(cc(ii,6),'%.1f')  ...
        ', $\tau=$' num2str(cc(ii,7),'%.1f') ', $\delta=$' num2str(cc(ii,8),'%.1f') ...
        ', $\theta=$' num2str(cc(ii,9),'%.1f') ', $D=$' dd(ii,:)];
end
% ------------------------------------------------------------------------- %

colore={[0 0.45 0.74],[0.64 0.08 0.18],[0.93 0.69 0.13],[0.85 0.33 0.1],... % colors blue, red, yellow orange, 
    [0.47 0.67 0.19],[0.49 0.18 0.56],[0 .81 .82] [.2 .2 .2],[.5 0 .93]};  % green, purple, dark turquoise, dark gray, 

grey=[.8 .8 .8];                                                            % colors for the grid
ntick=length(NT);

% Reorder DGp and eliminate 75-75 (75-75 is the first row)
idgp=[1 4 2 3 5 6 9 10];                                                    % Reorder DGP and eliminate Student-t and Laplace
Legenda=Legenda(idgp);                                                      % ----------------------------------------------
ZZ1=ZZ1(2:ntick,idgp,:);                                                    % Eliminate 75-75 (75-75 is the first row)
WW1=WW1(2:ntick,idgp,:);
ZZ2=ZZ2(2:ntick,idgp,:);
WW2=WW2(2:ntick,idgp,:);
ZZ3=ZZ3(2:ntick,idgp,:);
WW3=WW3(2:ntick,idgp,:);
ZZ4=ZZ4(2:ntick,idgp,:,:);
ZZ5=ZZ5(2:ntick,idgp,:);
WW5=WW5(2:ntick,idgp,:);

for nt=1:ntick          % ------------------------------------------------- % xticks with panel size
    xtl{nt,1}= ['$n=' num2str(NT(nt,1)) '$, $T=' num2str(NT(nt,2)) '$'];     %
end                     % ------------------------------------------------- %

xtl(1)=[]; ntick=ntick-1;

lqr=[num2str(QR(qr,1)) num2str(QR(qr,2))];
    
% common component
for qr=1    
    axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
    if qr==1
        for ii=3:size(ZZ1,2)
            pl1(ii)=plot(ZZ1(:,ii,qr),'-.','color',colore{ii},'linewidth',1.5);
        end
        a= ML_round(ML_min(ZZ1(:,:,qr)),1,1); b=ML_round(ML_min(ZZ1(:,:,qr),2),1,2);
    else
        a=ML_round(ML_min(ZZ1(:,1:2,qr)),1,1);  b=ML_round(ML_min(ZZ1(:,1:2,qr),2),1,2);
    end
    for ii=1:2
        pl1(ii)=plot(ZZ1(:,ii,qr),'color',colore{ii},'linewidth',2);
    end
    hold off; box on;
    axis([0.5 ntick+.5 a b]);        
    ytick=get(gca,'ytick')';
    set(gca,'xtick',1:5,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
    gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
    legend(pl1, Legenda,'interpreter','latex','location','NE')    
    if stampa==1 
         exportgraphics(gcf,[NameG '_chiMSE_' lqr '.pdf'],'ContentType','vector');
    end
    title('Common components - MSEs')
    clear pl1
end


for qr=1    
    axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
    if qr==1
        for ii=3:size(WW1,2)
            pl1(ii)=plot(WW1(:,ii,qr),'-.','color',colore{ii},'linewidth',1.5);
        end
        a=ML_round(ML_min(WW1(:,:,qr)),2,1); b=ML_round(ML_min(WW1(:,:,qr),2),2,2);
    else
        a=ML_round(ML_min(WW1(:,1:2,qr)),2,1); b=ML_round(ML_min(WW1(:,1:2,qr),2),2,2);
    end
    for ii=1:2
        pl1(ii)=plot(WW1(:,ii,qr),'color',colore{ii},'linewidth',2);
    end
    hold off; box on;    
    axis([0.5 ntick+.5 a b]);
    
    ytick=get(gca,'ytick')';
    set(gca,'xtick',1:5,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.3f'))
    gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
%     legend(pl1, Legenda,'interpreter','latex','location','best')         
    if stampa==1
        exportgraphics(gcf,[NameG '_chiMSErel_' lqr '.pdf'],'ContentType','vector');
    end
    title('Common components - MSEs - Relative to PCA')
    clear pl1
end



% % % ------------- Common factors - Trace statistics ------------------------- % 
% % for qr=1:2
% %     lqr=[num2str(QR(qr,1)) num2str(QR(qr,2))];
% %     axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% %     for ii=1:2; pl1(ii)=plot(ZZ2(:,ii,qr),'color', colore{ii},'linewidth',2); end
% %     hold off; box on;
% %     axis([0.5 ntick+.5 ML_round(ML_min(ZZ2(:,1:2,:)),1,1) 1])
% %     ytick=get(gca,'ytick')';
% %     set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
% %     gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% %     legend(pl1, Legenda,'interpreter','latex','location','SE')    
% %     if stampa==1; print('-depsc','-painters','-r600',[NameG '_fTRACE_' lqr '.eps']); end
% %     title('Common factors - Trace statistics')
% %     clear pl1
% % end 

% --------- Common factors - Trace statistics - Relative to PCA ---------- % 
for qr=1
    axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
    plot(0:ntick+1, ones(ntick+2,1),'k','linewidth',2)
    for ii=1:2; pl1(ii)=plot(WW2(:,ii,qr),'color', colore{ii},'linewidth',2); end    
    hold off; box on;
    axis([0.5 ntick+.5 ML_round(ML_min(WW2(:,1:2,:)),2,1) ML_round(ML_min(WW2(:,1:2,:),2),3,2)])
    ytick=get(gca,'ytick')';
    set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.3f'))
    gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
    legend(pl1, Legenda,'interpreter','latex','location','NE')    
    if stampa==1         
        exportgraphics(gcf,[NameG '_fTRACErel_' lqr '.pdf'],'ContentType','vector');
    end
    title('Common factors - Trace statistics - Relative to PCA')    
    clear pl1
end 

% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% plot(0:ntick+1, ones(ntick+2,1),'k','linewidth',2)
% for ii=1:2; pl1(ii)=plot(WW2(:,ii,1),'color', colore{ii},'linewidth',2); end
% for ii=1:2; plot(WW2(:,ii,2),'--','color', colore{ii},'linewidth',2); end
% hold off; box on;
% axis([0.5 ntick+.5 ML_round(ML_min(WW2(:,1:2,:)),2,1) ML_round(ML_min(WW2(:,1:2,:),2),2,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.3f'))
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend(pl1, Legenda,'interpreter','latex','location','NE')
% if stampa==1; print('-depsc','-painters','-r600',[NameG '_fTRACErel_.eps']); end
% title('Common factors - Trace statistics - Relative to PCA')
% clear pl1


% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% plot(0:ntick+1, ones(ntick+2,1),'k','linewidth',2)
% for ii=1:2; pl1(ii)=plot(WW2(:,ii,1),'color', colore{ii},'linewidth',2); end
% hold off; box on;
% axis([0.5 ntick+.5 ML_round(ML_min(WW2(:,1:2,1)),3,1) ML_round(ML_min(WW2(:,1:2,1),2),3,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.4f'))
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend(pl1, Legenda,'interpreter','latex','location','SE')
% if stampa==1; print('-depsc','-painters','-r600',[NameG '_fTRACErel_v2024.eps']); end % version 2024 no q<r
% title('Common factors - Trace statistics - Relative to PCA')
% clear pl1

    
% % 
% % % -------------- common factor R2 ----------------------------------------- % 
% % for qr=1:2
% %     lqr=[num2str(QR(qr,1)) num2str(QR(qr,2))];
% %     axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% %     pl1=plot(ZZ4(:,:,1,qr),'linewidth',1.5);
% %     plot(ZZ4(:,:,2,qr),'--','linewidth',1.5);
% %     plot(ZZ4(:,:,3,qr),':','linewidth',1.5);
% %     plot(ZZ4(:,:,4,qr),'-.','linewidth',1.5);
% %     hold off; box on;
% %     axis([0.5 ntick+.5 ML_round(ML_min(ZZ4),2,1) 1])
% %     ytick=get(gca,'ytick')';
% %     set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
% %     gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% %     legend(pl1, Legenda,'interpreter','latex','location','best')
% %     title('Common factors - R2')
% %     if stampa==1; print('-depsc','-painters','-r600',[NameG '_fR2_' lqr '.eps']); end
% %     clear pl1
% % end 



% % % ------------- Factors loadings - Trace statistics ----------------------- % 
% % for qr=1:2
% %     lqr=[num2str(QR(qr,1)) num2str(QR(qr,2))];
% %     axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% %     pl1=plot(ZZ3(:,:,qr),'linewidth',1.5);
% %     hold off; box on;
% %     axis([0.5 ntick+.5 ML_round(ML_min(ZZ3),2,1) 1])
% %     ytick=get(gca,'ytick')';
% %     set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
% %     gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% %     legend(pl1, Legenda,'interpreter','latex','location','best')
% %     title('Factors loadings - Trace statistics')   
% %     if stampa==1; print('-depsc','-painters','-r600',[NameG '_lambdaTRACE_' lqr '.eps']); end
% %     clear pl1
% % end


% --------- Factor loadings - Trace statistics - Relative to PCA --------- % 
for qr=1    
    axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
    plot(0:ntick+1, ones(ntick+2,1),'k','linewidth',2)
    for ii=1:2; pl1(ii)=plot(WW3(:,ii,qr),'color', colore{ii},'linewidth',2); end   
    hold off; box on;
    axis([0.5 ntick+.5 ML_round(ML_min(WW3(:,1:2,:)),3,1) ML_round(ML_min(WW3(:,1:2,:),2),3,2)])
    ytick=get(gca,'ytick')';
    set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.4f'))
    gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
    legend(pl1, Legenda,'interpreter','latex','location','best')
    if stampa==1
        exportgraphics(gcf,[NameG '_lambdaTRACErel_' lqr '.pdf'],'ContentType','vector');
    end
    title('Factor loadings - Trace statistics - Relative to PCA')        
    clear pl1
end

% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% plot(0:ntick+1, ones(ntick+2,1),'k','linewidth',2)
% for ii=1:2; pl1(ii)=plot(WW3(:,ii,1),'color', colore{ii},'linewidth',2); end
% for ii=1:2; plot(WW3(:,ii,2),'--','color', colore{ii},'linewidth',2); end
% hold off; box on;
% axis([0.5 ntick+.5 ML_round(ML_min(WW3(:,1:2,:)),3,1) ML_round(ML_min(WW3(:,1:2,:),2),3,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.4f'))
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend(pl1, Legenda,'interpreter','latex','location','NE')
% if stampa==1; print('-depsc','-painters','-r600',[NameG '_lambdaTRACErel.eps']); end
% title('Factor loadings - Trace statistics - Relative to PCA')
% clear pl1
% 
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% plot(0:ntick+1, ones(ntick+2,1),'k','linewidth',2)
% for ii=1:2; pl1(ii)=plot(WW3(:,ii,1),'color', colore{ii},'linewidth',2); end
% hold off; box on;
% axis([0.5 ntick+.5 ML_round(ML_min(WW3(:,1:2,1)),3,1) ML_round(ML_min(WW3(:,1:2,1),2),3,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.4f'))
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend(pl1, Legenda,'interpreter','latex','location','NE')
% if stampa==1; print('-depsc','-painters','-r600',[NameG '_lambdaTRACErel_v2024.eps']); end % version 2024 no q<r
% title('Factor loadings - Trace statistics - Relative to PCA')
% clear pl1


% ------------- Identification criteria common factors -------------------- % 
for qr=1    
    axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
    for ii=1:2; pl1(ii)=plot(ZZ5(:,ii,qr),'color', colore{ii},'linewidth',2); end     
    hold off; box on;
    axis([0.5 ntick+.5 0 ML_round(ML_min(ZZ5(:,1:2,qr),2),3,2)])
    ytick=get(gca,'ytick')';
    set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.3f'))
    gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
    legend(pl1, Legenda,'interpreter','latex','location','NE')      
    if stampa==1
        exportgraphics(gcf,[NameG '_ICf_' lqr '.pdf'],'ContentType','vector');
    end
    title('Identification criteria common factors') 
    clear pl1
end

% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% for ii=1:2; pl1(ii)=plot(ZZ5(:,ii,1),'color', colore{ii},'linewidth',2); end
% for ii=1:2; plot(ZZ5(:,ii,2),'--','color', colore{ii},'linewidth',2); end
% hold off; box on;
% axis([0.5 ntick+.5 0 ML_round(ML_min(ZZ5(:,1:2,:),2),2,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.3f'))
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend(pl1, Legenda,'interpreter','latex','location','NE')
% if stampa==1; print('-depsc','-painters','-r600',NameG '_ICf.eps'); end
% title('Identification criteria common factors')
% clear pl1
% 
% 
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% for ii=1:2; pl1(ii)=plot(ZZ5(:,ii,1),'color', colore{ii},'linewidth',2); end
% hold off; box on;
% axis([0.5 ntick+.5 0 ML_round(ML_min(ZZ5(:,1:2,1),2),3,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.3f'))
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend(pl1, Legenda,'interpreter','latex','location','NE')
% if stampa==1; print('-depsc','-painters','-r600',NameG '_ICf_2024.eps'); end % version 2024 no q<r
% title('Identification criteria common factors')
% clear pl1

% ------------- Identification criteria factor loadings-------------------- % 
for qr=1
    axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
    ww5=log(WW5(:,:,qr));
    for ii=1:2; pl1(ii)=plot(ww5(:,ii),'color', colore{ii},'linewidth',2); end 
    hold off; box on;
    axis([0.5 ntick+.5 ML_round(ML_min(ww5),0,1) ML_round(ML_min(ww5,2),0,2)])
    ytick=get(gca,'ytick')';
    set(gca,'xtick',1:ntick,'xticklabel',xtl)
    gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
    legend(pl1, Legenda,'interpreter','latex','location','NE')      
    if stampa==1
        exportgraphics(gcf,[NameG '_ICl_' lqr '.pdf'],'ContentType','vector');
    end
    title('Identification criteria factor loadings')  
end

% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% ww5=log([WW5(:,1:2,1) WW5(:,1:2,2)]);
% for ii=1:2; pl1(ii)=plot(ww5(:,ii),'color', colore{ii},'linewidth',2); end
% for ii=3:4; plot(ww5(:,ii),'--','color', colore{ii-2},'linewidth',2); end
% hold off; box on;
% axis([0.5 ntick+.5 ML_round(ML_min(ww5),0,1) ML_round(ML_min(ww5,2),0,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl)
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend(pl1, Legenda,'interpreter','latex','location','NE')
% if stampa==1; print('-depsc','-painters','-r600',NameG '_ICl.eps'); end
% title('Identification criteria factor loadings')
% 
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% ww5=log(WW5(:,1:2,1));
% for ii=1:2; pl1(ii)=plot(ww5(:,ii),'color', colore{ii},'linewidth',2); end
% hold off; box on;
% axis([0.5 ntick+.5 ML_round(ML_min(ww5),0,1) ML_round(ML_min(ww5,2),0,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl)
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend(pl1, Legenda,'interpreter','latex','location','NE')
% if stampa==1; print('-depsc','-painters','-r600',NameG '_ICl_2024.eps'); end    % version 2024 no q<r
% title('Identification criteria factor loadings')

stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filename='SimulationResults2021C.xls';
Foglio{1,1}='Chi';
Foglio{1,2}='Trace for F & L'; 
Foglio{1,3}='R2 for F & L'; 
Foglio{1,4}='Identification criteria';

Top1={'','q','r','T','N','mu','tau','delta','theta','D'};
Top2={'RMSE','RMSE', 'RTMSE','RTMSE', 'MAD','MAD', 'MAX','MAX'...
    'MSE','MSE', 'TMSE','TMSE',};
Top3={'PCA','MLE', 'PCA','MLE', 'PCA','MLE', 'PCA','MLE', 'PCA','MLE', 'PCA','MLE'};

xlswrite(filename, Top1, Foglio{1}, 'A1')
xlswrite(filename, PARAM, Foglio{1}, 'A3')
xlswrite(filename, Top2, Foglio{1}, 'K1')
xlswrite(filename, Top3, Foglio{1}, 'K2')
xlswrite(filename, TabellaChi, Foglio{1}, 'K3')

Top4={'RMSEtn','RMSEtn', 'RTMSEtn','RTMSEtn', 'MSEtn','MSEtn', 'TMSEtn','TMSEtn',};
xlswrite(filename, Top4, Foglio{1}, 'W1')
xlswrite(filename, {Top3{1:8}}, Foglio{1}, 'W2')
xlswrite(filename, TabellaChiTN, Foglio{1}, 'W3')


Top5={'F','F','F','F','L','L','L','L'};
Top6={'mean', 'mean', 'TM','TM','mean', 'mean', 'TM','TM'};

xlswrite(filename, Top1, Foglio{2}, 'A1')
xlswrite(filename, PARAM, Foglio{2}, 'A4')
xlswrite(filename, Top5, Foglio{2}, 'K1')
xlswrite(filename, Top6, Foglio{2}, 'K2')
xlswrite(filename, {Top3{1:8}}, Foglio{2}, 'K3')
xlswrite(filename, TabellaTR, Foglio{2}, 'K4')


Top7={'F','F','F','F','F','F','F','F','L','L','L','L','L','L','L'};
Top8={'MLE','MLE','MLE','MLE','PCA','PCA','PCA','PCA',...
    'MLE','MLE','MLE','MLE','PCA','PCA','PCA','PCA'};
    

xlswrite(filename, Top1, Foglio{3}, 'A1')
xlswrite(filename, PARAM, Foglio{3}, 'A3')
xlswrite(filename, Top7, Foglio{3}, 'K1')
xlswrite(filename, Top8, Foglio{3}, 'K2')
xlswrite(filename, [TabellaR2f TabellaR2l], Foglio{3}, 'K3')


xlswrite(filename, Top1,    Foglio{4}, 'A1')
xlswrite(filename, PARAM,   Foglio{4}, 'A2')
xlswrite(filename, {'F','L'},    Foglio{4}, 'K1')
xlswrite(filename, TabellaIdentifica, Foglio{4}, 'K2')


stop




% 
% for mu=.5
%     for tau=[0 0.5]; kk=kk+1; ii=0;
%         for delta=[0 2.3 2.5]; ii=ii+1; jj=0;                                  % AR idiosyncratic
%             for theta=[.5 1]; jj=jj+1;                    
%                 for stT=[0]; 
%                     J(:,ii,jj,kk)=find(PARAM(:,6)==mu & PARAM(:,7)==tau & PARAM(:,8)==delta & ...
%                         PARAM(:,9)==theta & PARAM(:,10)==stT);
%                     ZZ1(:,ii,jj,kk)=TabellaChi(J(:,ii,jj,kk),1);
%                     WW1(:,ii,jj,kk)=TabellaChi(J(:,ii,jj,kk),1)./TabellaChi(J(:,ii,jj,kk),2);
%                     
%                     ZZ2(:,ii,jj,kk)=TabellaTR(J(:,ii,jj,kk),1);
%                     WW2(:,ii,jj,kk)=TabellaTR(J(:,ii,jj,kk),1)./TabellaTR(J(:,ii,jj,kk),2);
% 
%                     ZZ3(:,ii,jj,kk)=TabellaTR(J(:,ii,jj,kk),5);
%                     WW3(:,ii,jj,kk)=TabellaTR(J(:,ii,jj,kk),5)./TabellaTR(J(:,ii,jj,kk),6);
%                     Legenda{ii,jj,kk}=['$\mu=$' num2str(mu)  ...
%                         ', $\tau=$' num2str(tau) ', $\theta=$' num2str(theta) ...
%                         ', $\delta=$' num2str(delta2(ii))];
%                 end
%             end
%         end
%     end
% end

% 
% % common component
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% for ii=1:3
%     pl1(ii)=plot(ZZ1(:,ii,1,1),'-','color',colore{ii},'linewidth',1.5);
%     pl2(ii)=plot(ZZ1(:,ii,1,2),'--','color',colore{ii},'linewidth',1.5);
%     pl3(ii)=plot(ZZ1(:,ii,2,1),':','color',colore{ii},'linewidth',1.5);
% end
% hold off; box on;
% axis([0.5 7.5 0 ML_round(ML_min(ZZ1,2),1,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
% xtickangle(15)
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend([pl1(1:3) pl2(1) pl3(1)],[Legenda(1:3,1,1); Legenda{1,1,2} ; Legenda{1,2,1}],...
%     'interpreter','latex','location','NE')
% title('Common components - MSEs')
% print('-depsc','-painters','-r600',NameG '_chiMSE.eps')
% 
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% plot(1:7, ones(7,1),'k','linewidth',1.5)
% for ii=1:3
%     pl1(ii)=plot(WW1(:,ii,1,1),'-','color',colore{ii},'linewidth',1.5);
%     pl2(ii)=plot(WW1(:,ii,1,2),'--','color',colore{ii},'linewidth',1.5);
%     pl3(ii)=plot(WW1(:,ii,2,1),':','color',colore{ii},'linewidth',1.5);
% end
% hold off; box on;
% axis([0.5 7.5 ML_round(ML_min(WW1),1,1) ML_round(ML_min(WW1,2),1,2)])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
% xtickangle(15)
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend([pl1(1:3) pl2(1) pl3(1)],[Legenda(1:3,1,1); Legenda{1,1,2} ; Legenda{1,2,1}],...
%     'interpreter','latex','location','SE')
% title('Common components - MSEs - Relative to PCA')
% print('-depsc','-painters','-r600',NameG '_chiMSErel.eps')
% 
% 
% % common factor
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% for ii=1:3
%     pl1(ii)=plot(ZZ2(:,ii,1,1),'-','color',colore{ii},'linewidth',1.5);
%     pl2(ii)=plot(ZZ2(:,ii,1,2),'--','color',colore{ii},'linewidth',1.5);
%     pl3(ii)=plot(ZZ2(:,ii,2,1),':','color',colore{ii},'linewidth',1.5);
% end
% hold off; box on;
% axis([0.5 7.5 ML_round(ML_min(ZZ2),1,1) 1])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:ntick,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
% xtickangle(15)
% gridxy(1:ntick,ytick,'color',grey,'linewidth',1);     % grid
% legend([pl1(1:3) pl2(1) pl3(1)],[Legenda(1:3,1,1); Legenda{1,1,2} ; Legenda{1,2,1}],...
%     'interpreter','latex','location','SE')
% title('Common factors - Trace statistics')
% print('-depsc','-painters','-r600',NameG '_fTRACE.eps')
% 
% 
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% plot(1:7, ones(7,1),'k','linewidth',1.5)
% for ii=1:3
%     pl1(ii)=plot(WW2(:,ii,1,1),'-','color',colore{ii},'linewidth',1.5);
%     pl2(ii)=plot(WW2(:,ii,1,2),'--','color',colore{ii},'linewidth',1.5);
%     pl3(ii)=plot(WW2(:,ii,2,1),':','color',colore{ii},'linewidth',1.5);
% end
% hold off; box on;
% axis([0.5 7.5 .99 1.05])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:7,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
% xtickangle(15)
% gridxy(1:7,ytick,'color',grey,'linewidth',1);     % grid
% legend([pl1(1:3) pl2(1) pl3(1)],[Legenda(1:3,1,1); Legenda{1,1,2} ; Legenda{1,2,1}],...
%     'interpreter','latex','location','NE')
% title('Common factors - Trace statistics - Relative to PCA')
% print('-depsc','-painters','-r600',NameG '_fTRACErel.eps')
% 
% 
% 
% % factor loadings
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% for ii=1:3
%     pl1(ii)=plot(ZZ3(:,ii,1,1),'-','color',colore{ii},'linewidth',1.5);
%     pl2(ii)=plot(ZZ3(:,ii,1,2),'--','color',colore{ii},'linewidth',1.5);
%     pl3(ii)=plot(ZZ3(:,ii,2,1),':','color',colore{ii},'linewidth',1.5);
% end
% hold off; box on;
% axis([0.5 7.5 ML_round(ML_min(ZZ3),1,1) 1])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:7,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.2f'))
% xtickangle(15)
% gridxy(1:7,ytick,'color',grey,'linewidth',1);     % grid
% legend([pl1(1:3) pl2(1) pl3(1)],[Legenda(1:3,1,1); Legenda{1,1,2} ; Legenda{1,2,1}],...
%     'interpreter','latex','location','SE')
% title('Factors loadings - Trace statistics')
% print('-depsc','-painters','-r600',NameG '_lambdaTRACE.eps')
% 
% 
% axes('Parent',figure,'FontSize',12,'TickLabelInterpreter','latex'); ML_FigureSize,hold on;
% plot(1:7, ones(7,1),'k')
% for ii=1:3
%     pl1(ii)=plot(WW3(:,ii,1,1),'-','color',colore{ii},'linewidth',1.5);
%     pl2(ii)=plot(WW3(:,ii,1,2),'--','color',colore{ii},'linewidth',1.5);
%     pl3(ii)=plot(WW3(:,ii,2,1),':','color',colore{ii},'linewidth',1.5);
% end
% hold off; box on;
% axis([0.5 7.5 .99 1.01])
% ytick=get(gca,'ytick')';
% set(gca,'xtick',1:7,'xticklabel',xtl,'yticklabel',num2str(ytick,'%.3f'))
% xtickangle(15)
% gridxy(1:7,ytick,'color',grey,'linewidth',1);     % grid
% legend([pl1(1:3) pl2(1) pl3(1)],[Legenda(1:3,1,1); Legenda{1,1,2} ; Legenda{1,2,1}],...
%     'interpreter','latex','location','NE')
% title('Factor loadings - Trace statistics - Relative to PCA')
% print('-depsc','-painters','-r600',NameG '_lambdaTRACErel.eps')
