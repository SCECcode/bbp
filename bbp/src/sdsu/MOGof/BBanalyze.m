function BBanalyze(typ);
%
% script to quantitatively and statistically compare
% the BB-simulated ground motions with respect to the
% data; reads in pre-computed *mat file with lots of
% parameters for all stations, and computes stuff from
% there, both station-wise as well as enssemble averages.
%
% EVASCERATED SCRIPT FOR KIM
%
% this generates black-and-white plots
%
% Martin Mai, JUNE 2005
% ---------------------


%% the large cell array contains structures for each site
if     typ == 'old'; load BBscatRun.mat;
elseif typ == 'ORG'; load BBscatRunORG_fQ1.mat;  BBscat = BBscatORG;
elseif typ == 'WLD'; load BBscatRunWLD_mg05.mat; BBscat = BBscatWLD;
elseif typ == 'MIX'; load BBscatRunMIX_fQ1.mat;  BBscat = BBscatMIX;
elseif typ == 'PSD'; load BBscatRunPSD_fQpar15.mat; BBscat = BBscatPSD;
elseif typ == 'CMB'; load BBscatRunCMB_mg05.mat; BBscat = BBscatCMB;
end  

%%% number of stations
%%%Ns = 30;
Ns = length(BBscat);
%% read PGA / PGV values for all stations, compute geometric mean
%% read also spectral accelerations etc
for kk = 1:Ns;
  
  sinfo = BBscat{kk};
  
  dPGAfn(kk) = sinfo.dataFN.pga;
  dPGAfp(kk) = sinfo.dataFP.pga;
  dPGVfn(kk) = sinfo.dataFN.pgv;
  dPGVfp(kk) = sinfo.dataFP.pgv;
  
  sPGAfn(kk) = sinfo.synFN.pga;
  sPGAfp(kk) = sinfo.synFP.pga;
  sPGVfn(kk) = sinfo.synFN.pgv;
  sPGVfp(kk) = sinfo.synFP.pgv; 
  
  %Repi(kk)   = sinfo.BBsyn.Repi;
  %Rhyp(kk)   = sinfo.BBsyn.Rhyp;

  
  %%% allocate arrays for later use
  if kk == 1; 
    T = sinfo.dataFN.T'; 
    a = zeros(length(T),Ns);
    dfnSA = a; dfpSA = a; dgmSA = a;
    sfnSA = a; sfpSA = a; sgmSA = a;
    rfnSA = a; rfpSA = a; resSA = a;
  end
  
  
  %%% read pre-computed response spectra in allocated arrays
  dfnSA(:,kk) = sinfo.dataFN.Sa(:);
  dfpSA(:,kk) = sinfo.dataFP.Sa(:);
  sfnSA(:,kk) = sinfo.synFN.Sa(:);
  sfpSA(:,kk) = sinfo.synFP.Sa(:);  
  
  
  %%% compute geometric means
  dgmSA(:,kk) = sqrt(dfnSA(:,kk).*dfpSA(:,kk));
  sgmSA(:,kk) = sqrt(sfnSA(:,kk).*sfpSA(:,kk));
  
  
  %%% compute residuals using the natural log of the ratio
  %%% of observed/simulated
  resSA(:,kk) = log(dgmSA(:,kk)./sgmSA(:,kk));
  rfnSA(:,kk) = log(dfnSA(:,kk)./sfnSA(:,kk));
  rfpSA(:,kk) = log(dfpSA(:,kk)./sfpSA(:,kk));

  
end


%%% compute the assemble averages for SA residuals
resMD = mean(resSA')';   %% this is what Rob calls model bias
rfnMD = mean(rfnSA')';
rfpMD = mean(rfpSA')';


%%% here I use Rob's approach for computing the standard error
a1 = resSA.^2; b1 = sum(a1')'; 
a2 = rfnSA.^2; b2 = sum(a2')'; 
a3 = rfpSA.^2; b3 = sum(a3')'; 
iNs = 1/(Ns-1);
resSD = sqrt(iNs*(b1-Ns.*resMD.^2));
rfnSD = sqrt(iNs*(b2-Ns.*rfnMD.^2));
rfpSD = sqrt(iNs*(b3-Ns.*rfpMD.^2));


res90 = 1.699*sqrt(1/Ns)*resSD;   %%% 90% confidence intervals of the mean
rfn90 = 1.699*sqrt(1/Ns)*rfnSD;   %%% these are SMALLER than the SE !!!!
rfp90 = 1.699*sqrt(1/Ns)*rfpSD; 



%%% plot the residual (model bias) and the estimated
%%% standard-deviation and 90% confidence limit as
%%% filled patched areas .....


figure

s1 = subplot(311);
a = resMD+resSD; a2 = resMD+res90;  %% a2 = resMD+2*resSD;
b = resMD-resSD; b2 = resMD-res90;  %% b2 = resMD-2*resSD;


pp=patch([(T)' fliplr((T'))],[a' fliplr(b')],[0.8 0.8 0.8]);
pp2 = get(pp,'parent'); set(pp2,'Xscale','log');
hold on; box on;
patch([(T)' fliplr((T'))],[a2' fliplr(b2')],[0.5 0.5 0.5]);
plot((T),resMD,'k-','LineW',2);
ll=legend('1-\sigma interval','90 % C.I.','median');
set(ll,'FontS',10,'box','off');
plot((T),zeros(size(T)),'k:','LineW',2);
set(s1,'LineW',2,'XTick',[0.01 0.1 0.2 0.3 0.5 1 2 3 5 8],...
       'XTicklabel',{'0.01' '0.1' '0.2' '0.3' '0.5' '1' '2'  '3' '5' '8'},...
       'YTick',[-2 -1.5 -1 -0.5 0 0.5 1 1.5 2],...
       'YTickLabel',{'-2' '-1.5' '-1' '-0.5' '0' '0.5' '1' '1.5' '2'},...
       'FontS',11);
ylabel('Model Bias Log (obs/sim)','FontS',10,'FontW','bo');
text(0.07,1.3,'Geometric Mean','FontS',12,'FontW','bo');
axis([0.05 max(T) -1.75 1.75]);
grid on
hold off



s2 = subplot(312);
a = rfnMD+rfnSD; a2 = rfnMD+rfn90;   %%a2 = rfnMD+2*rfnSD;
b = rfnMD-rfnSD; b2 = rfnMD-rfn90;   %%b2 = rfnMD-2*rfnSD;


pp=patch([(T)' fliplr((T'))],[a' fliplr(b')],[0.8 0.8 0.8]);
pp2 = get(pp,'parent'); set(pp2,'Xscale','log');
hold on; box on;
patch([(T)' fliplr((T'))],[a2' fliplr(b2')],[0.5 0.5 0.5]);
plot((T),rfnMD,'k-','LineW',2);
plot((T),zeros(size(T)),'k:','LineW',2);
set(s2,'LineW',2,'XTick',[0.01 0.1 0.2 0.3 0.5 1 2 3 5 8],...
       'XTicklabel',{'0.01' '0.1' '0.2' '0.3' '0.5' '1' '2'  '3' '5' '8'},...
       'YTick',[-2 -1.5 -1 -0.5 0 0.5 1 1.5 2],...
       'YTickLabel',{'-2' '-1.5' '-1' '-0.5' '0' '0.5' '1' '1.5' '2'},...
       'FontS',11);
ylabel('Model Bias Log (obs/sim)','FontS',10,'FontW','bo');
text(0.07,1.3,'Fault-Normal','FontS',12,'FontW','bo');
axis([0.05 max(T) -1.75 1.75]);
grid on
hold off



s3 = subplot(313);
a = rfpMD+rfpSD; a2 = rfpMD+rfp90;   %% a2 = rfpMD+2*rfpSD;
b = rfpMD-rfpSD; b2 = rfpMD-rfp90;   %% b2 = rfpMD-2*rfpSD;


plot((T),rfpMD,'r-','LineW',2); hold on
pp=patch([(T)' fliplr((T'))],[a' fliplr(b')],[0.8 0.8 0.8]);
pp2 = get(pp,'parent'); set(pp2,'Xscale','log');
hold on; box on;
patch([(T)' fliplr((T'))],[a2' fliplr(b2')],[0.5 0.5 0.5]);
plot((T),rfpMD,'k-','LineW',2);
plot((T),zeros(size(T)),'k:','LineW',2);
set(s3,'LineW',2,'XTick',[0.01 0.1 0.2 0.3 0.5 1 2 3 5 8],...
       'XTicklabel',{'0.01' '0.1' '0.2' '0.3' '0.5' '1' '2'  '3' '5' '8'},...
       'YTick',[-2 -1.5 -1 -0.5 0 0.5 1 1.5 2],...
       'YTickLabel',{'-2' '-1.5' '-1' '-0.5' '0' '0.5' '1' '1.5' '2'},...
       'FontS',11);
xlabel('Period (sec)','FontS',11,'FontW','bo');
ylabel('Model Bias Log (obs/sim)','FontS',10,'FontW','bo');
text(0.07,1.3,'Fault-Parallel','FontS',12,'FontW','bo');
axis([0.05 max(T) -1.75 1.75]);
grid on
hold off

set(s1,'pos',[0.12 0.72 0.78 0.22]);
set(s2,'pos',[0.12 0.44 0.78 0.22]);
set(s3,'pos',[0.12 0.16 0.78 0.22]);
