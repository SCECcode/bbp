 clear;
 warning off all
 load GOF.list
 load GOF_PGA.list
 load GOF_PGV.list
 load GOF_PGD.list
 load GOF_PSA.list
 load GOF_SPECFIT.list
 load GOF_ENERGYFIT.list
 load GOF_InElEl.list
 load GOF_SAFIT.list
 load GOF_SPECDUR.list
 load GOF_CROSSCOR.list
 load GOF_DUR.list
 load GOF_FS.list
 load InElS.19.16.NumSta.out
 load InElR.19.16.NumSta.out
 load PGA.out
 load PGV.out
 
 load FileR.out
 load FileS.out
 
 specPer=[.05,.075,.1,.15,.2,.25,.3,.4,.5,.75,1,1.5,2,3,4,5,7.5,10];
 
 load SpecDurSX.out
 SDSX = SpecDurSX;
 load SpecDurSY.out
 SDSY = SpecDurSY;
 load SpecDurRX.out
 SDRX = SpecDurRX;
 load SpecDurRY.out
 SDRY = SpecDurRY;
 
 
 load SpecOutSX.out
 SSX = SpecOutSX;
 load SpecOutSY.out
 SSY = SpecOutSY;
 load SpecOutRX.out
 SRX = SpecOutRX;
 load SpecOutRY.out
 SRY = SpecOutRY;
 
 load SpecOutHS.out
 SHS = SpecOutHS;
 load SpecOutHR.out
 SHR = SpecOutHR;
 
 for j = 1:901
    perList(j) = 0.001*real(j-1)+0.1;
 end
 for j = 902:991
    perList(j) = 0.1*real(j-901)+1.0;
 end
 
 nt =        8192 ;
 dt =  2.4417043E-02 ;
 CutL =  0.1000000     ;
 CutH =   10.00000     ;
 fmax = 1/CutL;
 fmin = 1/CutH;
 numSta =          10 ;
 wind = floor(0.1*  2.4417043E-02 *        8192 +1);
 InElRX = reshape(InElR(:,1), 19, 16, numSta);
 InElRY = reshape(InElR(:,2), 19, 16, numSta);
 InElRZ = reshape(InElR(:,3), 19, 16, numSta);
 InElSX = reshape(InElS(:,1), 19, 16, numSta);
 InElSY = reshape(InElS(:,2), 19, 16, numSta);
 InElSZ = reshape(InElS(:,3), 19, 16, numSta);
 RR = 1:0.5:10;
 PLTQ=questdlg ...
 ('Save plot results?', ...
 'Plot Question', ...
 'Yes','No','Yes');
 INELQ=questdlg ...
 ('Display Inelastic/Elastic curves?', ...
 'Inel/El Question', ...
 'Yes','No','Yes');
 FSQ=questdlg ...
 ('Fourier spectrum display?', ...
 'FS Question', ...
 'Smoothed','Unsmoothed','Smoothed');
 nms={'Output Directory','Comparison Name'};
 defans={'Plots', 'Results'};
 fields = {'dfilename','flname'};
 info = inputdlg(nms, 'input', 1, defans);
 info = cell2struct(info,fields);
 mkdir(info.dfilename);
 dirname = info.dfilename;
 flname = info.flname;
 DFname = sprintf('%s/%s',dirname,flname);
 
 figure(1)
 clf;
  bins = [5 15 25 35 45 55 65 75 85 95];
 hist(GOF(:,2:4), bins)
 title('GOF (XYZ)')
 pause(2)
 if(PLTQ(1,1) == 'Y')
 orient tall
 print('-dpdf',sprintf('%sHisto.pdf',DFname))
 end
 
 if(INELQ(1,1) == 'Y')
 for ii = 1:numSta
 figure(1)
 clf;
 title(['X Comp, Station = ', int2str(ii)])
 figure(2)
 clf;
 title(['Y Comp, Station = ', int2str(ii)])
 
 for jj = 1:16
 figure(1)
    subplot(4,4,jj)
    plot(RR,InElRX(:,jj,ii))
    hold on
    plot(RR,InElSX(:,jj,ii), '-r')
    legend('Set 1','Set 2')
 title(['Period = ',num2str(specPer(jj+2))])
 
 figure(2)
    subplot(4,4,jj)
    plot(RR,InElRY(:,jj,ii))
    hold on
    plot(RR,InElSY(:,jj,ii), '-r')
    legend('Set 1','Set 2')
 title(['Period = ',num2str(specPer(jj+2))])
 end
 if(PLTQ(1,1) == 'Y')
 orient landscape
 print('-dpdf',sprintf('%sINEL%d.pdf',DFname,ii))
 end
 end
 pause(2)
 end
 
 
 GOFC = cellstr([{'GOF'} {'PGA'} {'PGV'} ...
 {'PGD'} {'CumEn'} ...
 {'SA'} {'PSA'} {'Ratio'} ...
 {'SA16'} {'SpD'} ...
 {'CC'} {'DUR'} {'FSp'}]);
 
 ff(1:((nt/2)+1)) = (0:(nt/2))./(nt*dt);
 for j = ((nt/2)+2):nt
     ff(j) = (-(nt+1)+j)/(nt*dt);
 end
 
 BBscat = cell(numSta, 1);
 
 for ii = 1:numSta
    p1 = (ii-1)*nt+1;
    p2 = (ii)*nt;
    p3 = (ii-1)*991+1;
    p4 = (ii)*991;
 
    mSRRef(ii,:) = sqrt(SRX(p3:p4).*SRY(p3:p4));
    mSRSyn(ii,:) = sqrt(SSX(p3:p4).*SSY(p3:p4));
    figure(2)
    clf;
    title(['Station # = ', int2str(ii)])
 
    subplot(3,2,5)
    loglog(perList,SRX(p3:p4))
    hold on
    loglog(perList,SSX(p3:p4), 'r-')
 
 
    legend('Set 1','Set 2',3)
    axis([CutL CutH 0 max(max(SRX(p3:p4)),max(SSX(p3:p4)))])
 xlabel('Period (sec)')
 
    subplot(3,2,6)
    loglog(perList,SRY(p3:p4))
    hold on
    loglog(perList,SSY(p3:p4), 'r-')
 
 
 
    legend('Set 1','Set 2',3)
    axis([CutL CutH 0 max(max(SRY(p3:p4)),max(SSY(p3:p4)))])
 xlabel('Period (sec)')
 
    subplot(3,2,1)
    plot(FileR(p1:p2,1),FileR(p1:p2,2))
    hold on
    plot(FileS(p1:p2,1),FileS(p1:p2,2), 'r-')
    legend('Set 1','Set 2')
    title(['X component, Fit Value = ', int2str(round(GOF(ii,2)))])
 xlabel('Time (s)')
 
    subplot(3,2,2)
    plot(FileR(p1:p2,1),FileR(p1:p2,3))
    hold on
    plot(FileS(p1:p2,1),FileS(p1:p2,3), 'r-')
    legend('Set 1','Set 2')
    title(['Y component, Fit Value = ', int2str(round(GOF(ii,3)))])
 xlabel('Time (s)')
 
 XR = abs(fft(FileR(p1:p2,2)));
 XS = abs(fft(FileS(p1:p2,2)));
 YR = abs(fft(FileR(p1:p2,3)));
 YS = abs(fft(FileS(p1:p2,3)));
 IN1 = XR;
 IN2 = XS;
 IN3 = YR;
 IN4 = YS;
 if(FSQ(1,1) == 'S')
 for fsi = wind+1:(nt/2)-wind
       tmp1(fsi) = 0.0;
       tmp2(fsi) = 0.0;
       tmp3(fsi) = 0.0;
       tmp4(fsi) = 0.0;
       for fsj = -1*wind:wind
            tmp1(fsi) = tmp1(fsi) + (XR(fsi+fsj));
            tmp2(fsi) = tmp2(fsi) + (XS(fsi+fsj));
            tmp3(fsi) = tmp3(fsi) + (YR(fsi+fsj));
            tmp4(fsi) = tmp4(fsi) + (YS(fsi+fsj));
       end
       IN1(fsi) = tmp1(fsi)/(wind*2.0+1.0);
       IN2(fsi) = tmp2(fsi)/(wind*2.0+1.0);
       IN3(fsi) = tmp3(fsi)/(wind*2.0+1.0);
       IN4(fsi) = tmp4(fsi)/(wind*2.0+1.0);
 end
 XR = IN1;
 XS = IN2;
 YR = IN3;
 YS = IN4;
 end
    subplot(3,2,3)
    loglog(ff,XR)
    hold on
    loglog(ff,XS, 'r-')
    legend('Set 1','Set 2','location','southwest')
 axis([fmin fmax 0 max(max(XR),max(XS))])
 xlabel('Freq (Hz)')
 
    subplot(3,2,4)
    loglog(ff,YR)
    hold on
    loglog(ff,YS, 'r-')
    legend('Set 1','Set 2','location','southwest')
 axis([fmin fmax 0 max(max(YR),max(YS))])
 xlabel('Freq (Hz)')
 if(PLTQ(1,1) == 'Y')
 orient landscape
    print('-dpdf',sprintf('%sPlot%d.pdf',DFname,ii))
 end
 
 GOFXlist = [GOF(ii,2) GOF_PGA(ii,2) GOF_PGV(ii,2) ...
 GOF_PGD(ii,2) GOF_ENERGYFIT(ii,2)  ...
 GOF_SPECFIT(ii,2) GOF_PSA(ii,2) GOF_InElEl(ii,2) ...
 GOF_SAFIT(ii,2) GOF_SPECDUR(ii,2) ...
 GOF_CROSSCOR(ii,2) GOF_DUR(ii,2) GOF_FS(ii,2)];
 
 GOFYlist = [GOF(ii,3) GOF_PGA(ii,3) GOF_PGV(ii,3) ...
 GOF_PGD(ii,3) GOF_ENERGYFIT(ii,3)  ...
 GOF_SPECFIT(ii,3) GOF_PSA(ii,3) GOF_InElEl(ii,3) ...
 GOF_SAFIT(ii,3)  GOF_SPECDUR(ii,3) ...
 GOF_CROSSCOR(ii,3) GOF_DUR(ii,3) GOF_FS(ii,3)];
 
 GOFZlist = [GOF(ii,4) GOF_PGA(ii,4) GOF_PGV(ii,4) ...
 GOF_PGD(ii,4) GOF_ENERGYFIT(ii,4)  ...
 GOF_SPECFIT(ii,4) GOF_PSA(ii,4) GOF_InElEl(ii,4) ...
 GOF_SAFIT(ii,4)  GOF_SPECDUR(ii,4) ...
 GOF_CROSSCOR(ii,4) GOF_DUR(ii,4) GOF_FS(ii,4)];
 
 figure(3)
 clf;
 subplot(3,1,1)
  hold on
 for jj = 1:13
  plot(jj, GOFXlist(jj),'o')
  text(jj,15,char(GOFC(jj)),...
  'VerticalAlignment', 'middle',...
  'HorizontalAlignment','center',...
  'FontSize',8)
 end
  ylabel('X GOF values')
  axis([0 14 0 100])
 grid on
 
 subplot(3,1,2)
  hold on
 for jj = 1:13
  plot(jj, GOFYlist(jj),'o')
  text(jj,15,char(GOFC(jj)),...
  'VerticalAlignment', 'middle',...
  'HorizontalAlignment','center',...
  'FontSize',8)
 end
  ylabel('Y GOF values')
  axis([0 14 0 100])
 grid on
 
 subplot(3,1,3)
  hold on
 for jj = 1:13
  plot(jj, GOFZlist(jj),'o')
  text(jj,15,char(GOFC(jj)),...
  'VerticalAlignment', 'middle',...
  'HorizontalAlignment','center',...
  'FontSize',8)
 end
  ylabel('Z GOF values')
  axis([0 14 0 100])
 grid on
 orient landscape
 if(PLTQ(1,1) == 'Y')
    print('-dpdf',sprintf('%sGOF%d.pdf',DFname,ii))
 end
 sinfo.dataFN.pga = PGA(ii,1);
 sinfo.dataFP.pga = PGA(ii,2);
 sinfo.synFN.pga = PGA(ii,3);
 sinfo.synFP.pga = PGA(ii,4);
 sinfo.dataFN.pgv = PGV(ii,1);
 sinfo.dataFP.pgv = PGV(ii,2);
 sinfo.synFN.pgv = PGV(ii,3);
 sinfo.synFP.pgv = PGV(ii,4);
 sinfo.dataFN.T = perList;
 sinfo.dataFN.Sa = SpecOutRX(p3:p4);
 sinfo.dataFP.Sa = SpecOutRY(p3:p4);
 sinfo.synFN.Sa = SpecOutSX(p3:p4);
 sinfo.synFP.Sa = SpecOutSY(p3:p4);
 BBscat{ii} = sinfo;
 end
 save BBscatRun.mat BBscat
 BBanalyze('old')
 if(PLTQ(1,1) == 'Y')
 print('-dpdf',sprintf('%sBIAS.pdf',DFname))
 end
 
 ppRD = [0.0 0.1 0.15 0.2 0.25 0.3 0.4 ...
 0.5 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10 0.0];
 
 
 GOFlist = cellstr([{'GOF'};{'SpecR'}; ...
 {'PGA'};{'PGV'};{'PGD'};{'CC'}; ...
 {'Dur'};{'FFT'}; ...
 {'InEl'};{'Energy'};{'Spec Dur'}]);
 GOFlistX = cellstr([{'GOF X'}; ...
 {'SpecR X'}; {'PGA X'};{'PGV X'}; ...
 {'PGD X'};{'CC X'};{'Dur X'}; ...
 {'FFT X'};{'InEl X'};{'Energy X'}; ...
 {'Spec Dur X'}]);
 GOFlistY = cellstr([{'GOF Y'}; ...
 {'SpecR Y'}; {'PGA Y'};{'PGV Y'}; ...
 {'PGD Y'};{'CC Y'};{'Dur Y'}; ...
 {'FFT Y'};{'InEl Y'};{'Energy Y'}; ...
 {'Spec Dur Y'}]);
 GOFlistZ = cellstr([{'GOF Z'}; ...
 {'SpecR Z'}; {'PGA Z'};{'PGV Z'}; ...
 {'PGD Z'};{'CC Z'};{'Dur Z'}; 
 {'FFT Z'};{'InEl Z'};{'Energy Z'}; ...
 {'Spec Dur Z'}]);
  xlswrite(DFname, GOFlist', 'GOF', 'B1')
  xlswrite(DFname, GOFlistX', 'GOFX', 'B1')
  xlswrite(DFname, GOFlistY', 'GOFY', 'B1')
  xlswrite(DFname, GOFlistZ', 'GOFZ', 'B1')
  xlswrite(DFname, GOF(:,1), 'GOF', 'B2')
  xlswrite(DFname, GOF(:,2), 'GOFX', 'B2')
  xlswrite(DFname, GOF(:,3), 'GOFY', 'B2')
  xlswrite(DFname, GOF(:,4), 'GOFZ', 'B2')
  xlswrite(DFname, GOF_SPECFIT(:,1), 'GOF', 'C2')
  xlswrite(DFname, GOF_SPECFIT(:,2), 'GOFX', 'C2')
  xlswrite(DFname, GOF_SPECFIT(:,3), 'GOFY', 'C2')
  xlswrite(DFname, GOF_SPECFIT(:,4), 'GOFZ', 'C2')
  xlswrite(DFname, GOF_PGA(:,1), 'GOF', 'D2')
  xlswrite(DFname, GOF_PGA(:,2), 'GOFX', 'D2')
  xlswrite(DFname, GOF_PGA(:,3), 'GOFY', 'D2')
  xlswrite(DFname, GOF_PGA(:,4), 'GOFZ', 'D2')
  xlswrite(DFname, GOF_PGV(:,1), 'GOF', 'E2')
  xlswrite(DFname, GOF_PGV(:,2), 'GOFX', 'E2')
  xlswrite(DFname, GOF_PGV(:,3), 'GOFY', 'E2')
  xlswrite(DFname, GOF_PGV(:,4), 'GOFZ', 'E2')
  xlswrite(DFname, GOF_PGD(:,1), 'GOF', 'F2')
  xlswrite(DFname, GOF_PGD(:,2), 'GOFX', 'F2')
  xlswrite(DFname, GOF_PGD(:,3), 'GOFY', 'F2')
  xlswrite(DFname, GOF_PGD(:,4), 'GOFZ', 'F2')
  xlswrite(DFname, GOF_CROSSCOR(:,1), 'GOF', 'G2')
  xlswrite(DFname, GOF_CROSSCOR(:,2), 'GOFX', 'G2')
  xlswrite(DFname, GOF_CROSSCOR(:,3), 'GOFY', 'G2')
  xlswrite(DFname, GOF_CROSSCOR(:,4), 'GOFZ', 'G2')
  xlswrite(DFname, GOF_DUR(:,1), 'GOF', 'H2')
  xlswrite(DFname, GOF_DUR(:,2), 'GOFX', 'H2')
  xlswrite(DFname, GOF_DUR(:,3), 'GOFY', 'H2')
  xlswrite(DFname, GOF_DUR(:,4), 'GOFZ', 'H2')
  xlswrite(DFname, GOF_FS(:,1), 'GOF', 'I2')
  xlswrite(DFname, GOF_FS(:,2), 'GOFX', 'I2')
  xlswrite(DFname, GOF_FS(:,3), 'GOFY', 'I2')
  xlswrite(DFname, GOF_FS(:,4), 'GOFZ', 'I2')
  xlswrite(DFname, GOF_InElEl(:,1), 'GOF', 'J2')
  xlswrite(DFname, GOF_InElEl(:,2), 'GOFX', 'J2')
  xlswrite(DFname, GOF_InElEl(:,3), 'GOFY', 'J2')
  xlswrite(DFname, GOF_InElEl(:,4), 'GOFZ', 'J2')
  xlswrite(DFname, GOF_ENERGYFIT(:,1), 'GOF', 'K2')
  xlswrite(DFname, GOF_ENERGYFIT(:,2), 'GOFX', 'K2')
  xlswrite(DFname, GOF_ENERGYFIT(:,3), 'GOFY', 'K2')
  xlswrite(DFname, GOF_ENERGYFIT(:,4), 'GOFZ', 'K2')
  xlswrite(DFname, GOF_SPECDUR(:,1), 'GOF', 'L2')
  xlswrite(DFname, GOF_SPECDUR(:,2), 'GOFX', 'L2')
  xlswrite(DFname, GOF_SPECDUR(:,3), 'GOFY', 'L2')
  xlswrite(DFname, GOF_SPECDUR(:,4), 'GOFZ', 'L2')
  output.GOF.GOF = GOF;
  output.GOF.PGA = GOF_PGA;
  output.GOF.PGV = GOF_PGV;
  output.GOF.PGD = GOF_PGD;
  output.GOF.PSA = GOF_PSA;
  output.GOF.SPEC = GOF_SPECFIT;
  output.GOF.ENER = GOF_ENERGYFIT;
  output.GOF.InElEl = GOF_InElEl;
  output.GOF.SAFIT = GOF_SAFIT;
  output.GOF.SPDUR = GOF_SPECDUR;
  output.GOF.XC = GOF_CROSSCOR;
  output.GOF.DUR = GOF_DUR;
  output.GOF.FS = GOF_FS;
  output.DATA.PGA = PGA;
  output.DATA.PGV = PGV;
  output.DATA.InElS = InElS;
  output.DATA.InElR = InElR;
  output.DATA.SDSX = SDSX;
  output.DATA.SDSY = SDSY;
  output.DATA.SDRX = SDRX;
  output.DATA.SDRY = SDRY;
  output.DATA.SSX = SSX;
  output.DATA.SSY = SSY;
  output.DATA.SRX = SRX;
  output.DATA.SRY = SRY;
  output.DATA.SHS = SHS;
  output.DATA.SHR = SHR;
  output.DATA.DT = dt;
  output.DATA.NT = nt;
  output.DATA.FMAX = fmax;
  output.DATA.FMIN = fmin;
  output.DATA.MAXPER = CutH;
  output.DATA.MINPER = CutL;
  output.DATA.NUMSTA = numSta;
  output.DATA.PERLIST = perList;
  output.DATA.NGAPER = specPer(3:18);
  output.DATA.RR = RR;
 load GOF_SA16X_DATA.dat
 load GOF_SA16Y_DATA.dat
 load GOF_SA16Z_DATA.dat
 output.DATA.PER16 = reshape(GOF_SA16X_DATA(:,1),16,numSta);
 output.DATA.GOF16X = reshape(GOF_SA16X_DATA(:,2),16,numSta)*100;
 output.DATA.GOF16Y = reshape(GOF_SA16Y_DATA(:,2),16,numSta)*100;
 output.DATA.GOF16Z = reshape(GOF_SA16Z_DATA(:,2),16,numSta)*100;
 output.DATA.REF16X = reshape(GOF_SA16X_DATA(:,3),16,numSta);
 output.DATA.REF16Y = reshape(GOF_SA16Y_DATA(:,3),16,numSta);
 output.DATA.REF16Z = reshape(GOF_SA16Z_DATA(:,3),16,numSta);
 output.DATA.SYN16X = reshape(GOF_SA16X_DATA(:,4),16,numSta);
 output.DATA.SYN16Y = reshape(GOF_SA16Y_DATA(:,4),16,numSta);
 output.DATA.SYN16Z = reshape(GOF_SA16Z_DATA(:,4),16,numSta);
 output.DATA.SEIS_TIMESERIES = FileS(1:nt,1);
 output.DATA.SEIS_REFX = reshape(FileR(:,2),nt,numSta);
 output.DATA.SEIS_REFY = reshape(FileR(:,3),nt,numSta);
 output.DATA.SEIS_REFZ = reshape(FileR(:,4),nt,numSta);
 output.DATA.SEIS_SYNX = reshape(FileS(:,2),nt,numSta);
 output.DATA.SEIS_SYNY = reshape(FileS(:,3),nt,numSta);
 output.DATA.SEIS_SYNZ = reshape(FileS(:,4),nt,numSta);
 xlswrite(DFname, GOF_PSA(:,1), 'GOF', 'M2')
 xlswrite(DFname, GOF_PSA(:,2), 'GOFX', 'M2')
 xlswrite(DFname, GOF_PSA(:,3), 'GOFY', 'M2')
 xlswrite(DFname, GOF_PSA(:,4), 'GOFZ', 'M2')
 xlswrite(DFname, GOF_SAFIT(:,1), 'GOF', 'N2')
 xlswrite(DFname, GOF_SAFIT(:,2), 'GOFX', 'N2')
 xlswrite(DFname, GOF_SAFIT(:,3), 'GOFY', 'N2')
 xlswrite(DFname, GOF_SAFIT(:,4), 'GOFZ', 'N2')
  xlswrite(DFname, {'PSA'}, 'GOF', 'M1')
 xlswrite(DFname, {'PSAX'}, 'GOFX', 'M1')
 xlswrite(DFname, {'PSAY'}, 'GOFY', 'M1')
 xlswrite(DFname, {'PSAZ'}, 'GOFZ', 'M1')
 xlswrite(DFname, {'SA16'}, 'GOF', 'N1')
 xlswrite(DFname, {'SA16X'}, 'GOFX', 'N1')
 xlswrite(DFname, {'SA16Y'}, 'GOFY', 'N1')
 xlswrite(DFname, {'SA16Z'}, 'GOFZ', 'N1')
 output.DATA.GOF16 = ( ...
 output.DATA.GOF16X+output.DATA.GOF16Y ...
 +output.DATA.GOF16Z)./3;
 xlswrite(DFname, output.DATA.GOF16', 'SA_GOF', 'B2')
 xlswrite(DFname, specPer(3:18), 'SA_GOF', 'B1')
 load GOF_InElEl_All.list
 numPerS = length(GOF_InElEl_All(:,1))/numSta;
 output.DATA.InEl_GOFX = reshape(GOF_InElEl_All(:,2),numPerS,numSta);
 xlswrite(DFname, GOF_InElEl_All(1:numPerS,1)', 'GOF_InElX', 'B1')
 xlswrite(DFname, output.DATA.InEl_GOFX', 'GOF_InElX', 'B2')
 output.DATA.InEl_GOFY = reshape(GOF_InElEl_All(:,3),numPerS,numSta);
 xlswrite(DFname, GOF_InElEl_All(1:numPerS,1)', 'GOF_InElY', 'B1')
 xlswrite(DFname, output.DATA.InEl_GOFY', 'GOF_InElY', 'B2')
 output.DATA.InEl_GOFZ = reshape(GOF_InElEl_All(:,4),numPerS,numSta);
 xlswrite(DFname, GOF_InElEl_All(1:numPerS,1)', 'GOF_InElZ', 'B1')
 xlswrite(DFname, output.DATA.InEl_GOFZ', 'GOF_InElZ', 'B2')
 save(sprintf('%s_DATA_STRUCT.mat',DFname),'output')
