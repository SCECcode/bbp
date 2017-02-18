C --------------------------GetParam----------------------------------
C Subroutine to read in the input parameter file
C --------------------------------------------------------------------
      subroutine GetParam(inputS, inputR, ntS, ntR, leng, 
     +M, del, lam, Ztor, TAG, numSta, VAL, numHEAD, 
     +CutL, CutH, inputNGA)

      implicit none
      integer                :: nts, ntR, numSta, i, numHEAD
      real                   :: leng, CutL, CutH
      real                   :: M, del, lam, Ztor
      character*80           :: inputR, inputS, inputNGA
      character*80           :: TAG
      real, dimension(20)    :: VAL

      open(1, file='PARAM.dat', status='old')
      VAL(20) = 0.0

C Read in the 1st variable
      read(1,*) inputR
      read(1,*) inputS
      read(1,*) numHEAD
      read(1,*) numSta
      read(1,*) ntR
      read(1,*) ntS         
      read(1,*) TAG     
      read(1,*) leng

      read(1,*) M
      read(1,*) del
      read(1,*) lam
      read(1,*) Ztor
      read(1,*) inputNGA

      read(1,*) CutL
      read(1,*) CutH
      read(1,*) VAL(1)
      read(1,*) VAL(2)
      read(1,*) VAL(3)
      read(1,*) VAL(4)
      read(1,*) VAL(5)
      read(1,*) VAL(6)
      read(1,*) VAL(7)
      read(1,*) VAL(8)
      read(1,*) VAL(9)
      read(1,*) VAL(10)
      read(1,*) VAL(11)
      read(1,*) VAL(12)

      VAL(13) = 0.0
      VAL(14) = 0.0
      VAL(15) = 0.0
      VAL(16) = 0.0
      VAL(17) = 0.0
      VAL(18) = 0.0
      VAL(19) = 0.0

      do i = 1,19
         VAL(20) = VAL(20) + VAL(i)
      enddo

      close(1)
      end subroutine GetParam
C --------------------------------------------------------------------
C --------------------------MFile-------------------------------------
C Subroutine to read in the output M-file
C --------------------------------------------------------------------
      subroutine MFile(numSta, dt2, nt2, CutL, CutH)
      implicit none
      integer    :: nt2, numSta
      real       :: dt2, CutL, CutH
C---------------------------------------------------------------------
      open(9999, file='out/OUTPUT.m', status='unknown')

      write(9999,*) 'clear;'
      write(9999,*) 'warning off all'
      write(9999,*) 'load GOF.list'
      write(9999,*) 'load GOF_PGA.list'
      write(9999,*) 'load GOF_PGV.list'
      write(9999,*) 'load GOF_PGD.list'
      write(9999,*) 'load GOF_PSA.list'
      write(9999,*) 'load GOF_SPECFIT.list'
      write(9999,*) 'load GOF_ENERGYFIT.list'
      write(9999,*) 'load GOF_InElEl.list'
      write(9999,*) 'load GOF_SAFIT.list'
      write(9999,*) 'load GOF_SPECDUR.list'
      write(9999,*) 'load GOF_CROSSCOR.list'
      write(9999,*) 'load GOF_DUR.list'
      write(9999,*) 'load GOF_FS.list'

      write(9999,*) 'load InElS.19.16.NumSta.out'
      write(9999,*) 'load InElR.19.16.NumSta.out'

      write(9999,*) 'load PGA.out'
      write(9999,*) 'load PGV.out'


      write(9999,*) ''
      write(9999,*) 'load FileR.out'
      write(9999,*) 'load FileS.out'
      write(9999,*) ''
      write(9999,*) 'specPer=[.05,.075,.1,.15,.2,.25,.3,.4,.5,.75,
     +1,1.5,2,3,4,5,7.5,10];'
      write(9999,*) ''
      write(9999,*) 'load SpecDurSX.out'
      write(9999,*) 'SDSX = SpecDurSX;'
      write(9999,*) 'load SpecDurSY.out'
      write(9999,*) 'SDSY = SpecDurSY;'
      write(9999,*) 'load SpecDurRX.out'
      write(9999,*) 'SDRX = SpecDurRX;'
      write(9999,*) 'load SpecDurRY.out'
      write(9999,*) 'SDRY = SpecDurRY;'
      write(9999,*) ''
      write(9999,*) 'load NGAplotCB.out'
      write(9999,*) 'NGACB = NGAplotCB;'
      write(9999,*) 'load NGAplotBA.out'
      write(9999,*) 'NGABA = NGAplotBA;'
      write(9999,*) ''
      write(9999,*) 'load SpecOutSX.out'
      write(9999,*) 'SSX = SpecOutSX;'
      write(9999,*) 'load SpecOutSY.out'
      write(9999,*) 'SSY = SpecOutSY;'
      write(9999,*) 'load SpecOutRX.out'
      write(9999,*) 'SRX = SpecOutRX;'
      write(9999,*) 'load SpecOutRY.out'
      write(9999,*) 'SRY = SpecOutRY;'
      write(9999,*) ''
      write(9999,*) 'load SpecOutHS.out'
      write(9999,*) 'SHS = SpecOutHS;'
      write(9999,*) 'load SpecOutHR.out'
      write(9999,*) 'SHR = SpecOutHR;'
      write(9999,*) ''
      write(9999,*) 'for j = 1:901'
      write(9999,*) '   perList(j) = 0.001*real(j-1)+0.1;'
      write(9999,*) 'end'
      write(9999,*) 'for j = 902:991'
      write(9999,*) '   perList(j) = 0.1*real(j-901)+1.0;'
      write(9999,*) 'end'
      write(9999,*) ''
      write(9999,*) 'nt =',nt2,';'
      write(9999,*) 'dt =',dt2,';'
      write(9999,*) 'CutL =',CutL,';'
      write(9999,*) 'CutH =',CutH,';'
      write(9999,*) 'fmax = 1/CutL;'
      write(9999,*) 'fmin = 1/CutH;'
      write(9999,*) 'numSta =',numSta,';'
      write(9999,*) 'wind = floor(0.1*',dt2,'*',nt2,'+1);'

      write(9999,*) 'InElRX = reshape(InElR(:,1), 19, 16, numSta);'
      write(9999,*) 'InElRY = reshape(InElR(:,2), 19, 16, numSta);'
      write(9999,*) 'InElRZ = reshape(InElR(:,3), 19, 16, numSta);'

      write(9999,*) 'InElSX = reshape(InElS(:,1), 19, 16, numSta);'
      write(9999,*) 'InElSY = reshape(InElS(:,2), 19, 16, numSta);'
      write(9999,*) 'InElSZ = reshape(InElS(:,3), 19, 16, numSta);'

      write(9999,*) 'RR = 1:0.5:10;'

      write(9999,*) 'PLTQ=questdlg ...'
      write(9999,*) '(''Save plot results?'', ...'
      write(9999,*) '''Plot Question'', ...'
      write(9999,*) '''Yes'',''No'',''Yes'');'
      
      write(9999,*) 'INELQ=questdlg ...'
      write(9999,*) '(''Display Inelastic/Elastic curves?'', ...'
      write(9999,*) '''Inel/El Question'', ...'
      write(9999,*) '''Yes'',''No'',''Yes'');'

      write(9999,*) 'NGAQ=questdlg ...'
      write(9999,*) '(''Which NGA comparison?'', ...'
      write(9999,*) '''NGA Question'', ...'
      write(9999,*) '''C&B08'',''B&A08'',''C&B08'');'

      write(9999,*) 'FSQ=questdlg ...'
      write(9999,*) '(''Fourier spectrum display?'', ...'
      write(9999,*) '''FS Question'', ...'
      write(9999,*) '''Smoothed'',''Unsmoothed'',''Smoothed'');'


      write(9999,*)'nms={''Output Directory'',''Comparison Name''};'
      write(9999,*)'defans={''Plots'', ''Results''};'
      write(9999,*)'fields = {''dfilename'',''flname''};'
      write(9999,*)'info = inputdlg(nms, ''input'', 1, defans);'
      write(9999,*)'info = cell2struct(info,fields);'
      write(9999,*)'mkdir(info.dfilename);'
      write(9999,*)'dirname = info.dfilename;'
      write(9999,*)'flname = info.flname;'

      write(9999,*) 'DFname = sprintf(''%s/%s'',dirname,flname);'

C      Get directory name
C      Show Inel/El plots
C      Smoothed or unsmoothed
C      C&B vs B&A
      
      write(9999,*) ''
      write(9999,*) 'figure(1)'
      write(9999,*) 'clf;'
      write(9999,*) 'hist(GOF(:,2:4))'
      write(9999,*) 'title(''GOF (XYZ)'')'
      write(9999,*) 'pause(2)'
      write(9999,*) 'if(PLTQ(1,1) == ''Y'')'
      write(9999,*) 'orient tall'
      write(9999,*) 'print(''-dpdf'',sprintf(''%sHisto.pdf'',DFname))'
      write(9999,*) 'end'
      
      write(9999,*) ''
      write(9999,*) 'if(INELQ(1,1) == ''Y'')'
      write(9999,*) 'for ii = 1:numSta'
      write(9999,*) 'figure(1)'
      write(9999,*) 'clf;'
      write(9999,*) 'title([''X Comp, Station = '', int2str(ii)])'
      write(9999,*) 'figure(2)'
      write(9999,*) 'clf;'
      write(9999,*) 'title([''Y Comp, Station = '', int2str(ii)])'
      write(9999,*) ''
      write(9999,*) 'for jj = 1:16'
      write(9999,*) 'figure(1)'
      write(9999,*) '   subplot(4,4,jj)'
      write(9999,*) '   plot(RR,InElRX(:,jj,ii))'
      write(9999,*) '   hold on'
      write(9999,*) '   plot(RR,InElSX(:,jj,ii), ''-r'')'
      write(9999,*) '   legend(''Set 1'',''Set 2'')'
      write(9999,*) 'title([''Period = '',num2str(specPer(jj+2))])'
      write(9999,*) ''
      write(9999,*) 'figure(2)'
      write(9999,*) '   subplot(4,4,jj)'
      write(9999,*) '   plot(RR,InElRY(:,jj,ii))'
      write(9999,*) '   hold on'
      write(9999,*) '   plot(RR,InElSY(:,jj,ii), ''-r'')'
      write(9999,*) '   legend(''Set 1'',''Set 2'')'
      write(9999,*) 'title([''Period = '',num2str(specPer(jj+2))])'
      write(9999,*) 'end'
      write(9999,*) 'if(PLTQ(1,1) == ''Y'')'
      write(9999,*) 'orient landscape'
      write(9999,*) '   print(''-dpdf'',sprintf(''%sINEL%d.pdf''
     +,DFname,ii))'
      write(9999,*) 'end'
      write(9999,*) 'end'
      write(9999,*) 'pause(2)'
      write(9999,*) 'end'
      write(9999,*) ''
      
      write(9999,*) ''
      write(9999,*)'GOFC = cellstr([{''GOF''} {''PGA''} {''PGV''} ...'
      write(9999,*)'{''PGD''} {''CumEn''} ...'
      write(9999,*)'{''SA''} {''PSA''} {''Ratio''} ...'
      write(9999,*)'{''SA16''} {''SpD''} ...'
      write(9999,*)'{''CC''} {''DUR''} {''FSp''}]);'
      write(9999,*) ''
      write(9999,*) 'ff(1:((nt/2)+1)) = (0:(nt/2))./(nt*dt);'
      write(9999,*) 'for j = ((nt/2)+2):nt'
      write(9999,*) '    ff(j) = (-(nt+1)+j)/(nt*dt);'
      write(9999,*) 'end'
      write(9999,*) ''
      write(9999,*) 'BBscat = cell(numSta, 1);'
      write(9999,*) ''
      
      
      
      
      
      
      write(9999,*) 'for ii = 1:numSta'
      write(9999,*) '   p1 = (ii-1)*nt+1;'
      write(9999,*) '   p2 = (ii)*nt;'
      write(9999,*) '   p3 = (ii-1)*991+1;'
      write(9999,*) '   p4 = (ii)*991;'
      write(9999,*) ''
      write(9999,*) '   mSRRef(ii,:) = sqrt(SRX(p3:p4).*SRY(p3:p4));'
      write(9999,*) '   mSRSyn(ii,:) = sqrt(SSX(p3:p4).*SSY(p3:p4));'
      write(9999,*) '   figure(2)'
      write(9999,*) '   clf;'
      write(9999,*) '   title([''Station # = '', int2str(ii)])'
      write(9999,*) ''
      write(9999,*) '   subplot(3,2,5)'
      write(9999,*) '   loglog(perList,SRX(p3:p4))'
      write(9999,*) '   hold on'
      write(9999,*) '   loglog(perList,SSX(p3:p4), ''r-'')'
      write(9999,*) ''
      
      
      
      
      
      write(9999,*) 'if(NGAQ(1,1) == ''C'')'
      write(9999,*) '   loglog(NGACB(((ii-1)*17+2):ii*17,1),
     +NGACB(((ii-1)*17+2):ii*17,2), ''g-'')'
      write(9999,*) 'end'
      
      write(9999,*) 'if(NGAQ(1,1) == ''B'')'
      write(9999,*) '   loglog(NGABA(((ii-1)*17+2):ii*17,1),
     +NGABA(((ii-1)*17+2):ii*17,2), ''y-'')'
      write(9999,*) 'end'
      
      
      write(9999,*) ''
      write(9999,*) 'if(NGAQ(1,1) == ''C'')'
      write(9999,*) '   loglog(NGACB(((ii-1)*17+2):ii*17,1),
     +NGACB(((ii-1)*17+2):ii*17,3), ''g.'')'
      write(9999,*) '   loglog(NGACB(((ii-1)*17+2):ii*17,1),
     +NGACB(((ii-1)*17+2):ii*17,4), ''g.'')'
      write(9999,*) 'end'
      write(9999,*) ''
      write(9999,*) 'if(NGAQ(1,1) == ''B'')'
      write(9999,*) '   loglog(NGABA(((ii-1)*17+2):ii*17,1),
     +NGABA(((ii-1)*17+2):ii*17,3), ''y.'')'
      write(9999,*) '   loglog(NGABA(((ii-1)*17+2):ii*17,1),
     +NGABA(((ii-1)*17+2):ii*17,4), ''y.'')'
      write(9999,*) 'end'
      write(9999,*) '   legend(''Set 1'',''Set 2'',NGAQ,3)'
      write(9999,*) ''
      write(9999,*) '   axis([CutL CutH 0 max(max(SRX(p3:p4)),
     +max(SSX(p3:p4)))])'
      write(9999,*) 'xlabel(''Period (sec)'')'
      write(9999,*) ''
      
      
      
      
      
      write(9999,*) '   subplot(3,2,6)'
      write(9999,*) '   loglog(perList,SRY(p3:p4))'
      write(9999,*) '   hold on'
      write(9999,*) '   loglog(perList,SSY(p3:p4), ''r-'')'
      write(9999,*) ''
      
      
      
      write(9999,*) 'if(NGAQ(1,1) == ''C'')'
      write(9999,*) '   loglog(NGACB(((ii-1)*17+2):ii*17,1),
     +NGACB(((ii-1)*17+2):ii*17,2), ''g-'')'
      write(9999,*) 'end'
     
      write(9999,*) 'if(NGAQ(1,1) == ''B'')'
      write(9999,*) '   loglog(NGABA(((ii-1)*17+2):ii*17,1),
     +NGABA(((ii-1)*17+2):ii*17,2), ''k-'')'
      write(9999,*) 'end'
      write(9999,*) ''
      
      write(9999,*) 'if(NGAQ(1,1) == ''C'')'
      write(9999,*) '   loglog(NGACB(((ii-1)*17+2):ii*17,1),
     +NGACB(((ii-1)*17+2):ii*17,3), ''g.'')'
      write(9999,*) '   loglog(NGACB(((ii-1)*17+2):ii*17,1),
     +NGACB(((ii-1)*17+2):ii*17,4), ''g.'')'
      write(9999,*) 'end'
      write(9999,*) ''
      
      write(9999,*) 'if(NGAQ(1,1) == ''B'')'
      write(9999,*) '   loglog(NGABA(((ii-1)*17+2):ii*17,1),
     +NGABA(((ii-1)*17+2):ii*17,3), ''k.'')'
      write(9999,*) '   loglog(NGABA(((ii-1)*17+2):ii*17,1),
     +NGABA(((ii-1)*17+2):ii*17,4), ''k.'')'
      write(9999,*) 'end'
     

      write(9999,*) ''
      write(9999,*) '   legend(''Set 1'',''Set 2'',NGAQ,3)'
      write(9999,*) '   axis([CutL CutH 0 max(max(SRY(p3:p4)),
     +max(SSY(p3:p4)))])'
      write(9999,*) 'xlabel(''Period (sec)'')'
      write(9999,*) ''
      write(9999,*) '   subplot(3,2,1)'
      write(9999,*) '   plot(FileR(p1:p2,1),FileR(p1:p2,2))'
      write(9999,*) '   hold on'
      write(9999,*) '   plot(FileS(p1:p2,1),FileS(p1:p2,2), ''r-'')'
      write(9999,*) '   legend(''Set 1'',''Set 2'')'
      write(9999,*) '   title([''X component, Fit Value = '', 
     +int2str(round(GOF(ii,2)))])'
      write(9999,*) 'xlabel(''Time (s)'')'
      write(9999,*) ''
      write(9999,*) '   subplot(3,2,2)'
      write(9999,*) '   plot(FileR(p1:p2,1),FileR(p1:p2,3))'
      write(9999,*) '   hold on'
      write(9999,*) '   plot(FileS(p1:p2,1),FileS(p1:p2,3), ''r-'')'
      write(9999,*) '   legend(''Set 1'',''Set 2'')'
      write(9999,*) '   title([''Y component, Fit Value = '', 
     +int2str(round(GOF(ii,3)))])'
      write(9999,*) 'xlabel(''Time (s)'')'
      write(9999,*) ''
      
      
      write(9999,*) 'XR = abs(fft(FileR(p1:p2,2)));'
      write(9999,*) 'XS = abs(fft(FileS(p1:p2,2)));'
      write(9999,*) 'YR = abs(fft(FileR(p1:p2,3)));'
      write(9999,*) 'YS = abs(fft(FileS(p1:p2,3)));'
      write(9999,*) 'IN1 = XR;'
      write(9999,*) 'IN2 = XS;'
      write(9999,*) 'IN3 = YR;'
      write(9999,*) 'IN4 = YS;'
      
      write(9999,*) 'if(FSQ(1,1) == ''S'')'
      write(9999,*) 'for fsi = wind+1:(nt/2)-wind'
      write(9999,*) '      tmp1(fsi) = 0.0;'
      write(9999,*) '      tmp2(fsi) = 0.0;'
      write(9999,*) '      tmp3(fsi) = 0.0;'
      write(9999,*) '      tmp4(fsi) = 0.0;'
      write(9999,*) '      for fsj = -1*wind:wind'
      write(9999,*) '           tmp1(fsi) = tmp1(fsi) + (XR(fsi+fsj));'
      write(9999,*) '           tmp2(fsi) = tmp2(fsi) + (XS(fsi+fsj));'
      write(9999,*) '           tmp3(fsi) = tmp3(fsi) + (YR(fsi+fsj));'
      write(9999,*) '           tmp4(fsi) = tmp4(fsi) + (YS(fsi+fsj));'
      write(9999,*) '      end'
      write(9999,*) '      IN1(fsi) = tmp1(fsi)/(wind*2.0+1.0);'
      write(9999,*) '      IN2(fsi) = tmp2(fsi)/(wind*2.0+1.0);'
      write(9999,*) '      IN3(fsi) = tmp3(fsi)/(wind*2.0+1.0);'
      write(9999,*) '      IN4(fsi) = tmp4(fsi)/(wind*2.0+1.0);'
      write(9999,*) 'end'
      write(9999,*) 'XR = IN1;'
      write(9999,*) 'XS = IN2;'
      write(9999,*) 'YR = IN3;'
      write(9999,*) 'YS = IN4;'
      write(9999,*) 'end'
      
      
      
      write(9999,*) '   subplot(3,2,3)'
      write(9999,*) '   loglog(ff,XR)'
      write(9999,*) '   hold on'
      write(9999,*) '   loglog(ff,XS, ''r-'')'
      write(9999,*) '   legend(''Set 1'',''Set 2'',
     +''location'',''southwest'')'
      write(9999,*) 'axis([fmin fmax 0 max(max(XR),max(XS))])'
      write(9999,*) 'xlabel(''Freq (Hz)'')'
      write(9999,*) ''
      write(9999,*) '   subplot(3,2,4)'
      write(9999,*) '   loglog(ff,YR)'
      write(9999,*) '   hold on'
      write(9999,*) '   loglog(ff,YS, ''r-'')'
      write(9999,*) '   legend(''Set 1'',''Set 2'',
     +''location'',''southwest'')'
      write(9999,*) 'axis([fmin fmax 0 max(max(YR),max(YS))])'
      write(9999,*) 'xlabel(''Freq (Hz)'')'
      
      
      
      write(9999,*) 'if(PLTQ(1,1) == ''Y'')'
      write(9999,*) 'orient landscape'
      write(9999,*) '   print(''-dpdf'',sprintf(''%sPlot%d.pdf''
     +,DFname,ii))'
      write(9999,*) 'end'

      write(9999,*) ''
      write(9999,*) 'GOFXlist = [GOF(ii,2) GOF_PGA(ii,2) 
     +GOF_PGV(ii,2) ...'
      write(9999,*) 'GOF_PGD(ii,2) GOF_ENERGYFIT(ii,2) 
     + ...'
      write(9999,*) 'GOF_SPECFIT(ii,2) GOF_PSA(ii,2) 
     +GOF_InElEl(ii,2) ...'
      write(9999,*) 'GOF_SAFIT(ii,2) 
     +GOF_SPECDUR(ii,2) ...'
      write(9999,*) 'GOF_CROSSCOR(ii,2) GOF_DUR(ii,2) 
     +GOF_FS(ii,2)];'
      write(9999,*) ''
      write(9999,*) 'GOFYlist = [GOF(ii,3) GOF_PGA(ii,3) 
     +GOF_PGV(ii,3) ...'
      write(9999,*) 'GOF_PGD(ii,3) GOF_ENERGYFIT(ii,3) 
     + ...'
      write(9999,*) 'GOF_SPECFIT(ii,3) GOF_PSA(ii,3) 
     +GOF_InElEl(ii,3) ...'
      write(9999,*) 'GOF_SAFIT(ii,3)  
     +GOF_SPECDUR(ii,3) ...'
      write(9999,*) 'GOF_CROSSCOR(ii,3) GOF_DUR(ii,3) 
     +GOF_FS(ii,3)];'
      write(9999,*) ''
      write(9999,*) 'GOFZlist = [GOF(ii,4) GOF_PGA(ii,4) 
     +GOF_PGV(ii,4) ...'
      write(9999,*) 'GOF_PGD(ii,4) GOF_ENERGYFIT(ii,4) 
     + ...'
      write(9999,*) 'GOF_SPECFIT(ii,4) GOF_PSA(ii,4) 
     +GOF_InElEl(ii,4) ...'
      write(9999,*) 'GOF_SAFIT(ii,4)  
     +GOF_SPECDUR(ii,4) ...'
      write(9999,*) 'GOF_CROSSCOR(ii,4) GOF_DUR(ii,4) 
     +GOF_FS(ii,4)];'
      write(9999,*) ''
      write(9999,*) 'figure(3)'
      write(9999,*) 'clf;'
      write(9999,*) 'subplot(3,1,1)'
      write(9999,*) ' hold on'
      write(9999,*) 'for jj = 1:13'
      write(9999,*) ' plot(jj, GOFXlist(jj),''o'')'
      write(9999,*) ' text(jj,15,char(GOFC(jj)),...'
      write(9999,*) ' ''VerticalAlignment'', ''middle'',...'
      write(9999,*) ' ''HorizontalAlignment'',''center'',...'
      write(9999,*) ' ''FontSize'',8)'
      write(9999,*) 'end'
      write(9999,*) ' ylabel(''X GOF values'')'
      write(9999,*) ' axis([0 14 0 100])'
      write(9999,*) 'grid on'
      write(9999,*) ''
      write(9999,*) 'subplot(3,1,2)'
      write(9999,*) ' hold on'
      write(9999,*) 'for jj = 1:13'
      write(9999,*) ' plot(jj, GOFYlist(jj),''o'')'
      write(9999,*) ' text(jj,15,char(GOFC(jj)),...'
      write(9999,*) ' ''VerticalAlignment'', ''middle'',...'
      write(9999,*) ' ''HorizontalAlignment'',''center'',...'
      write(9999,*) ' ''FontSize'',8)'
      write(9999,*) 'end'
      write(9999,*) ' ylabel(''Y GOF values'')'
      write(9999,*) ' axis([0 14 0 100])'
      write(9999,*) 'grid on'
      write(9999,*) ''
      write(9999,*) 'subplot(3,1,3)'
      write(9999,*) ' hold on'
      write(9999,*) 'for jj = 1:13'
      write(9999,*) ' plot(jj, GOFZlist(jj),''o'')'
      write(9999,*) ' text(jj,15,char(GOFC(jj)),...'
      write(9999,*) ' ''VerticalAlignment'', ''middle'',...'
      write(9999,*) ' ''HorizontalAlignment'',''center'',...'
      write(9999,*) ' ''FontSize'',8)'
      write(9999,*) 'end'
      write(9999,*) ' ylabel(''Z GOF values'')'
      write(9999,*) ' axis([0 14 0 100])'
      write(9999,*) 'grid on'
      write(9999,*) 'if(PLTQ(1,1) == ''Y'')'
      write(9999,*) 'orient landscape'
      write(9999,*) '   print(''-dpdf'',sprintf(''%sGOF%d.pdf''
     +,DFname,ii))'
      write(9999,*) 'end'
      write(9999,*) 'sinfo.dataFN.pga = PGA(ii,1);'
      write(9999,*) 'sinfo.dataFP.pga = PGA(ii,2);'
      write(9999,*) 'sinfo.synFN.pga = PGA(ii,3);'
      write(9999,*) 'sinfo.synFP.pga = PGA(ii,4);'

      write(9999,*) 'sinfo.dataFN.pgv = PGV(ii,1);'
      write(9999,*) 'sinfo.dataFP.pgv = PGV(ii,2);'
      write(9999,*) 'sinfo.synFN.pgv = PGV(ii,3);'
      write(9999,*) 'sinfo.synFP.pgv = PGV(ii,4);'

      write(9999,*) 'sinfo.dataFN.T = perList;'

      write(9999,*) 'sinfo.dataFN.Sa = SpecOutRX(p3:p4);'
      write(9999,*) 'sinfo.dataFP.Sa = SpecOutRY(p3:p4);'
      write(9999,*) 'sinfo.synFN.Sa = SpecOutSX(p3:p4);'
      write(9999,*) 'sinfo.synFP.Sa = SpecOutSY(p3:p4);'
      write(9999,*) 'BBscat{ii} = sinfo;'
      write(9999,*) 'end'
      write(9999,*) 'save BBscatRun.mat BBscat'
      write(9999,*) 'BBanalyze(''old'')'
      write(9999,*) 'if(PLTQ(1,1) == ''Y'')'
      write(9999,*) 'print(''-dpdf'',sprintf(''%sBIAS.pdf'',DFname))'
      write(9999,*) 'end'




 
      write(9999,*) ' load NGA_GMRotI50.out'
 
      write(9999,*) ' rrr = reshape(NGA_GMRotI50(:,1),17,numSta);'
      write(9999,*) ' sss = reshape(NGA_GMRotI50(:,2),17,numSta);'
 
      write(9999,*) ' medCB = reshape(NGAplotCB(:,2),17,numSta);'
      write(9999,*) ' stdhCB = reshape(NGAplotCB(:,3),17,numSta);'
      write(9999,*) ' stdlCB = reshape(NGAplotCB(:,4),17,numSta);'
 
      write(9999,*) ' medBA = reshape(NGAplotBA(:,2),17,numSta);'
      write(9999,*) ' stdhBA = reshape(NGAplotBA(:,3),17,numSta);'
      write(9999,*) ' stdlBA = reshape(NGAplotBA(:,4),17,numSta);'
 
      write(9999,*) 'ppRD = [0.0 0.1 0.15 0.2 0.25 0.3 0.4 ...'
      write(9999,*) '0.5 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10];'
 
      write(9999,*) ' xlswrite(DFname, ppRD, ''NGA_GMRotI50_1''
     + ,''B1'')'
      write(9999,*) ' xlswrite(DFname, {''PGA''}, 
     +''NGA_GMRotI50_1'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, rrr'', 
     +''NGA_GMRotI50_1'', ''B2'')'
      write(9999,*) ' xlswrite(DFname, ppRD, 
     +''NGA_GMRotI50_2'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, {''PGA''}, 
     +''NGA_GMRotI50_2'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, sss'', 
     +''NGA_GMRotI50_2'', ''B2'')'
 
      write(9999,*) ' xlswrite(DFname, ppRD, 
     +''NGA_Med_CB'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, {''PGA''}, 
     +''NGA_Med_CB'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, ppRD, 
     +''NGA_StDH_CB'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, {''PGA''}, 
     +''NGA_StDH_CB'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, ppRD, 
     +''NGA_StDL_CB'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, {''PGA''}, 
     +''NGA_StDL_CB'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, medCB'', 
     +''NGA_Med_CB'', ''B2'')'
      write(9999,*) ' xlswrite(DFname, stdhCB'', 
     +''NGA_StDH_CB'', ''B2'')'
      write(9999,*) ' xlswrite(DFname, stdlCB'', 
     +''NGA_StDL_CB'', ''B2'')'
 
      write(9999,*) ' xlswrite(DFname, ppRD, 
     +''NGA_Med_BA'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, {''PGA''}, 
     +''NGA_Med_BA'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, ppRD, 
     +''NGA_StDH_BA'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, {''PGA''}, 
     +''NGA_StDH_BA'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, ppRD, 
     +''NGA_StDL_BA'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, {''PGA''}, 
     +''NGA_StDL_BA'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, medBA'', 
     +''NGA_Med_BA'', ''B2'')'
      write(9999,*) ' xlswrite(DFname, stdhBA'', 
     +''NGA_StDH_BA'', ''B2'')'
      write(9999,*) ' xlswrite(DFname, stdlBA'', 
     +''NGA_StDL_BA'', ''B2'')'


      write(9999,*) 'GOFlist = cellstr([{''GOF''};{''SpecR''}; ...'
      write(9999,*) '{''PGA''};{''PGV''};{''PGD''};{''CC''}; ...'
      write(9999,*) '{''Dur''};{''FFT''}; ...'
      write(9999,*) '{''InEl''};{''Energy''};{''Spec Dur''}]);'
      
      write(9999,*) 'GOFlistX = cellstr([{''GOF X''}; ...'
      write(9999,*) '{''SpecR X''}; {''PGA X''};{''PGV X''}; ...'
      write(9999,*) '{''PGD X''};{''CC X''};{''Dur X''}; ...'
      write(9999,*) '{''FFT X''};{''InEl X''};{''Energy X''}; ...'
      write(9999,*) '{''Spec Dur X''}]);'
      write(9999,*) 'GOFlistY = cellstr([{''GOF Y''}; ...'
      write(9999,*) '{''SpecR Y''}; {''PGA Y''};{''PGV Y''}; ...'
      write(9999,*) '{''PGD Y''};{''CC Y''};{''Dur Y''}; ...'
      write(9999,*) '{''FFT Y''};{''InEl Y''};{''Energy Y''}; ...'
      write(9999,*) '{''Spec Dur Y''}]);'
      write(9999,*) 'GOFlistZ = cellstr([{''GOF Z''}; ...'
      write(9999,*) '{''SpecR Z''}; {''PGA Z''};{''PGV Z''}; ...'
      write(9999,*) '{''PGD Z''};{''CC Z''};{''Dur Z''}; '
      write(9999,*) '{''FFT Z''};{''InEl Z''};{''Energy Z''}; ...'
      write(9999,*) '{''Spec Dur Z''}]);'
 
      write(9999,*) ' xlswrite(DFname, GOFlist'',
     + ''GOF'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, GOFlistX'', 
     +''GOFX'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, GOFlistY'', 
     +''GOFY'', ''B1'')'
      write(9999,*) ' xlswrite(DFname, GOFlistZ'', 
     +''GOFZ'', ''B1'')'
 
      write(9999,*) ' xlswrite(DFname, GOF(:,1), 
     +''GOF'', ''B2'')'
      write(9999,*) ' xlswrite(DFname, GOF(:,2), 
     +''GOFX'', ''B2'')'
      write(9999,*) ' xlswrite(DFname, GOF(:,3), 
     +''GOFY'', ''B2'')'
      write(9999,*) ' xlswrite(DFname, GOF(:,4), 
     +''GOFZ'', ''B2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_SPECFIT(:,1), 
     +''GOF'', ''C2'')'
      write(9999,*) ' xlswrite(DFname, GOF_SPECFIT(:,2), 
     +''GOFX'', ''C2'')'
      write(9999,*) ' xlswrite(DFname, GOF_SPECFIT(:,3), 
     +''GOFY'', ''C2'')'
      write(9999,*) ' xlswrite(DFname, GOF_SPECFIT(:,4), 
     +''GOFZ'', ''C2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_PGA(:,1), 
     +''GOF'', ''D2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGA(:,2), 
     +''GOFX'', ''D2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGA(:,3), 
     +''GOFY'', ''D2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGA(:,4), 
     +''GOFZ'', ''D2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_PGV(:,1), 
     +''GOF'', ''E2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGV(:,2), 
     +''GOFX'', ''E2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGV(:,3), 
     +''GOFY'', ''E2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGV(:,4), 
     +''GOFZ'', ''E2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_PGD(:,1), 
     +''GOF'', ''F2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGD(:,2), 
     +''GOFX'', ''F2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGD(:,3), 
     +''GOFY'', ''F2'')'
      write(9999,*) ' xlswrite(DFname, GOF_PGD(:,4), 
     +''GOFZ'', ''F2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_CROSSCOR(:,1), 
     +''GOF'', ''G2'')'
      write(9999,*) ' xlswrite(DFname, GOF_CROSSCOR(:,2), 
     +''GOFX'', ''G2'')'
      write(9999,*) ' xlswrite(DFname, GOF_CROSSCOR(:,3), 
     +''GOFY'', ''G2'')'
      write(9999,*) ' xlswrite(DFname, GOF_CROSSCOR(:,4), 
     +''GOFZ'', ''G2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_DUR(:,1), 
     +''GOF'', ''H2'')'
      write(9999,*) ' xlswrite(DFname, GOF_DUR(:,2), 
     +''GOFX'', ''H2'')'
      write(9999,*) ' xlswrite(DFname, GOF_DUR(:,3), 
     +''GOFY'', ''H2'')'
      write(9999,*) ' xlswrite(DFname, GOF_DUR(:,4), 
     +''GOFZ'', ''H2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_FS(:,1), 
     +''GOF'', ''I2'')'
      write(9999,*) ' xlswrite(DFname, GOF_FS(:,2), 
     +''GOFX'', ''I2'')'
      write(9999,*) ' xlswrite(DFname, GOF_FS(:,3), 
     +''GOFY'', ''I2'')'
      write(9999,*) ' xlswrite(DFname, GOF_FS(:,4), 
     +''GOFZ'', ''I2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_InElEl(:,1), 
     +''GOF'', ''J2'')'
      write(9999,*) ' xlswrite(DFname, GOF_InElEl(:,2), 
     +''GOFX'', ''J2'')'
      write(9999,*) ' xlswrite(DFname, GOF_InElEl(:,3), 
     +''GOFY'', ''J2'')'
      write(9999,*) ' xlswrite(DFname, GOF_InElEl(:,4), 
     +''GOFZ'', ''J2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_ENERGYFIT(:,1), 
     +''GOF'', ''K2'')'
      write(9999,*) ' xlswrite(DFname, GOF_ENERGYFIT(:,2), 
     +''GOFX'', ''K2'')'
      write(9999,*) ' xlswrite(DFname, GOF_ENERGYFIT(:,3), 
     +''GOFY'', ''K2'')'
      write(9999,*) ' xlswrite(DFname, GOF_ENERGYFIT(:,4), 
     +''GOFZ'', ''K2'')'
 
      write(9999,*) ' xlswrite(DFname, GOF_SPECDUR(:,1), 
     +''GOF'', ''L2'')'
      write(9999,*) ' xlswrite(DFname, GOF_SPECDUR(:,2), 
     +''GOFX'', ''L2'')'
      write(9999,*) ' xlswrite(DFname, GOF_SPECDUR(:,3), 
     +''GOFY'', ''L2'')'
      write(9999,*) ' xlswrite(DFname, GOF_SPECDUR(:,4), 
     +''GOFZ'', ''L2'')'
 
      write(9999,*) ' output.GOF.GOF = GOF;'
      write(9999,*) ' output.GOF.PGA = GOF_PGA;'
      write(9999,*) ' output.GOF.PGV = GOF_PGV;'
      write(9999,*) ' output.GOF.PGD = GOF_PGD;'
      write(9999,*) ' output.GOF.PSA = GOF_PSA;'
      write(9999,*) ' output.GOF.SPEC = GOF_SPECFIT;'
      write(9999,*) ' output.GOF.ENER = GOF_ENERGYFIT;'
      write(9999,*) ' output.GOF.InElEl = GOF_InElEl;'
      write(9999,*) ' output.GOF.SAFIT = GOF_SAFIT;'
      write(9999,*) ' output.GOF.SPDUR = GOF_SPECDUR;'
      write(9999,*) ' output.GOF.XC = GOF_CROSSCOR;'
      write(9999,*) ' output.GOF.DUR = GOF_DUR;'
      write(9999,*) ' output.GOF.FS = GOF_FS;'
      write(9999,*) ' output.DATA.PGA = PGA;'
      write(9999,*) ' output.DATA.PGV = PGV;'
      write(9999,*) ' output.DATA.InElS = InElS;'
      write(9999,*) ' output.DATA.InElR = InElR;'
 
      write(9999,*) ' output.DATA.SDSX = SDSX;'
      write(9999,*) ' output.DATA.SDSY = SDSY;'
      write(9999,*) ' output.DATA.SDRX = SDRX;'
      write(9999,*) ' output.DATA.SDRY = SDRY;'
      write(9999,*) ' output.DATA.SSX = SSX;'
      write(9999,*) ' output.DATA.SSY = SSY;'
      write(9999,*) ' output.DATA.SRX = SRX;'
      write(9999,*) ' output.DATA.SRY = SRY;'
      write(9999,*) ' output.DATA.SHS = SHS;'
      write(9999,*) ' output.DATA.SHR = SHR;'
 
      write(9999,*) ' output.DATA.DT = dt;'
      write(9999,*) ' output.DATA.NT = nt;'
      
      write(9999,*) ' output.DATA.FMAX = fmax;'
      write(9999,*) ' output.DATA.FMIN = fmin;'
      
      write(9999,*) ' output.DATA.MAXPER = CutH;'
      write(9999,*) ' output.DATA.MINPER = CutL;'
      
      write(9999,*) ' output.DATA.NUMSTA = numSta;'
      write(9999,*) ' output.DATA.PERLIST = perList;'
      write(9999,*) ' output.DATA.NGAPER = specPer(3:18);'
      write(9999,*) ' output.DATA.RR = RR;'

      write(9999,*) ' output.DATA.NGA_GMRotD50R = rrr;'
      write(9999,*) ' output.DATA.NGA_GMRotD50S = sss;'
      write(9999,*) ' output.DATA.NGAplotBA = NGAplotBA;'
      write(9999,*) ' output.DATA.NGAplotCB = NGAplotCB;'

      write(9999,*) 'load GOF_SA16X_DATA.dat'
      write(9999,*) 'load GOF_SA16Y_DATA.dat'
      write(9999,*) 'load GOF_SA16Z_DATA.dat'

C  '%Period - GOF - Set 1 - Set 2'

      write(9999,*) 'output.DATA.PER16 = 
     +reshape(GOF_SA16X_DATA(:,1),16,numSta);'

      write(9999,*) 'output.DATA.GOF16X = 
     +reshape(GOF_SA16X_DATA(:,2),16,numSta)*100;'
      write(9999,*) 'output.DATA.GOF16Y = 
     +reshape(GOF_SA16Y_DATA(:,2),16,numSta)*100;'
      write(9999,*) 'output.DATA.GOF16Z = 
     +reshape(GOF_SA16Z_DATA(:,2),16,numSta)*100;'

      write(9999,*) 'output.DATA.REF16X = 
     +reshape(GOF_SA16X_DATA(:,3),16,numSta);'
      write(9999,*) 'output.DATA.REF16Y = 
     +reshape(GOF_SA16Y_DATA(:,3),16,numSta);'
      write(9999,*) 'output.DATA.REF16Z = 
     +reshape(GOF_SA16Z_DATA(:,3),16,numSta);'

      write(9999,*) 'output.DATA.SYN16X = 
     +reshape(GOF_SA16X_DATA(:,4),16,numSta);'
      write(9999,*) 'output.DATA.SYN16Y = 
     +reshape(GOF_SA16Y_DATA(:,4),16,numSta);'
      write(9999,*) 'output.DATA.SYN16Z = 
     +reshape(GOF_SA16Z_DATA(:,4),16,numSta);'
     
      write(9999,*) 'output.DATA.SEIS_TIMESERIES = 
     +FileS(1:nt,1);'
     
      write(9999,*) 'output.DATA.SEIS_REFX = 
     +reshape(FileR(:,2),nt,numSta);'
      write(9999,*) 'output.DATA.SEIS_REFY = 
     +reshape(FileR(:,3),nt,numSta);'
      write(9999,*) 'output.DATA.SEIS_REFZ = 
     +reshape(FileR(:,4),nt,numSta);'

      write(9999,*) 'output.DATA.SEIS_SYNX = 
     +reshape(FileS(:,2),nt,numSta);'
      write(9999,*) 'output.DATA.SEIS_SYNY = 
     +reshape(FileS(:,3),nt,numSta);'
      write(9999,*) 'output.DATA.SEIS_SYNZ = 
     +reshape(FileS(:,4),nt,numSta);'
 


      write(9999,*) 'xlswrite(DFname, GOF_PSA(:,1), ''GOF'', ''M2'')'
      write(9999,*) 'xlswrite(DFname, GOF_PSA(:,2), ''GOFX'', ''M2'')'
      write(9999,*) 'xlswrite(DFname, GOF_PSA(:,3), ''GOFY'', ''M2'')'
      write(9999,*) 'xlswrite(DFname, GOF_PSA(:,4), ''GOFZ'', ''M2'')'
  
      write(9999,*) 'xlswrite(DFname, GOF_SAFIT(:,1), ''GOF'', ''N2'')'
      write(9999,*) 'xlswrite(DFname, GOF_SAFIT(:,2), ''GOFX'', ''N2'')'
      write(9999,*) 'xlswrite(DFname, GOF_SAFIT(:,3), ''GOFY'', ''N2'')'
      write(9999,*) 'xlswrite(DFname, GOF_SAFIT(:,4), ''GOFZ'', ''N2'')'

      write(9999,*) ' xlswrite(DFname, {''PSA''}, ''GOF'', ''M1'')'
      write(9999,*) 'xlswrite(DFname, {''PSAX''}, ''GOFX'', ''M1'')'
      write(9999,*) 'xlswrite(DFname, {''PSAY''}, ''GOFY'', ''M1'')'
      write(9999,*) 'xlswrite(DFname, {''PSAZ''}, ''GOFZ'', ''M1'')'


      write(9999,*) 'xlswrite(DFname, {''SA16''}, ''GOF'', ''N1'')'
      write(9999,*) 'xlswrite(DFname, {''SA16X''}, ''GOFX'', ''N1'')'
      write(9999,*) 'xlswrite(DFname, {''SA16Y''}, ''GOFY'', ''N1'')'
      write(9999,*) 'xlswrite(DFname, {''SA16Z''}, ''GOFZ'', ''N1'')'
      write(9999,*) 'output.DATA.GOF16 = ( ...'
      write(9999,*) 'output.DATA.GOF16X+output.DATA.GOF16Y ...'
      write(9999,*) '+output.DATA.GOF16Z)./3;'
      write(9999,*) 'xlswrite(DFname, output.DATA.GOF16'', 
     +''SA_GOF'', ''B2'')'
      write(9999,*) 'xlswrite(DFname, specPer(3:18), 
     +''SA_GOF'', ''B1'')'

      write(9999,*) 'load GOF_InElEl_All.list'
      write(9999,*) 'numPerS = length(GOF_InElEl_All(:,1))/numSta;'

      write(9999,*) 'output.DATA.InEl_GOFX = 
     +reshape(GOF_InElEl_All(:,2),numPerS,numSta);'
      write(9999,*) 'xlswrite(DFname, GOF_InElEl_All(1:numPerS,1)'', 
     +''GOF_InElX'', ''B1'')'
      write(9999,*) 'xlswrite(DFname, output.DATA.InEl_GOFX'', 
     +''GOF_InElX'', ''B2'')'

      write(9999,*) 'output.DATA.InEl_GOFY = 
     +reshape(GOF_InElEl_All(:,3),numPerS,numSta);'
      write(9999,*) 'xlswrite(DFname, GOF_InElEl_All(1:numPerS,1)'', 
     +''GOF_InElY'', ''B1'')'
      write(9999,*) 'xlswrite(DFname, output.DATA.InEl_GOFY'', 
     +''GOF_InElY'', ''B2'')'


      write(9999,*) 'output.DATA.InEl_GOFZ = 
     +reshape(GOF_InElEl_All(:,4),numPerS,numSta);'
      write(9999,*) 'xlswrite(DFname, GOF_InElEl_All(1:numPerS,1)'', 
     +''GOF_InElZ'', ''B1'')'
      write(9999,*) 'xlswrite(DFname, output.DATA.InEl_GOFZ'', 
     +''GOF_InElZ'', ''B2'')'


 


      write(9999,*)'save(sprintf(''%s_DATA_STRUCT.mat'',DFname)
     +,''output'')'
C----------------------------------------------------------------------
      close(9999)
C----------------------------------------------------------------------
      end subroutine MFile
