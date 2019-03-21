      SUBROUTINE match(filt_flag)
      !----------------------------------------------------------
      ! Description:
      !
      !   matched filterin
      !   from Matched Filters by Robert Graves (Dec. 5, 2016)
      !   based on match.py in BBPlatform.
      !
      !   input: velocity [m/s]
      !          filt_flag:
      !                   0=LFs (lf)
      !                   1=HFs (conv_seis) 
      !
      !   based from "program band_pass_filter" in filt.f
      !   calling XAPIIR written by Dave Harris (September 12, 1990)
      !
      ! Authors: R. Takedatsu
      !
      ! Modified: February 2019 (v2.0)
      !
      use def_kind
      use scattering, only: npts
      use waveform

      implicit none

      integer(kind=i_single),intent(in)             :: filt_flag
      integer(kind=i_single)                        :: i,nt
      real(kind=r_single)                           :: dt
      real(kind=r_single),allocatable,dimension(:,:):: pri_filt_seis
      real(kind=r_single),allocatable,dimension(:,:):: pri_filt_lf
      real(kind=r_single),allocatable,dimension(:,:):: tmp_seis
      real(kind=r_single),allocatable,dimension(:,:):: acc_seis
      
      ! for HFs =====================================================================
      if (filt_flag .eq. 1) then

         nt=npts
         dt=lf_len/(v_npts-1) !  dt=lf_len/(nt-1)

         ! save original conv_seis in pri_filt_seis and clear conv_seis
         if (.not. allocated(pri_filt_seis)) then
            allocate(pri_filt_seis(nt,3))
         endif
         pri_filt_seis=conv_seis
         conv_seis(1:nt,1:3)=0.0

         if (.not. allocated(acc_seis)) then
            allocate(acc_seis(nt,3),tmp_seis(nt,3))
         endif

         ! vel2acc
         do i=1,3
            call vel2acc(pri_filt_seis(:,i),dt,nt,acc_seis(:,i))
         enddo

         ! HP filtering
         do i=1,3
            call bp_filter(acc_seis(:,i),nt,dt,   &
                           filt_flag,tmp_seis(:,i))
         enddo

         ! acc2velc
         do i=1,3
            call acc2velc(pri_filt_seis(1,i),tmp_seis(:,i),   &
                          dt,npts,conv_seis(:,i))
         enddo

         ! deallocate temporary arrays
         deallocate(acc_seis,tmp_seis) 
         deallocate(pri_filt_seis)

      endif

      ! for LFs =====================================================================
      if (filt_flag .eq. 0) then

         nt=lf_npts
         dt=lf_len/lf_npts

         ! save original lf_seis in pri_filt_lf and clear lf_seis
         if (.not. allocated(pri_filt_lf)) then
            allocate(pri_filt_lf(nt,3))
         endif
         pri_filt_lf=lf_seis
         lf_seis(1:nt,1:3)=0.0

         if (.not. allocated(acc_seis)) then
            allocate(acc_seis(nt,3),tmp_seis(nt,3))
         endif

         ! vel2acc
         do i=1,3
            call vel2acc(pri_filt_lf(:,i),dt,nt,acc_seis(:,i))
         enddo

         ! LP filtering
         do i=1,3
            call bp_filter(acc_seis(:,i),nt,dt,   &
                           filt_flag,tmp_seis(:,i))
         enddo

         ! acc2velc
         do i=1,3
            call acc2velc(pri_filt_lf(1,i),tmp_seis(:,i),   &
                          dt,nt,lf_seis(:,i))
         enddo

         ! deallocate temporary arrays
         deallocate(acc_seis,tmp_seis) 
         deallocate(pri_filt_lf)

      endif


      END SUBROUTINE match
!================================================================

      SUBROUTINE bp_filter(in_data,nt,dt,filt_flag,out_data)
      !----------------------------------------------------------
      ! Description:
      !
      !   Filtering LFs or HFs.
      !   Originally from filt_bbp.f (2017)
      !   Use for Matched Filters by Robert Graves (Dec. 5, 2016)
      !   based on match.py in BBPlatform.
      !
      ! Input:
      !   in_data:
      !   order: IORD = 4
      !   fh:   should be fl < fh here, in RG's filtering, fh < fl
      !   fl:
      !   nt:
      !   dt:
      !   phase: PASSES = 2
      !
      ! =========================================================
      !
      ! program band_pass_filter
      ! implicit none
      ! integer nt, ns, i, j, k, m, nrec
      ! real*4 fl, fh, dt, a, b, c, t
      ! real*4 x(600000), x1(200000,100), t1(200000)
      ! character s*128, sfl*32, sfh*32, snt*32, sdt*32

      use def_kind
      use matching, only: targ_fr ! add v163

      implicit none
      ! input
      integer(kind=i_single),intent(in)             :: nt,filt_flag
      real(kind=r_single),intent(in)                :: dt
      real(kind=r_single),dimension(nt),intent(in)  :: in_data
      ! output
      real(kind=r_single),dimension(nt),intent(out) :: out_data

      ! order (IORD) and phase (PASSES)
      integer(kind=r_single)                         :: ord,ph
      ! fl < fh here
      real(kind=r_single)                            :: fl,fh
      real(kind=r_single)                            :: val
      ! tmp array
      real(kind=r_single),dimension(nt)             :: tmp_data
      
      !----------------------------------------------------------

      if (filt_flag .eq. 0) then
         ord=4
         fl=0.0
         fh=targ_fr
         ph=2
         val=-1.0/(2.0*ord)
         fh=(fh*exp(val*log(sqrt(2.0)-1.0)))
         tmp_data = in_data
         CALL XAPIIR( tmp_data, nt,'BU', 0.0, 0.0, ord,'LP', &
                      fl, fh, dt, ph )

      elseif (filt_flag .eq. 1) then 
         ord=4
         fl=targ_fr
         ph=2
         val=1.0/(2.0*ord)
         fl=(fl*exp(val*log(sqrt(2.0)-1.0)))
         tmp_data = in_data
         fh=1.0e+15
         CALL XAPIIR( tmp_data, nt,'BU', 0.0, 0.0, ord,'HP', &
                   fl, fh, dt, ph )
      endif

      out_data = tmp_data

      END SUBROUTINE bp_filter 

! XAPIIR -- SUBROUTINE:   IIR FILTER DESIGN AND IMPLEMENTATION                  

!  AUTHOR:  Dave Harris                                                         
!  LAST MODIFIED:  September 12, 1990                                           
!  ARGUMENTS:                                                                   
!  ----------                                                                   
!    DATA           REAL ARRAY CONTAINING SEQUENCE TO BE FILTERED               
!                     ORIGINAL DATA DESTROYED, REPLACED BY FILTERED DATA        
!
!    NSAMPS         NUMBER OF SAMPLES IN DATA                                   
!
!    APROTO         CHARACTER*8 VARIABLE, CONTAINS TYPE OF ANALOG               
!                     PROTOTYPE FILTER                                          
!                     '(BU)TTER  ' -- BUTTERWORTH FILTER                        
!                     '(BE)SSEL  ' -- BESSEL FILTER                             
!                     'C1      ' -- CHEBYSHEV TYPE I                            
!                     'C2      ' -- CHEBYSHEV TYPE II                           
!
!    TRBNDW         TRANSITION BANDWIDTH AS FRACTION OF LOWPASS                 
!                   PROTOTYPE FILTER CUTOFF FREQUENCY.  USED                    
!                   ONLY BY CHEBYSHEV FILTERS.                                  
!
!    A              ATTENUATION FACTOR.  EQUALS AMPLITUDE                       
!                   REACHED AT STOPBAND EDGE.  USED ONLY BY                     
!                   CHEBYSHEV FILTERS.                                          
!
!    IORD           ORDER (#POLES) OF ANALOG PROTOTYPE                          
!                   NOT TO EXCEED 10 IN THIS CONFIGURATION.  4 - 5              
!                   SHOULD BE AMPLE.                                            
!
!    TYPE           CHARACTER*8 VARIABLE CONTAINING FILTER TYPE                 
!                     'LP' -- LOW PASS                                          
!                     'HP' -- HIGH PASS                                         
!                     'BP' -- BAND PASS                                         
!                     'BR' -- BAND REJECT                                       
!
!    FLO            LOW FREQUENCY CUTOFF OF FILTER (HERTZ)                      
!                   IGNORED IF TYPE = 'LP'                                      
!
!    FHI            HIGH FREQUENCY CUTOFF OF FILTER (HERTZ)                     
!                   IGNORED IF TYPE = 'HP'                                      
!
!    TS             SAMPLING INTERVAL (SECONDS)                                 
!
!    PASSES           INTEGER VARIABLE CONTAINING THE NUMBER OF PASSES          
!                   1 -- FORWARD FILTERING ONLY                                 
!                   2 -- FORWARD AND REVERSE (I.E. ZERO PHASE) FILTERING        



!  SUBPROGRAMS REFERENCED:  BILIN2, BUROOTS, WARP, CUTOFFS, LPTHP, LPTBP,       
!    LP, LPTBR, BEROOTS, C1ROOTS, C2ROOTS, CHEBPARM, DESIGN, APPLY              

      SUBROUTINE XAPIIR( DATA, NSAMPS, APROTO, TRBNDW, A, IORD, TYPE, &
                         FLO, FHI, TS, PASSES)

        DIMENSION DATA(1)                                               
        !CHARACTER*8 TYPE, APROTO                                        
        CHARACTER*2 TYPE, APROTO                                        
        INTEGER NSAMPS, PASSES, IORD                                    
        REAL*4 TRBNDW, A, FLO, FHI, TS, SN(30), SD(30)                  
        LOGICAL ZP                                                      

!-------------------                                                                               
!  Filter designed                                                              

!-------------------                                                                               

        CALL DESIGN( IORD, TYPE(1:2), APROTO(1:2), A, TRBNDW, &        
                     FLO, FHI, TS, SN, SD, NSECTS )                     

!-------------------                                                                               
!  Filter data                                                                  
                                                                       
!-------------------                                                                               

        IF (   PASSES .EQ. 1 ) THEN                                     
          ZP = .FALSE.                                                  
        ELSE                                                            
          ZP = .TRUE.                                                   
        END IF                                                          

        CALL APPLY( DATA, NSAMPS, ZP, SN, SD, NSECTS )                  

      RETURN                                                            

      END                                                                       

!                                                               APPLY           
!
!  Subroutine to apply an iir filter to a data sequence.                        
!    The filter is assumed to be stored as second order sections.               
!    Filtering is in-place.                                                     
!    Zero-phase (forward and reverse) is an option.                             
!
!  Input Arguments:                                                             
!  ----------------                                                             
!
!    DATA                           Array containing data                       
!
!    NSAMPS                         Number of data samples                      
!
!    ZP                             Logical variable, true for                  
!                                     zero phase filtering, false               
!                                     for single pass filtering                 
!
!    SN                             Numerator polynomials for second            
!                                     order sections.                           
!
!    SD                             Denominator polynomials for second          
!                                     order sections.                           
!
!    NSECTS                         Number of second-order sections             
!
!  Output Arguments:                                                            
!  -----------------                                                            
!
!    DATA                          Data array (same as input)                   

      SUBROUTINE APPLY( DATA, NSAMPS, ZP, SN, SD, NSECTS )              

        !REAL*4 SN(1), SD(1), DATA(1)                                    
        REAL*4 SN(*), SD(*), DATA(*)                                    
        REAL*4 OUTPUT                                                   
        LOGICAL ZP                                                      

        JPTR = 1                                                        

        DO    1 J = 1, NSECTS                                           
          X1 = 0.0                                                      
          X2 = 0.0                                                      
          Y1 = 0.0                                                      
          Y2 = 0.0                                                      
          B0 = SN(JPTR)                                                 
          B1 = SN(JPTR+1)                                               
          B2 = SN(JPTR+2)                                               
          A1 = SD(JPTR+1)                                               
          A2 = SD(JPTR+2)                                               

          DO    2 I = 1, NSAMPS                                         
            OUTPUT = B0*DATA(I) + B1*X1 + B2*X2                         
            OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )                         
            Y2 = Y1                                                     
            Y1 = OUTPUT                                                 
            X2 = X1                                                     
            X1 = DATA(I)                                                
            DATA(I) = OUTPUT                                            

    2     CONTINUE                                                      

          JPTR = JPTR + 3                                               

    1   CONTINUE                                                        

        IF (   ZP ) THEN                                                

          JPTR = 1                                                      

          DO    3 J = 1, NSECTS                                         
            X1 = 0.0                                                    
            X2 = 0.0                                                    
            Y1 = 0.0                                                    
            Y2 = 0.0                                                    
            B0 = SN(JPTR)                                               
            B1 = SN(JPTR+1)                                             
            B2 = SN(JPTR+2)                                             
            A1 = SD(JPTR+1)                                             
            A2 = SD(JPTR+2)                                             

            DO    4 I = NSAMPS, 1, -1                                   
              OUTPUT = B0*DATA(I) + B1*X1 + B2*X2                       
              OUTPUT = OUTPUT - ( A1*Y1 + A2*Y2 )                       
              Y2 = Y1                                                   
              Y1 = OUTPUT                                               
              X2 = X1                                                   
              X1 = DATA(I)                                              
              DATA(I) = OUTPUT                                          

    4       CONTINUE                                                    

            JPTR = JPTR + 3                                             

    3     CONTINUE                                                      

        END IF                                                          

      RETURN                                                            

      END                                                               


!  DESIGN -- Subroutine to design IIR digital filters from analog              
!    prototypes.
!
!
!  Input Arguments:                                                             
!  ----------------                                                             
!
!    IORD                Filter order (10 MAXIMUM)                              
!
!    TYPE                Character*2 variable containing filter type            
!                          LOWPASS (LP)                                         
!                          HIGHPASS (HP)                                        
!                          BANDPASS (BP)                                        
!                          BANDREJECT (BR)                                      
!
!    APROTO              Character*2 variable designating analog prototype      
!                          Butterworth (BU)                                     
!                          Bessel (BE)                                          
!                          Chebyshev Type I (C1)                                
!                          Chebyshev Type II (C2)                               
!
!    A                   Chebyshev stopband attenuation factor                  
!
!    TRBNDW              Chebyshev transition bandwidth (fraction of            
!                          lowpass prototype passband width)                    
!
!    FL                  Low-frequency cutoff                                   
!
!    FH                  High-frequency cutoff                                  
!
!    TS                  Sampling interval (in seconds)                         
!
!
!  Output Arguments:                                                            
!  -----------------                                                            
!
!    SN                  Array containing numerator coefficients of             
!                        second-order sections packed head-to-tail.             
!
!    SD                  Array containing denominator coefficients              
!                        of second-order sections packed head-to-tail.          
!
!    NSECTS              Number of second-order sections.                       


      SUBROUTINE DESIGN( IORD, TYPE, APROTO, A, TRBNDW, &
                        FL, FH, TS, SN, SD, NSECTS )                   

        COMPLEX P(10), Z(10)                                            
        CHARACTER*2 TYPE, APROTO                                        
        CHARACTER*3 STYPE(10)                                           
        !REAL*4 SN(1), SD(1)                                             
        REAL*4 SN(*), SD(*)                                             

!  Analog prototype selection                                                   

        IF (     APROTO .EQ. 'BU' ) THEN                                
          CALL BUROOTS( P, STYPE, DCVALUE, NSECTS, IORD )               
        ELSE IF (    APROTO .EQ. 'BE' ) THEN                            
          CALL BEROOTS( P, STYPE, DCVALUE, NSECTS, IORD )               
        ELSE IF (    APROTO .EQ. 'C1' ) THEN                            
          CALL CHEBPARM( A, TRBNDW, IORD, EPS, RIPPLE )                 
          CALL C1ROOTS( P, STYPE, DCVALUE, NSECTS, IORD, EPS )          
        ELSE IF (    APROTO .EQ. 'C2' ) THEN                            
          OMEGAR = 1. + TRBNDW                                          
          CALL C2ROOTS( P, Z, STYPE, DCVALUE, NSECTS, IORD, A, OMEGAR ) 
        END IF                                                          

!  Analog mapping selection                                                     

        IF (     TYPE .EQ. 'BP' ) THEN                                  
          FLW = WARP( FL*TS/2., 2. )                                    
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LPTBP( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )  
        ELSE IF (   TYPE .EQ. 'BR' ) THEN                               
          FLW = WARP( FL*TS/2., 2. )                                    
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LPTBR( P, Z, STYPE, DCVALUE, NSECTS, FLW, FHW, SN, SD )  
        ELSE IF (    TYPE .EQ. 'LP' ) THEN                              
          FHW = WARP( FH*TS/2., 2. )                                    
          CALL LP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )               
          CALL CUTOFFS( SN, SD, NSECTS, FHW )                           
        ELSE IF (    TYPE .EQ. 'HP' ) THEN                              
          FLW = WARP( FL*TS/2., 2. )                                    
          CALL LPTHP( P, Z, STYPE, DCVALUE, NSECTS, SN, SD )            
          CALL CUTOFFS( SN, SD, NSECTS, FLW )                           
        END IF                                                          

!  Bilinear analog to digital transformation                                    

        CALL BILIN2( SN, SD, NSECTS )                                   

      RETURN                                                            

      END                                                               


! BUROOTS -- SUBROUTINE TO COMPUTE BUTTERWORTH POLES FOR                        
!   NORMALIZED LOWPASS FILTER                                                   
!
!
! LAST MODIFIED:  SEPTEMBER 7, 1990                                             
!
!  OUTPUT ARGUMENTS:                                                            
!  -----------------                                                            
!
!      P              COMPLEX ARRAY CONTAINING POLES                            
!                       CONTAINS ONLY ONE FROM EACH                             
!                       COMPLEX CONJUGATE PAIR, AND                             
!                       ALL REAL POLES                                          
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
!                       TYPE:                                                   
!                         (SP)  SINGLE REAL POLE                                
!                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
!
!      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
!
!
!  INPUT ARGUMENTS:                                                             
!  ----------------                                                             
!
!      IORD           DESIRED FILTER ORDER                                      


      SUBROUTINE BUROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )             

        !COMPLEX P(1)                                                    
        COMPLEX P(*)                                                    
        INTEGER HALF                                                    
        !CHARACTER*3 RTYPE(1)                                            
        CHARACTER*3 RTYPE(*)                                            

        PI=3.14159265                                                   

        HALF = IORD/2                                                   

! TEST FOR ODD ORDER, AND ADD POLE AT -1                                        

        NSECTS = 0                                                      
        IF (    2*HALF .LT. IORD ) THEN                                 
          P(1) = CMPLX( -1., 0. )                                       
          RTYPE(1) = 'SP'                                               
          NSECTS = 1                                                    
        END IF                                                          

        DO    1  K = 1, HALF                                            
          ANGLE = PI * ( .5 + FLOAT(2*K-1) / FLOAT(2*IORD) )            
          NSECTS = NSECTS + 1
          P(NSECTS) = CMPLX( COS(ANGLE), SIN(ANGLE) )
          RTYPE(NSECTS) = 'CP'                                          
    1   CONTINUE                                                        

        DCVALUE = 1.0                                                   

      RETURN                                                            

      END  


! BEROOTS -- SUBROUTINE TO RETURN BESSEL POLES FOR                              
!   NORMALIZED LOWPASS FILTER                                                   
!
!
! LAST MODIFIED:  April 15, 1992. Changed P and RTYPE to adjustable 
!                 array by using an "*" rather than a "1".     
!
!  OUTPUT ARGUMENTS:                                                            
!  -----------------                                                            
!
!      P              COMPLEX ARRAY CONTAINING POLES                            
!                       CONTAINS ONLY ONE FROM EACH                             
!                       COMPLEX CONJUGATE PAIR, AND                             
!                       ALL REAL POLES                                          
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
!                       TYPE:                                                   
!                         (SP)  SINGLE REAL POLE                                
!                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
!
!      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
!
!  INPUT ARGUMENTS:                                                             
!  ----------------                                                             
!
!      IORD           DESIRED FILTER ORDER                                      



      SUBROUTINE BEROOTS( P, RTYPE, DCVALUE, NSECTS, IORD )             

        COMPLEX P(*)                                                    
        INTEGER NSECTS, IORD                                            
        CHARACTER*3 RTYPE(*)                                            

        IF (   IORD .EQ. 1 ) THEN                                       
          P(1) = CMPLX( -1.0, 0.0 )                                     
          RTYPE(1) = 'SP'                                               
        ELSE IF (  IORD .EQ. 2 ) THEN                                   
          P(1) = CMPLX( -1.1016013,  0.6360098 )                        
          RTYPE(1) = 'CP'                                               
        ELSE IF (  IORD .EQ. 3 ) THEN                                   
          P(1) = CMPLX( -1.0474091, 0.9992645 )                         
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3226758, 0.0 )                               
          RTYPE(2) = 'SP'                                               
        ELSE IF (  IORD .EQ. 4 ) THEN                                   
          P(1) = CMPLX( -0.9952088,  1.2571058 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3700679, 0.4102497 )                         
          RTYPE(2) = 'CP'                                               
        ELSE IF (  IORD .EQ. 5 ) THEN                                   
          P(1) = CMPLX( -0.9576766,  1.4711244 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3808774,  0.7179096 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.5023160, 0.0 )                               
          RTYPE(3) = 'SP'                                               
        ELSE IF (  IORD .EQ. 6 ) THEN                                   
          P(1) = CMPLX( -0.9306565,  1.6618633 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3818581,  0.9714719 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.5714904,  0.3208964 )                        
          RTYPE(3) = 'CP'                                               
        ELSE IF (  IORD .EQ. 7 ) THEN                                   
          P(1) = CMPLX( -0.9098678,  1.8364514 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3789032,  1.1915667 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.6120388,  0.5892445 )                        
          RTYPE(3) = 'CP'                                               
          P(4) = CMPLX( -1.6843682, 0.0 )                               
          RTYPE(4) = 'SP'                                               
        ELSE IF (  IORD .EQ. 8 ) THEN                                   
          P(1) = CMPLX( -0.8928710,  1.9983286 )                        
          RTYPE(1) = 'CP'                                               
          P(2) = CMPLX( -1.3738431,  1.3883585 )                        
          RTYPE(2) = 'CP'                                               
          P(3) = CMPLX( -1.6369417,  0.8227968 )                        
          RTYPE(3) = 'CP'                                               
          P(4) = CMPLX( -1.7574108,  0.2728679 )                        
          RTYPE(4) = 'CP'                                               
        END IF                                                          

        NSECTS = IORD - IORD/2                                          
        DCVALUE = 1.0                                                   

!  DONE                                                                         

      RETURN                                                            

      END                                                               

!  CHEBPARM - Calculates Chebyshev type I and II design parameters              
!
!  INPUT ARGUMENTS                                                              
!  ---------------                                                              
!
!       A                Desired stopband attenuation                           
!                          i.e. max stopband amplitude is 1/ATTEN               
!
!       TRBNDW           Transition bandwidth between stop and passbands        
!                          as a fraction of the passband width                  
!
!       IORD             Filter order (number of poles)                         
!
!  OUTPUT ARGUMENTS                                                             
!  ----------------                                                             
!
!       EPS              Chebyshev passband parameter                           
!
!       RIPPLE           Passband ripple                                        


      SUBROUTINE CHEBPARM( A, TRBNDW, IORD, EPS, RIPPLE )               

          OMEGAR  =  1. + TRBNDW                                        
          ALPHA = ( OMEGAR + SQRT( OMEGAR**2 - 1. ) ) ** IORD           
          G = ( ALPHA**2 + 1. ) / (2.*ALPHA)                            
          EPS = SQRT( A**2 - 1. ) / G                                   
          RIPPLE = 1. / SQRT( 1. + EPS**2 )                             

      RETURN                                                            

      END                                                               


! C1ROOTS -- SUBROUTINE TO COMPUTE CHEBYSHEV TYPE I POLES FOR                   
!   NORMALIZED LOWPASS FILTER                                                   
!
!
! LAST MODIFIED:  SEPTEMBER 7, 1990                                             
!
!
!  OUTPUT ARGUMENTS:                                                            
!  -----------------                                                            
!
!      P              COMPLEX ARRAY CONTAINING POLES                            
!                       CONTAINS ONLY ONE FROM EACH                             
!                       COMPLEX CONJUGATE PAIR, AND                             
!                       ALL REAL POLES                                          
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
!                       TYPE:                                                   
!                         (SP)  SINGLE REAL POLE                                
!                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
!
!      DCVALUE        RESPONSE OF FILTER AT ZERO FREQUENCY                      
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
!
!  INPUT ARGUMENTS:                                                             
!  ----------------                                                             
!
!      IORD           DESIRED FILTER ORDER                                      
!
!      EPS            CHEBYSHEV PARAMETER RELATED TO PASSBAND RIPPLE            

      SUBROUTINE C1ROOTS( P, RTYPE, DCVALUE, NSECTS, IORD, EPS )        

        COMPLEX P(1)                                                    
        INTEGER HALF                                                    
        CHARACTER*3 RTYPE(1)                                            

        PI = 3.14159265                                                 
        HALF = IORD/2                                                   

!  INTERMEDIATE DESIGN PARAMETERS                                               

        GAMMA = ( 1. + SQRT( 1. + EPS*EPS ) ) / EPS                     
        GAMMA = ALOG(GAMMA) / FLOAT(IORD)                               
        GAMMA = EXP(GAMMA)                                              
        S = .5 * ( GAMMA - 1./GAMMA )                                   
        C = .5 * ( GAMMA + 1./GAMMA )                                   

!  CALCULATE POLES                                                              

        NSECTS = 0                                                      

        DO    1  I = 1 ,  HALF                                          
          RTYPE(I) = 'CP'                                               
          ANGLE = FLOAT(2*I-1) * PI/FLOAT(2*IORD)                       
          SIGMA = -S * SIN(ANGLE)                                       
          OMEGA =  C * COS(ANGLE)                                       
          P(I) = CMPLX( SIGMA, OMEGA )                                  
          NSECTS = NSECTS + 1                                           
    1   CONTINUE                                                        
        IF (   2*HALF .LT. IORD ) THEN                                  
          RTYPE( HALF + 1 ) = 'SP'                                      
          P(HALF+1) = CMPLX( -S, 0.0 )                                  
          NSECTS = NSECTS + 1                                           
          DCVALUE = 1.0                                                 
        ELSE                                                            
          DCVALUE = 1./SQRT( 1 + EPS**2 )                               
        END IF                                                          

!  DONE                                                                         

      RETURN                                                            

      END                                                               


! C2ROOTS -- SUBROUTINE TO COMPUTE ROOTS FOR NORMALIZED LOWPASS                 
!   CHEBYSHEV TYPE 2 FILTER                                                     
!
!
! LAST MODIFIED:  SEPTEMBER 7, 1990                                             
!
!
!  OUTPUT ARGUMENTS:                                                            
!  -----------------                                                            
!
!      P              COMPLEX ARRAY CONTAINING POLES                            
!                       CONTAINS ONLY ONE FROM EACH                             
!                       COMPLEX CONJUGATE PAIR, AND                             
!                       ALL REAL POLES                                          
!
!      Z              COMPLEX ARRAY CONTAINING ZEROS                            
!                       CONTAINS ONLY ONE FROM EACH                             
!                       COMPLEX CONJUGATE PAIR, AND                             
!                       ALL REAL ZEROS                                          
!
!      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION              
!                       TYPE:                                                   
!                         (SP)  SINGLE REAL POLE                                
!                         (CP)  COMPLEX CONJUGATE POLE PAIR                     
!                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS               
!
!      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY                     
!
!      NSECTS         NUMBER OF SECOND ORDER SECTIONS                           
!
!
!  INPUT ARGUMENTS:                                                             
!  ----------------                                                             
!
!      IORD           DESIRED FILTER ORDER                                      
!
!      A              STOPBAND ATTENUATION FACTOR                               
!
!      OMEGAR         CUTOFF FREQUENCY OF STOPBAND                              
!                     PASSBAND CUTOFF IS AT 1.0 HERTZ                           


      SUBROUTINE C2ROOTS( P, Z, RTYPE, DCVALUE, NSECTS, IORD, A, OMEGAR)

        COMPLEX P(1), Z(1)                                              
        INTEGER HALF                                                    
        CHARACTER*3 RTYPE(1)                                            

        PI = 3.14159265                                                 

        HALF = IORD/2                                                   

!  INTERMEDIATE DESIGN PARAMETERS                                               

        GAMMA = (A+SQRT(A*A-1.))                                        
        GAMMA = ALOG(GAMMA)/FLOAT(IORD)                                 
        GAMMA = EXP(GAMMA)                                              
        S = .5*(GAMMA-1./GAMMA)                                         
        C = .5*(GAMMA+1./GAMMA)                                         

        NSECTS = 0                                                      

        DO    1 I = 1, HALF                                             

!  CALCULATE POLES                                                              

          RTYPE(I) = 'CPZ'                                              

          ANGLE = FLOAT(2*I-1) * PI/FLOAT(2*IORD)                       
          ALPHA = -S*SIN(ANGLE)                                         
          BETA = C*COS(ANGLE)                                           
          DENOM = ALPHA*ALPHA + BETA*BETA                               
          SIGMA = OMEGAR*ALPHA/DENOM                                    
          OMEGA = -OMEGAR*BETA/DENOM                                    
          P(I) = CMPLX( SIGMA, OMEGA )                                  

!  CALCULATE ZEROS                                                              

          OMEGA = OMEGAR/COS(ANGLE)                                     
          Z(I) = CMPLX( 0.0, OMEGA )                                    
          NSECTS = NSECTS + 1                                           

    1   CONTINUE                                                        

!  ODD-ORDER FILTERS                                                            

        IF (  2*HALF .LT. IORD ) THEN                                   
          RTYPE(HALF+1) = 'SP'                                          
          P(HALF+1) = CMPLX( -OMEGAR/S, 0.0 )                           
          NSECTS = NSECTS + 1                                           
        END IF                                                          

!  DC VALUE                                                                     

        DCVALUE = 1.0                                                   

!  DONE                                                                         

      RETURN                                                            

      END                                                               


! WARP -- FUNCTION, APPLIES TANGENT FREQUENCY WARPING TO COMPENSATE             
!         FOR BILINEAR ANALOG -> DIGITAL TRANSFORMATION                         
!
! ARGUMENTS:                                                                    
! ----------                                                                    
!
!      F       ORIGINAL DESIGN FREQUENCY SPECIFICATION (HERTZ)                  
!
!      TS      SAMPLING INTERVAL (SECONDS)                                      
!
!  LAST MODIFIED:  SEPTEMBER 20, 1990                                           


      REAL FUNCTION WARP( F , TS )                                      

        TWOPI = 6.2831853                                               
        ANGLE = TWOPI*F*TS/2.                                           
        WARP = 2.*TAN(ANGLE)/TS                                         
        WARP = WARP/TWOPI                                               

      RETURN                                                            

      END      

!  Subroutine to generate second order section parameterization                 
!    from an pole-zero description for lowpass filters.                         
!
!
!  Input Arguments:                                                             
!  ----------------                                                             
!
!    P                       Array containing poles                             
!
!    Z                       Array containing zeros                             
!
!    RTYPE                   Character array containing root type information   
!                              (SP)  Single real pole or                        
!                              (CP)  Complex conjugate pole pair                
!                              (CPZ) Complex conjugate pole and zero pairs      
!
!    DCVALUE                 Zero-frequency value of prototype filter           
!
!    NSECTS                  Number of second-order sections                    
!
!  Output Arguments:                                                            
!  -----------------                                                            
!
!    SN                      Numerator polynomials for second order             
!                              sections.                                        
!
!    SD                      Denominator polynomials for second order           
!                              sections.                                        


      SUBROUTINE LP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )             

        COMPLEX P(*), Z(*)                                              
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    

        IPTR = 1                                                        

        DO    1 I = 1, NSECTS                                           
          IF (   RTYPE(I) .EQ. 'CPZ' ) THEN                             
            SCALE = REAL( P(I) * CONJG( P(I) ) ) &                
                  / REAL( Z(I) * CONJG( Z(I) ) )                        
            SN( IPTR )     = REAL( Z(I) * CONJG( Z(I) ) ) * SCALE       
            SN( IPTR + 1 ) = -2. * REAL( Z(I) ) * SCALE                 
            SN( IPTR + 2 ) = 1. * SCALE                                 
            SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )               
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             

          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         

            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
            SN( IPTR )     = SCALE                                      
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P(I) * CONJG( P(I) ) )               
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             

          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          

            SCALE = -REAL( P(I) )                                       
            SN( IPTR )     = SCALE                                      
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = -REAL( P(I) )                              
            SD( IPTR + 1 ) = 1.                                         
            SD( IPTR + 2 ) = 0.                                         
            IPTR = IPTR + 3                                             

          END IF                                                        

    1   CONTINUE                                                        

        SN(1) = DCVALUE * SN(1)                                         
        SN(2) = DCVALUE * SN(2)                                         
        SN(3) = DCVALUE * SN(3)                                         

      RETURN                                                            

      END

!                                                    LPTBP                      
!
!                                                                               
!
!  Subroutine to convert an prototype lowpass filter to a bandpass filter via   
!    the analog polynomial transformation.  The lowpass filter is               
!    described in terms of its poles and zeros (as input to this routine).      
!    The output consists of the parameters for second order sections.           
!
!  Input Arguments:                                                             
!  ----------------                                                             
!
!    P                       Array containing poles                             
!
!    Z                       Array containing zeros                             
!
!    RTYPE                   Character array containing type information        
!                              (SP) single real pole  or                        
!                              (CP) complex conjugate pole pair  or             
!                              (CPZ) complex conjugate pole/zero pairs          
!
!    DCVALUE                 Zero frequency value of filter                     
!
!    NSECTS                  Number of second-order sections upon input         
!
!    FL                      Low-frequency cutoff                               
!
!    FH                      High-frequency cutoff                              
!
!  Output Arguments:                                                            
!  -----------------                                                            
!
!    SN                      Numerator polynomials for second order             
!                              sections.                                        
!
!    SD                      Denominator polynomials for second order           
!                              sections.                                        
!
!    NSECTS                  Number of second order sections upon output        
!                              This subroutine doubles the number of            
!                              sections.                                        

      SUBROUTINE LPTBP( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )  

        COMPLEX P(*), Z(*), CTEMP, P1, P2, Z1, Z2, S, H                 
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    

        PI = 3.14159265                                                 
        TWOPI = 2.*PI                                                   
        A = TWOPI*TWOPI*FL*FH                                           
        B = TWOPI*( FH - FL )                                           
        N = NSECTS                                                      
        NSECTS = 0                                                      
        IPTR = 1                                                        

        DO    1 I = 1, N                                                

          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN                            

            CTEMP = ( B*Z(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            Z1 = 0.5*( B*Z(I) + CTEMP )                                 
            Z2 = 0.5*( B*Z(I) - CTEMP )                                 
            CTEMP = ( B*P(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*P(I) + CTEMP )                                 
            P2 = 0.5*( B*P(I) - CTEMP )                                 
            SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z1 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z2 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            NSECTS = NSECTS + 2                                         

          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         

            CTEMP = ( B*P(I) )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*P(I) + CTEMP )                                 
            P2 = 0.5*( B*P(I) - CTEMP )                                 
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            NSECTS = NSECTS + 2                                         

          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          

            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = B                                          
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = A                                          
            SD( IPTR + 1 ) = -B*REAL( P(I) )                            
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            NSECTS = NSECTS + 1                                         

          END IF                                                        

    1   CONTINUE                                                        

!  Scaling - use the fact that the bandpass filter amplitude at sqrt( omega_l * 
!            equals the amplitude of the lowpass prototype at d.c.              

        S = CMPLX( 0., SQRT(A) )                                        
        H = CMPLX( 1., 0. )                                             
        IPTR = 1                                                        

        DO    2 I = 1, NSECTS                                           
          H = H * ( ( SN(IPTR+2)*S + SN(IPTR+1) )*S + SN(IPTR) ) &     
             / ( ( SD(IPTR+2)*S + SD(IPTR+1) )*S + SD(IPTR) )            
          IPTR = IPTR + 3                                               

    2   CONTINUE                                                        

        SCALE = DCVALUE / SQRT( REAL( H ) * CONJG( H ) )                

        SN(1) = SN(1) * SCALE                                           
        SN(2) = SN(2) * SCALE                                           
        SN(3) = SN(3) * SCALE                                           

      RETURN                                                            

      END                                                               

!                                                    LPTBR                      
!
!                                                                               
!
!  Subroutine to convert a lowpass filter to a band reject filter               
!    via an analog polynomial transformation.  The lowpass filter is            
!    described in terms of its poles and zeros (as input to this routine).      
!    The output consists of the parameters for second order sections.           
!
!  Input Arguments:                                                             
!  ----------------                                                             
!
!    P                       Array containing poles                             
!
!    Z                       Array containing zeros                             
!
!    RTYPE                   Character array containing type information        
!                              (SP)  single real pole or                        
!                              (CP)  complex conjugate pole pair                
!                              (CPZ) complex conjugate pole/zero pairs          
!
!    DCVALUE                 Zero-frequency value of prototype filter           
!
!    NSECTS                  Number of second-order sections                    
!                              prior to transformation                          
!
!    FL                      Low-frequency cutoff                               
!
!    FH                      High-frequency cutoff                              
!
!  Output Arguments:                                                            
!  -----------------                                                            
!
!    SN                      Numerator polynomials for second order             
!                              sections.                                        
!
!    SD                      Denominator polynomials for second order           
!                              sections.                                        
!
!    NSECTS                  Number of second order sections following          
!                              transformation.  The number is doubled.          


      SUBROUTINE LPTBR( P, Z, RTYPE, DCVALUE, NSECTS, FL, FH, SN, SD )  

        COMPLEX P(*), Z(*), CINV, CTEMP, P1, P2, Z1, Z2                 
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*)                                             

        PI = 3.14159265                                                 
        TWOPI = 2.*PI                                                   
        A = TWOPI*TWOPI*FL*FH                                           
        B = TWOPI*( FH - FL )                                           
        N = NSECTS                                                      
        NSECTS = 0                                                      
        IPTR = 1                                                        

        DO    1 I = 1, N                                                

          IF (    RTYPE(I) .EQ. 'CPZ' ) THEN                            

            CINV = 1./Z(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            Z1 = 0.5*( B*CINV + CTEMP )                                 
            Z2 = 0.5*( B*CINV - CTEMP )                                 
            CINV = 1./P(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*CINV + CTEMP )                                 
            P2 = 0.5*( B*CINV - CTEMP )                                 
            SN( IPTR )     = REAL( Z1 * CONJG( Z1 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z1 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = REAL( Z2 * CONJG( Z2 ) )                   
            SN( IPTR + 1 ) = -2. * REAL( Z2 )                           
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            NSECTS = NSECTS + 2                                         

          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         

            CINV = 1./P(I)                                              
            CTEMP = ( B*CINV )**2 - 4.*A                                
            CTEMP = CSQRT( CTEMP )                                      
            P1 = 0.5*( B*CINV + CTEMP )                                 
            P2 = 0.5*( B*CINV - CTEMP )                                 
            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P1 * CONJG( P1 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P1 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = REAL( P2 * CONJG( P2 ) )                   
            SD( IPTR + 1 ) = -2. * REAL( P2 )                           
            SD( IPTR + 2 ) = 1.                                         
            IPTR = IPTR + 3                                             
            NSECTS = NSECTS + 2                                         

          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          

            SN( IPTR )     = A                                          
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = 1.                                         
            SD( IPTR )     = -A*REAL( P(I) )                            
            SD( IPTR + 1 ) = B                                          
            SD( IPTR + 2 ) = -REAL( P(I) )                              
            IPTR = IPTR + 3                                             
            NSECTS = NSECTS + 1                                         

          END IF                                                        

    1   CONTINUE                                                        


!  Scaling - use the fact that the bandreject filter amplitude  at d.c.         
!            equals the lowpass prototype amplitude at d.c.                     

        H = 1.0                                                         
        IPTR = 1                                                        

        DO    2 I = 1, NSECTS                                           
          H = H * SN(IPTR) / SD(IPTR)                                   
          IPTR = IPTR + 3                                               
    2   CONTINUE                                                        
        SCALE = DCVALUE / ABS(H)                                        
        SN(1) = SN(1) * SCALE                                           
        SN(2) = SN(2) * SCALE                                           
        SN(3) = SN(3) * SCALE                                           

      RETURN                                                            

      END                                                               

!                                                    LPTHP                      
!
!                                                                               
!
!  Subroutine to convert a lowpass filter to a highpass filter via              
!    an analog polynomial transformation.  The lowpass filter is                
!    described in terms of its poles and zeroes (as input to this routine).     
!    The output consists of the parameters for second order sections.           
!
!  Input Arguments:                                                             
!  ----------------                                                             
!
!    P                       Array containing poles                             
!
!    Z                       Array containing zeroes                            
!
!    RTYPE                   Character array containing root type information   
!                              (SP) single real pole or                         
!                              (CP)  complex conjugate pair                     
!                              (CPZ) complex pole/zero pairs                    
!
!    DCVALUE                 Zero-frequency value of prototype filter           
!
!    NSECTS                  Number of second-order sections                    
!
!  Output Arguments:                                                            
!  -----------------                                                            
!
!    SN                      Numerator polynomials for second order             
!                              sections.                                        
!
!    SD                      Denominator polynomials for second order           
!                              sections.                                        

      SUBROUTINE LPTHP( P, Z, RTYPE, DCVALUE, NSECTS, SN, SD )          

        COMPLEX P(*), Z(*)                                              
        CHARACTER*3 RTYPE(*)                                            
        REAL*4 SN(*), SD(*), DCVALUE                                    

        IPTR = 1                                                        

        DO    1 I = 1, NSECTS                                           

          IF (     RTYPE(I) .EQ. 'CPZ' ) THEN                           

            SCALE = REAL( P(I) * CONJG( P(I) ) ) &                     
                  / REAL( Z(I) * CONJG( Z(I) ) )                        

            SN( IPTR )     = 1.  *  SCALE                               
            SN( IPTR + 1 ) = -2. * REAL( Z(I) )  *  SCALE               
            SN( IPTR + 2 ) = REAL( Z(I) * CONJG( Z(I) ) )  *  SCALE     
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )               
            IPTR = IPTR + 3                                             

          ELSE IF (   RTYPE(I) .EQ. 'CP' ) THEN                         

            SCALE = REAL( P(I) * CONJG( P(I) ) )                        
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = 0.                                         
            SN( IPTR + 2 ) = SCALE                                      
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -2. * REAL( P(I) )                         
            SD( IPTR + 2 ) = REAL( P(I) * CONJG( P(I) ) )               
            IPTR = IPTR + 3                                             

          ELSE IF (  RTYPE(I) .EQ. 'SP' ) THEN                          

            SCALE = -REAL( P(I) )                                       
            SN( IPTR )     = 0.                                         
            SN( IPTR + 1 ) = SCALE                                      
            SN( IPTR + 2 ) = 0.                                         
            SD( IPTR )     = 1.                                         
            SD( IPTR + 1 ) = -REAL( P(I) )                              
            SD( IPTR + 2 ) = 0.                                         
            IPTR = IPTR + 3                                             

          END IF                                                        

    1   CONTINUE                                                        

        SN(1) = SN(1) * DCVALUE                                         
        SN(2) = SN(2) * DCVALUE                                         
        SN(3) = SN(3) * DCVALUE                                         

      RETURN                                                            

      END             

!                                                    CUTOFFS                    
!
!
!  Subroutine to alter the cutoff of a filter.  Assumes that the                
!    filter is structured as second order sections.  Changes                    
!    the cutoffs of normalized lowpass and highpass filters through             
!    a simple polynomial transformation.                                        
!
!  Input Arguments:                                                             
!  ----------------                                                             
!
!    F                       New cutoff frequency                               
!
!  Input/Output Arguments:                                                      
!  -----------------------                                                      
!
!    SN                      Numerator polynomials for second order             
!                              sections.                                        
!
!    SD                      Denominator polynomials for second order           
!                              sections.                                        
!
!    NSECTS                  Number of second order sectionsects                
!                                                                               

      SUBROUTINE CUTOFFS( SN, SD, NSECTS, F )                           

        !REAL*4 SN(1), SD(1)                                             
        REAL*4 SN(*), SD(*)                                             

        SCALE = 2.*3.14159265*F                                         

        IPTR = 1                                                        

        DO    1 I = 1, NSECTS                                           

          SN( IPTR + 1 ) = SN( IPTR + 1 ) / SCALE                       
          SN( IPTR + 2 ) = SN( IPTR + 2 ) / (SCALE*SCALE)               
          SD( IPTR + 1 ) = SD( IPTR + 1 ) / SCALE                       
          SD( IPTR + 2 ) = SD( IPTR + 2 ) / (SCALE*SCALE)               
          IPTR = IPTR + 3                                               

    1   CONTINUE                                                        

      RETURN                                                            

      END                                                               


!  Transforms an analog filter to a digital filter via the bilinear transformati
!    Assumes both are stored as second order sections.  The transformation is   
!    done in-place.                                                             
!
!                                                                               
!
!  Input Arguments:                                                             
!  ----------------                                                             
!
!    SN                   Array containing numerator polynomial coefficients for
!                           second order sections.  Packed head-to-tail.        
!
!    SD                   Array containing denominator polynomial coefficients f
!                           second order sections.  Packed head-to-tail.        
!
!    NSECTS               Number of second order sections.                      
!                                                                               

      SUBROUTINE BILIN2( SN, SD, NSECTS )                               

        !REAL*4 SN(1), SD(1)                                             
        REAL*4 SN(*), SD(*)                                             

        IPTR = 1                                                        

        DO    1 I = 1, NSECTS                                           

          A0 = SD(IPTR)                                                 
          A1 = SD(IPTR+1)                                               
          A2 = SD(IPTR+2)                                               

          SCALE = A2 + A1 + A0                                          
          SD(IPTR)   = 1.                                               
          SD(IPTR+1) = (2.*(A0 - A2)) / SCALE                           
          SD(IPTR+2) = (A2 - A1 + A0) / SCALE                           

          A0 = SN(IPTR)                                                 
          A1 = SN(IPTR+1)                                               
          A2 = SN(IPTR+2)                                               

          SN(IPTR)   = (A2 + A1 + A0) / SCALE                           
          SN(IPTR+1) = (2.*(A0 - A2)) / SCALE                           
          SN(IPTR+2) = (A2 - A1 + A0) / SCALE                           

          IPTR = IPTR + 3                                               

    1   CONTINUE                                                        

      RETURN                                                            

      END                                                               

