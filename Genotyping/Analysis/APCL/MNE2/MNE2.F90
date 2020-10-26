!-----------------------------------------------------------------------!
!                                                                       !
!                               M N e                                   !
!                                                                       !
!-----------------------------------------------------------------------!
!       MNE uses maximum likelihood and moment methods to estimate      !
!   the effective size and migration rate of a metapopulation from      !
!   temporal data on gene frequencies. As a special case, it can also   !
!   estimate the effective size for a single isolated population only.  !
!   For more details, see WANG & WHITLOCK (2003, Genetics) and          !
!   Wang (2001, Genet. Res).                                            !
!                                                                       !
!   Notes:                                                              !
!   1. The maximum Ne allowed can be increased to 2 folds if the        !
!     variables TRMX,TRMX2,COL_LOG in the module is defined "SP"        !
!   2. The codes are compiled with OPENMP in Intel's compiler. For other!
!     compilers without openmp, the few openmp lines must be deleted    !
!                                                                       !
!   Update Records:                                                     !
!   Date        Changes                                                 !
!   13/10/04    Discarding alleles with no copies in all focal samples  !
!   28/01/05    Fixed a bug in calculating PBDIS, L703-709              !
!   28/04/05    Refixed the bug in calculating PBDIS, L697-721          !
!                                                                       !
!       If you detect any bugs, please report to:                       !
!   Dr. Jinliang Wang                                                   !
!   Institute of Zoology                                                !
!   Regent's Park, London NW1 4RY, UK                                   !
!   Tel: +44 20 74496620  Email: jinliang.wang@ioz.ac.uk                !
!-----------------------------------------------------------------------!
      MODULE DATA
      INTEGER,PARAMETER:: SP=KIND(1.0)
      INTEGER,PARAMETER:: DP=SELECTED_REAL_KIND(2*PRECISION(1._SP))
      INTEGER:: NLOCI,NTLOCI,NSAMB,MAXLC,NE20,M_ESTIMATE,  &
        NPARA(0:10000,2),NFUNC,NINTPA,INIT_POINT,MONITOR,NumThread
      REAL(DP):: PREC,PREC_LOG,FPARA(0:10000)
      CHARACTER*10:: FNAME    
      INTEGER,ALLOCATABLE::NSAMP(:,:),NGEN(:),NS2(:,:),IBPPAA(:,:),KALLELE(:),IB(:,:)
      REAL(DP),ALLOCATABLE::TRMX(:),TRMX2(:),COL_LOG(:,:)
      REAL(DP),ALLOCATABLE::FACTOR(:),PPAA(:,:),FAC(:),YFC_TEM(:,:),PA0(:),PT(:,:,:)
      REAL(DP),PRIVATE:: XY(2)
      CHARACTER*4,PRIVATE:: ANSWER
      CONTAINS
        SUBROUTINE READ_IN_DATA
          OPEN(UNIT=10,FILE='MNE_DATA')
          READ(10,*,IOSTAT=IOS,ERR=999,End=999)M_ESTIMATE
          READ(10,*,IOSTAT=IOS,ERR=999,End=999)MaximumNe
          NE20=2*MaximumNe
          READ(10,*,IOSTAT=IOS,ERR=999,End=999)Monitor
          READ(10,*,IOSTAT=IOS,ERR=999,End=999)NumThread          
          READ(10,*,IOSTAT=IOS,ERR=999,End=999)NLOCI
          ALLOCATE(KALLELE(1:NLOCI))
          READ(10,*,IOSTAT=IOS,ERR=999,End=999)(KALLELE(JJ),JJ=1,NLOCI)
          MAXLC=SUM(KALLELE)
          READ(10,*,IOSTAT=IOS,ERR=999,End=999)NSAMB
          ALLOCATE(NGEN(NSAMB),NS2(MAXLC,0:NSAMB),NSAMP(MAXLC,0:NSAMB))
          READ(10,*,IOSTAT=IOS,ERR=999,End=999)(NGEN(I),I=1,NSAMB)
          DO I=1,NSAMB
            N=0
            DO J=1,NLOCI
              READ(10,*,IOSTAT=IOS,ERR=999,End=999) &
                (NSAMP(JJ,I),JJ=N+1,N+KALLELE(J))
              N2=SUM(NSAMP(N+1:N+KALLELE(J),I))
              NS2(N+1:N+KALLELE(J),I)=N2
              N=N+KALLELE(J)
            ENDDO
          ENDDO
          IF(M_ESTIMATE==1)THEN
            N=0
            DO J=1,NLOCI
              READ(10,*,IOSTAT=IOS,ERR=999,End=999) &
                (NSAMP(JJ,0),JJ=N+1,N+KALLELE(J))
              N2=SUM(NSAMP(N+1:N+KALLELE(J),0))
              NS2(N+1:N+KALLELE(J),0)=N2
              N=N+KALLELE(J)
            ENDDO
            READ(10,*,IOSTAT=IOS,ERR=999,End=999) INIT_POINT
          ELSE
            INIT_POINT=1
          ENDIF
          KK=SUM(KALLELE(1:NLOCI))
          JJ=NLOCI
          N=KK+1
          DO J=NLOCI,1,-1
            DO K=KALLELE(J),1,-1
              N=N-1
              IF(ALL(NSAMP(N,1-M_ESTIMATE:NSAMB)==0).OR.  &
                 ALL(NSAMP(N,1-M_ESTIMATE:NSAMB)==NS2(N,1-M_ESTIMATE:NSAMB)))THEN
                DO I=N+1,KK
                  NSAMP(I-1,1-M_ESTIMATE:NSAMB)=NSAMP(I,1-M_ESTIMATE:NSAMB)
                  NS2(I-1,1-M_ESTIMATE:NSAMB)=NS2(I,1-M_ESTIMATE:NSAMB)
                ENDDO
                KALLELE(J)=KALLELE(J)-1
              ENDIF
            ENDDO
            IF(KALLELE(J)==0)THEN
              DO K=J+1,JJ
                KALLELE(K-1)=KALLELE(K)
              ENDDO
              NLOCI=NLOCI-1
            ENDIF
          ENDDO
          FNAME='MNE_OUT'
          PREC=1.E-30_DP
          PREC_LOG=LOG(PREC)
          NINTPA=100
          CLOSE(10) 
          J=Max(Maxval(NS2(:MAXLC,1-M_ESTIMATE:NSAMB)),NE20)        
          ALLOCATE(fac(0:J))                  
          fac(0:1)=0._DP
          DO i=2,J
            fac(i)=fac(i-1)+LOG(DBLE(i))
          END DO          
          XY(1)=1._DP/NE20
          IF(M_ESTIMATE==1)THEN
            XY(2)=0.1_DP
          ELSE
            XY(2)=0._DP
          ENDIF         
          CALL TrmxDim(NE20,XY(1),XY(2),MAXDIM)                      
          IF(MAXDIM>0.AND.MAXDIM<200000000)THEN
            ALLOCATE(TRMX(0:MAXDIM),STAT=KK)
          ELSE
            KK=1
          ENDIF
          IF(KK/=0)THEN
            ISTEP=-MIN(MAX(NINT(NE20*.02),2),500)
            DO I=NE20,2,ISTEP
              IF(ALLOCATED(TRMX))DEALLOCATE(TRMX)
              XY(1)=1._DP/I
              CALL TrmxDim(I,XY(1),XY(2),MAXDIM)      
              IF(MAXDIM>0.AND.MAXDIM<200000000)ALLOCATE(TRMX(0:MAXDIM),STAT=KK)                            
              IF(KK==0)EXIT
            ENDDO
            WRITE(*,'(a,I8/a/A,I8)') &
              'The maximum value of Ne set in your data is:',NE20/2, &
              ' which is too large that your computer has not enough RAM to run the program. ',&
              ' The suggested maximum value of Ne is:', NINT(I*.5)
              WRITE(*,'(/A/A/A)')'The program will run with this maximum value of Ne. However,', &
                ' if memory problem still exists in the run, please set maximum value of Ne', &
                ' in your data input to 90% of the value suggested, and re-run the program'              
            WRITE(*,'(/A)')'Do you accept this new maximum value of Ne? (Y/N)'
            READ(*,*)ANSWER
            IF(ANSWER=='Y'.OR.ANSWER=='y')THEN
              NE20=I
              MaximumNe=NINT(I*.5)
            ELSE
              STOP
            ENDIF
          ENDIF
          IF(ALLOCATED(TRMX))DEALLOCATE(TRMX)                    
          GOTO 888
999       WRITE(*,*)'Errors in DATA INPUT. Insufficient data or incorrect', &
            '  format. Please check your DATA in file MNE_DATA and re-run'
          STOP
888     END SUBROUTINE READ_IN_DATA
      END MODULE DATA
!
      PROGRAM MAIN
      USE DATA
      CHARACTER*3:: PARANAME(2)
      CHARACTER*12:: DATE,TIME       
      REAL(SP):: MLNE2,MLM,MTNE,MTM
      REAL(DP):: X1,Y1,DY,fmax,fmax2,xmax(2),xmax2(2),VarStep(2),XII(2,2), &
        XI(2,2),XIII(2,2),FUNC,MAX_LD,CI95(-100:100,2,2),RANGE(2,2),XA(3),YA(3)
      INTEGER:: BOUND(2,2)
      DATA  XII/8._DP,2*0._DP,0.05_DP/,PARANAME/' Ne','  m'/, &
        XIII/4._DP,2*0._DP,0.001_DP/        
      CALL DATE_AND_TIME(DATE,TIME)        
      CALL READ_IN_DATA      
      ALLOCATE(FACTOR(MAXLC),IBPPAA(MAXLC,2),YFC_TEM(MAXLC,NSAMB),PA0(MAXLC))
      NPARA(0,1)=-1
      NTLOCI=SUM(KALLELE(1:NLOCI))
      CALL MTEST(MTNE,MTM)
      IF(MONITOR>0)THEN
        IF(M_ESTIMATE==1)THEN
          write(*,'(A,F12.2,A,F12.4)') &
            'Moment estimates: Ne=',MTNE,'  m=',MTM
        ELSE
          write(*,'(A,F12.2)')'Moment estimates: Ne=',MTNE
        ENDIF
      ENDIF
      IF(MTNE<0.)WRITE(*,*)'Ne<0 indicates an infinitely large estimate of Ne'
      NTLOCI=0
      DO I=1,NLOCI
        IF(KALLELE(I)==2)THEN
          X1=1._DP
        ELSE
          X1=1._DP-1._DP/KALLELE(I)
        ENDIF
        DO J=1,KALLELE(I)
          NTLOCI=NTLOCI+1
          IF(KALLELE(I)==2.AND.J==2)THEN
            FACTOR(NTLOCI)=-1._DP
          ELSE
            FACTOR(NTLOCI)=X1
          ENDIF
        ENDDO
      ENDDO
      K=0
      DO I=1,NTLOCI
        IF(FACTOR(I)<0._DP.AND.I<=NTLOCI-K)THEN
          DO J=I+1,NTLOCI
            FACTOR(J-1)=FACTOR(J)
            NSAMP(J-1,0:NSAMB)=NSAMP(J,0:NSAMB)
            NS2(J-1,0:NSAMB)=NS2(J,0:NSAMB)
          END DO
          K=K+1
        ENDIF
      END DO
      NTLOCI=NTLOCI-K
      DO I=1,NTLOCI
        IF(NSAMP(I,0)>NS2(I,0)*0.5) &
          NSAMP(I,0:NSAMB)=NS2(I,0:NSAMB)-NSAMP(I,0:NSAMB)
      ENDDO
      IDUM=-1
      IF(M_ESTIMATE==1)THEN
        NINTPA=NINTPA/2
        PA0(1:NTLOCI)=NSAMP(1:NTLOCI,0)/DBLE(NS2(1:NTLOCI,0))
        DO
          NINTPA=NINTPA*2
          ALLOCATE(PPAA(NTLOCI,0:NINTPA))
          CALL GTDIS(NINTPA_FLAG)
          IF(NINTPA_FLAG==0)EXIT
          DEALLOCATE(PPAA)
        ENDDO
      ENDIF
      NFUNC=0
      LFLAG=0
      IF(MTNE>0.5.AND.MTNE<=NE20*.5)THEN
        XMAX(1)=MTNE*2._DP  
      ELSEIF(MTNE>NE20*.5.AND.MTNE<NE20)THEN
        XMAX(1)=MTNE
      ELSE
        XMAX(1)=RAND1(IDUM)*NE20*.1+1._DP
        LFLAG=1
      ENDIF
      IF(MTM<=0.001)THEN
        XMAX(2)=0.001
      ELSEIF(MTM>0.999)THEN
        XMAX(2)=0.999
      ELSE
        XMAX(2)=MTM
      ENDIF
      XMAX2=XMAX
      X1=FUNC(XMAX,2,NINT(XMAX(1)))
      IF(X1>1.E99_DP)THEN
        DO
          XMAX(1)=RAND1(IDUM)*NE20*.1+1._DP
          IF(LFLAG>0)THEN
            LFLAG=LFLAG+1
            XMAX(2)=XMAX2(2)
          ELSE
            XMAX(2)=RAND1(IDUM)
          ENDIF
          FMAX=FUNC(XMAX,2,NINT(XMAX(1)))
          IF(FMAX<1.E99_DP)EXIT
          IF(LFLAG>1000)LFLAG=0
        ENDDO
      ENDIF
      XI=XII
      CALL POWELL(xmax,XI,ITER,MAX_LD,M_ESTIMATE,0,MONITOR)
      MLNE2=XMAX(1)
      MLM=XMAX(2)
      DO II=1,INIT_POINT-1   !Trying different random starting points
        DO
          XMAX(1)=RAND1(IDUM)*NE20*.1+1._DP
          XMAX(2)=RAND1(IDUM)
          FMAX=FUNC(XMAX,2,NINT(XMAX(1)))
          IF(FMAX>1.E99_DP)CYCLE
          IF(MONITOR==1)WRITE(*,'(A,I4,F12.2,F9.4)') &
            'ML estimate with random starting point', &
            II,(XMAX(J),J=1,1+M_ESTIMATE)
          XI=XII
          CALL POWELL(xmax,XI,ITER,FMAX,M_ESTIMATE,0,MONITOR)
          IF(FMAX<MAX_LD)THEN
            MAX_LD=FMAX
            MLNE2=XMAX(1)
            MLM=XMAX(2)
          ENDIF
          IF(MONITOR==1)THEN
            IF(M_ESTIMATE==1)THEN
              write(*,'(A,F12.2,2(A,F12.4)/)')'  Estimate: Ne=',Xmax(1)*.5, &
                '    m=',XMAX(2),'   Log(ML)=',-FMAX
            ELSE
              write(*,'(A,F12.2,A,F12.4/)')'  Estimate: Ne=',Xmax(1)*.5, &
                '   Log(ML)=',-FMAX
            ENDIF
          ENDIF
          EXIT
        ENDDO
      ENDDO
      IF(MONITOR>0)THEN
        IF(M_ESTIMATE==1)THEN
          WRITE(*,'(A,A)')'          ML_Ne           ML_m', &
            '        Log(ML)          MT_Ne           MT_m'
          WRITE(*,'(2(F15.2,2F15.4)/)')MLNE2/2,MLM,-MAX_LD,MTNE,MTM
        ELSE
          WRITE(*,'(A)')'          ML_Ne        Log(ML)          MT_Ne'
          WRITE(*,'(F15.2,2F15.4/)')MLNE2/2,-MAX_LD,MTNE
        ENDIF
        IF(MTNE<0.)WRITE(*,*)'Ne<0 indicates an infinitely large estimate of Ne'
      ENDIF
      RANGE(1,1)=2._DP
      RANGE(2,1)=NE20
      RANGE(:,2)=(/0.001_DP,1.0_DP/)
      K2=0
111   IF(K2/=0)THEN
        IF(M_ESTIMATE==1)THEN
          XI=XIII
          CALL POWELL(xmax,XI,ITER,MAX_LD,M_ESTIMATE,0,MONITOR)
        ELSE
          MAX_LD=FMAX
        ENDIF
        MLNE2=XMAX(1)
        MLM=XMAX(2)
      ENDIF
      CI95=0._DP
      CI95(0,1,1)=MLNE2
      CI95(0,1,2)=MLM
      CI95(0,2,:)=MAX_LD
      VarStep(1)=MAX(MLNE2*.01_DP,2._DP)
      VarStep(2)=MAX(MLM*0.05_DP,0.01_DP)
      BOUND=0
! Finding 95% confidence interval for each parameter
      DO K=1,1+M_ESTIMATE
        DO K2=2,1,-1
          IF(K2==1)THEN
            WRITE(*,'(/A,A)') &
              'Finding +95% confidence interval for parameter',PARANAME(K)
            XSTEP=VarStep(K)
            KSTEP=1
          ELSE
            WRITE(*,'(/A,A)') &
              'Finding -95% confidence interval for parameter',PARANAME(K)
            XSTEP=-VarStep(K)
            KSTEP=-1
          ENDIF
          XMAX2(:)=CI95(0,1,:)
          FMAX2=CI95(0,2,K)
          KK=0
          YSTEP=0.005
          DO
            XMAX=XMAX2
            XSTEP2=XSTEP
            IF(YSTEP<=0.08)then
              YSTEP2=YSTEP
              YSTEP=YSTEP+YSTEP
            ENDIF
            DO
              XMAX(K)=XMAX(K)+XSTEP2
              IF(XMAX(K)<RANGE(1,K).OR.XMAX(K)>RANGE(2,K))EXIT
              FMAX=FUNC(XMAX,2,NINT(XMAX(1)))
              IF(FMAX+0.01_DP<CI95(0,2,K)) GOTO 111
              X1=ABS(FMAX-FMAX2)
              IF(X1>YSTEP+YSTEP)THEN
                IF(ABS(XSTEP2)>XIII(K,K)+1.E-5)THEN
                  XMAX(K)=XMAX(K)-XSTEP2
                  XSTEP2=MAX(ABS(XSTEP2*.5_DP),XIII(K,K))
                  XSTEP2=SIGN(XSTEP2,XSTEP)
                ELSE
                  EXIT
                ENDIF
              ELSEIF(X1>YSTEP)THEN
                XSTEP2=XSTEP2*.5
                EXIT
              ELSE IF(X1<YSTEP2)THEN
                XTEM=MIN(XMAX(K)*.1,DBLE(ABS(XSTEP2+XSTEP2)))
                XTEM=MAX(DBLE(XTEM),XIII(K,K))
                XSTEP2=SIGN(XTEM,XSTEP2)
              ENDIF
            ENDDO
            XSTEP=XSTEP2
            IF((XMAX(K)<RANGE(1,K).and.CI95(0,1,k)<1.1*RANGE(1,K)).OR. &
              (XMAX(K)>RANGE(2,K).and.CI95(0,1,k)>.99*RANGE(2,K))) EXIT
            IF(XMAX(K)<RANGE(1,K))THEN
              XMAX(K)=RANGE(1,K)
              KFLAG=1
            ELSEIF(XMAX(K)>RANGE(2,K))THEN
              XMAX(K)=RANGE(2,K)
              KFLAG=1
            ELSE
              KFLAG=-1
            ENDIF
            IF(M_ESTIMATE==1)THEN
              XI=XIII
              CALL POWELL(XMAX,XI,ITER,FMAX,M_ESTIMATE,K,MONITOR)
              IF(FMAX+0.01_DP<CI95(0,2,K)) GOTO 111
            ELSE
              IF(KFLAG==1)FMAX=FUNC(XMAX,2,NINT(XMAX(1)))
            ENDIF
            KK=KK+KSTEP
            CI95(KK,2,K)=FMAX
            CI95(KK,1,K)=XMAX(K)
            X1=FMAX-2._DP-CI95(0,2,K)
            IF(X1>-0.001_DP.AND.X1<0.001_DP)THEN
              BOUND(K2,K)=KK
              EXIT
            ELSEIF(X1>=0.001_DP.OR.X1>-0.1)THEN
              IF(ABS(KK)>2)THEN
                NPL=3
              ELSE
                NPL=2
              ENDIF
              XA(1:NPL)=CI95(KK:KK-(NPL-1)*KSTEP:-KSTEP,2,K)
              YA(1:NPL)=CI95(KK:KK-(NPL-1)*KSTEP:-KSTEP,1,K)
              X1=CI95(0,2,K)+2._DP
              CALL polint(xa,ya,NPL,x1,y1,dy)
              CI95(KK,1,K)=Y1
              CI95(KK,2,K)=X1
              BOUND(K2,K)=KK
              EXIT
            ENDIF
            IF(KFLAG==1)EXIT
            XMAX2=XMAX
            FMAX2=FMAX
          ENDDO
        ENDDO
      ENDDO
!
      OPEN(UNIT=10,FILE=fname)
      WRITE(10,'(A/)')'                   Output from the program'
      IF(M_ESTIMATE==0)THEN
        WRITE(10,*) &
          'Estimating Ne only, assuming a single isolated population (m==0)'
      ELSE
        WRITE(10,*) &
          'Estimating Ne and m jointly, assuming an infinite source'
      ENDIF
      WRITE(10,'(/A/A)') &
        '**************** Likelihood estimates ***************', &
        '   PARAMETER      POINT_EST         +95%CI         -95%CI'
      IF(ALL(BOUND(:,1)/=0))THEN
        WRITE(10,'(A,3ES15.4)')'          Ne',CI95(0,1,1)/2, &
          CI95(BOUND(1,1),1,1)/2,CI95(BOUND(2,1),1,1)/2
      ELSEIF(BOUND(1,1)/=0)THEN
        WRITE(10,'(A,2ES15.4,A,ES12.4)')'          Ne',CI95(0,1,1)/2, &
          CI95(BOUND(1,1),1,1)/2,'  <',RANGE(1,1)/2
      ELSEIF(BOUND(2,1)/=0)THEN
        WRITE(10,'(A,ES15.4,A,ES12.4,ES15.4)')'          Ne',CI95(0,1,1)/2, &
          '  >',RANGE(2,1)/2,CI95(BOUND(2,1),1,1)/2
      ELSE
        WRITE(10,'(A,ES15.4,A,ES12.4,A,ES12.4)')'          Ne', &
          CI95(0,1,1),'  >',RANGE(2,1)/2,'  <',RANGE(1,1)/2
      ENDIF
      IF(M_ESTIMATE==1) THEN
        IF(ALL(BOUND(:,2)/=0))THEN
          WRITE(10,'(A,3ES15.4)')'           m',CI95(0,1,2), &
            CI95(BOUND(1,2),1,2),CI95(BOUND(2,2),1,2)
        ELSEIF(BOUND(1,2)/=0)THEN
          WRITE(10,'(A,2ES15.4,A,ES12.4)')'           m',CI95(0,1,2), &
            CI95(BOUND(1,2),1,2),'  <',RANGE(1,2)
        ELSEIF(BOUND(2,2)/=0)THEN
          WRITE(10,'(A,ES15.4,A,ES12.4,F15.4)')'           m',CI95(0,1,2), &
            '  >',RANGE(2,2),CI95(BOUND(2,2),1,2)
        ELSE
          WRITE(10,'(A,ES15.4,A,ES12.4,A,ES12.4)')'           m', &
          CI95(0,1,2),'  >',RANGE(2,2),'  <',RANGE(1,2)
        ENDIF
      ENDIF
      WRITE(10,'(/A/A)') &
        'Relative profile log-likelihood as a function of Ne', &
        '              Ne          Log(L)'
      DO I=-100,100
        IF(ABS(CI95(I,1,1))>1.E-50_DP) &
          WRITE(10,'(F16.2,F16.4)')CI95(I,1,1)/2,CI95(0,2,1)-CI95(I,2,1)
      ENDDO
      IF(M_ESTIMATE==1)THEN
        WRITE(10,'(/A/A)') &
          'Relative profile log-likelihood as a function of  m', &
          '               m          Log(L)'
        DO I=-100,100
          IF(ABS(CI95(I,1,2))>1.E-50_DP) &
            WRITE(10,'(F16.4,F16.4)')CI95(I,1,2),CI95(0,2,2)-CI95(I,2,2)
        ENDDO
      ENDIF
      WRITE(10,'(/A/)')'**************** Moment     estimates ***************'
      IF(M_ESTIMATE==1)THEN
        IF(MTNE>0.AND.MTM>-1.)THEN
          WRITE(10,'(A,F12.2,A,F12.4)')'MT_Ne=',MTNE,'    MT_m=',MTM
        ELSEIF(MTNE>0.)THEN
          WRITE(10,'(A,F12.2,A)')'MT_Ne=',MTNE,'    NO valid MT estimate of m'
        ELSEIF(MTM>-1.)THEN
          WRITE(10,'(A,F12.2)')' NO valid MT estimate of Ne, MT_m=',MTM
        ELSE
          WRITE(10,'(A)')' NO valid MT estimates of Ne and m'
        ENDIF
      ELSE
        IF(MTNE>0)THEN
          WRITE(10,'(A,F12.2,A)')'MT_Ne=',MTNE
        ELSE
          WRITE(10,'(A)')' Estimate of Ne INFINITE'
        ENDIF
      ENDIF
      IF(MTNE<0.)WRITE(10,*)'Ne<0 indicates an infinitely large estimate of Ne'
      WRITE(10,'(//)')
      WRITE(10,*)'Program started:  Date(d/m/y)=',date(7:8),'/',date(5:6), &
        '/',date(1:4),',  Time(s/m/hr)=',time(5:6),'/',time(3:4),'/',time(1:2)
      CALL DATE_AND_TIME(DATE,TIME)
      WRITE(10,*)'Program finished: Date(d/m/y)=',date(7:8),'/',date(5:6), &
        '/',date(1:4),',  Time(s/m/hr)=',time(5:6),'/',time(3:4),'/',time(1:2)        
      CLOSE(10)
      WRITE(*,*)'Output to file:',FNAME
      END PROGRAM MAIN
!
      SUBROUTINE MTEST(MTNE,MTM)
      USE DATA
      REAL(SP),PARAMETER:: TINY=1.E-8,BIG=1.-TINY
      REAL(SP):: MTNE,MTM,MTNE2,MTM2
      MTNE=0.
      MTM=0.
      NN=0
      DO KK=1,NSAMB-1
        KK1=KK+1
        IF(M_ESTIMATE==1)THEN
          DBT=0.
          WNE=0.
          DAB0=0.
          WW=0.
          WF=0.
        ELSE
          FC=0.
          NLOCI_FC=0
        ENDIF
        DO I=1,NTLOCI
          B0=NSAMP(I,KK)/REAL(NS2(I,KK))
          BT=NSAMP(I,KK1)/REAL(NS2(I,KK1))
          IF(M_ESTIMATE==1)THEN
            A0=NSAMP(I,0)/REAL(NS2(I,0))
            XT=abs(A0-B0)
            DAB0=DAB0-XT*XT
            IF(XT>TINY)THEN
              WI=(A0-B0)**2
              WW=WW+WI
              FI=(BT-B0)/(A0-B0)
              WF=WF+WI*FI
            ENDIF
            PQ=(B0*(1.-B0)+BT*(1.-BT))*0.5
            WNE=WNE+PQ
            DBT=DBT+(BT-B0)**2-B0*(1-B0)/NS2(I,KK)-BT*(1.-BT)/NS2(I,KK1)
          ELSE
            XT=BT+B0
            IF(XT>=TINY.AND.XT<2.)THEN
              FC=FC+2.*(BT-B0)**2*(1./XT+1./(2.-XT))-1./NS2(I,KK)-1./NS2(I,KK1)
              NLOCI_FC=NLOCI_FC+1
            ENDIF
          ENDIF
        END DO
        NT=NGEN(KK1)-NGEN(KK)
        IF(M_ESTIMATE==1)THEN
          IF(WW<=0.)CYCLE
          NN=NN+1
          MTM2=1.-ABS(1.-WF/WW)**(1./nt)
          MTM=MTM+MTM2
          IF(ABS(MTM2)>TINY)THEN
            MTNE2=(DBT+DAB0*(1.-(1.-MTM2)**NT)**2)/(WNE*  &
                 (1.-(1.-MTM2)**(2*NT))/(MTM2*(2.-MTM2)))
          ELSE
            MTNE2=DBT/(WNE*NT)
          ENDIF
        ELSE
          MTNE2=FC/NLOCI_FC/REAL(NT)
          NN=NN+1
        ENDIF
        IF(MTNE2>0.)MTNE=MTNE+MTNE2
      ENDDO
      IF(NN==0)THEN
        MTNE=-1.
        MTM=-10.
      ELSE
        IF(MTNE>TINY)THEN
          MTNE=NN/MTNE*.5
        ELSE
          MTNE=-1.
        ENDIF
        IF(M_ESTIMATE==1) MTM=MTM/NN
      ENDIF
      END SUBROUTINE MTEST
!
      SUBROUTINE GTDIS(NINTPA_FLAG)
      USE DATA
      REAL(DP):: X2,X3,xmax,xmax2,YLOG
      INTEGER:: LSTEP(2),LBOUND(2)
      DATA LSTEP/-1,1/,LBOUND/0,1/
      LBOUND(2)=NINTPA
      DO L=1,NTLOCI
        NSPA=NS2(L,0)
        N1=NSAMP(L,0)
        N2=NSPA-N1
        X2=FAC(NSPA)-NSPA*(FAC(NINTPA)-FAC(NINTPA-1))-FAC(N1)-FAC(N2)
        XMAX=-1._DP
        NI=NINT(N1*NINTPA/DBLE(NSPA))
        IF(NI<=0)THEN
          NI=1
        ELSEIF(NI>=NINTPA)THEN
          NI=NINTPA-1
        ENDIF
        DO JJ=1,2
          IBPPAA(L,JJ)=LBOUND(JJ)
          DO K=NI+JJ-1,LBOUND(JJ),LSTEP(JJ)
            CALL BINOM(X3,YLOG,X2,K,NINTPA,N1,N2)
            IF(K==0.AND.N1==0)THEN
              IF(ANY(NSAMP(L,1:NSAMB)>0))THEN
                X3=0._DP
                YLOG=-1.E300_DP
              ENDIF
            ELSE IF(K==NINTPA.AND.N1==NSPA)THEN
              IF(ANY(NSAMP(L,1:NSAMB)<NS2(L,1:NSAMB)))THEN
                X3=0._DP
                YLOG=-1.E300_DP
              ENDIF
            ENDIF
            PPAA(L,K)=X3
            IF(X3>XMAX)THEN
              XMAX=X3
              XMAX2=YLOG+PREC_LOG
              JMAX=K
            ELSEIF(YLOG<XMAX2)THEN
              IBPPAA(L,JJ)=K
              EXIT
            ENDIF
          END DO 
        ENDDO
        XMAX2=XMAX*PREC
        DO JJ=JMAX,LBOUND(1),-1
          IF(PPAA(L,JJ)<XMAX2)THEN
            IBPPAA(L,1)=JJ
            EXIT
          ENDIF
        END DO 
      ENDDO
      NINTPA_FLAG=0
      DO L=1,NTLOCI
        X3=SUM(PPAA(L,IBPPAA(L,1):IBPPAA(L,2)))
         IF(X3>1.E-30_DP)THEN
          X3=1._DP/X3
          DO i=IBPPAA(L,1),IBPPAA(L,2)
            PPAA(L,I)=PPAA(L,I)*X3
          ENDDO
        ELSE
          NINTPA_FLAG=1
          EXIT
        ENDIF
      END DO
      END SUBROUTINE GTDIS
!
      FUNCTION func(XMIN,NP,IDIM)
      USE DATA
      REAL(DP):: X0,X1,X2,X3,XT,Y,XMAX,XMAX2,XMIN(NP),XDIM, &
        FUNCL(MAXLC),FUNC,YLOG,PBDIS(0:IDIM)
      INTEGER:: IA(2),LSTEP(2),LBOUNDS(2),IFLAG(MAXLC),MAXDIM
      INTEGER:: OMP_GET_NUM_PROCS
      INTEGER, PARAMETER:: CHUNK=1
      SAVE XDIM,MAXDIM
      DATA LSTEP/-1,1/,LBOUNDS/0,1/                  
      IF(M_ESTIMATE==0)XMIN(2)=0._DP
      IF(XMIN(1)<2.OR.XMIN(1)>NE20.OR.XMIN(2)<0._DP  &
          .OR.XMIN(2)>1._DP)THEN
        FUNC=1.E100_DP
        RETURN
      ENDIF
      IMM=NINT(XMIN(2)*1000)
      DO I=1,NFUNC
        IF(IDIM==NPARA(I,1).AND.IMM==NPARA(I,2))THEN
          FUNC=FPARA(I)
          RETURN
        ENDIF
      END DO
      LBOUNDS(2)=IDIM
      IF(NPARA(NFUNC,1)/=IDIM)THEN
        XDIM=1._DP/IDIM
        DO I=1,NSAMB
          DO J=1,NTLOCI
            NNN=NS2(J,I)
            YFC_TEM(J,I)=FAC(NNN)-(FAC(IDIM)-FAC(IDIM-1))*NNN
          ENDDO
        ENDDO
        CALL TrmxDim(IDIM,XDIM,XMIN(2),MAXDIM)
      ENDIF
      ALLOCATE(TRMX(0:MAXDIM),IB(3,0:IDIM),STAT=KK)
      IF(KK/=0)THEN
        WRITE(*,'(a/a/A)') &
         'The maximum value of Ne set in your data is too large', &
          ' that your computer has not enough RAM to run the program. ',&
          ' Try a smaller maximum value of Ne and re-run the program'
        STOP
      ENDIF
      ALLOCATE(PT(NTLOCI,0:1,0:IDIM))
      FUNCL(1:NTLOCI)=0._DP 
      IF(M_ESTIMATE==0)THEN
        PBDIS=1._DP/(IDIM+1)
      ELSE 
        X0=(IDIM+IDIM)*XMIN(2)
        X3=GAMMLN(X0)
      ENDIF
      DO IPA=0,NINTPA   
        IF(M_ESTIMATE==1)THEN   
          IF(IPA>0.AND.IPA<NINTPA)THEN
            XT=DBLE(IPA)/NINTPA
            X1=X0*XT-1._DP
            X2=X0*(1._DP-XT)-1._DP 
            PBDIS(0)=EXP(X3-GAMMLN(X1+2.D0)-GAMMLN(X2+1.D0))*XDIM**(X1+1.D0)
            PBDIS(IDIM)=EXP(X3-GAMMLN(X1+1.D0)-GAMMLN(X2+2.D0))*XDIM**(X2+1.D0)            
            Y=PBDIS(0)+PBDIS(IDIM)  
            XT=EXP(X3-GAMMLN(X1+1.D0)-GAMMLN(X2+1.D0))*XDIM
            DO I=1,IDIM-1
              PBDIS(I)=XT*(I*XDIM)**X1*(1.D0-I*XDIM)**X2
              Y=Y+PBDIS(I)
            ENDDO          
            Y=1._DP/Y
            DO I=0,IDIM
              PBDIS(I)=PBDIS(I)*Y
            ENDDO
          ELSEIF(IPA==0)THEN
            PBDIS(0)=1.D0
            PBDIS(1:)=0.D0
          ELSE
            PBDIS(IDIM)=1.D0
            PBDIS(0:IDIM-1)=0.D0
          ENDIF
        ENDIF                    
        IF(XMIN(2)>0._DP)THEN
          JJ=0
          DO I=1,NTLOCI
            IF(IPA>=IBPPAA(I,1).AND.IPA<=IBPPAA(I,2))THEN
              IFLAG(I)=1
              JJ=1
            ELSE
              IFLAG(I)=0
            ENDIF
          END DO
          IF(JJ==0)CYCLE
        ENDIF                 
        CALL GTRANS(IDIM,XDIM,XMIN(2),IPA,MAXDIM)        
! Open MP parallel for multiple processors
        IF(NumThread<=0)then
          call OMP_SET_NUM_THREADS(OMP_GET_NUM_PROCS())
        ELSE
          call OMP_SET_NUM_THREADS(NumThread)
        ENDIF
!$OMP PARALLEL DEFAULT(SHARED) &
 !$OMP private(NA1,NA2,JM,XMAX,X1,X2,K,IA,J,YLOG,JMAX,XMAX2,KFLAG,KK,XT,MARK,NG)
!$OMP DO SCHEDULE(dynamic,CHUNK)         
        DO I=1,NTLOCI     
          IF(XMIN(2)>0._DP.AND.IFLAG(I)==0) CYCLE
          NA1=NSAMP(I,1)
          NA2=NS2(I,1)-NA1
          JM=NINT(NA1*IDIM/DBLE(NS2(I,1)))         
          IF(JM<=0)THEN
            JM=1
          ELSEIF(JM>=IDIM)THEN
            JM=IDIM-1
          ENDIF          
          XMAX=-1._DP
          X2=YFC_TEM(I,1)-FAC(NA1)-FAC(NA2)                           
          DO K=1,2        
            IA(K)=LBOUNDS(K)
            MARK=0                       
            DO J=JM+K-1,LBOUNDS(K),LSTEP(K)
              IF(MARK==1)CYCLE
              CALL BINOM(X1,YLOG,X2,J,IDIM,NA1,NA2)                         
              X1=X1*PBDIS(J)              
              PT(I,0,J)=X1                  
              IF(X1>XMAX)THEN
                XMAX=X1
                JMAX=J
                XMAX2=X1*PREC
              ELSEIF(X1<XMAX2)THEN
                IA(K)=J
                MARK=1
              ENDIF
            END DO             
          ENDDO                            
          DO NG=1,NGEN(NSAMB)                                     
            CALL DISTR(IDIM,JMAX,IA,I)                         
            KFLAG=0
            DO KK=2,NSAMB
              IF(NGEN(KK)/=NG)CYCLE
              KFLAG=1
              JM=JMAX
              IF(JM<=0)THEN
                JM=1
              ELSEIF(JM>=IDIM)THEN
                JM=IDIM-1
              ENDIF
              NA1=Nsamp(I,KK)
              NA2=NS2(I,KK)-NA1
              XT=YFC_TEM(I,KK)-FAC(NA1)-FAC(NA2)
              XMAX=-1._DP              
              DO K=1,2
                MARK=0
                DO J=JM+K-1,IA(K),LSTEP(K)
                  IF(MARK==0)then
                    CALL BINOM(X1,YLOG,XT,J,IDIM,NA1,NA2)
                    X2=PT(I,1,J)*X1
                    PT(I,0,J)=X2
                    IF(X2>XMAX)THEN
                      XMAX=X2
                      JMAX=J
                      XMAX2=X2*PREC
                    ELSEIF(X2<XMAX2)THEN
                      IA(K)=J
                      MARK=1
                    ENDIF
                  ENDIF
                END DO
              ENDDO               
            ENDDO            
            IF(KFLAG==0) PT(I,0,IA(1):IA(2))=PT(I,1,IA(1):IA(2))
          END DO
          X2=SUM(PT(I,0,IA(1):IA(2)))
          IF(XMIN(2)>0._DP)THEN
            FUNCL(I)=FUNCL(I)+X2*PPAA(I,IPA)
          ELSE
            FUNCL(I)=X2
          ENDIF
        ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL          
        IF(XMIN(2)<=0._DP)EXIT
      ENDDO
      DEALLOCATE(TRMX,PT,IB)
      FUNC=0._DP
      DO I=1,NTLOCI
        IF(FUNCL(I)>1.E-300_DP)THEN
          FUNC=FUNC-LOG(FUNCL(I))*FACTOR(I)
        ELSE
          FUNC=1.E100_DP
          EXIT
        ENDIF
      END DO
      NFUNC=NFUNC+1
      NPARA(NFUNC,1)=IDIM
      NPARA(NFUNC,2)=IMM
      FPARA(NFUNC)=FUNC
      IF(MONITOR>1.AND.IDIM>5000.or.MONITOR>2)THEN
        IF(M_ESTIMATE==1)THEN
          WRITE(*,'(A,I10,A,F8.4,A,ES16.8)') &
            'Ne=',IDIM/2,'  m=',XMIN(2),'  Log(L)=',-FUNC
        ELSE
          WRITE(*,'(A,I10,A,ES16.8)') &
            'Ne=',IDIM/2,'     Log(L)=',-FUNC
        ENDIF
      ENDIF
      END FUNCTION func
!
      SUBROUTINE TrmxDim(IDIM,XDIM,Rmig,MAXDIM)
      USE DATA
      REAL(DP):: X1,XDIM,Rmig
      INTEGER::  IDIM   
        IF(IDIM>500)THEN
          X1=PREC_LOG*2._DP*XDIM+2._DP*(FAC(IDIM/2)-FAC(IDIM/2-1))
          DO J=IDIM/2+1,IDIM-1
            IF(FAC(J)-FAC(J-1)+FAC(IDIM-J)-FAC(IDIM-J-1)<=X1)THEN
              X1=(J*2._DP-IDIM)*IDIM
              EXIT
            ENDIF
          ENDDO
          IF(X1<2147483646)THEN
            MAXDIM=NINT(X1)
            IF(M_ESTIMATE/=1)THEN
              IF(IDIM<1000)THEN
                MAXDIM=INT(MAXDIM*0.89)
              ELSE
                MAXDIM=INT(MAXDIM*0.8)
              ENDIF
            ELSE
              IF(Rmig<=0.1) MAXDIM=INT(MAXDIM*0.89)
            ENDIF
          ELSE
            MAXDIM=-1
          ENDIF
        ELSE
          MAXDIM=(IDIM+1)*(IDIM+1)
        ENDIF
      END SUBROUTINE TrmxDim
!
      SUBROUTINE DISTR(NR,JMAX,IA,Locus)
      USE DATA
      REAL(DP):: X1,XMAX,XMAX2
      INTEGER:: IA(2),LSTEP(2),LBOUNDS(2)
      DATA LSTEP/-1,1/,LBOUNDS/0,1/           
      LBOUNDS(2)=NR
      JJM=JMAX
      k3=IA(1)
      k4=IA(2)
      XMAX=-1._DP      
      DO KI=1,2     
        IA(KI)=LBOUNDS(KI)
        DO i=JJM+KI-1,LBOUNDS(KI),LSTEP(KI)
          L1=IB(1,I)         
          L=IB(2,I)-L1+1
          K1=IB(3,I)
          K2=K1+L-1
          J1=MAX(K1,k3)
          J2=MIN(K2,k4)
          X1=0._DP
          KL=L1-K1
          DO J=J1,J2
            X1=X1+TRMX(J+KL)*Pt(Locus,0,J)
          ENDDO
          Pt(Locus,1,i)=X1         
          IF(X1>XMAX)THEN
            XMAX=X1
            JMAX=I
            XMAX2=X1*1.E-15_DP
          ELSEIF(X1<XMAX2)THEN
            IA(KI)=I
            EXIT
          ENDIF
        ENDDO       
      ENDDO
      NN=IA(1)
      DO I=JMAX-1,NN,-1
        IF(PT(Locus,1,I)<XMAX2)THEN
          IA(1)=I
          EXIT
        ENDIF
      ENDDO
      END SUBROUTINE DISTR
!
      SUBROUTINE POWELL(P,XI,ITER,FRET,M_ESTIMATE,IVAR,MONITOR)
      INTEGER,PARAMETER:: DP=SELECTED_REAL_KIND(2*PRECISION(1.))
      INTEGER,PARAMETER:: N=2,ITMAX=200
      REAL(DP):: DEL,FP,FPTT,T,PT(N),PTT(N),XIT(N),FRET,P(N),XI(N,N),FUNC
      REAL(DP),PARAMETER:: FTOL=1.E-6_DP
      EXTERNAL FUNC 
      FRET=FUNC(P,N,NINT(P(1)))
      PT=P
      ITER=0
      DO
        ITER=ITER+1
        FP=FRET
        IBIG=0
        DEL=0._DP
        DO I=1,N
          IF(I/=IVAR)THEN
            XIT(1:N)=XI(1:N,I)
            FPTT=FRET
            CALL LINMIN(P,XIT,N,FRET,FTOL)
            IF(ABS(FPTT-FRET)>DEL)THEN
              DEL=ABS(FPTT-FRET)
              IBIG=I
            ENDIF
            IF(MONITOR>0)WRITE(*,'(1x,2(a,i3),a,f16.6)') &
             '    Iter=',ITER,' Line minimization=',I,'   Log(L)=',-FRET
            IF(M_ESTIMATE==0) RETURN
          ENDIF
        END DO
        IF(IVAR>0)RETURN
        IF(2._DP*ABS(FP-FRET)<=FTOL*(ABS(FP)+ABS(FRET)))RETURN
        IF(ITER==ITMAX)WRITE(*,*)'POWELL EXCEEDING MAXIMUM ITERATIONS'
        PTT(:)=2._DP*P(:)-PT(:)
        XIT(:)=P(:)-PT(:)
        PT(:)=P(:)
        FPTT=FUNC(PTT,N,NINT(PTT(1)))
        IF(FPTT>=FP)CYCLE
        T=2._DP*(FP-2._DP*FRET+FPTT)*(FP-FRET-DEL)**2-DEL*(FP-FPTT)**2
        IF(T>=0._DP)CYCLE
        CALL LINMIN(P,XIT,N,FRET,FTOL)
        IF(IVAR==N)then
          NN=N-1
        ELSE
          NN=N
        ENDIF
        XI(:,IBIG)=XI(:,NN)
        XI(:,NN)=XIT(:)
      ENDDO
      END SUBROUTINE POWELL
!
      SUBROUTINE LINMIN(P,XI,N,FRET,TOL)
      INTEGER,PARAMETER:: DP=SELECTED_REAL_KIND(2*PRECISION(1.))
      REAL(DP):: FRET,P(N),XI(N),TOL,F1DIM
      REAL(DP):: AX,BX,FA,FB,FX,XMIN,XX,PCOM(2),XICOM(2),BRENT
      COMMON/F1COM/PCOM,XICOM,NCOM
      EXTERNAL F1DIM
      NCOM=2
      PCOM=P
      XICOM=XI
      AX=1._DP
      XX=0._DP
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
      XI(1:N)=XMIN*XI(1:N)
      P(1:N)=P(1:N)+XI(1:N)
      END SUBROUTINE LINMIN
!
      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUN)
      INTEGER,PARAMETER:: DP=SELECTED_REAL_KIND(2*PRECISION(1.))
      REAL(DP):: AX,BX,CX,FA,FB,FC,FUN,GOLD,GLIMIT,TINY,DUM,FU,Q,R,U,ULIM
      EXTERNAL FUN
      PARAMETER(GOLD=1.618034_DP,GLIMIT=100._DP,TINY=1.E-20_DP)
      FA=FUN(AX)
      FB=FUN(BX)
      IF(FB>FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUN(CX)
      DO
        IF(FB<FC)RETURN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2._DP*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX)>0._DP)THEN
          FU=FUN(U)
          IF(FU<FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            RETURN
          ELSEIF(FU>FB)THEN
            CX=U
            FC=FU
            RETURN
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUN(U)
        ELSE IF((CX-U)*(U-ULIM)>0._DP)THEN
          FU=FUN(U)
          IF(FU<FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUN(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX)>=0._DP)THEN
          U=ULIM
          FU=FUN(U)
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUN(U)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
      ENDDO
      END SUBROUTINE MNBRAK
!
      FUNCTION F1DIM(X)
      INTEGER,PARAMETER:: DP=SELECTED_REAL_KIND(2*PRECISION(1.))
      REAL(DP):: F1DIM,FUNC,X,PCOM(2),XICOM(2),XT(2)
      COMMON/F1COM/PCOM,XICOM,NCOM
      XT(1:NCOM)=PCOM(1:NCOM)+X*XICOM(1:NCOM)
      F1DIM=FUNC(XT,2,NINT(XT(1)))
      END FUNCTION F1DIM
!
      FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
      INTEGER,PARAMETER:: DP=SELECTED_REAL_KIND(2*PRECISION(1.))
      REAL(DP):: BRENT,AX,BX,CX,TOL,XMIN,F,CGOLD,ZEPS, &
        A,B,D,E,ETEMP,FU,FV,FW,FX,P,Q,R,TOL1,TOL2,U,V,W,X,XXMM,TINY
      PARAMETER(ITMAX=100,CGOLD=0.381966_DP,ZEPS=1.0E-10_DP,TINY=1.E-100_DP)
      EXTERNAL F
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0._DP
      FX=F(X)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XXMM=0.5_DP*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2._DP*TOL1
        IF(ABS(X-XXMM)<=(TOL2-0.5_DP*(B-A)))GOTO 3
        IF(ABS(E)>TOL1)THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2._DP*(Q-R)
          IF(Q>0._DP)P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P)>=ABS(0.5_DP*Q*ETEMP).OR.P<=Q*(A-X).OR. &
            P>=Q*(B-X))GOTO 1
          D=P/Q
          U=X+D
          IF(U-A<TOL2.OR.B-U<TOL2)D=SIGN(TOL1,XXMM-X)
          GOTO 2
        ENDIF
1       IF(X>=XXMM)THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D)>=TOL1)THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U)
        IF(FU<=FX)THEN
          IF(U>=X)THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U<X)THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU<=FW.OR.ABS(W-X)<TINY)THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSEIF(FU<=FV.OR.ABS(V-X)<TINY.OR.ABS(V-W)<TINY)THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      WRITE(*,*)'BRENT EXCEED MAXIMUM ITERATIONS'
3     XMIN=X
      BRENT=FX
      END FUNCTION BRENT
!
      SUBROUTINE BINOM(X,YLOG,XT,J,IDIM,NA1,NA2)
      USE DATA
      REAL(DP):: X,YLOG,XT
      X=0._DP
      YLOG=-1.E300_DP
      IF(J==0)THEN
        IF(NA1==0)THEN
          X=1._DP
          YLOG=0._DP
        ENDIF
      ELSEIF(J==IDIM)THEN
        IF(NA2==0)THEN
          X=1._DP
          YLOG=0._DP
        ENDIF
      ELSE
        YLOG=XT+NA1*(FAC(J)-FAC(J-1))+NA2*(FAC(IDIM-J)-FAC(IDIM-J-1))
        IF(YLOG>-100._DP)X=EXP(YLOG)
      ENDIF
      END SUBROUTINE BINOM
!
      SUBROUTINE GTRANS(NR,XR,XM,IPA,MAXDIM)
      USE DATA
      REAL(DP):: X1,X3,XMAX,XMAX2,XR,XM,XM1,XM2,XM3,XM4,Y1,YLOG,PTT(0:NR)
      INTEGER:: LSTEP(2),LBOUNDS(2),LB(2),KB(2),LB2,MAXDIM2
      PARAMETER(TINY=1.E-9,BIG=1.-TINY)
      integer OMP_GET_NUM_PROCS,CHUNK
      PARAMETER (CHUNK=1)      
      DATA LSTEP/-1,1/,LBOUNDS/0,1/,NR_OLD/-1/
      SAVE NR_OLD
      LBOUNDS(2)=NR
      XM1=(1._DP-XM)*XR
      XM2=IPA*XM/NINTPA
      XM3=XM2*NR
      XM4=1._DP-XM
      IF(ABS(XM4)<TINY)THEN
        XM4=1._DP
      ELSE
        XM4=1._DP/XM4
      ENDIF
      IF(NR/=NR_OLD)THEN
        IF(NR_OLD/=-1)DEALLOCATE(COL_LOG)
        ALLOCATE(COL_LOG(0:NR,0:3))
      ENDIF
      X1=XM2/PREC
      X3=-XM1
      DO J=0,NR
        X3=X3+XM1
        Y1=X3+XM2
        COL_LOG(J,0)=Y1
        IF(Y1>TINY.AND.Y1<BIG)THEN
          IF(NR/=NR_OLD.OR.X1>=X3)THEN
            COL_LOG(J,1)=LOG(Y1)
            COL_LOG(J,2)=LOG(1.-Y1)
            COL_LOG(J,3)=Y1/(1.-Y1)
          ENDIF
        ENDIF
      ENDDO
      M=0
! Open MP parallel for multiple processors
      IF(NumThread<=0)then
        call OMP_SET_NUM_THREADS(OMP_GET_NUM_PROCS())
      ELSE
        call OMP_SET_NUM_THREADS(NumThread)
      ENDIF
!$OMP PARALLEL DEFAULT(SHARED) &
 !$OMP private(XMAX,XMAX2,i2,X3,JJ,KB,LB,L,J,Y1,X1,YLOG,PTT,LB2,NN,JMAX)
!$OMP DO SCHEDULE(dynamic,CHUNK)        
      DO I=0,NR
        XMAX=-1._DP
        I2=NR-I
        x3=FAC(NR)-fac(i)-fac(I2)
        JJ=NINT((I-XM3)*XM4)
        IF(JJ<=0)THEN
          JJ=1
        ELSEIF(JJ>=NR)THEN
          JJ=NR-1
        ENDIF
        KB(1)=JJ-3
        IF(KB(1)<0)KB(1)=0
        KB(2)=JJ+3
        IF(KB(2)>NR)KB(2)=NR
        DO L=1,2
          LB(L)=LBOUNDS(L)
          LB2=0
          DO j=JJ+L-1,LBOUNDS(L),LSTEP(L)
            IF(LB2/=0)CYCLE
            Y1=COL_LOG(J,0)
            X1=0._DP
            YLOG=-1.E300_DP
            IF(Y1<=TINY)THEN
              IF(I==0)THEN
                X1=1._DP
                YLOG=0._DP
              ENDIF
            ELSEIF(Y1>=BIG)THEN
              IF(I==NR)THEN
                X1=1._DP
                YLOG=0._DP
              ENDIF
            ELSE
              YLOG=X3+I*COL_LOG(J,1)+I2*COL_LOG(J,2)
              IF(YLOG>-100)X1=EXP(YLOG)
            ENDIF
            PTT(J)=X1
            IF(X1>XMAX)THEN
              XMAX=X1
              XMAX2=YLOG+PREC_LOG
              JMAX=J
            ELSEIF(YLOG<XMAX2)THEN
              LB(L)=J
              LB2=1
            ENDIF
            IF(XMAX<TINY.AND.J==KB(L))THEN
              LB(L)=J
              LB2=1
            ENDIF
          ENDDO
        ENDDO
        IF(XMAX>=TINY)THEN
          NN=LB(1)
          XMAX2=XMAX*PREC
          LB2=0
          DO J=JMAX-1,NN,-1
            IF(LB2/=0)CYCLE
            IF(PTT(J)<XMAX2)THEN
              LB(1)=J
              LB2=1
            ENDIF
          ENDDO
        ELSE
          LB=JMAX
        ENDIF
!$OMP CRITICAL 
        IB(1,I)=M+1        
        DO J=LB(1),LB(2)
          IF(M+1>MAXDIM)THEN
            MAXDIM2=INT(MAXDIM*1.1)
            ALLOCATE(TRMX2(MAXDIM2),STAT=KK1)
            IF(KK1==0)THEN
              TRMX2(1:M)=TRMX(1:M)
              DEALLOCATE(TRMX)
              MAXDIM=MAXDIM2
              ALLOCATE(TRMX(MAXDIM),STAT=KK2)
              IF(KK2==0)THEN
                TRMX(1:M)=TRMX2(1:M)
                DEALLOCATE(TRMX2)
              ENDIF
            ENDIF
            IF(KK1/=0.OR.KK2/=0)THEN
              WRITE(*,'(a/a/A)') &
                'The maximum value of Ne set in your data is too large', &
                ' that your computer has not enough RAM to run the program. ',&
                ' Try a smaller maximum value of Ne and re-run the program'
              STOP
            ENDIF
          endif      
          M=M+1             
          TRMX(M)=PTT(J)
        ENDDO
        IB(2,I)=M
        IB(3,I)=LB(1)
!$OMP END CRITICAL                     
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL       
      NR_OLD=NR
      END SUBROUTINE GTRANS
!
      FUNCTION RAND1(IDUM)
      INTEGER,PARAMETER:: SP=KIND(1.0)
      INTEGER IDUM,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL(SP):: RAND1,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211, &
        IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER IDUM2,J,K,IV(NTAB),IY
      SAVE IV,IY,IDUM2
      DATA IDUM2/123456789/,IV/NTAB*0/,IY/0/
      IF(IDUM<=0)THEN
        IDUM=MAX(-IDUM,1)
        IDUM2=IDUM
        DO J=NTAB+8,1,-1
          K=IDUM/IQ1
          IDUM=IA1*(IDUM-K*IQ1)-K*IR1
          IF(IDUM<0)IDUM=IDUM+IM1
          IF(J<=NTAB)IV(J)=IDUM
        ENDDO
        IY=IV(1)
      ENDIF
      K=IDUM/IQ1
      IDUM=IA1*(IDUM-K*IQ1)-K*IR1
      IF(IDUM<0)IDUM=IDUM+IM1
      K=IDUM2/IQ2
      IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2
      IF(IDUM2<0)IDUM2=IDUM2+IM2
      J=1+IY/NDIV
      IY=IV(J)-IDUM2
      IV(J)=IDUM
      IF(IY<1)IY=IY+IMM1
      RAND1=MIN(AM*IY,RNMX)
      END FUNCTION RAND1
!
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER,PARAMETER:: DP=SELECTED_REAL_KIND(2*PRECISION(1.))
      REAL(DP):: dy,x,y,xa(n),ya(n),TINY
      PARAMETER (NMAX=10,TINY=1.E-50_DP)
      REAL(DP):: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do i=1,n
        dift=abs(x-xa(i))
        if (dift<dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(abs(den)<=tiny)then
            write(*,*)'Failure in polint, program stopped.'
            STOP
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if (2*ns<n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo
      END SUBROUTINE polint

      FUNCTION gammln(xx)
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6),gammln,xx
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,  &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,   &
       -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      DO j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
      ENDDO
      gammln=tmp+log(stp*ser/x)
      END FUNCTION gammln          