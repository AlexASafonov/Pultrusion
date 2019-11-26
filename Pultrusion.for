C----      Used to save data  
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TIME(2)
C
      DOUBLE PRECISION :: k_min_al_1, k_max_T_1, k_max_T_2,
     1 k_Strength, k_smax(22), k_Str(22), k_smaxT(22) 

      common k_min_al_1, k_min_al_2, k_max_T_1, k_max_T_2, k_max_S_1, 
     1 k_max_S_2, k_Strength, k_smax, k_Str, k_smaxT   
      
      character*256 OUTDIR
      
        CALL GETOUTDIR( OUTDIR, LENOUTDIR )        
        OUTDIR = trim(OUTDIR) // trim('\Results.txt')      
      
      
      IF (LOP.EQ.0) THEN

        open(unit=105, file=OUTDIR)
        
        k_min_al_1=1.0
        k_max_T_1=0.0
        k_max_T_2=0.0
        k_Strength=1.0
     
      END IF

      IF (LOP.EQ.1) THEN

C----  Determination of the maximum temperature in profile    
      
      IF (k_max_T_1.GT.k_max_T_2) THEN
    
            k_max_T_2=k_max_T_1
            
      END IF   
      
      k_max_T_1=0.0
 
     
C----  Strength test: 1-satisfied, 0-not satisfied,

       IF (TIME(2).GT.10.0) THEN

             DO i=1,22
       
               IF (k_Str(i).LT.k_smax(i)) THEN  
            
                 k_Strength=0.0
            
               END IF    
       
             END DO

       END IF
      
      END IF
      
C----  Save data to Results.txt      
      
      IF (LOP.EQ.3) THEN

C----   Results of Strength test
        write(105,'(F12.2)') k_Strength
        
C----   Minimum degree of cure in cross-section by the end of the process
        write(105,'(F12.6)') k_min_al_1

C----   Maximum temperature in the profile         
        write(105,'(F12.6)') k_max_T_2
        
C----   Maximum temperature in cross-section by the end of the process        
        write(105,'(F12.6)') k_max_T_1   
                
C----   Distribution of maximum transversal stresses over the cross-section
        DO i=1,22
        
          write(105,'(F12.2)') k_smaxT(i)  
        
        ENDDO        
        
C----   Time when the calculation was completed
        write(105,'(F12.6)') TIME(2) 
       
        close(105)
        
      END IF

      RETURN
      END  


C----  Setting temperature boundary conditions
      SUBROUTINE FILM(H,SINK,TEMP,KSTEP,KINC,TIME,NOEL,NPT,
     1 COORDS,JLTYP,FIELD,NFIELD,SNAME,NODE,AREA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION H(2),TIME(2),COORDS(3), FIELD(NFIELD)

      CHARACTER*80 SNAME
      DOUBLE PRECISION :: z
      DOUBLE PRECISION, PARAMETER :: Hmax=5000
      DOUBLE PRECISION, PARAMETER :: Hair=10.0
C     Vel - pulling speed
      DOUBLE PRECISION, PARAMETER :: Vel=0.000833333
C     kx1 - entrance of Zone 1    
      DOUBLE PRECISION, PARAMETER :: kx1=0.3
C     kx2 - end of Zone 1 
      DOUBLE PRECISION, PARAMETER :: kx2=0.5
C     kx3 - entrance of Zone 2      
      DOUBLE PRECISION, PARAMETER :: kx3=0.6
C     kx4 - end of Zone 2      
      DOUBLE PRECISION, PARAMETER :: kx4=0.8
C     kx5 - end of Die      
      DOUBLE PRECISION, PARAMETER :: kx5=1.0
C     kT1 - Temp at the entrance of Die       
      DOUBLE PRECISION, PARAMETER :: kT0=75.0
C     kT1 - Temp of Zone 1 
      DOUBLE PRECISION, PARAMETER :: kT1=150.0
C     kT2 - Temp of Zone 2      
      DOUBLE PRECISION, PARAMETER :: kT2=190.0
C     kT3 - Temp at the end of Die
      DOUBLE PRECISION, PARAMETER :: kT3=160.0
      
      z=TIME(2)*Vel
      H(2)=0
      
      if (z.LE.0) then
         H(1)=0
         goto 111
      endif  
      if ((z.GT.0).and.(z.LE.kx1)) then
         H(1)=Hmax
         SINK=(kT1-kT0)*z/kx1+kT0
         goto 111
      endif  
      if ((z.GT.kx1).and.(z.LE.kx2)) then
         H(1)=Hmax
         SINK=kT1
         goto 111
      endif  
      if ((z.GT.kx2).and.(z.LE.kx3)) then
         H(1)=Hmax
         SINK=kT1*(z-kx3)/(kx2-kx3)+kT2*(z-kx2)/(kx3-kx2)
         goto 111
      endif  
      if ((z.GT.kx3).and.(z.LE.kx4)) then
         H(1)=Hmax
         SINK=kT2
         goto 111
      endif  
      if ((z.GT.kx4).and.(z.LE.kx5)) then
         H(1)=Hmax
         SINK=kT2*(z-kx5)/(kx4-kx5)+kT3*(z-kx4)/(kx5-kx4)
         goto 111
      endif  
      
      H(1)=Hair
      SINK=25
111   RETURN
      END
      
C---- Define field variable (degree of cure) 
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
     
	  FIELD(1)=STATEV(1)
	  
      RETURN
      END
      
C---- Define a heat flux due to internal heat generation	  
      SUBROUTINE HETVAL(CMNAME,TEMP,TIME,DTIME,STATEV,FLUX,
     1 PREDEF,DPRED)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION TEMP(2),STATEV(*),PREDEF(*),TIME(2),FLUX(2),
     1 DPRED(*)
     
       FLUX(1)=STATEV(3)
       
      RETURN
      END
      
C---- Define incremental thermal and chemical shrinkage strains
      SUBROUTINE UEXPAN(EXPAN,DEXPANDT,TEMP,TIME,DTIME,PREDEF,
     1 DPRED,STATEV,CMNAME,NSTATV,NOEL)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION EXPAN(*),DEXPANDT(*),TEMP(2),TIME(2),PREDEF(*),
     1 DPRED(*),STATEV(NSTATV)
      EXPAN(1)=-STATEV(5)+STATEV(8)
      EXPAN(2)=-STATEV(6)+STATEV(9)
      EXPAN(3)=-STATEV(7)+STATEV(10)
      DEXPANDT(1)=STATEV(11)
      DEXPANDT(2)=STATEV(12)
      DEXPANDT(3)=STATEV(13)
      RETURN
      END
      
C---- Define the mechanical constitutive behavior of a material
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 DSTRESS(NTENS)
      DOUBLE PRECISION :: S(3,3), DD(3,3), mtmp(3), D(6,6),
     1 Em, Ef, num, nuf, Gf, Gm, Vf, Vm, X, dX_dt, dX_dt1, X1, dX,
     2 CTEf, CTEm, kf, km, kt, tmp, tmp1, tmp2, tmp3, tmp4,
     3 tmp7, tmp8, E1, E2, nu12, nu23, G12, G23, Emix, CTE1, CTE2,
     4 CSC1, CSC2, EPSCm, dEC1, dEC2, dET1, dET2, TTg, kmm
      integer :: arm_dir, ks_number(22)
      DOUBLE PRECISION, parameter :: volStrinkage=0.07
 
      DOUBLE PRECISION :: k_min_al_1, k_max_T_1, k_max_T_2,
     1 k_Strength, k_smax(22), k_Str(22), k_smaxT(22) 

      common k_min_al_1, k_min_al_2, k_max_T_1, k_max_T_2, k_max_S_1, 
     1 k_max_S_2, k_Strength, k_smax, k_Str, k_smaxT   
             
C     	1 - X (degree of cure)
C	  2 - dX_dt
C	  3 - q
C       4 - VR - chemical volumetric strain of thermoreactive matrix
C       5 - dEc1 - chemical strain increment dir 1
C       6 - dEc2 - chemical strain increment dir 2
C       7 - dEc3 - chemical strain increment dir 3
C       8 - dEt1 - thermal strain increment dir 1
C       9 - dEt2 - thermal strain increment dir 2
C      10 - dEt3 - thermal strain increment dir 3
C      11 - CTE1 - thermal coefficient dir 1
C      12 - CTE2 - thermal coefficient dir 2
C      13 - CTE3 - thermal coefficient dir 3
C
C
C PARAMS:
C      1 - rho*Htot*(1-Vf)
C      2 - Ef
C	 3 - nuf
C      4 - CTEf
C      5 - CTEm
C      6 - Vf
C      7 - arm dir (1,2,3)

       interface
       function inverse(MM) result(N)
        DOUBLE PRECISION, dimension(3,3), intent(in):: MM
        DOUBLE PRECISION, dimension(3,3)            :: N
       end function inverse

       function STRENGTH(TEMP, CURE)
        DOUBLE PRECISION, intent(in) :: TEMP, CURE
        DOUBLE PRECISION             :: STRENGTH
       end function STRENGTH

       function MATRIX_YOUNG(TEMP, CURE)
        DOUBLE PRECISION, intent(in) :: TEMP, CURE
        DOUBLE PRECISION             :: MATRIX_YOUNG
       end function MATRIX_YOUNG
       
       function MATRIX_BULK(TEMP, CURE)
        DOUBLE PRECISION, intent(in) :: TEMP, CURE
        DOUBLE PRECISION             :: MATRIX_BULK
       end function MATRIX_BULK

       function Tg(CURE)
        DOUBLE PRECISION, intent(in) :: CURE
        DOUBLE PRECISION             :: Tg
       end function Tg

       function CURE_RATE(CURE, TEMP)
        DOUBLE PRECISION, intent(in) :: CURE, TEMP
        DOUBLE PRECISION             :: CURE_RATE
       end function CURE_RATE

       function VOLUME_STRINKAGE(CURE)
        DOUBLE PRECISION, intent(in) :: CURE
        DOUBLE PRECISION             :: VOLUME_STRINKAGE
       end function VOLUME_STRINKAGE

       function CTEmSCALE(TEMP, CURE)
        DOUBLE PRECISION, intent(in) :: TEMP, CURE
        DOUBLE PRECISION             :: CTEmSCALE
       end function CTEmSCALE
      end interface 
 
 
C---- List of finite elements along radius 
      ks_number=(/ 220, 198, 176, 154, 132, 110, 88, 66, 44, 22, 263, 
     1 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274 /)            
     
C---- Find degree of cure    
            
      X = STATEV(1)

      dX_dt=CURE_RATE(X, TEMP+273)
      dX=0
      IF (dX_dt==0) GOTO 10
      dX = dX_dt*DTIME
      X=X+dX
      
      IF (X.GT.1.0) THEN
      
      X=1.0
      
      END IF
	
C---- calculating q for HETVAL
10    STATEV(3)=PROPS(1)*STATEV(2)

C---- calculating effective properties
      Ef=PROPS(2)
      nuf=PROPS(3)
      CTEf=PROPS(4)
      Vf=PROPS(6)
      arm_dir=PROPS(7)
      Gf=0.5*Ef/(1+nuf)
      Vm=1-Vf
      TTg=Tg(X)
      
      STATEV(21)=TTg
      
      Em=MATRIX_YOUNG(TTg-TEMP, X)
      STATEV(20)=Em
      
      kmm=MATRIX_BULK(TTg-TEMP, X)
      num=(3.0*kmm-Em)/(6.0*kmm)
      STATEV(23)=num
       
      Gm=0.5*Em/(1+num)
      kf=0.5*Ef/(1-nuf-2*nuf*nuf)
      km=0.5*Em/(1-num-2*num*num)
      tmp=(kf+Gm)*km+(kf-km)*Gm*Vf
      tmp2=Gm*Vm*Vf
      E1=Ef*Vf+Em*Vm+4*(num-nuf*nuf)*kf*km*tmp2/tmp
      tmp7=kf+Gm
      tmp8=(kf-km)*Vf
      kt=(tmp7*km+tmp8*Gm)/(tmp7-tmp8)
      tmp3=Gf+Gm
      tmp4=(Gf-Gm)*Vf
      G12=Gm*(tmp3+tmp4)/(tmp3-tmp4)
      G23=Gm*(km*tmp3+2*Gf*Gm+km*tmp4)/
     1    (km*tmp3+2*Gf*Gm-(km+2*Gm)*tmp4)
      nu12=nuf*Vf+num*Vm+(num-nuf)*(km-kf)*tmp2/tmp
      E2=1/(0.25/kt+0.25/G23+nu12*nu12/E1)
      nu23=0.5*(2*E1*kt-E1*E2-4*nu12*nu12*kt*E2)/E1/kt

C---- calculating thermal and chemical expansion
      CTEm=PROPS(5)*CTEmSCALE(TTg-TEMP, X)
      Emix=Ef*Vf+Em*Vm
      CTE1=(CTEf*Ef*Vf+CTEm*Em*Vm)/Emix
      CTE2=(CTEf+nuf*CTEf)*Vf+CTEm*(1+num)*Vm
     1 -(nuf*Vf+num*Vm)*CTE1

      CSC1=Em*Vm/Emix
      CSC2=(1+num)*Vm-(nuf*Vf+num*Vm)*CSC1
      
      EPSCm=volStrinkage/3.0*dX
      
      STATEV(27)=nu12
      STATEV(28)=nu23
      STATEV(29)=G12
      STATEV(30)=G23
      STATEV(31)=CTE1
      STATEV(32)=CTE2
      STATEV(33)=CSC1
      STATEV(34)=CSC2
    
      dEC1=CSC1*EPSCm
      dEC2=CSC2*EPSCm
      STATEV(4)=EPSCm

      dET1=CTE1*DTEMP
      dET2=CTE2*DTEMP
      STATEV(5)=dEC2
      STATEV(6)=dEC2
      STATEV(7)=dEC2
      STATEV(4+arm_dir)=dEC1

      STATEV(8)=dET2
      STATEV(9)=dET2
      STATEV(10)=dET2
      STATEV(7+arm_dir)=dET1

      STATEV(11)=CTE2
      STATEV(12)=CTE2
      STATEV(13)=CTE2
      STATEV(10+arm_dir)=CTE1
      STATEV(14)=STATEV(14)+STATEV(5)
      STATEV(15)=STATEV(15)+STATEV(6)
      STATEV(16)=STATEV(16)+STATEV(7)
      STATEV(17)=STATEV(17)+STATEV(8)
      STATEV(18)=STATEV(18)+STATEV(9)
      STATEV(19)=STATEV(19)+STATEV(10)

C---- calculating Jacobian
      S(1,1)=1.0/E1
      S(1,2)=-nu12/E1
      S(2,1)=S(1,2)
      S(2,2)=1.0/E2
      S(1,3)=-nu12/E1
      S(3,1)=S(1,3)
      S(3,3)=S(2,2)
      S(3,2)=-nu23/E2
      S(2,3)=S(3,2)
      D=0.0
      DD=inverse(S)
      IF (arm_dir.EQ.1) GOTO 5
      mtmp=DD(1,:)
      DD(1,:)=DD(arm_dir,:)
      DD(arm_dir,:)=mtmp
      mtmp=DD(:,1)
      DD(:,1)=DD(:,arm_dir)
      DD(:,arm_dir)=mtmp
5     D(1:3,1:3)=DD
      D(4,4)=G12
      D(5,5)=G12
      D(6,6)=G12
      D(7-arm_dir, 7-arm_dir)=G23
      DDSDDE(1:NTENS,1:NTENS)=D(1:NTENS,1:NTENS)

C---- calculating stresses
      DO I=1,NDI
       DSTRESS(I)=0
       DO J=1,NDI
        DSTRESS(I)=DSTRESS(I)+DDSDDE(I,J)*DSTRAN(J)
       ENDDO
      ENDDO
      DO I=1,NSHR
        DSTRESS(I+3)=DDSDDE(I+3,I+3)*DSTRAN(I+3)
      ENDDO

C---- calculating Energy
       EN=0
       DO I=1,NTENS
        EN=EN+(STRESS(I)+0.5*DSTRESS(I))*DSTRAN(I)
       ENDDO
       SSE=EN

C---- applying data
      STATEV(1)=X
      STATEV(2)=dX_dt
      DO I=1,NTENS
       STRESS(I)=STRESS(I)+DSTRESS(I)
      ENDDO
      IF (dX.GE.0.1) THEN
       PNEWDT=0.5
      ENDIF
      IF (dX.LT.0.005) THEN
       PNEWDT=2
      ENDIF
 
C---- Find minimum degree of cure in cross-section     
      IF (NOEL.EQ.1) THEN
      
      k_min_al_1=STATEV(1)
      
      END IF  
      
      IF (STATEV(1).LT.k_min_al_1) THEN
      
      k_min_al_1=STATEV(1)
      
      END IF            

C---- Find maximum temperature in cross-section
 
      IF (NOEL.EQ.1) THEN
      
      k_max_T_1=TEMP
      
      END IF 
      
      IF (TEMP.GT.k_max_T_1) THEN
      
      k_max_T_1=TEMP
      
      END IF  
      
C----  Store the transversal stresses distribution
      
      DO i=1,22
      
            IF (NOEL.EQ.ks_number(i)) THEN      
      
                STATEV(35+i)=STRESS(1)
      
            END IF      
      
      ENDDO 
      
 
       DO i=1,22
      
            IF (STATEV(35+i).GT.STATEV(57+i)) THEN      
      
                STATEV(57+i)=STATEV(35+i)
      
            END IF      
      
      ENDDO 
      
 
       DO i=1,22
       
            IF (NOEL.EQ.ks_number(i)) THEN
        
               k_smax(i)=STATEV(57+i) 
               
                 IF (TEMP.GE.101.5) THEN
                    
                    k_smaxT(i)=STATEV(57+i) 
               
                 END IF
               
               k_Str(i)=STRENGTH(TEMP, STATEV(1))
             
            END IF    
        
       ENDDO 
       
       
C----  Store the strength distribution      
 
      STATEV(25)=STRENGTH(TEMP, STATEV(1)) 

      RETURN
      END
 
 
C---- Function for calculation of matrix strength   
        DOUBLE PRECISION FUNCTION STRENGTH(TEMP, CURE)
        
        DOUBLE PRECISION, intent(in) :: TEMP, CURE
       DOUBLE PRECISION, PARAMETER :: S20=51.14e6
       DOUBLE PRECISION, PARAMETER :: S40=44.87e6
       DOUBLE PRECISION, PARAMETER :: S60=38.82e6
       DOUBLE PRECISION, PARAMETER :: S80=30.03e6
       DOUBLE PRECISION, PARAMETER :: S100=1.03e6

        STRENGTH=0
        
        IF (TEMP.LT.20) STRENGTH=S20
        IF ((TEMP.GE.20).AND.(TEMP.LT.40)) THEN
             STRENGTH=S20-(S20-S40)*(TEMP-20)/20
        ENDIF
        IF ((TEMP.GE.40).AND.(TEMP.LT.60)) THEN
             STRENGTH=S40-(S40-S60)*(TEMP-40)/20
        ENDIF
        IF ((TEMP.GE.60).AND.(TEMP.LT.80)) THEN
             STRENGTH=S60-(S60-S80)*(TEMP-60)/20
        ENDIF
        IF ((TEMP.GE.80).AND.(TEMP.LT.101.5)) THEN
             STRENGTH=S80-(S80-S100)*(TEMP-80)/21.5
        ENDIF
        
         IF (TEMP.GE.101.5) THEN
             STRENGTH=S100
        ENDIF
        
        END FUNCTION STRENGTH      

C---- Function for calculation of curing degree rate       
       DOUBLE PRECISION FUNCTION CURE_RATE(CURE, TEMP)
       
        DOUBLE PRECISION, INTENT(IN):: CURE, TEMP
        DOUBLE PRECISION            :: RT	  
        DOUBLE PRECISION, PARAMETER :: K0=31.4
        DOUBLE PRECISION, PARAMETER :: E=127.2e3
        DOUBLE PRECISION, PARAMETER :: n=1.8
        RT=8.31434*TEMP
        CURE_RATE=dexp(K0-E/RT)*(1-CURE)**n
       END FUNCTION CURE_RATE
       
C---- Function for calculation of matrix Young modulus
       DOUBLE PRECISION FUNCTION MATRIX_YOUNG(TEMP, CURE)

       DOUBLE PRECISION, INTENT(IN):: TEMP, CURE
       DOUBLE PRECISION, PARAMETER :: E101=3.769e7
       DOUBLE PRECISION, PARAMETER :: E80=2.321e9
       DOUBLE PRECISION, PARAMETER :: E60=3.025e9
       DOUBLE PRECISION, PARAMETER :: E40=3.277e9
       DOUBLE PRECISION, PARAMETER :: E20=3.475e9              
       DOUBLE PRECISION, PARAMETER :: T101=0.0
       DOUBLE PRECISION, PARAMETER :: T80=21.5
       DOUBLE PRECISION, PARAMETER :: T60=41.5
       DOUBLE PRECISION, PARAMETER :: T40=61.5
       DOUBLE PRECISION, PARAMETER :: T20=81.5  
       DOUBLE PRECISION, PARAMETER :: C1=0.6
       
       MATRIX_YOUNG=E101
       
       IF (CURE.LE.C1) GOTO 200
       
       IF ((TEMP.GT.T101).AND.(TEMP.LE.T80)) THEN
       MATRIX_YOUNG=E101+(TEMP-Tc1)/(T80-T101)*(E80-E101)
       GOTO 200
       ENDIF
  
       IF ((TEMP.GT.T80).AND.(TEMP.LE.T60)) THEN
       MATRIX_YOUNG=E80+(TEMP-Tc1)/(T60-T80)*(E60-E80)
       GOTO 200
       ENDIF
 
       IF ((TEMP.GT.T60).AND.(TEMP.LE.T40)) THEN
       MATRIX_YOUNG=E60+(TEMP-Tc1)/(T40-T60)*(E40-E60)
       GOTO 200
       ENDIF 
       
       IF ((TEMP.GT.T40).AND.(TEMP.LE.T20)) THEN
       MATRIX_YOUNG=E40+(TEMP-Tc1)/(T20-T40)*(E20-E40)
       GOTO 200
       ENDIF 
             
       IF (TEMP.GE.T20) MATRIX_YOUNG=E20
       
200    END FUNCTION MATRIX_YOUNG

C---- Function for calculation of matrix Bulk modulus
       DOUBLE PRECISION FUNCTION MATRIX_BULK(TEMP, CURE)
 
       DOUBLE PRECISION, INTENT(IN):: TEMP, CURE
       DOUBLE PRECISION, PARAMETER :: E0=1.931e9
       DOUBLE PRECISION, PARAMETER :: E8=4.826e9
       DOUBLE PRECISION, PARAMETER :: Tc1=0.0
       DOUBLE PRECISION, PARAMETER :: Tc2=21.5
       DOUBLE PRECISION, PARAMETER :: C1=0.6
       MATRIX_BULK=E0
       IF (CURE.LE.C1) GOTO 210
       IF ((TEMP.GT.Tc1).AND.(TEMP.LT.Tc2)) THEN
       MATRIX_BULK=E0+(TEMP-Tc1)/(Tc2-Tc1)*(E8-E0)
       GOTO 210
       ENDIF
       IF (TEMP.GE.Tc2) MATRIX_BULK=E8
210    END FUNCTION MATRIX_BULK

C---- Function for calculation of Tg
       DOUBLE PRECISION FUNCTION Tg(CURE)

       DOUBLE PRECISION, INTENT(IN):: CURE
       DOUBLE PRECISION, PARAMETER :: Tg0=-41.0
       DOUBLE PRECISION, PARAMETER :: Tg8=101.5
       DOUBLE PRECISION, PARAMETER :: lambda=0.44
       Tg=Tg0+lambda*CURE*(Tg8-Tg0)/(1-(1-lambda)*CURE)
       END FUNCTION Tg
                                                                       
C---- Function for CTEm scaling     
       DOUBLE PRECISION FUNCTION CTEmSCALE(TEMP, CURE)
       
       DOUBLE PRECISION, INTENT(IN):: TEMP, CURE
       DOUBLE PRECISION, PARAMETER :: dT=1.0
       DOUBLE PRECISION, PARAMETER :: C1=0.6      
       CTEmSCALE=2.5
       IF (CURE.LE.C1) GOTO 222      
       CTEmSCALE=2.5
       IF ((TEMP.GT.-dT).AND.(TEMP.LT.dT)) THEN
       CTEmSCALE=2.5-0.75*(TEMP+dT)/dT
       GOTO 222
       ENDIF
       IF (TEMP.GT.dT) CTEmSCALE=1
222    END FUNCTION CTEmSCALE    
              
C---- Function for inverting 3x3 matrix 
      function inverse(MM) result(N)
      
      INCLUDE 'ABA_PARAM.INC'
       DOUBLE PRECISION, dimension(3,3), intent(in) :: MM
       DOUBLE PRECISION, dimension(3,3)             :: N, M
       DOUBLE PRECISION                             :: tmp
       integer :: i,j
       N = 0.0
       do i=1,3
        N(i,i) = 1.0
       enddo
       M=MM
       do i=1,3
        tmp=M(i,i)
        M(i,:)=M(i,:)/tmp
        N(i,:)=N(i,:)/tmp
        do j=i+1,3
         tmp=M(j,i)
         M(j,:)=M(j,:)-M(i,:)*tmp
         N(j,:)=N(j,:)-N(i,:)*tmp
        enddo
       enddo
       do i=3,2,-1
        do j=i-1,1,-1
         tmp=M(j,i)
         M(j,:)=M(j,:)-M(i,:)*tmp
         N(j,:)=N(j,:)-N(i,:)*tmp
        enddo
       enddo
       return
      end function inverse
