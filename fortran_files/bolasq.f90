!		Q-balls -  novembro de 2021 - Felipe Console

!		Resolve a equação para Q-balls com potential 
!		U =  1.1 * phi^2 - 2* phi^4 + phi^6
!		usando a frequencia do campo escalar, omega, como variável dinâmica.
!		Input ao programa colsys:
! 			Equações:
!				ff1 = 0 (omega é constante)
!				ff2 = ... (equaçao para o campo escalar phi'' = ...)
!
!			z1 = omega
!			z2 = omega'
!			z3 = phi
!			z4 = phi'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		condições de contorno (mstar = 4)
!			omega'(0) = 0 
!			phi(0) = phi_0 ( = cf)	
!			phi'(0) = 0 
!			phi(infty) = 0
!	
	  implicit real*8 (a-h,o-z)              
      
	  common /parb/ xl,ipr
      double precision zeta(4),fspace(1000000) ,dx,fator(120000)
      double precision tol(4),z(4),dmval(2),ff(2)
      double precision xx(1000),omega,cf,cargan,tomega
	  double precision phi2(1002),delta(1002),h(1002),soma,xl
	  integer :: i,ip,ipp, np	  
      integer ispace(74000),m(3),ipar(11),ltol(4),npas(4)  
      double precision inte(5001) 
	  common /pard/ cf, omega, oomega,tomega
      external colsys                            
      external fsub,dfsub,gsub,dgsub,solutn
	  

		open(unit = 3, file='solu.new')
		open(unit = 4, file='solu.old')
		open(unit = 5, file='solu.dat')
		! open(unit = 75, file='sopri.dat')
		!open(unit = 105, file='xl_415_tol_d_6_bc2_168pts.dat',position='append')
		open(unit = 95, file='parametros.dat')
		open(unit = 155, file='data.dat',position='append')
!
!			PARAMETERS
!
      pi=4.d0*datan(1.d0)                
      NPRI=1
	!   cf = 1.d0
	  read(95,*)cf
	  ! write(*,*)cf
	  ipr = 415
	
	                            
!     INPUT FOR SUBROUTINE COLSYS    
!            
		  NCOMP = 2 ! número de equações
			
						! ordem da equação i .... lembre que mstar = m(1) + m(2) + ... + m(ncomp)
						! neste caso, mstar = m(1) = 2
			M(1) = 2	! equação para omega
			M(2) = 2    ! equação para phi
			
		  ALEFT = 0.d0 			! valor da fronteira esquerda
		  ARIGHT = 1200.d0		! valor da fronteira direita
		  xl = aright
		  xh = 0.d0
		  
			ZETA(1) = ALEFT         ! onde a primeira equação de contorno deve ser satisfeita 
			ZETA(2) = ALEFT         ! onde a segunda equação de contorno deve ser satisfeita  
			ZETA(3) = ALEFT
			ZETA(4) = ARIGHT
			
		  tolin = 1.d-6
			
		  TOL(1)=tolin			! tolerâncias
		  TOL(2)=tolin    		! tolerâncias
		  TOL(3)=tolin    		! tolerâncias
		  TOL(4)=tolin    		! tolerâncias
			  
		  LTOL(1)=1      		!!!!!!!!!!
		  LTOL(2)=2     		!!!!!!!!!!
		  LTOL(3)=3     		!!!!!!!!!!
		  LTOL(4)=4     		!!!!!!!!!!
		  
		  IPAR(1)=1                       ! = 0  se as equações são lineares e  = 1 se são não-lineares 
		  IPAR(2)=3                       !número de pontos de colocação por subintervalo
		  IPAR(3)=ipr+1                     !número de subintervalos na mesh inicial
		  IPAR(4)=4                       ! devemos ter """" 0 .lt. ntol .le. mstar """"
		  IPAR(5)=1000000                  !dimensão de fspace 
		  IPAR(6)=740000                   !dimensão de ispace
		  IPAR(7)=1                       !controle de output
		  
		  IPAR(8)=0		                     !se = 0 implica o colsys gerar uma rede inicial uniforme
										  !se = 1 implica que o usuário tem que fornecer uma rede
		  
		  
		  IPAR(9)=1                       ! = 0 (não há chute inicial p/ solução)
		  IPAR(10)=0                      ! = 0 assumimos que o problema é regular
		  IPAR(11)=0      				  ! número de pontos fixos além dos extremos do intervalo aleft e aright.

		
			ipp = ipr + 1
			xdelta = 0.d0


		!!!	Lire la meche	

				do i = 1, ipr+1
				read(3,*) x,xxx,yyy
				fspace(i) = x + xdelta
				end do
				rewind 3
				
				call intpol
				
			
      write(6,189)
  189 format(" now colsys start ")
      CALL COLSYS(NCOMP,M,ALEFT,ARIGHT,ZETA,IPAR,LTOL,TOL,FIXPNT,ISPACE,FSPACE,IFLAG,FSUB,DFSUB,GSUB,DGSUB,SOLUTN)                  
                                                          

	WRITE(6,6)IFLAG                                    
    6 FORMAT('IFLAG= ',I2) 
		
		
  	if(iflag.ne.1) stop 
	goto 300
      
 300   continue
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
		do i = 1,ipr+1 
		x = fspace(i)
		call appsln(x,z,fspace,ispace)
		omega = z(1)
		!write(55,27)x, omega, 0.0
		end do
	27 format(e16.10,2x,e16.10,2x,e16.10)
		
		do i = 1,ipr+1 	
		x = fspace(i)
		call appsln(x,z,fspace,ispace)
		phi = z(3)
		phip = z(4)
		!write(55,27)x, phi, phip
		end do
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
		
		DO i = 1,ipr+1    
		x = fspace(i)
		CALL APPSLN(x,Z,FSPACE,ISPACE)
		omega = z(1)
		phi = z(3)
		phip = z(4)
		xpotential = 1.1d0*phi**2 -2.d0*phi**4 + phi**6
		xenergydensity = omega**2*phi**2 + phip**2 + xpotential  
		write(5,26) x,omega,phi,phip,x**2*xenergydensity, xenergydensity
		end do
				
 26		FORMAT(E16.10,2X,E16.10,2X,E16.10,2x,e16.10,2x,e16.10,2x,e16.10) 		
 

		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		
		! cálculo da  carga de noether		
		! regra de simpson
		
		somaQ =0.d0
		somaE =0.d0
		np = 12000
		xnp = dfloat(np+1)
		dx = ( xl - xh ) / xnp
		
		fator(1) = 1.d0/3.d0
		fator(np+1) = 1.d0/3.d0
		fator(np) = 4.d0/3.d0
		
		do i =2,np-2,2
		ip =i+1
		fator(i) = 4.d0/3.d0
		fator(ip) = 2.d0/3.d0
		end do

		x=xh
		
		do j=1,np+1	
		call appsln(x,z,fspace,ispace)
		
		omega = z(1)
		phi = z(3)
		phip = z(4)
		
		somaQ = somaQ + x**2*phi**2*dx*fator(j)
		xpotential = 1.1d0*phi**2 -2.d0*phi**4 + phi**6
		somaE = somaE + x**2*(omega**2*phi**2 + phip**2  + xpotential)*dx*fator(j)
		
		x = x+dx
		
		end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
		
		x = xh
		call appsln(xm,z,fspace,ispace)
		omega=z(1)
		phi = z(3)
		cargaN = 8.d0*pi*omega*somaQ
		xenergy = 4.d0*pi*somaE
		write(155,36)cf,omega,cargaN,xenergy
36		format(e16.10,2x,e16.10,2x,e16.10,2x,e16.10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

		
      STOP                                
      END  

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCE

      SUBROUTINE FSUB(X,Z,FF)            
      IMPLICIT real*8 (A-H,O-Z) 
	  
      DOUBLE PRECISION Z(4),FF(2)    
	  
		z1=z(1) !omega
		z2=z(2)	!omega'
		z3=z(3)	!phi
		z4=z(4)	!phip

		ff1 = 0.d0
		ff(1) = ff1

		
		ff2 = (1.1d0 - z1**2.d0)*z3 - 4.d0*z3**3.d0 + 3.d0*z3**5.d0 - (2.d0*z4)/x
		ff(2) = ff2
		
		return	

      END      
	  
		SUBROUTINE DFSUB(X,Z,DFF)         
		
		IMPLICIT real*8 (A-H,O-Z)                               
		DOUBLE PRECISION Z(4),DFF(2,4),omega,CF,OOMEGA,XL
		common /pard/ cf, omega, oomega,tomega
		
		z1=z(1)
		z2=z(2)
		z3=z(3)
		z4=z(4)
		
		ff1_1 = 0.d0
	 	ff1_2 = 0.d0
		ff1_3 = 0.d0
		ff1_4 = 0.d0
		
		ff2_1 = -2.d0*z1*z3
		ff2_2 = 0.d0
		ff2_3 = 1.1d0 - z1**2.d0 - 12.d0*z3**2.d0 + 15.d0*z3**4.d0
		ff2_4 = -2.d0/x
			
		dff(1,1) = ff1_1 
		dff(1,2) = ff1_2 
		dff(1,3) = ff1_3 
		dff(1,4) = ff1_4 
	
		dff(2,1) = ff2_1 
		dff(2,2) = ff2_2 
		dff(2,3) = ff2_3 
		dff(2,4) = ff2_4 
		
		  
		return

      END 

                                                                   
      SUBROUTINE GSUB(I,Z,G)                                            
      IMPLICIT real*8 (A-H,O-Z)                                     
	  real*8 Z(4), G,cf,XL,OOMEGA
	  common /pard/ cf, omega, oomega,tomega
			
		z1=z(1)
		
		oomega = dsqrt(1.1d0-z(1)**2)
		
	 
		  GOTO(1,2,3,4),I
		

1		 G = Z(3)-cf
		 RETURN                                                                                                                    
2		 G = Z(2)-0.d0
		 RETURN                                                                                                                       
3		 G = z(4)-0.d0
		 RETURN
! 4	    G = xl*Z(4) + z(3)*xl*oomega+z(3)	
4		G = z(3)-0.d0
		 RETURN                                       
      END
!c				
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c

      SUBROUTINE DGSUB(I,Z,DG)                                          
      IMPLICIT real*8 (A-H,O-Z)                               
      double precision z(4), DG(4),xl,omega,oomega	  
	  common /pard/ cf, omega, oomega,tomega
			
			
			z1 = z(1)
			z2 = z(2)
			z3 = z(3)
			z4 = z(4)
	        
			oomega = dsqrt(1.1d0-z(1)**2)      
				  
		  GOTO(1,2,3,4),I
1		 DG(1) =  0.d0
		 DG(2) =  0.d0	
		 DG(3) =  1.d0	
		 DG(4) =  0.d0	
		  RETURN                                                            
2		 DG(1) =  0.d0
		 DG(2) =  1.d0	
		 DG(3) =  0.d0	
		 DG(4) =  0.d0	
		  RETURN                                                            
3		 DG(1) =  0.d0
		 DG(2) =  0.d0	
		 DG(3) =  0.d0
		 DG(4) =  1.d0		
		 RETURN
! 4		DG(1) =  -z1*z3*xl/oomega	
4		DG(1) = 0.d0	
		 DG(2) =  0.d0
		 DG(3) = 1.d0 
		 ! DG(3) =  xl*oomega+1.d0
		 DG(4) = 0.d0
		 ! DG(4) =  xl		
		 RETURN       
		 
	
		END
! !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      
      SUBROUTINE INTPOL                                                 
      IMPLICIT real*8 (A-H,O-Z)     
      DIMENSION Z(4), DMVAL(2), X1(1002)                            
         COMMON /FIELDS/ XI(1002),V1(1002),V2(1002),V3(1002),V4(1002),V5(1002),V6(1002),&
		 W1(1002),W2(1002),W3(1002),W4(1002)                         
      COMMON /PARB/ xl,IPR
      COMMON / PARc / delx(1002)
    
      IMAX = IPR + 2
      imm = imax-1
      inn = imax-2
	  
		DO 19 I=1,imm
		READ(3,*)XI(I),V1(I),V2(I)
   19 	WRITE(4,25)XI(I),V1(I),V2(I)
	
		DO 29 I=1,imm
		READ(3,*)XI(I),V3(I),V4(I)
   29 	WRITE(4,25)XI(I),V3(I),V4(i)
		

   25 FORMAT(E16.10,2X,E16.10,2X,E16.10) 

		REWIND 3

      XMAX=XI(imm)                                                       
!C                                                                       
      DXI = XI(2)- XI(1)                                                   
      Delx(1) = XI(2) - XI(1)                                                   
      DXIo = XI(imm) - XI(inn)                                            
      Delx(imm)=XI(imm)-XI(inn)                                                   
      DXI2=2.D0*DXI                                                     
!C                                                                       
      W1(1)=(V2(2)-V2(1))/DXI                                                                              
      W2(1)=(V4(2)-V4(1))/DXI                                                                                 
                                         
	  W1(imm)=(V2(imm)-V2(inn))/DXIo                                                                                                                                                       
      W2(imm)=(V4(imm)-V4(inn))/DXIo                                                                                
                                      
      DO 40 I=2,inn                                                    
      IP=I+1                                                            
      IM=I-1                                                            
      dxi2=xi(ip)-xi(im)
      delx(i)=xi(ip)-xi(im)
      W1(I)=(V2(IP)-V2(IM))/DXI2                                                                               
      W2(I)=(V4(IP)-V4(IM))/DXI2                                        
      
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
                                       


		! SUBROUTINE SOLUTN(X,Z,DMVAL)                                      
		! IMPLICIT real*8 (A-H,O-Z)     
		! DIMENSION Z(4), DMVAL(2), X1(1002)                            
		! COMMON /FIELDS/ XI(1002),V1(1002),V2(1002),V3(1002),V4(1002),W1(1002),W2(1002)                         
		! COMMON /PARB/ xl,IPR
		! COMMON / PARc / delx(1002)
    
      ! xmax=xl
      ! imax=ipr+1

      ! IF(X.EQ.0.D0) GO TO 1                                                
        ! IF(X.GE.XMAX) GO TO 2                                                
        ! imm=imax-1
        ! do 3 i=1,imm
      ! 3 if(x.ge.xi(i).and.x.lt.xi(i+1)) go to 5
    
	  ! 5 idx=i
        ! IDXP=IDX+1                                                        

        ! del=(x-xi(idx))/(xi(idxp)-xi(idx))
	  
        ! Z(1)=V1(IDX)+(V1(IDXP)-V1(IDX))*del
        ! Z(2)=V2(IDX)+(V2(IDXP)-V2(IDX))*del
        ! Z(3)=V3(IDX)+(V3(IDXP)-V3(IDX))*del
        ! Z(4)=V4(IDX)+(V4(IDXP)-V4(IDX))*del

        ! DMVAL(1)=W1(IDX)+(W1(IDXP)-W1(IDX))*del
        ! DMVAL(2)=W2(IDX)+(W2(IDXP)-W2(IDX))*del

        ! RETURN                                                            

      ! 1 CONTINUE
        ! Z(1)=V1(1)                                                        
        ! Z(2)=V2(1)                                                        
        ! Z(3)=V3(1)                                                        
        ! Z(4)=V4(1)                                                        

        ! DMVAL(1)=W1(1)                                                    
        ! DMVAL(2)=W2(1)                                                    
	                                              
        ! RETURN                                                            
      ! 2 CONTINUE
        ! Z(1)=V1(IMAX)                                                       
        ! Z(2)=V2(IMAX)
        ! Z(3)=V3(IMAX)                                                       
        ! Z(4)=V4(IMAX)                                                       
	                                                 
        ! DMVAL(1)=W1(IMAX)                                                   
        ! DMVAL(2)=W2(IMAX)                                                   

        ! RETURN                                                            
        ! END  

									   
       SUBROUTINE SOLUTN(X,Z,DMVAL) !!!!!solução inicial phi gaussiano                                      
       IMPLICIT real*8 (A-H,O-Z)     
       DIMENSION Z(4),DMVAL(2)                                           
	   COMMON /FIELDS/ XI(1002),V1(1002),V2(1002),V3(1002),V4(1002),V5(1002),V6(1002),&
		 W1(1002),W2(1002),W3(1002),W4(1002)                         
       COMMON / PARc / delx(1002)
	   DIMENSION X1(1002),X2(1002),X3(1002),X4(1002)
	   COMMON /PARB/ xl,IPR
	   common /pard/ cf, omega, oomega,tomega
	   double precision xl,omega, oomega, cf
	  
		
		 tomega = 0.99d0
		 phinot = 1.0d0
		 aux = dsqrt(1.1d0-tomega**2)
		 phiseed = phinot*exp(-aux*x**2)
		
		
		 z(1) = tomega
		 z(2) = 0.d0
		 z(3) = phiseed
		 z(4) = -2.d0*aux*x*phiseed
		 
		
		 dmval(1) = 0.d0
		 dmval(2) = -2.d0*aux*phiseed+4.d0*x**2*aux**2*phiseed
		
		
	   RETURN                                                            
       END                                                               

                                                                                                                
	   