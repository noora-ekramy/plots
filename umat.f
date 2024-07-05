c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +  Abaqus Umat, Abaqus Umat, Abaqus Umat, Abaqus Umat, Abaqus Umat        +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine umat(stress,statev,ddsdde,
     &                sse,spd,scd,rpl,
     &                ddsddt,drplde,drpldt,stran,dstran, 
     &                time,dtime,temp,dtemp,
     &                predef,dpred,cmname,ndi,nshr,ntens, 
     &                nstatv,props,nprops,
     &                coords,drot,pnewdt,celent, 
     &                dfgrd0,dfgrd1,noel,npt,layer,
     &                kspt,kstep,kinc)
c                       
c      include "aba_param.inc"
c
      implicit none
      character*8 cmname
      integer noel,npt,layer,kspt,kstep,kinc
      integer ntens,nstatv,nprops,ndi,nshr
      real(8) stress(ntens), statev(nstatv), 
     &        ddsdde(ntens, ntens),
     &        ddsddt(ntens), drplde(ntens), 
     &        stran(ntens), dstran(ntens),
     &        predef(1), dpred(1), props(nprops), 
     &        coords(3), drot(3,3),
     &        dfgrd0(3,3), dfgrd1(3,3)
      real(8) dtime,time(2),sse,spd,scd,rpl
      real(8) pnewdt,celent,temp,dtemp,drpldt
C
      integer Ising
      integer IB1(9),IB2(9)
      integer i, j, nBK, iLoop, converged, ind_alpha
      integer  nLoop
      real(8) XI33(3,3)
      real(8) elastic_modulus, Q_inf, b, D_inf, a,
     &        shear_modulus, bulk_modulus, poission_ratio, 
     &        mu2, lame_first
      real(8) R0, ep_eq, ep0_eq, a_temp,
     &        RQ, RD, hard_iso_total, a_dot_n,
     &        plastic_mult, p_mult_numer, p_mult_denom, 
     &        yield_function, hardI, hardK,
     &        stress_relM, strain_trace, alpha_trace, e_k,
     &        ID2_out_ID2, n_out_n, stress_hydro, sigma_vm,
     &        Lam, n33_check, alpha_out_n, 
     &        beta, theta_1, theta_2, theta_3,srn2
      real(8), dimension(:, :), ALLOCATABLE :: alpha_k, alpha0_k
      real(8), dimension(:),    ALLOCATABLE :: c_bk, gm_bk
      real(8), dimension(6, 6) :: ID4, c_mat
      real(8) strain_tens(6), strain_plastic(6),
     &        stress_relN(6), alpha0(6), alpha(6), 
     &        strain_trial(6), stress_rel(6),
     &        stress_dev(6), ID2(6), stress_tri(6), 
     &        check(6), dstran_tens(6), alpha_diff(6),
     &        alpha_upd(6), dpe(6)
      real(8) TOL, SQRT23
      real(8) propsX(100)
      real(8) Fg0(3,3),IFg0(3,3),Fg(3,3)
      real(8) mx33_1(3,3),mx33_2(3,3),mx33_3(3,3)
      real(8) x1,x2,x3
C-------------------------------------------------
      pnewdt=1
      IB1=[1,2,3, 1,1,2 ,2,3,3]
      IB2=[1,2,3, 2,3,3 ,1,1,2]
      XI33(1,:)=[1.d0, 0.d0, 0.d0]
      XI33(2,:)=[0.d0, 1.d0, 0.d0]
      XI33(3,:)=[0.d0, 0.d0, 1.d0]
      ID2=0.d0 
      ID4=0.d0 
      DO i=1,3
         ID2(i  )=1.d0
         ID4(i,i)=1.d0  
      ENDDO
      DO i=4,6
         ID2(i  )=0.d0
         ID4(i,i)=0.5d0
      END DO
      SQRT23=dSQRT(2.0D0/3.0D0)
      TOL=1.0D-10
      nLoop=1000
C------------------------------------------
C     model parameters parameters
C------------------------------------------
      elastic_modulus = props(1)
      poission_ratio  = props(2)
      R0              = props(3)
      Q_inf           = props(4)
      b               = props(5)
      D_inf           = props(6)
      a               = props(7)
      nBK             = int(props(8))
      ALLOCATE(c_bk(nBK))
      ALLOCATE(gm_bk(nBK))
      ALLOCATE(alpha_k(nBK, 6))
      ALLOCATE(alpha0_k(nBK, 6))
      do i=1,nBK
         c_bk(i)      = props(8+(i-1)*2+1)
         gm_bk(i)     = props(8+(i-1)*2+2)
      enddo
      shear_modulus=elastic_modulus/(2*(1+poission_ratio))
      bulk_modulus=elastic_modulus/(3*(1-2*poission_ratio))
      mu2=2*shear_modulus

C------------------------------------------
C     internal variables: plasticity and backstress
c       1            Equivalent plastic strain
c       2 - 7        Plastic strain (ep11,ep22,ep33,2*ep12,2*ep13,2*ep23)
c       8 - 13       Backstress-1 (s11,s22,s33,s12,s13,s23)
c      14 - 19       Backstress-2 (s11,s22,s33,s12,s13,s23)
c      ...
C------------------------------------------
      ep_eq = statev(1)
      ep0_eq = ep_eq
      CALL ROTSIG(statev(2), drot, strain_plastic, 2, ndi, nshr)
      alpha=0.d0
      DO i = 1, nBK
         CALL ROTSIG(statev(8+(i-1)*6),drot,alpha_k(i,:),1,ndi,nshr)
         alpha = alpha + alpha_k(i, :)
      END DO

C ----------------------------------------------------------------------C
C
C     Elastic trial step
C
C ----------------------------------------------------------------------C
      strain_tens = stran + dstran
      DO j=1,6
      DO i=1,6
         c_mat(i,j)=ID2(i)*ID2(j)*bulk_modulus 
     &             +mu2*(ID4(i,j)-ID2(i)*ID2(j)/3.d0)
      ENDDO
      ENDDO
      stress_tri = stress + MATMUL(c_mat, dstran)
      x1=SUM(stress_tri(1:3))/3.d0
      stress_rel = stress_tri - x1*[1,1,1,0,0,0] - alpha
      x1=sum(stress_rel(1:3)**2)+2*sum(stress_rel(1:3)**2) 
      stress_relM = dsqrt(x1)                     !--magnitude of relative stress
      stress_relN = stress_rel/(TOL+stress_relM)  !--normal of relative stress
     
C ----------------------------------------------------------------------C
C
C     First check of yield condition
C
C ----------------------------------------------------------------------C
      RQ = Q_inf * (1-dEXP(-b*ep_eq))
      RD = D_inf * (1-dEXP(-a*ep_eq))
      yield_function = stress_relM - SQRT23*(R0+RQ-RD)
      IF (yield_function > TOL) THEN
         converged = 0
      ELSE
         converged = 1
      ENDIF

C ----------------------------------------------------------------------C
C
C     Radial return mapping if plastic loading
C
C ----------------------------------------------------------------------C
      iLoop = 0
      plastic_mult=0.d0
      alpha0=alpha
      alpha0_k=alpha_k
      DO WHILE (converged==0 .AND. iLoop<nLoop)
        iLoop = iLoop + 1
C
        !--Calculate the isotropic hardening parameters
        RQ = Q_inf*(1-dEXP(-b*ep_eq))
        RD = D_inf*(1-dEXP(-a*ep_eq))
        hardI = b*(Q_inf-RQ)-a*(D_inf-RD)

        !--Calculate the kinematic hardening parameters
        hardK=0.d0
        DO i = 1, nBK
           x1=dexp(-gm_bk(i)*(ep_eq-ep0_eq))
           x2=dot_product(stress_relN(1:3), alpha_k(i,1:3))
           x3=dot_product(stress_relN(4:6), alpha_k(i,4:6))*2
           hardK=hardK+x1*(c_bk(i)-dSQRT(3.d0/2)*gm_bk(i)*(x2+x3))
        END DO

        alpha=0.d0
        DO i = 1, nBK
          x1=dEXP(-gm_bk(i)*(ep_eq-ep0_eq))
          x2=(1-x1)
          alpha_k(i,:)=x1*alpha0_k(i,:) 
     &                +x2*SQRT23*stress_relN*c_bk(i)/gm_bk(i)  
          alpha=alpha+alpha_k(i, :)
        END DO
        x1=dot_product(alpha(1:3)-alpha(1:3), stress_relN(1:3))
        x2=dot_product(alpha(4:6)-alpha(4:6), stress_relN(4:6))*2
        p_mult_numer = stress_relM -
     1     ((x1+x2) + SQRT23*(R0+RQ-RD) + mu2*plastic_mult)
        p_mult_denom = -mu2*(1+(hardK+hardI)/(3*shear_modulus))
C
        ! Update variables
        plastic_mult = plastic_mult - p_mult_numer / p_mult_denom
        ep_eq = ep0_eq + SQRT23 * plastic_mult
C
        IF (ABS(p_mult_numer) < TOL) THEN
          converged = 1
        END IF
      END DO

C ----------------------------------------------------------------------C
C
C     calculate stress, plastic strain and backstress
C
C ----------------------------------------------------------------------C
      IF (iLoop==0) THEN  ! Elastic loading
         stress = stress_tri
      ELSE  ! Plastic loading
         dpe=plastic_mult*stress_relN
         dpe(4:6)=dpe(4:6)*2
         strain_plastic = strain_plastic + dpe
         stress = stress_tri - MATMUL(c_mat, dpe)
      ENDIF

C ----------------------------------------------------------------------C
C
C     calculate stiffness
C
C ----------------------------------------------------------------------C
      IF (iLoop==0) THEN  ! Elastic loading
         ddsdde= c_mat
      ELSE  ! Plastic loading
         alpha_diff=alpha-alpha0
         x1=dot_product(stress_relN(1:3),alpha_diff(1:3))
         x2=dot_product(stress_relN(4:6),alpha_diff(4:6))*2
         beta    = 1.d0+(hardK+hardI)/(3*shear_modulus)
         theta_1 = 1.d0 - mu2 * plastic_mult / stress_relM
         theta_3 = 1.d0/(beta*stress_relM)
         theta_2 = 1.d0/beta +(x1+x2)*theta_3-(1.d0-theta_1)
         DO j=1,6
         DO i=1,6
            ddsdde(i, j) = bulk_modulus * ID2(i)*ID2(j)
     1      +mu2*theta_1*(ID4(i,j)-1.d0/3*ID2(i)*ID2(j))  !--elastic influence
     2      -mu2*theta_2*stress_relN(i)*stress_relN(j)    !--plastic influence
     3      +mu2*theta_3*alpha_diff(i)*stress_relN(j)     !--back stress influence
         ENDDO
         ENDDO
      !   ddsdde = (TRANSPOSE(ddsdde) + ddsdde)/2
      ENDIF
C ----------------------------------------------------------------------C
C
      ! Update the state variables
C
C ----------------------------------------------------------------------C
      statev(1)=ep_eq
      DO i=1,6
        statev(i+1)=strain_plastic(i)
      ENDDO
      DO i=1,nBK
      DO j=1,6
         statev(7+(i-1)*6+j)=alpha_k(i,j)
      ENDDO
      ENDDO
C ----------------------------------------------------------------------C
C
C     Reduce time increment if did not converge
C
C ----------------------------------------------------------------------C
      IF (iLoop>=nLoop) PNEWDT=0.25
c
      RETURN
      END      


c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   Uexternaldb...                                                        +
c +       lop = 0 beginning of analysis                                     +
c +             1 start of increment                                        +
c +             2 end of increment                                          +
c +             3 end of analysis                                           +
c +             4 beginning of restart                                      +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
      implicit none
      integer lop,lrestart,kstep,kinc
      real(8) time(2),dtime
      real(8) x1,x2,x3
      real(8) v3_1(3),v3_2(3),v3_3(3)
      integer i,j,ix,ip
      integer ix1,ix2,ix3,ix4,ix5
      integer iv5_1(5)
      integer np_x,np_y,np_z
      character(len=255) :: Fpathx
      character(len=255) :: PathFilex
      integer get_thread_id,getnumthreads,I_thread,N_thread
c
      I_thread=0
      N_thread=1
c      I_thread=get_thread_id()
c      N_thread=getnumthreads()
c      
      if(lop==0)then   !ini
         print '("First call",2I3,I5,10E15.5)' 
c     &         ,I_thread,N_thread,kinc,dtime,time
      endif
      if(lop==1)then
         print '("dt-start",  2I3,I5,10E15.5)' 
c     &         ,I_thread,N_thread,kinc,dtime,time
      endif
      if(lop==2)then
         print '("dt-end",    2I3,I5,10E15.5)' 
c     &         ,I_thread,N_thread,kinc,dtime,time
      endif
     
      return
      end

      
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c +                                                                         +
c +   useful subroutines                                                    +
c +                                                                         +
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine ROTSIG(v6_1, Rm, v6_2, ix, ndi, nshr)
      implicit none
      integer i,j,k,ndi,nshr,ix
      real(8) v6_1(6),v6_2(6),Rm(3,3)
      real(8) x1,x2,x3
      v6_2 = v6_1
      return
      end      
c-----------------------------------------------------      
      subroutine inverseMatrix(A, A_inv, Ising)
      implicit none
      real(8) A(3,3), A_inv(3,3)
      real(8) det
      integer i, j, Ising
      Ising=0
      det = +A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))  
     &      -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))  
     &      +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      if (dabs(det)>1.d-10) then
         A_inv(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / det
         A_inv(1,2) = (A(1,3)*A(3,2) - A(1,2)*A(3,3)) / det
         A_inv(1,3) = (A(1,2)*A(2,3) - A(1,3)*A(2,2)) / det
         A_inv(2,1) = (A(2,3)*A(3,1) - A(2,1)*A(3,3)) / det
         A_inv(2,2) = (A(1,1)*A(3,3) - A(1,3)*A(3,1)) / det
         A_inv(2,3) = (A(1,3)*A(2,1) - A(1,1)*A(2,3)) / det
         A_inv(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1)) / det
         A_inv(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2)) / det
         A_inv(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / det
      else
         A_inv = 0.0d0
         Ising=1
      end if
      return
      end

      

