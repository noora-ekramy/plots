! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                               +
! +   Strain control CP deformation modelling                     +
! +                                                               +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      include 'umat.f'
      include 'other_code.f'
c
      implicit none
      integer,parameter:: ndi=3
      integer,parameter:: nshr=3   
      integer,parameter:: ntens=6  
      integer,parameter:: nstatv=200
      integer,parameter:: nprops=20 
      integer,parameter:: layer=1
      integer,parameter:: kspt=1
      integer,parameter:: kstep=1
      integer,parameter:: noel=1
      integer,parameter:: npt=1
      integer,parameter:: kinc=1
c---------------------------------------------------------------------------
      character*80 cmname          !user defined material name
      real(8) drpldt               !jacobian drpl_dt
      real(8) dtime                !time increment dt
      real(8) temp                 !temperature at t0
      real(8) dtemp                !increment of temperature.
      real(8) celent               !characteristic element length
      real(8) sse                  !specific elastic strain energy
      real(8) spd                  !specific plastic dissipation
      real(8) scd                  !specific creep dissipation
      real(8) rpl                  !volumetric heat generation per unit time
      real(8) pnewdt               !dt_next/dt_now
c---------------------------------------------------------------------------
      real(8) stress(ntens)        !stress tensor
      real(8) ddsdde(ntens,ntens)  !jacobian ds_de
      real(8) ddsddt(ntens)        !jacobian ds_dt
      real(8) drplde(ntens)        !jacobian drpl_de
      real(8) stran (ntens)        !strains at t0
      real(8) dstran(ntens)        !strain increments
c---------------------------------------------------------------------------
      real(8) statev(nstatv)       !state variables
      real(8) props (nprops)       !material constants 
      real(8) propsX(nprops)       !material constants 
      real(8) dfgrd0(3,3)          !deformation gradient at t0
      real(8) dfgrd1(3,3)          !deformation gradient at t0+dt
      real(8) drot  (3,3)          !rotation increment matrix
      real(8) coords(3)            !coordinates of this point
      real(8) time  (2)            !1:step time; 2:total time, At t0
      real(8) predef(1)            !predefined field variables at t0
      real(8) dpred (1)            !incr of predefined field vrbs
c---------------------------------------------------------------------------
      integer lop,lrestart
      integer ising,icut,istep,ipct                                     
      integer icol,iit_gnd,ix
      integer i,j,k,ii,jj,is,icomp,ip,ipx
      real(8) x0,x1,x2,x3,dzeit,rand
      real(8) phi1,phi,phi2,Qm(3,3),XI33(3,3)
      real(8) Lg(3,3),IFg0(3,3),Fg0(3,3),Fg1(3,3)
      real(8) Em(3,3),Sm(3,3),Ev(6),Sv(6)
c---------------------------------------------------------------------------
      real(8) mx33_0(3,3),mx33_1(3,3),mx33_2(3,3),mx66_1(6,6)
      real(8) vx6_1(6),vx6_2(6)
      real(8) av1(3),av2(3),av3(3),cv(3),a,c
      real(8) vct4(4),vct3(3)
      real(8) strain_max,strain_det0,strain_det,strainABS_crt
      real(8) strain,strainABS
c
      integer,parameter:: ngr=1
c
      integer igr,npoint,npointO
      real(8) Fg0_RVE(ngr,3,3)
      real(8) Fe0_RVE(ngr,3,3)
      real(8) Fp0_RVE(ngr,3,3)
      real(8) Fg1_RVE(ngr,3,3)
      real(8) Fe1_RVE(ngr,3,3)
      real(8) Fp1_RVE(ngr,3,3)
      real(8) statev0_RVE(ngr,nstatv)
      real(8) statev1_RVE(ngr,nstatv)
      real(8) props_RVE(ngr,nprops) 
      real(8) sv_RVE(ngr,6) 
      real(8) ev_RVE(ngr,6)
      real(8) eRate,timQu
c
      integer ip1,ip2

      character(5) cnum
c
      real(8) evAvg(6),svAvg(6)
      real(8) dsdeAvg(6,6)
      real(8) IdsdeAvg(6,6)
      real(8) dvAvg0(6),dvAvg(6),detDv(6)
      real(8) dmAvg0(3,3),dmAvg(3,3),vfAvg 
      integer IB1(9),IB2(9)
      integer iloop,nloop,ix1,ix2,ix3,Iload
c-----------------------------------
      dzeit = 1.d-2      ! max time step
      nloop = 10         ! stress relaxation number 
      eRate = 5.d-3      ! uniaxial strain rate
      timQu = 19.53      ! 1/4 of period
      strain_det0=0.0001

c-----------------------------------
      XI33(1,:)=[1,0,0]
      XI33(2,:)=[0,1,0]
      XI33(3,:)=[0,0,1]
      IB1=[1,2,3, 1,1,2, 2,3,3]
      IB2=[1,2,3, 2,3,3, 1,1,2]
c
      open(300,file='parameters.dat')
      read(300,*) propsX(1)  !Young's module [MPa]   
      read(300,*) propsX(2)  !Poisson ratio 
      read(300,*) propsX(3)  !Initial flow stress R0 [MPa] 
      read(300,*) propsX(4)  !Isotropic hardening coefficient Q_inf [MPa] 
      read(300,*) propsX(5)  !Isotropic hardening parameter b 
      read(300,*) propsX(6)  !Isotropic recovery coefficient D_inf [MPa] 
      read(300,*) propsX(7)  !Isotropic recovery parameter a 
      read(300,*) propsX(8)  !total number of backstress nBK 
      do i=1,int(propsX(8))
         read(300,*) propsX(8+(i-1)*2+1)  !#i backstress parameter c(i)      
         read(300,*) propsX(8+(i-1)*2+2)  !#i backstress parameter c(i)      
      enddo
      close(300)

      do igr=1,ngr
         props_RVE(igr,:)=propsX
      enddo

c-----------------------------------
      do Iload=3,3
         write(cnum,'(I1.1)') Iload
         open(100,file='se'//trim(cnum)//'.dat')
      
         dvAvg0=0.0 
         if(Iload==1)then
            dvAvg0(1)= eRate/2 
            dvAvg0(2)=-eRate/2 
            dvAvg0(3)=-eRate/2
         elseif(Iload==2)then 
            dvAvg0(1)=-eRate/2 
            dvAvg0(2)= eRate 
            dvAvg0(3)=-eRate/2
         elseif(Iload==3)then 
            dvAvg0(1)=-eRate/2 
            dvAvg0(2)=-eRate/2 
            dvAvg0(3)= eRate
         endif
         call icams_conv6to33(dvAvg0,ib1,ib2,dmAvg0)
         do igr=1,ngr
            Fg0_RVE(igr,:,:)=XI33
         enddo
      
c-----------------------------------
         dtime=0.d0; time=0.d0
         lop=0; lrestart=0
         call uexternaldb(lop,0,time,dtime,kstep,kinc)

         evAvg=0 
         dtime=dzeit
         strain=0
         strainABS=0
         strainABS_crt=0
         npoint=0
         npointO=0

         do istep=1,100000
103         continue
            lop=1; lrestart=0
            call uexternaldb(lop,lrestart,time,dtime,kstep,kinc)

            if(    time(2)<timQu*1)then
               dvAvg      = +dvAvg0
               strain_det = +strain_det0
            elseif(time(2)<timQu*2)then
               dvAvg      = -dvAvg0
               strain_det = -strain_det0

            elseif(time(2)<timQu*3)then
               dvAvg      = -dvAvg0
               strain_det = -strain_det0
            elseif(time(2)<timQu*4)then
               dvAvg      = +dvAvg0
               strain_det = +strain_det0

            elseif(time(2)<timQu*5)then
               dvAvg      = +dvAvg0
               strain_det = +strain_det0
            elseif(time(2)<timQu*6)then
               dvAvg      = -dvAvg0
               strain_det = -strain_det0

            elseif(time(2)<timQu*7)then
               dvAvg      = -dvAvg0
               strain_det = -strain_det0
            elseif(time(2)<timQu*8)then
               dvAvg      = +dvAvg0
               strain_det = +strain_det0

            elseif(time(2)<timQu*9)then
               dvAvg      = +dvAvg0
               strain_det = +strain_det0
            elseif(time(2)<timQu*10)then
               dvAvg      = -dvAvg0
               strain_det = -strain_det0

            elseif(time(2)<timQu*11)then
               dvAvg      = -dvAvg0
               strain_det = -strain_det0
            elseif(time(2)<timQu*12)then
               dvAvg      = +dvAvg0
               strain_det = +strain_det0

c            elseif(time(2)<timQu*4.3)then
c               dvAvg = +dvAvg0
            else
c               stop
               goto 104
            endif         

            svAvg=0
            IdsdeAvg=0 
            do iloop=1,nloop 

               vx6_1=svAvg
               vx6_1(Iload)=0

               detDv=matmul(IdsdeAvg,vx6_1)/dtime

               x1=sum(dabs(detDv))/6 
               x2=sum(dabs(dvAvg0))/6 
               if( iloop>1 .and. x1<1.d-3*x2 ) goto 112

               dvAvg = dvAvg - detDv
               call icams_conv6to33(dvAvg,ib1,ib2,dmAvg)

               vfAvg=0
               svAvg=0
               dsdeAvg=0
               do igr=1,ngr
                  coords(1:3)=[0.,0.,0.]
                  props=props_RVE(igr,:)
                  statev=statev0_RVE(igr,:) 
                  Fg0=Fg0_RVE(igr,:,:) 
                  Fg1=Fg0+matmul(dmAvg,Fg0)*dtime
                  dfgrd0=Fg0 
                  dfgrd1=Fg1 


                  call inverseMatrix(Fg0,IFg0,Ising)
                  if(Ising/=0)then
                     pnewdt=0.5
                     return
                  endif
                  mx33_1=matmul(Fg1,IFg0)
                  mx33_2=(mx33_1+transpose(mx33_1))/2-XI33
                  do i=1,6
                     dstran(i)=mx33_2(IB1(i),IB2(i))
                  enddo
                  dstran(4:6)=dstran(4:6)*2


                  call umat(stress,statev,ddsdde,
     &                      sse,spd,scd,rpl,
     &                      ddsddt,drplde,drpldt,
     &                      stran,dstran,
     &                      time,dtime,
     &                      temp,dtemp,
     &                      predef,dpred,
     &                      cmname,
     &                      ndi,nshr,ntens,nstatv,
     &                      props,nprops,coords,drot,
     &                      pnewdt,
     &                      celent,
     &                      dfgrd0,dfgrd1,
     &                      noel,npt,layer,kspt,kstep,kinc)

                 ! call pm(dfgrd1,3,3)
                 ! call pv(stress,6)
                 ! call pm(ddsdde,6,6)
                 ! print*,pnewdt
                 ! read*

                  if(pnewdt<1.d0)then
                     dtime=dtime*pnewdt
                     goto 103
                  endif

                  vfAvg=vfAvg+sum(statev(76+19:76+24))/0.129
                  svAvg=svAvg+stress
                  dsdeAvg=dsdeAvg+ddsdde

                  Fg1_RVE(igr,:,:)=Fg1
                  sv_RVE(igr,:)=stress
                  statev1_RVE(igr,:)=statev

                  
               enddo
               vfAvg=vfAvg/ngr
               svAvg=svAvg/ngr
               dsdeAvg=dsdeAvg/ngr
               call gaussj(dsdeAvg,6,IdsdeAvg,Ising)
               if(Ising/=0)then
                  print*,'dsdeAvg is non-invertible'
                  stop               
               endif

            enddo

112         continue
            Fg0_RVE=Fg1_RVE
            statev0_RVE=statev1_RVE

            print*,'statev0 is updated'

            evAvg=evAvg+dvAvg*dtime
            npoint=npoint+1

            strain=strain+dvAvg(Iload)*dtime
            strainABS=strainABS+dabs(dvAvg(Iload))*dtime

            write(*,'(I5,50e12.3)') npoint,strain*100,svAvg,
     &       statev0_RVE(1,3)
            if(strainABS>strainABS_crt)then
               npointO=npointO+1
               strainABS_crt = strainABS_crt + dabs(strain_det)
               write(100,'(I5,50e15.5)') npointO,strain*100,svAvg,
     &                                   statev0_RVE(1,1:11)


            endif

            lop=2; lrestart=0
            call uexternaldb(lop,lrestart,time,dtime,kstep,kinc)
            time(2)=time(2)+dtime

            if(pnewdt>=1.d0)then
               x1=dtime*1.2
               dtime=min(x1,dzeit)
            endif

         enddo
104      continue
         close(100)
      enddo
      end

c----------------------------------------------      
      subroutine xit
      return
      end  





