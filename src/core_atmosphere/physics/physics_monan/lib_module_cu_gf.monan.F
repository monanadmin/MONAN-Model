   !---------------------------------------------------------------------------------------------------
   subroutine set_Tq_pertub (use_excess,its,ite,itf,xlandi,ztexec,zqexec,cum_ztexec,cum_zqexec)
         implicit none
         integer, intent(in) :: its, ite, itf, use_excess
         real, dimension(:), intent(in)  :: ztexec,zqexec,xlandi
         real, dimension(:), intent(out) :: cum_ztexec,cum_zqexec
         integer :: i
         !
         !-- set minimum/max for excess of T and Q
         if(use_excess == 0) then
            cum_ztexec(:)= 0.
            cum_zqexec(:)= 0.
         elseif (use_excess == 1) then
            cum_ztexec(:)= ztexec(:)
            cum_zqexec(:)= zqexec(:)
         elseif (use_excess == 2) then
            do i=its,itf
               cum_zqexec(i)=min(5.e-4, max(1.e-4,zqexec(i)))! kg kg^-1
               cum_ztexec(i)=min(0.5,   max(0.2  ,ztexec(i)))! Kelvin
            enddo
         else
            do i=its,itf
               if(xlandi(i) > 0.98) then ! ocean
                  cum_zqexec(i)=min(8.e-4, max(5.e-4,zqexec(i)))! kg kg^-1
                  cum_ztexec(i)=min(1.,    max(0.5  ,ztexec(i)))! Kelvin
               else                      ! land
                  cum_ztexec(i)= ztexec(i)
                  cum_zqexec(i)= zqexec(i)
               endif
            enddo
         endif

   end subroutine set_Tq_pertub
   !---------------------------------------------------------------------------------------------------
   subroutine get_capmax(cumulus,itf,ktf,its,ite,kts,kte,cap_max_inc &
                        ,cap_max_increment,cap_max,cap_maxs,MOIST_TRIGGER)
      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      integer  ,intent (in )                 :: MOIST_TRIGGER
      real     ,intent (out)                 :: cap_max_inc,cap_maxs
      real     ,intent (out), dimension(:)   :: cap_max_increment,cap_max

      !-local vars
      integer :: i,k
      !
      !--- maximum depth (mb) of capping inversion (larger cap = no convection)
      if( MOIST_TRIGGER==0) then
         if(cumulus == 'deep'   ) then
            cap_max_inc=20.
         endif 
         if(cumulus == 'mid'    ) then
            cap_max_inc=10.
         endif 
         if(cumulus == 'shallow') then
            cap_max_inc=25.
         endif
      else
         if(cumulus == 'deep'   ) then
            cap_max_inc=90.
         endif
         if(cumulus == 'mid'    ) then
            cap_max_inc=90.
         endif 
         if(cumulus == 'shallow') then
            cap_max_inc=10.
         endif 
      endif
      cap_max_increment(:)= cap_max_inc
      cap_max          (:)= cap_maxs

   end subroutine get_capmax
   !---------------------------------------------------------------------------------------------------
   subroutine get_lambdaU(cumulus,itf,ktf,its,ite,kts,kte,lambau_dp,lambau_dn &
                         ,lambau_deep,lambau_shdn,pgcon)

      implicit none
      character*(*), intent (in)             :: cumulus
      integer  ,intent (in )                 :: itf,ktf, its,ite, kts,kte
      real     ,intent (in )                 :: lambau_deep,lambau_shdn,pgcon
      real     ,intent (out), dimension(:)   :: lambau_dp, lambau_dn
      !-local vars
      integer :: i,k

      !--- lambda_U parameter for momentum transport
      !
      if(cumulus == 'deep'   ) then
         lambau_dp (:) = lambau_deep
         lambau_dn (:) = lambau_shdn
      endif
      if(cumulus == 'mid'    ) then
         lambau_dp (:) = lambau_shdn
         lambau_dn (:) = lambau_shdn
      endif
      if(cumulus == 'shallow') then
         lambau_dp (:) = lambau_shdn
         lambau_dn (:) = lambau_shdn
      endif

      if(pgcon .ne. 0.) then
         lambau_dp (:) = 0.
         lambau_dn (:) = 0.
      endif

   end subroutine get_lambdaU
   !------------------------------------------------------------------------------------
   subroutine get_interp(q_old,t_old,po_cup,q_new,t_new)
      implicit none
      real    ,intent (in   ) :: po_cup ! original
      real    ,intent (inout) :: q_old,t_old,q_new,t_new ! extrapolated

      !---locals
      real ::  zqp, zcond1, zcor, zqsat
      real ::  psp, pt , pq, ptare
      real ::  foealfcu, foeewmcu,foedemcu,foeldcpmcu
      
      !real, parameter :: &
      !  c_rd=287.06                             &
      ! ,c_rv=461.52                             &
      ! ,c_rcpd=1004.71                          &
      ! ,c_rtt=273.16                            &
      ! ,c_rhoh2o=1000.                          &
      ! ,c_rlvtt=2.5008e+6                       &
      ! ,c_rlstt=2.8345e+6                       &
      ! ,c_retv = c_rv/c_rd-1.0                  &
      ! ,c_rlmlt= c_rlstt-c_rlvtt                &
      ! ,c_rcpv=4.*c_rv                          &
      ! ,c_r2es=611.21*c_rd/c_rv                 &
      ! ,c_r3les=17.502                          &
      ! ,c_r3ies=22.587                          &
      ! ,c_r4les=32.19                           &
      ! ,c_r4ies=-0.7                            &
      ! ,c_r5les= c_r3les*(c_rtt-c_r4les)        &
      ! ,c_r5ies= c_r3ies*(c_rtt-c_r4ies)        &
      ! ,c_r5alvcp= c_r5les*c_rlvtt/c_rcpd       &
      ! ,c_r5alscp= c_r5ies*c_rlstt/c_rcpd       &
      ! ,c_ralvdcp= c_rlvtt/c_rcpd               &
      ! ,c_ralsdcp= c_rlstt/c_rcpd               &
      ! ,c_ralfdcp= c_rlmlt/c_rcpd               &
      ! ,c_rtwat= c_rtt                          &
      ! ,c_rtber= c_rtt-5.                       &
      ! ,c_rtbercu= c_rtt-5.0                    &
      ! ,c_rtice= c_rtt-23.                      &
      ! ,c_rticecu= c_rtt-23.                    &
      ! ,c_rtwat_rtice_r=1./(c_rtwat-c_rtice)    &
      ! ,c_rtwat_rticecu_r=1./(c_rtwat-c_rticecu)&
      ! ,c_rvtmp2= c_rcpv/c_rcpd-1.              &
      ! ,c_zqmax=0.5

      integer :: i

      pt  = t_old       ! k
      pq  = q_old       ! kg/kg
      psp = po_cup*100. ! hpa

      !-- for testing
      !              psp                   temp        q                     zcond1
      ! input    27940.0000000000        236.604976804749       3.220181796223121e-004
      ! output   27940.0000000000        236.361132108860       4.084506812610067e-004
      !  pt  = 236.604976804749      ! k
      !  pq  = 3.220181796223121e-004       ! kg/kg
      !  psp = 27940. ! hpa
      !----------------------
      !print*,"1",psp,pt,pq

      zqp   =1.0/psp
      do i=1,2
         ptare = pt

         foealfcu = min(1.0,((max(c_rticecu,min(c_rtwat,ptare))-c_rticecu)*c_rtwat_rticecu_r)**2)
         foeewmcu = c_r2es *(foealfcu *exp(c_r3les*(ptare-c_rtt)/(ptare-c_r4les))+&
                  (1.0-foealfcu)*exp(c_r3ies*(ptare-c_rtt)/(ptare-c_r4ies)))
         zqsat=foeewmcu*zqp

         !    if(1.0-retv  *zqsat == 0.) then
         !
         !      print*,"zqsat=",zqp,foeewmcu,q_old,t_old,po_cup,q_new,t_new
         !3.5491847e-02   46.36052      0.5000000       249.8219
         !  0.2817549      0.5000000       249.8219
         !      call flush(6)
         !      stop 3333
         !    endif

         zcor=1.0/(1.0-c_retv  *zqsat)
         zqsat=zqsat*zcor

         foedemcu =  foealfcu     *c_r5alvcp*(1.0/(ptare-c_r4les)**2)+&
                    (1.0-foealfcu)*c_r5alscp*(1.0/(ptare-c_r4ies)**2)


         zcond1=(pq-zqsat)/(1.0+zqsat*zcor*foedemcu)

         foeldcpmcu= foealfcu*c_ralvdcp+(1.0-foealfcu)*c_ralsdcp
         pt=pt+foeldcpmcu*zcond1
         pq=pq-zcond1
      enddo
      !-- FINAL --------------------------
      q_new =  PQ
      t_new =  PT
      !print*,"2",PSP,PT,PQ
      !print*,"E",100*(PT-236.361132108860)/236.361132108860,100*(PQ-4.084506812610067E-004)/4.084506812610067E-004
   end subroutine get_interp
   !------------------------------------------------------------------------------------
   subroutine get_wetbulb(jmin,qo_cup,t_cup,po_cup ,q_wetbulb,t_wetbulb)

      implicit none
      integer ,intent (in   ) :: jmin
      real    ,intent (in   ) :: qo_cup,t_cup,po_cup
      real    ,intent (inout) :: q_wetbulb,t_wetbulb

      !---locals
      real ::  zqp, zcond, zcond1, zcor, zqsat
      real :: psp, pt , pq
      real :: z3es,   z4es, z5alcp, zaldcp
      real :: ptare, evap
      real :: foedelta,foeewmcu,foealfcu,foedemcu,foeldcpmcu

      !-- for testing
      !              PSP                   TEMP        Q                     ZCOND1
      ! input   85090.0000000000        289.140030372766     1.105078557441815E-002
      ! output  85090.0000000000        287.230570412846     1.181792062536557E-002 -2.761256206705639E-005
      ! PT  = 289.140030372766
      ! PQ  = 1.105078557441815E-002
      ! PSP = 85090.
      !----------------------

      !-- environmental values
      PT  = t_cup       ! K
      PQ  = qo_cup      ! kg/kg
      psp = po_cup*100. ! hpa

      if (pt > c_rtt) then
         z3es=c_r3les
         z4es=c_r4les
         z5alcp=c_r5alvcp
         zaldcp=c_ralvdcp
      else
         z3es=c_r3ies
         z4es=c_r4ies
         z5alcp=c_r5alscp
         zaldcp=c_ralsdcp
      endif

      !--- get wet bulb thermo properties --------------------------
      ptare = pt
      zqp    =1.0/psp

      foealfcu = min(1.0,((max(c_rticecu,min(c_rtwat,ptare))-c_rticecu)*c_rtwat_rticecu_r)**2)
      foeewmcu = c_r2es *(foealfcu *exp(c_r3les*(ptare-c_rtt)/(ptare-c_r4les))+&
                     (1.0-foealfcu)*exp(c_r3ies*(ptare-c_rtt)/(ptare-c_r4ies)))
      zqsat=foeewmcu*zqp

      zqsat=min(p_max_qsat,zqsat)
      zcor=1.0/(1.0-c_retv  *zqsat)
      zqsat=zqsat*zcor

      foedemcu =  foealfcu *c_r5alvcp*(1.0/(ptare-c_r4les)**2)+&
             (1.0-foealfcu)*c_r5alscp*(1.0/(ptare-c_r4ies)**2)

      zcond=(pq-zqsat)/(1.0+zqsat*zcor*foedemcu)

      zcond=min(zcond,0.0)

      foeldcpmcu= foealfcu*c_ralvdcp+(1.0-foealfcu)*c_ralsdcp
      pt=pt+foeldcpmcu*zcond

      pq=pq-zcond

      !--update ptare
      ptare = pt

      foealfcu = min(1.0,((max(c_rticecu,min(c_rtwat,ptare))-c_rticecu)*c_rtwat_rticecu_r)**2)
      foeewmcu = c_r2es *(foealfcu *exp(c_r3les*(ptare-c_rtt)/(ptare-c_r4les))+&
                     (1.0-foealfcu)*exp(c_r3ies*(ptare-c_rtt)/(ptare-c_r4ies)))
      zqsat=foeewmcu*zqp

      zqsat=min(0.5,zqsat)
      zcor=1.0/(1.0-c_retv  *zqsat)
      zqsat=zqsat*zcor

      foedemcu =  foealfcu *c_r5alvcp*(1.0/(ptare-c_r4les)**2)+&
             (1.0-foealfcu)*c_r5alscp*(1.0/(ptare-c_r4ies)**2)
      zcond1=(pq-zqsat)/(1.0+zqsat*zcor*foedemcu)

      if(zcond == 0.0)zcond1=min(zcond1,0.0)
      foeldcpmcu= foealfcu*c_ralvdcp+(1.0-foealfcu)*c_ralsdcp
      pt=pt+foeldcpmcu*zcond1
      pq=pq-zcond1

      !-- set output --------------------------
      q_wetbulb =  pq
      t_wetbulb =  pt
      evap      = -ZCOND1 != q_wetbulb-qo_cup, source for water vapor
   end subroutine get_wetbulb
   !-----------------------------------------------------------------------------------------
   subroutine fct1d3 (ktop,n,dt,z,tracr,massflx,trflx_in,del_out)

      ! --- modify a 1-D array of tracer fluxes for the purpose of maintaining
      ! --- monotonicity (including positive-definiteness) in the tracer field
      ! --- during tracer transport.

      ! --- the underlying transport equation is   (d tracr/dt) = - (d trflx/dz)
      ! --- where  dz = |z(k+1)-z(k)| (k=1,...,n) and  trflx = massflx * tracr

      ! --- note: tracr is carried in grid cells while z and fluxes are carried on
      ! --- interfaces. interface variables at index k are at grid location k-1/2.
      ! --- sign convention: mass fluxes are considered positive in +k direction.

      ! --- massflx and trflx_in  must be provided independently to allow the
      ! --- algorithm to generate an auxiliary low-order (diffusive) tracer flux
      ! --- as a stepping stone toward the final product trflx_out.

      implicit none
      integer,intent(in) :: n,ktop                        ! number of grid cells
      real   ,intent(in) :: dt                            ! transport time step
      real   ,intent(in) :: z(n+0)                        ! location of cell interfaces
      real   ,intent(in) :: tracr(n)                      ! the transported variable
      real   ,intent(in) :: massflx  (n+0)                ! mass flux across interfaces
      real   ,intent(in) :: trflx_in (n+0)                ! original tracer flux
      real   ,intent(out):: del_out  (n+0)                ! modified tracr flux
      real               :: trflx_out(n+0)                ! modified tracr flux
      integer k,km1,kp1
      logical :: NaN, error=.false., vrbos=.false.
      real dtovdz(n),trmax(n),trmin(n),flx_lo(n+0),antifx(n+0),clipped(n+0),  &
         soln_hi(n),totlin(n),totlout(n),soln_lo(n),clipin(n),clipout(n),arg
      real,parameter :: epsil=1.e-22           ! prevent division by zero
      real,parameter :: damp=1.                ! damper of antidff flux (1=no damping)

      logical, parameter :: hi_order = .false.

      NaN(arg) = .not. (arg.ge.0. .or. arg.le.0.)        ! NaN detector
      soln_lo(:)=0.
      antifx (:)=0.
      clipout(:)=0.
      flx_lo (:)=0.

      do k=1,ktop
         dtovdz(k)=.01*dt/abs(z(k+1)-z(k))                ! time step / grid spacing
      !     if (z(k).eq.z(k+1)) error=.true.
      end do
      if (vrbos .or. error) print '(a/(8es10.3))','(fct1d) dtovdz =',dtovdz(1:ktop)

      do k=2,ktop
         if (massflx(k) > 0.) then
            flx_lo(k)=massflx(k)*tracr(k-1)              ! low-order flux, upstream
         else
            flx_lo(k)=massflx(k)*tracr(k)                ! low-order flux, upstream
         endif
         antifx(k)=trflx_in(k)-flx_lo(k)                ! antidiffusive flux
      end do
      flx_lo(  1)   =trflx_in(  1)
      flx_lo(ktop+1)=trflx_in(ktop+1)
      antifx(  1)   =0.
      antifx(ktop+1)=0.
      ! --- clip low-ord fluxes to make sure they don't violate positive-definiteness
      do k=1,ktop
         totlout(k)=max(0.,flx_lo(k+1))-min(0.,flx_lo(k  ))         ! total flux out
         clipout(k)=min(1.,tracr(k)/max(epsil,totlout(k))/ (1.0001*dtovdz(k)))
      end do

      do k=2,ktop
         if (massflx(k).ge.0.)  then
            flx_lo(k)=flx_lo(k)*clipout(k-1)
         else
            flx_lo(k)=flx_lo(k)*clipout(k)
         endif
      end do
      if (massflx(1)     .lt.0.) flx_lo(1)     =flx_lo(1)     *clipout(1)
      if (massflx(ktop+1).gt.0.) flx_lo(ktop+1)=flx_lo(ktop+1)*clipout(ktop)

      ! --- a positive-definite low-order (diffusive) solution can now be  constructed

      do k=1,ktop
         soln_lo  (k)=tracr(k)-(flx_lo(k+1)-flx_lo(k))*dtovdz(k)        ! low-ord solutn
         del_out  (k)=-c_grav*(flx_lo(k+1)-flx_lo(k))*dtovdz(k)/dt
      end do

      if(.not. hi_order) return

      soln_hi  (:)=0.
      clipin   (:)=0.
      trmin    (:)=0.
      trmax    (:)=0.
      clipped  (:)=0.
      trflx_out(:)=0.


      do k=1,ktop
         km1=max(1,k-1)
         kp1=min(n,k+1)
         trmax(k)=       max(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
            tracr  (km1),tracr  (k),tracr  (kp1))        ! upper bound
         trmin(k)=max(0.,min(soln_lo(km1),soln_lo(k),soln_lo(kp1),        &
            tracr  (km1),tracr  (k),tracr  (kp1)))       ! lower bound
      end do

      do k=1,ktop
         totlin (k)=max(0.,antifx(k  ))-min(0.,antifx(k+1))                ! total flux in
         totlout(k)=max(0.,antifx(k+1))-min(0.,antifx(k  ))                ! total flux out

         clipin (k)=min(damp,(trmax(k)-soln_lo(k))/max(epsil,totlin (k)) / (1.0001*dtovdz(k)))
         clipout(k)=min(damp,(soln_lo(k)-trmin(k))/max(epsil,totlout(k)) / (1.0001*dtovdz(k)))

         if (NaN(clipin (k))) print *,'(fct1d) error: clipin is NaN,  k=',k
         if (NaN(clipout(k))) print *,'(fct1d) error: clipout is NaN,  k=',k

         if (clipin(k).lt.0.) then
            print 100,'(fct1d) error: clipin < 0 at k =',k,                        &
               'clipin',clipin(k),'trmax',trmax(k),'soln_lo',soln_lo(k),        &
               'totlin',totlin(k),'dt/dz',dtovdz(k)
            error=.true.
         end if
         if (clipout(k).lt.0.) then
            print 100,'(fct1d) error: clipout < 0 at k =',k,                        &
               'clipout',clipout(k),'trmin',trmin(k),'soln_lo',soln_lo(k),        &
               'totlout',totlout(k),'dt/dz',dtovdz(k)
            error=.true.
         end if
100      format (a,i3/(4(a10,"=",es9.2)))
      end do

      do k=2,ktop
         if (antifx(k).gt.0.)  then
            clipped(k)=antifx(k)*min(clipout(k-1),clipin(k))
         else
            clipped(k)=antifx(k)*min(clipout(k),clipin(k-1))
         end if
         trflx_out(k)=flx_lo(k)+clipped(k)
         if (NaN(trflx_out(k)))  then
            print *,'(fct1d) error: trflx_out is NaN,  k=',k
            error=.true.
         end if
      end do

      trflx_out(     1)=trflx_in(     1)
      trflx_out(ktop+1)=trflx_in(ktop+1)
      do k=1,ktop
         soln_hi(k)=tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
         del_out(k) =     -c_grav*(trflx_out(k+1)-trflx_out(k))*dtovdz(k)/dt
        !write(32,*)'3',k,soln_lo(k),soln_hi(k)
      end do

      if (vrbos .or. error) then
         do k=2,ktop
            write(32,99)k,                   &
               'tracr(k)', tracr(k),            &
               'flx_in(k)', trflx_in(k),        &
               'flx_in(k+1)', trflx_in(k+1),    &
               'flx_lo(k)', flx_lo(k),          &
               'flx_lo(k+1)', flx_lo(k+1),      &
               'soln_lo(k)', soln_lo(k),        &
               'trmin(k)', trmin(k),            &
               'trmax(k)', trmax(k),            &
               'totlin(k)', totlin(k),          &
               'totlout(k)', totlout(k),        &
               'clipin(k-1)', clipin(k-1),      &
               'clipin(k)', clipin(k),          &
               'clipout(k-1)', clipout(k-1),    &
               'clipout(k)', clipout(k),        &
               'antifx(k)', antifx(k),          &
               'antifx(k+1)', antifx(k+1),      &
               'clipped(k)', clipped(k),        &
               'clipped(k+1)', clipped(k+1),    &
               'flx_out(k)', trflx_out(k),      &
               'flx_out(k+1)', trflx_out(k+1),  &
               'dt/dz(k)', dtovdz(k),           &
               'final', tracr(k)-(trflx_out(k+1)-trflx_out(k))*dtovdz(k)
99          format ('(trc1d)   k =',i4/(3(a13,'=',es13.6)))
         end do
         if (error) stop '(fct1d error)'
      end if
   end subroutine fct1d3
   !---------------------------------------------------------------------------------------------------
   subroutine tridiag (m,a,b,c,f)
      !-- this routine solves the problem: aa*f(k-1,t+1) + bb*f(k,t+1) + cc*f(k+1,t+1) = dd
      !-- an updated "f" at time t+1 is the output
      implicit none
      integer, intent(in) :: m
      real, dimension(m), intent(inout) :: a,b,c
      real, dimension(m), intent(inout) :: f
      !--locals
      real, dimension(m) :: q
      integer :: k
      real :: p

      c(m)=0.
      q(1)=-c(1)/b(1)
      f(1)= f(1)/b(1)
      do k=2,m
         p  = 1./( b(k)+a(k)*q(k-1) )
         q(k) = -c(k)*p
         f(k) = p*(f(k) - a(k)*f(k-1))
      enddo
      do k=m-1,1,-1
         f(k) = f(k) +q(k)*f(k+1)
      enddo
   end subroutine tridiag

  ! ---------------------------------------------------------------------------------------------------
   subroutine setGradsVar(i_in, k_in, nvar, f_in, name1, name2, name3)
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'setGradsVar' ! Subroutine Name
   
      !Variables (input, output, inout)
      integer, intent(in)    :: i_in, k_in
      
      real, intent(in) :: f_in
      
      character(len=*), intent(in) :: name1
      character(len=*), intent(in) :: name2
      character(len=*), intent(in) :: name3

      integer, intent(inout) :: nvar

      cupout(nvar)%varp(i_in, k_in) = f_in
      cupout(nvar)%varn(1) = name1
      cupout(nvar)%varn(2) = name2
      cupout(nvar)%varn(3) = name3
      nvar = nvar + 1
      if (nvar > p_nvar_grads) stop 'nvar>nvar_grads'

   end subroutine setGradsVar
   ! ------------------------------------------------------------------------------------
   subroutine wrtBinCtl(n, mzp, p2d, cumulus,ntimes)
      implicit none
      !Parameters:
      character(len=*), parameter :: procedureName = 'wrtBinCtl' ! Subroutine Name

      real, parameter :: p_undef = -9.99e33
   
      !Variables (input, output, inout)
      integer, intent(in):: n, mzp, ntimes

      real, intent(in) :: p2d(:)

      character(len=*), intent(in) :: cumulus
      
      !Local variables:
      integer:: nvartotal, klevgrads(200), jk, int_byte_size, nvar, maxklevgrads
      real   :: real_byte_size
      integer :: nrec, rec_size

      nrec = 0
      maxklevgrads = min(60, mzp)
      runname = '15geos5_'//cumulus
      runlabel = runname

      print *, "writing grads control file:',trim(runname)//'.ctl", ntimes
      !call flush (6)
      !
      !number of variables to be written
      nvartotal = 0
      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") nvartotal = nvartotal + 1
         if (cupout(nvar)%varn(3) == "3d") klevgrads(nvar) = maxklevgrads
         if (cupout(nvar)%varn(3) == "2d") klevgrads(nvar) = 1
      end do

      !- binary file
      inquire (iolength=int_byte_size) real_byte_size  ! inquire by output list

      print *, 'opening grads file:', trim(runname)//'.gra'
      rec_size = size(cupout(nvar)%varp, 1)*real_byte_size
      if (ntimes == 1) then
         open (19, file=trim(runname)//'.gra', form='unformatted', &
               access='direct', status='replace', recl=rec_size)
      else
         open (19, file=trim(runname)//'.gra', form='unformatted', &
               access='direct', status='old', recl=rec_size)
      end if

      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") then
            do jk = 1, klevgrads(nvar)
               nrec = nrec + 1
               !write(19)          real((cupout(nvar)%varp(:,jk)),4)
               write (19, rec=nrec) real((cupout(nvar)%varp(:, jk)), 4)
            end do
         end if
      end do
      close (19)
      !-setting vertical dimension '0' for 2d var
      where (klevgrads == 1) klevgrads = 0
      !- ctl file
      open (20, file=trim(runname)//'.ctl', status='unknown')
      write (20, 2001) '^'//trim(runname)//'.gra'
      write (20, 2002) 'undef -9.99e33'
      write (20, 2002) 'options sequential byteswapped' ! zrev'
      write (20, 2002) 'title '//trim(runlabel)
      write (20, 2003) 1, 0., 1. ! units m/km
      write (20, 2004) n, 1., 1.
      write (20, 2005) maxklevgrads, (p2d(jk), jk=1, maxklevgrads)
      write (20, 2006) ntimes, '00:00Z24JAN1999', '10mn'
      write (20, 2007) nvartotal
      do nvar = 1, p_nvar_grads
         if (cupout(nvar)%varn(1) .ne. "xxxx") then
            !
            write (20, 2008) cupout(nvar)%varn(1) (1:len_trim(cupout(nvar)%varn(1))), klevgrads(nvar) &
               , cupout(nvar)%varn(2) (1:len_trim(cupout(nvar)%varn(2)))
         end if
      end do
      write (20, 2002) 'endvars'
      close (20)

2001  format('dset ', a)
2002  format(a)
2003  format('xdef ', i4, ' linear ', 2f15.3)
2004  format('ydef ', i4, ' linear ', 2f15.3)
2005  format('zdef ', i4, ' levels ', 60f8.3)
2006  format('tdef ', i4, ' linear ', 2a15)
2007  format('vars ', i4)
2008  format(a10, i4, ' 99 ', a40)!'[',a8,']')
!2055  format(60f7.0)
!133   format(1x, F7.0)

   end subroutine wrtBinCtl
   !----------------------------------------------------------------------
   subroutine gen_random(its,ite,use_random_num,random)
      implicit none
      integer, intent(in)  :: its,ite
      real,    intent(in)  :: use_random_num
      real,    intent(out) :: random(:)

      !-local vars
      integer   :: i
      integer(8) :: iran, ranseed = 0
      integer(8), parameter :: ix=2147483646
      call system_clock(ranseed)
      ranseed=mod(ranseed,ix)+1 !seed between 1 and 2^31-2
      iran = -ranseed

      !-- ran1 produces numbers between [ 0,1]
      !-- random        will be between [-1,1]
      !-- with use_random_num the interval will be [-use_random_num,+use_random_num]
      do i=its,ite
         random(i) = use_random_num * 2.0*(0.5-real(RAN1(IRAN),4))
        !print*,"ran=",i,random(i)
      enddo

      if(maxval(abs(random)) > use_random_num) stop "random > use_random_num"

   end subroutine 
   !----------------------------------------------------------------------
   real(8) function ran1(idum)

      ! This is contributed code standardized by Yong Wang
      ! Random number generator taken from Press et al.
      !
      ! Returns numbers in the range 0-->1
      !
      ! Their description...
      ! "Minimal" random number generator of Park and Miller with Bays-Durham
      ! shuffle and added safeguards. Returns a uniform deviate between 0.0 and 1.0
      ! (exclusive of the endpoint values). Call with idum a negative integer to
      ! initialize; thereafter, do not alter idum between successive calls in a
      ! sequence. RNMX should approximate the largest floating value that is less
      ! than 1.

      !use shr_kind_mod,    only: r8 => shr_kind_r8, i8 => shr_kind_i8
      implicit none
      integer(8), parameter:: ntab = 32,iq = 127773,ia = 16807,ir = 2836,ium=1
      integer(8), parameter:: im = 2147483647,ndiv = 1+(im-1)/ntab
      real(8)   , parameter:: am = 1.0/im,eps = 1.2e-7,rnmx = 1.0-eps
      
      integer(8), intent(inout):: idum

      integer(8):: iy
      integer(8), dimension(ntab):: iv
      !save iv,iy
      data iv /ntab*0/, iy /0/
      integer(8):: j,k

      !
      if (idum.le.0.or.iy.eq.0) then
         ! initialize
         idum = max(-idum,ium)
         do j = ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum = idum+im
            if (j.le.ntab) iv(j) = idum
         end do
         iy = iv(1)
      end if
      !
      k = idum/iq
      ! compute idum = mod(ia*idum,im) without overflows by schrage's method
      idum = ia*(idum-k*iq)-ir*k
      if (idum.lt.0) idum = idum+im
      ! j will be in the range 1-->ntab
      j = 1+iy/ndiv
      ! output previously stored value and refill the shuffle table
      iy = iv(j)
      iv(j) = idum
      ran1 = min(am*iy,rnmx)

   end function ran1
   !---------------------------------------------------------------------------------------------------
   subroutine calc_lcl(t0,pp0,r0,tlcl,plcl,dzlcl)
      implicit none
      real,intent(in ) :: t0,pp0,r0
      real,intent(out) :: tlcl,plcl,dzlcl
      real :: ttd

      ttd=td(pp0,r0)
      tlcl=ttd-(0.001296*ttd+0.1963)*(t0-ttd)
      plcl=pp0*(tlcl/t0)**c_cpor
      dzlcl=127*(t0-ttd)
      if(dzlcl.le.0.)dzlcl=-999.
   
   end subroutine calc_lcl
   !---------------------------------------------------------------------------------------------------
   real function td(p,rs)
      implicit none
      real :: rr,rs,es,esln,p
      rr=rs+1e-8
      es=p*rr/(.622+rr)
      esln=log(es)
      td=(35.86*esln-4947.2325)/(esln-23.6837)
      return
   end function td
   !---------------------------------------------------------------------------------------------------
   subroutine reset_1d(its,ite,ierrc,xland,xland1,aa0,aa1,aa2,aa3&
                   ,aa1_bl,aa1_fa,aa0_bl,q_adv,aa1_radpbl,aa1_adv,alpha_adv,cin1&
                   ,xk_x,edt,edto,tau_bl,q_wetbulb,t_wetbulb,tau_ecmwf,xf_dicycle&
                   ,x_add_buoy,xf_coldpool,wlpool_bcon,ke_gustfront,random,mbdt  &
                   ,ave_from_surface)
      implicit none       
      logical, intent (inout)                     :: ave_from_surface
      integer, intent (in)                        :: its,ite
      character*(*), intent (out), dimension(:)   :: ierrc
      real          ,intent (in ), dimension(:)   :: xland       
      real          ,intent (out), dimension(:)   ::   &       
         xland1,& 
         aa0   ,& 
         aa1   ,& 
         aa2   ,& 
         aa3   ,& 
         aa1_bl,& 
         aa1_fa,& 
         aa0_bl,& 
         q_adv ,& 
         aa1_radpbl,&
         aa1_adv   ,&
         alpha_adv ,&
         cin1      ,&
         xk_x      ,&
         edt       ,&
         edto      ,&
         tau_bl    ,&
         q_wetbulb ,&
         t_wetbulb ,&
         tau_ecmwf ,&
         xf_dicycle,&
         x_add_buoy  ,&
         xf_coldpool ,&
         wlpool_bcon ,&
         ke_gustfront,&
         random      ,&
         mbdt
         !---
         ierrc  (:) = "ierrtxt"
         
         xland1 (:) = xland(:) ! 1.
       
         aa0    (:) = 0.0
         aa1    (:) = 0.0
         aa2    (:) = 0.0
         aa3    (:) = 0.0
         aa1_bl (:) = 0.0
         aa1_fa (:) = 0.0
         aa0_bl (:) = 0.0
         q_adv  (:) = 0.0
         aa1_radpbl (:) = 0.0
         aa1_adv    (:) = 0.0
         alpha_adv  (:) = 0.0
         cin1       (:) = 0.0
         xk_x       (:) = 0.0
         edt        (:) = 0.0
         edto       (:) = 0.0
         tau_bl     (:) = 0.0
         q_wetbulb  (:) = 0.0
         t_wetbulb  (:) = 0.0
         tau_ecmwf  (:) = 0.0
         xf_dicycle (:) = 0.0
         x_add_buoy (:) = 0.0
         xf_coldpool(:) = 0.0
         wlpool_bcon(:) = 0.0
         ke_gustfront(:)= 0.0
         random      (:)= 0.0
         !
         !--- mbdt ~ xmb * timescale
         mbdt(:)= 0.1
        !mbdt(i)= 100.*(p_cup(kbcon(i),i)-p_cup(i,kbcon(i)+1))/(c_grav*dtime)
         !-- default flag for get_cloud_bc (only get_lcl uses it as 'true')
         ave_from_surface = .false. 

   end subroutine reset_1d
   !---------------------------------------------------------------------------------------------------
   subroutine reset_2d(its,ite,kts,kte,zo,z,xz,hcdo,cupclw,qrcdo&
                      ,hcot,xf_ens,pr_ens,evap_bcb,uc,vc,hc,hco,zuo,zdo,zenv,Tpert)
      implicit none       
      integer  ,intent (in )                   :: kts,kte,its,ite
      real     ,intent (in ), dimension(:,:)   ::   zo   
      real     ,intent (out), dimension(:,:)   ::   & 
           z     ,&
           xz    ,&
           hcdo  ,&
           cupclw,&
           qrcdo ,&
           hcot  ,&
           xf_ens,&
           pr_ens,&
           evap_bcb, &
           uc  ,&
           vc  ,&
           hc  ,&
           hco ,&
           zuo ,&
           zdo ,&
           zenv,&
           Tpert
 
      z     (:,:) = zo(:,:)
      xz    (:,:) = zo(:,:)
      hcdo  (:,:) = 0.0
      uc    (:,:) = 0.0
      vc    (:,:) = 0.0
      hc    (:,:) = 0.0
      hco   (:,:) = 0.0
      zuo   (:,:) = 0.0
      zdo   (:,:) = 0.0
      zenv  (:,:) = 0.0
      cupclw(:,:) = 0.0
      qrcdo (:,:) = 0.0
      hcot  (:,:) = 0.0
      xf_ens(:,:) = 0.0
      pr_ens(:,:) = 0.0
      Tpert (:,:) = 0.0
      evap_bcb(:,:) = 0.0

   end subroutine reset_2d
   !---------------------------------------------------------------------------------------------------
   
