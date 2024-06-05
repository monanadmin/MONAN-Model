module lib_module_cu_gf_monan
   
   use modGate, only: cupout, rundata, p_nvar_grads, jl, p_use_gate, ppres, ptemp, pq, pu &
                  ,   runlabel, runname, pv, pvervel, pgeo, zqr, zadvq, zadvt

   use modConstants, only: c_rgas, c_cp, c_rm, c_p00, c_tcrit, c_grav, c_cpor, c_alvl, c_pi &    
                        ,  c_akmin, c_ccnclean, c_t00, c_t_ice, c_xlf, p_max_qsat           &
                        ,  p_xmbmaxshal, p_mintracer, p_smaller_qv, c_t01, c_t100, c_temp0i &
                        ,  c_rgas_atm, c_hplus, c_r2es, c_r3les, c_r4ies, c_r4les, c_retv   &
                        ,  c_rticecu, c_rtwat_rticecu_r, c_r3ies, c_r5alscp, c_r5alvcp      &
                        ,  c_ralsdcp, c_ralvdcp, c_rtice, c_rtwat_rtice_r, p_tkmin, c_tice  &
                        ,  p_ccnclean,  c_rtt, c_rtwat
   !
   contains
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
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
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
      !! ## ???
      !!
      !! Author: Saulo Freitas [SRF] e Georg Grell [GAG]
      !!
      !! E-mail: <mailto:saulo.r.de.freitas@gmail.com>, <mailto:georg.a.grell@noaa.gov>
      !!
      !! Date: 2014
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! ???
      !!
      !! ** History**:
      !!
      !! - 
      !! ---
      !! <img src="https://www.gnu.org/graphics/gplv3-127x51.png" width="63">
      !! 
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
   subroutine bidiag (m,b,c,f)
      !-- this routine solves the problem:  bb*f(k,t+1) + cc*f(k+1,t+1) = dd
      !-- an updated "f" at time t+1 is the output
      implicit none
      integer, intent(in) :: m
      real, dimension(m), intent(inout) :: b,c
      real, dimension(m), intent(inout) :: f
      !--locals
      real, dimension(m) :: q
      integer :: k
      real :: p

      c(m)=0.
      q(1)=-c(1)/b(1)
      f(1)= f(1)/b(1)
      do k=2,m
         p  = 1./b(k)
         q(k) = -c(k)*p
         f(k) =  f(k)*p
      enddo
      do k=m-1,1,-1
         f(k) = f(k) +q(k)*f(k+1)
      enddo
   end subroutine bidiag

   !------------------------------------------------------------------------------------
   subroutine get_Qadv(cumulus,itf,ktf,its,ite,kts,kte,ierr,dt,q,qo,qo_adv,po,po_cup &
                      ,qeso, Q_adv,col_sat_adv,alpha_adv,tau_bl,zo_cup,kbcon,ktop)

      implicit none
      real    :: alpha_adv_tuning      = 0.8  != tuning parameter for the Becker et al (2021) closure
      real    :: col_sat_adv_threshold = 0.94 != suppress Qadv closure for col_sat_adv > col_sat_adv_threshold

      character *(*), intent (in)                      :: cumulus
      integer ,intent (in)                             :: itf,ktf, its,ite, kts,kte
      real    ,intent (in)                             :: dt
      integer ,intent (in) ,dimension(:)   :: ierr,kbcon,ktop

      real    ,intent (in) ,dimension(:,:) :: q,qo,qo_adv,po_cup,qeso,po,zo_cup
      real    ,intent (in) ,dimension(:)   :: tau_bl

      real    ,intent (inout) ,dimension(:)   :: Q_adv,col_sat_adv,alpha_adv

      !--locals
      integer :: i,k
      real :: dp, layer,H_cloud,dz
      real ,parameter :: ptop = 60.
      !-- get the advective moisture tendency scaled with the relative humidity
      !--  Q_adv = integral( q/q*  DQv/Dt_adv dp), see Eq 1 Becker et al(2021 QJRMS)
      !-- units here are "J m^-3" _or_  "J kg^-1"

      do i=its,itf
         col_sat_adv(i) = 0.   !check if it needs be inout, perhavps only local var

         if(ierr(i) /= 0) cycle

         alpha_adv  (i) = alpha_adv_tuning
         layer = 0.

         loopN: do k=kts,ktf
            if(po(k,i) < ptop) exit loopN

            !dp=100.*(po_cup(k+1,i)-po_cup(k,i)) ! dp < 0.
            dz=      zo_cup(k+1,i)-zo_cup(k,i)  ! dz > 0

            !-- integral over dp
            !Q_adv(i) = Q_adv(i) + dp*(qo(k,i)/qeso(k,i))*(qo_adv(k,i)-q(k,i))/dt

            !-- integral over dz
            Q_adv(i) = Q_adv(i) + dz*(qo(k,i)/qeso(k,i))*(qo_adv(k,i)-q(k,i))/dt

            col_sat_adv(i) = col_sat_adv(i) + dz* qo(k,i)/qeso(k,i)

            layer = layer + dz

         enddo loopN
         !-- get the column-average saturation fraction
         col_sat_adv(i) = col_sat_adv(i)/(1.e-8+layer)

         !--check if the col-ave saturation fraction is over the threshold
         if(col_sat_adv(i) > col_sat_adv_threshold) then

            alpha_adv(i) = 0.0
            cycle

         endif

         !-- check if cloud top _OR_cloud layer   !<<<< check this
         H_cloud =  zo_cup(ktop(i),i)- zo_cup(kbcon(i),i)

         !-- convert Q_adv to units as in Eq (1) => J m^-3
         !Q_adv(i) = - Q_adv(i) * tau_bl(i) * c_alvl / (g * H_cloud)

         !-- convert Q_adv to units as in cloud work function => J kg-1
         Q_adv(i) =  Q_adv(i) * tau_bl(i) * c_alvl / (H_cloud)

          !if(abs(q_adv(i))>1.) print*,"Qadv=",i,q_adv(i),Q_adv_dz(i)call flush(6)
      enddo

   end subroutine get_Qadv


   !------------------------------------------------------------------------------------
   subroutine get_condensation(q_old,t_old,po_cup,q_new,t_new)

      !-- calculate condensation and adjust t and q accordingly
      implicit none
      real    ,intent (in   ) :: po_cup,q_old,t_old ! before condensation
      real    ,intent (inout) ::        q_new,t_new ! after  condensation

      !---locals
      real ::  zqp, zcond, zcond1, zcor, zqsat,zi,zl,zf
      real :: psp, pt , pq
      real :: z3es,   z4es, z5alcp, zaldcp
      real :: ptare, cond
      real :: foeewmcu,foealfcu,foedemcu,foeldcpmcu

      ! real, parameter :: &
      !     RD=287.06                             &
      !    ,RV=461.52                             &
      !    ,RCPD=1004.71                          &
      !    ,RTT=273.16                            &
      !    ,RHOH2O=1000.                          &
      !    ,RLVTT=2.5008e+6                       &
      !    ,RLSTT=2.8345e+6                       &
      !    ,RETV=RV/RD-1.0                        &
      !    ,RLMLT=RLSTT-RLVTT                     &
      !    ,RCPV=4.*RV                            &
      !    ,R2ES=611.21*RD/RV                     &
      !    ,R3LES=17.502                          &
      !    ,R3IES=22.587                          &
      !    ,R4LES=32.19                           &
      !    ,R4IES=-0.7                            &
      !    ,R5LES=R3LES*(RTT-R4LES)               &
      !    ,R5IES=R3IES*(RTT-R4IES)               &
      !    ,R5ALVCP=R5LES*RLVTT/RCPD              &
      !    ,R5ALSCP=R5IES*RLSTT/RCPD              &
      !    ,RALVDCP=RLVTT/RCPD                    &
      !    ,RALSDCP=RLSTT/RCPD                    &
      !    ,RALFDCP=RLMLT/RCPD                    &
      !    ,RTWAT=RTT                             &
      !    ,RTBER=RTT-5.                          &
      !    ,RTBERCU=RTT-5.0                       &
      !    ,RTICE=RTT-23.                         &
      !    ,RTICECU=RTT-23.                       &
      !    ,RTWAT_RTICE_R=1./(RTWAT-RTICE)        &
      !    ,RTWAT_RTICECU_R=1./(RTWAT-RTICECU)    &
      !    ,RVTMP2=RCPV/RCPD-1.                   &
      !    ,ZQMAX=0.5

      !----------------------
      !-- for testing
      !----------------------
      !                          PSP (hPA)            TEMP(K)              Q(kg/kg)                ZCOND1(kg/kg)
      ! input          1   98020.0000000000        295.163188640152     1.745679200989956E-002
      ! output         1   98020.0000000000        295.513490789916     1.731605637618801E-002 -9.779916453243843E-007
      !----------------------
      ! input       157   85090.0000000000        288.089188935407     1.399404805052166E-002
      ! output      157   85090.0000000000        289.294751760460     1.350970717999820E-002 -1.146268822756454E-005
      !----------------------
      ! PT  = 288.089188935407
      ! PQ  = 1.399404805052166E-002
      ! PSP = 85090
      !----------------------

      !-- initial values
      PT  = t_old       ! K
      PQ  = q_old       ! kg/kg
      PSP = po_cup*100. ! hPa


      !--- get condensation in moist ascent --------------------------
      PTARE = PT
      ZQP    =1.0/PSP

      ZL=1.0/(PT-c_R4LES)
      ZI=1.0/(PT-c_R4IES)

      FOEALFCU = MIN(1.0,((MAX(c_RTICECU,MIN(c_RTWAT,PTARE))-c_RTICECU)*c_RTWAT_RTICECU_R)**2)
      ZQSAT=c_R2ES *(FOEALFCU *EXP(c_R3LES*(PTARE-c_RTT)*ZL)+&
                (1.0-FOEALFCU)*EXP(c_R3IES*(PTARE-c_RTT)*ZI))

      ZQSAT=ZQSAT*ZQP
      ZQSAT=MIN(0.5,ZQSAT)
      ZCOR=1.0-c_RETV*ZQSAT

      ZF=FOEALFCU*c_R5ALVCP*ZL**2 + (1.0-FOEALFCU)*c_R5ALSCP*ZI**2
      ZCOND=(PQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)

      if(ZCOND > 0.0)then
         FOELDCPMCU= FOEALFCU*c_RALVDCP+(1.0-FOEALFCU)*c_RALSDCP
         PT=PT+FOELDCPMCU*ZCOND
         PTARE = PT
         PQ=PQ-ZCOND

         ZL=1.0/(PT-c_R4LES)
         ZI=1.0/(PT-c_R4IES)


         FOEALFCU = MIN(1.0,((MAX(c_RTICECU,MIN(c_RTWAT,PTARE))-c_RTICECU)*c_RTWAT_RTICECU_R)**2)
         ZQSAT=c_R2ES *(FOEALFCU* EXP(c_R3LES*(PT-c_RTT)*ZL)+&
                   (1.0-FOEALFCU)*EXP(c_R3IES*(PT-c_RTT)*ZI))

         ZQSAT=ZQSAT*ZQP
         ZQSAT=ZQSAT-0.5*(ABS(0.5-ZQSAT)-(0.5-ZQSAT))


         ZCOR=1.0-c_RETV*ZQSAT
         ZF=FOEALFCU*c_R5ALVCP*ZL**2 + (1.0-FOEALFCU)*c_R5ALSCP*ZI**2

         ZCOND1=(PQ*ZCOR**2-ZQSAT*ZCOR)/(ZCOR**2+ZQSAT*ZF)
         if(ZCOND ==  0.0)ZCOND1=0.0
         FOELDCPMCU= FOEALFCU*c_RALVDCP+(1.0-FOEALFCU)*c_RALSDCP
         PT=PT+FOELDCPMCU*ZCOND1
         PQ=PQ-ZCOND1
      endif

      !-- FINAL --------------------------
      q_new =  PQ
      t_new =  PT
      cond  = -ZCOND1 != q_old-qnew, source for the liquid water
   end subroutine get_condensation


!********************
!********************
!********************
!********************
!******************** OBSOLETOS *********************************




   !--- get cloud fraction
   !
   ! do i=its,itf
   !    clfrac(i,:)=0.
   !    if(ierr(i) /= 0) cycle
   !    dummy1(kts:ktf) = xmb(i)* zuo(i,kts:ktf)
   !    dummy2(kts:ktf) = 100.*po_cup(i,kts:ktf)
   !    call get_cloud_fraction(ktf,kts,ktf                                                   &
   !     ,dummy2(kts:ktf),zo_cup(i,kts:ktf),tn_cup(i,kts:ktf),qo_cup(i,kts:ktf) &
   !     ,qco (i,kts:ktf),  qrco(i,kts:ktf),  dummy1(kts:ktf),clfrac(i,kts:ktf) )
   ! enddo
   !--------------------------------------------------------------------------------------------!
 
   subroutine get_cloud_fraction(  mzp, kts, ktf, &
      PPABS, PZZ, PT, PRV, QCO, QRCO, PMFLX, PCLDFR ,PRC, PRI )
      !!    PURPOSE
      !!    -------
      !!**  Routine to diagnose cloud fraction and liquid and ice condensate mixing ratios
      !!**  METHOD
      !!    ------
      !!    Based on the large-scale fields of temperature, water vapor, and possibly
      !!    liquid and solid condensate, the conserved quantities r_t and h_l are constructed
      !!    and then fractional cloudiness, liquid and solid condensate is diagnosed.
      !!
      !!    The total variance is parameterized as the sum of  stratiform/turbulent variance
      !!    and a convective variance.
      !!    The turbulent variance is parameterized as a function of first-order moments, and
      !!    the convective variance is modelled as a function of the convective mass flux (units kg/s m^2)
      !!    as provided by a mass flux convection scheme.
      !!
      !!    Nota: if the host model does not use prognostic values for liquid and solid condensate
      !!    or does not provide a convective mass flux, put all these values to zero.
      !!    Also, it is supposed that vertical model levels are numbered from
      !!    1 to MZP, where 1 is the first model level above the surface
      !!
      !!    ------------------
      !!    REFERENCE
      !!    ---------
      !!      Chaboureau J.P. and P. Bechtold (J. Atmos. Sci. 2002)
      !!      Chaboureau J.P. and P. Bechtold (JGR/AGU 2005)
      !!
      !!    AUTHOR
      !!    ------
      !!      P. BECHTOLD       * Laboratoire d'Aerologie *
      !!
      !!    MODIFICATIONS
      !!    -------------
      !!      Original    13/06/2001
      !!      modified    20/03/2002 : add convective Sigma_s and improve turbulent
      !!                               length-scale in boundary-layer and near tropopause
      !!      adapted     09/12/2016 : adapted to GEOS-5 by Saulo Freitas
      !-------------------------------------------------------------------------------
      !*       0.    DECLARATIONS
      !              ------------
      implicit none
      !
      !-------------------------------------------------------------------------------
      !
      !*       1.    Set the fundamental thermodynamical constants
      !              these have the same values (not names) as in ARPEGE IFS
      !              -------------------------------------------------------
      real, parameter :: XP00   = 1.e5        ! reference pressure
      real, parameter :: XPI    = 3.141592654 ! C_pi
      real, parameter ::  XG    = 9.80665     ! gravity constant
      real, parameter :: XMD    = 28.9644e-3  ! molecular weight of dry air
      real, parameter :: XMV    = 18.0153e-3  ! molecular weight of water vapor
      real, parameter :: XRD    = 287.05967   ! gaz constant for dry air
      real, parameter :: XRV    = 461.524993  ! gaz constant for water vapor
      real, parameter :: XCPD   = 1004.708845 ! specific heat of dry air
      real, parameter :: XCPV   = 1846.1      ! specific heat of water vapor
      real, parameter :: XRHOLW = 1000.       ! density of liquid water
      real, parameter :: XCL    = 4218.       ! specific heat of liquid water
      real, parameter :: XCI    = 2106.       ! specific heat of ice
      real, parameter :: XTT    = 273.16      ! triple point temperature
      real, parameter :: C_ALVLTT  = 2.5008e6    ! latent heat of vaporisation at XTT
      real, parameter :: XLSTT  = 2.8345e6    ! latent heat of sublimation at XTT
      real, parameter :: XLMTT  = 0.3337e6    ! latent heat of melting at XTT
      real, parameter :: XESTT  = 611.14      ! saturation pressure at XTT
      real, parameter :: XALPW  = 60.22416    ! constants in saturation pressure over liquid water
      real, parameter :: XBETAW = 6822.459384
      real, parameter :: XGAMW  = 5.13948
      real, parameter :: XALPI  = 32.62116    ! constants in saturation pressure over ice
      real, parameter :: XBETAI = 6295.421
      real, parameter :: XGAMI  = 0.56313
      logical, parameter :: LUSERI = .true. ! logical switch to compute both
                                                ! liquid and solid condensate (LUSERI=.TRUE.)
                                                ! or only liquid condensate (LUSERI=.FALSE.)
      !
      !*       0.1   Declarations of dummy arguments :
      !
      !
      integer,              intent(in)   :: mzp     ! vertical dimension
      integer,              intent(in)   :: kts     ! vertical  computations start at
      !                                             ! KTS that is at least 1
      integer,              intent(in)   :: ktf     ! vertical computations can be
                                                    ! limited to MZP + 1 - KTF
                                                    ! default=1
      real, dimension(mzp), intent(in)    :: PPABS  ! pressure (Pa)
      real, dimension(mzp), intent(in)    :: PZZ    ! height of model levels (m)
      real, dimension(mzp), intent(in)    :: PT     ! grid scale T  (K)
      real, dimension(mzp), intent(in)    :: PRV    ! grid scale water vapor mixing ratio (kg/kg)
      real, dimension(mzp), intent(in)    :: PMFLX  ! convective mass flux (kg/(s m^2))
      real, dimension(mzp), intent(in)    :: QRCO   ! sub-grid scale liq water mixing ratio (kg/kg)
      real, dimension(mzp), intent(in)    :: QCO    ! in-cloud water mixing ratio (kg/kg)
      real, dimension(mzp), intent(inout),optional :: PRC    ! grid scale r_c mixing ratio (kg/kg)
      real, dimension(mzp), intent(inout),optional :: PRI    ! grid scale r_i (kg/kg)
      real, dimension(mzp), intent(out)   :: PCLDFR ! fractional cloudiness (between 0 and 1)
      !
      !
      !*       0.2   Declarations of local variables :
      !
      integer  ::  JKT, JKP, JKM,K     ! loop index
      real, dimension(mzp) :: ZTLK, ZRT       ! work arrays for T_l, r_t
      real, dimension(mzp) :: ZL              ! length-scale
      integer   :: ITPL    ! top levels of tropopause/highest inversion
      real      :: ZTMIN   ! min Temp. related to ITPL
      real, dimension(mzp) :: LOC_PRC,LOC_PRI
      !
      real :: ZTEMP, ZLV, ZLS, ZTL, ZPV, ZQSL, ZPIV, ZQSI, ZFRAC, ZCOND, ZCPD ! thermodynamics
      real :: ZLL, DZZ, ZZZ ! length scales
      real :: ZAH, ZA, ZB, ZSBAR, ZQ1, ZSIGMA, ZDRW, ZDTL ! related to computation of Sig_s
      real :: ZSIG_CONV,  ZSIGMA_NOCONV,  ZQ1_NOCONV      ! convective part of Sig_s
      !
      !*       0.3  Definition of constants :
      !
      !-------------------------------------------------------------------------------
      !
      real :: ZL0     = 600.        ! tropospheric length scale
                                    ! changed to 600 m instead of 900 m to give a consistent
                                    ! value (linear increase) in general 500 m deep oceanic
                                    ! mixed layer - but could be put back to 900 m if wished
      real :: ZCSIGMA = 0.2         ! constant in sigma_s parameterization
      real :: ZCSIG_CONV = 0.30e-2  ! scaling factor for ZSIG_CONV as function of mass flux
      !
      !
      logical :: ONLY_CONVECTIVE_CLOUD_FRACTION=.true. ! set .false. for the total cloud fraction
      !-------------------------------------------------------------------------------
      !RETURN
      !
      if(PRESENT(PRC)) then
         LOC_PRC(:)=PRC(:)
      else
         LOC_PRC(:)=0.0
      endif
      if(PRESENT(PRI)) then
         LOC_PRI(:)=PRI(:)
      else
         LOC_PRI(:)=0.0
      endif

      PCLDFR(:) = 0. ! Initialize values
      !

      JKT = MZP+1-KTS
      !-will limit the model vertical column to 60 hPa
      do K=KTF,KTS,-1
         if(PPABS(k) > 60.*100.) then
            JKT = k
            !PRINT*,"JKT=",K,MZP+1-KTS ;CALL FLUSH(6)
            exit
         endif
      enddo

      do K=KTS,JKT
         ZTEMP  = PT(k)
          !latent heat of vaporisation/sublimation
         ZLV    = C_ALVLTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
         ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

         !store temperature at saturation and total water mixing ratio
         ZRT(k)   = PRV(k) + LOC_PRC(k) + LOC_PRI(k)
         ZCPD     = XCPD  + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
         ZTLK(k)  = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD
      end do

      !-------------------------------------------------------------------------------
      ! Determine tropopause/inversion  height from minimum temperature

      ITPL  = KTS+1
      ZTMIN = 400.
      do k = KTS+1,JKT-1
         if ( PT(k) < ZTMIN ) then
            ZTMIN = PT(k)
            ITPL  = K
         endif
      end do

      ! Set the mixing length scale - used for computing the "turbulent part" of Sigma_s

      ZL(:) = 20.
      do k = KTS+1,JKT

         ! free troposphere
         ZL(k) = ZL0
         JKP   = ITPL
         ZZZ   = PZZ(k) -  PZZ(KTS)
            ! approximate length for boundary-layer : linear increase
         if ( ZL0 > ZZZ )  ZL(k) = ZZZ
            ! gradual decrease of length-scale near and above tropopause/top inversion
         if ( ZZZ > 0.9*(PZZ(JKP)-PZZ(KTS)) ) &
            ZL(k) = .6 * ZL(K-1)
      end do
      !-------------------------------------------------------------------------------

      do k=KTS+1,JKT-1
         JKP=k+1
         JKM=k-1
         ZTEMP  = PT(k)

         !latent heat of vaporisation/sublimation
         ZLV    = C_ALVLTT + ( XCPV - XCL ) * ( ZTEMP - XTT )
         ZLS    = XLSTT + ( XCPV - XCI ) * ( ZTEMP - XTT )

         ZCPD   = XCPD + XCPV*PRV(k) + XCL*LOC_PRC(k) + XCI*LOC_PRI(k)
         !temperature at saturation
         ZTL    = ZTEMP - ZLV*LOC_PRC(k)/ZCPD - ZLS*LOC_PRI(k)/ZCPD

         !saturated water vapor mixing ratio over liquid water
         ZPV    = MIN(EXP( XALPW - XBETAW / ZTL - XGAMW * LOG( ZTL ) ),0.99*PPABS(k))
         ZQSL   = XRD / XRV * ZPV / ( PPABS(k) - ZPV )

          !saturated water vapor mixing ratio over ice
         ZPIV   = MIN(EXP( XALPI - XBETAI / ZTL - XGAMI * LOG( ZTL ) ),0.99*PPABS(k))
         ZQSI   = XRD / XRV * ZPIV / ( PPABS(k) - ZPIV )

         !interpolate between liquid and solid as function of temperature
         ! glaciation interval is specified here to 20 K
         ZFRAC = ( ZTL  - 250.16 ) / ( XTT - 250.16 )  ! liquid/solid fraction
         ZFRAC = MAX( 0., MIN(1., ZFRAC ) )

         if(.not. LUSERI) ZFRAC=1.
         ZQSL = ( 1. - ZFRAC ) * ZQSI + ZFRAC * ZQSL
         ZLV  = ( 1. - ZFRAC ) * ZLS  + ZFRAC * ZLV

         !coefficients a and b
         ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 ) * (XRV * ZQSL / XRD + 1.)
         !orig  ZAH  = ZLV * ZQSL / ( XRV * ZTL**2 )

         ZA   = 1. / ( 1. + ZLV/ZCPD * ZAH )
         ZB   = ZAH * ZA

         !- parameterize Sigma_s with first_order closure
         DZZ    =  PZZ (JKP)  - PZZ(JKM)
         ZDRW   =  ZRT (JKP)  - ZRT(JKM)
         ZDTL   =  ZTLK(JKP) - ZTLK(JKM) + XG/ZCPD * DZZ
         ZLL    =  ZL(k)

         !- standard deviation due to convection
         ZSIG_CONV = ZCSIG_CONV * PMFLX(k) / ZA

         !- turb + conv
         ZSIGMA = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
            ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
            + ZB*ZB*ZDTL*ZDTL                      ) &
            + ZSIG_CONV * ZSIG_CONV ) )

         !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
         ZSIGMA = MAX( ZSIGMA, 1.e-10 )

         !- normalized saturation deficit
         ZSBAR = ZA * ( ZRT (k) - ZQSL )
         !- "Q1" parameter
         ZQ1   = ZSBAR / ZSIGMA

         !- total cloud fraction
         PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.36*ATAN(1.55*ZQ1)) )

         if(ONLY_CONVECTIVE_CLOUD_FRACTION) then
            !- get cloud fraction associated with ONLY the sub-grid scale convective part
            !- this sigma does not include the sub-grid scale convective part
            ZSIGMA_NOCONV = SQRT( MAX( 1.e-25, ZCSIGMA*ZCSIGMA* ZLL*ZLL/(DZZ*DZZ) * ( &
               ZA*ZA*ZDRW*ZDRW - 2.*ZA*ZB*ZDRW*ZDTL   &
               + ZB*ZB*ZDTL*ZDTL  )))
            !- zsigma should be of order 4.e-4 in lowest 5 km of atmosphere
            ZSIGMA_NOCONV = MAX( ZSIGMA_NOCONV, 1.e-10 )
            ZQ1_NOCONV = ZSBAR / ZSIGMA_NOCONV

            !- cloud fraction associated with ONLY convective part ("total-turb")
            PCLDFR(k) = 0.36*(ATAN(1.55*ZQ1)-ATAN(1.55*ZQ1_NOCONV))

            PCLDFR(k) = MAX( 0., MIN(1.,PCLDFR(k)) )

         endif
         !- newer formulation, see GMD 2015
         !PCLDFR(k) = MAX( 0., MIN(1.,0.5+0.34*ATAN(1.85*ZQ1+2.33)) )
         !- this is area fraction of cloud cores
         !PCLDFR(k) = MAX( 0., MIN(1.,0.292/ZQ1**2) )

         cycle
         !total condensate diagnostic (not being used)
         if (ZQ1 > 0. .and. ZQ1 <= 2. ) then
            !orig   ZCOND =     EXP(-1.)+.66*ZQ1+.086*ZQ1*ZQ1
            ZCOND = MIN(EXP(-1.)+.66*ZQ1+.086*ZQ1**2, 2.) ! We use the MIN function for continuity
         else if (ZQ1 > 2.) then
            ZCOND = ZQ1
         else
            ZCOND = EXP( 1.2*ZQ1-1. )
         end if
         ZCOND = ZCOND * ZSIGMA

         if ( zcond < 1.e-12) then
            zcond = 0.
            pcldfr(k) = 0.
         end if
         if ( pcldfr(k) == 0.) then
            zcond = 0.
         end if

         LOC_PRC(k) = ZFRAC * ZCOND ! liquid condensate
         if (LUSERI) then
            LOC_PRI(k) = (1.-ZFRAC) * ZCOND   ! solid condensate
         end if

      !---
      ! compute s'rl'/Sigs^2
      ! used in w'rl'= w's' * s'rl'/Sigs^2
      !  PSIGRC(k) = PCLDFR(k)   ! Gaussian relation
      !
      ! s r_c/ sig_s^2
      !    PSIGRC(JI,JJ,JK) = PCLDFR(JI,JJ,JK)  ! use simple Gaussian relation
      !
      !    multiply PSRCS by the lambda3 coefficient
      !
      !      PSIGRC(JI,JJ,JK) = 2.*PCLDFR(JI,JJ,JK) * MIN( 3. , MAX(1.,1.-ZQ1) )
      ! in the 3D case lambda_3 = 1.
      !      INQ1 = MIN( MAX(-22,FLOOR(2*ZQ1) ), 10)
      !      ZINC = 2.*ZQ1 - INQ1
      !
      !      PSIGRC(k) =  MIN(1.,(1.-ZINC)*ZSRC_1D(INQ1)+ZINC*ZSRC_1D(INQ1+1))
      !
      !      PSIGRC(k) = PSIGRC(k)* MIN( 3. , MAX(1.,1.-ZQ1) )
      !---
      end do
     !
   end subroutine get_cloud_fraction

 !------------------------------------------------------------------------------------
   function auto_rk(n,step,aux,xexp,qrc1) result(PW)
      integer, intent(in) :: n
      real   , intent(in) :: step,aux,qrc1,xexp
      real                :: PW

      PW=step*qrc1*(1.0-exp(-aux**xexp))/float(n)

   end function auto_rk
!+---+-----------------------------------------------------------------+
   !DSM {
   pure function intfuncgamma(x, y) result(z)
      real :: z
      real, intent(in) :: x, y

      z = x**(y-1.0) * exp(-x)
   end function intfuncgamma

   function gammaBrams(a) result(g)
      real :: g
      real, intent(in) :: a

      real, parameter :: small = 1.0e-4
      integer, parameter :: points = 100000

      real :: infty, dx, p, sp(2, points), x
      integer :: i
      logical :: correction

      x = a

      correction = .false.
      ! value with x<1 gives \infty, so we use
      ! \Gamma(x+1) = x\Gamma(x)
      ! to avoid the problem
      if ( x < 1.0 ) then
         correction = .true.
         x = x + 1
      end if

      ! find a "reasonable" infinity...
      ! we compute this integral indeed
      ! \int_0^M dt t^{x-1} e^{-t}
      ! where M is such that M^{x-1} e^{-M} ≤ \epsilon
      infty = 1.0e4
      do while ( intfuncgamma(infty, x) > small )
         infty = infty * 10.0
      end do

      ! using simpson
      dx = infty/real(points)
      sp = 0.0
      forall(i=1:points/2-1) sp(1, 2*i) = intfuncgamma(2.0*(i)*dx, x)
      forall(i=1:points/2) sp(2, 2*i - 1) = intfuncgamma((2.0*(i)-1.0)*dx, x)
      g = (intfuncgamma(0.0, x) + 2.0*sum(sp(1,:)) + 4.0*sum(sp(2,:)) + &
         intfuncgamma(infty, x))*dx/3.0

      if ( correction ) g = g/a

   end function gammaBrams
   !DSM}   !------------------------------------------------------------------------------------
   real function coldPoolStart_orig(CNV_TR)
      implicit none
      real,intent(in)  :: CNV_TR
      real             :: f1
      real,parameter   :: width = 100. !orig 100.
      real    :: mx_buoy1       = (c_cp*5.0 + c_alvl*2.e-3)*0.025  !=   250.5 J/kg
      real    :: mx_buoy2       = (c_cp*10. + c_alvl*4.e-3)        != 20004.0 J/kg: temp exc=10 K, q deficit=4 g/kg (=> mx_buoy ~ 20 kJ/kg)

      f1= min (mx_buoy2,CNV_TR)
      !--- f1 > mx_buoy1 => coldPoolStart ---> 1 
      coldPoolStart_orig =  (1.35+atan( (f1-mx_buoy1)/width))/2.8
      coldPoolStart_orig =  max(0.00,min(coldPoolStart_orig,1.00))
   end function
   !------------------------------------------------------------------------------------
   subroutine fix_scale_dep(cumulus,kts,kte,ktf,zo_cup,w,rho,v_ratio)
      implicit none
      character *(*),intent (in) :: cumulus
      integer       ,intent (in) :: kts,kte,ktf
      real          ,intent (in) :: zo_cup(kts:kte),w(kts:kte),rho(kts:kte)
      
      real          ,intent (out) :: v_ratio
      !-- local vars
      integer :: k
      real :: vert_int_w,total_dz,dz
      real :: vert_int_w_threshold  = 1.0  != grid-scale vertical velocity threshold (m/s)
      real :: vert_height_threshold = 3000.!= max height to calculate the mean grid-scale vertical velocity (m)

      total_dz    = 0.
      vert_int_w  = 0.

      loopK: do k = kts,ktf 
      
         dz         = zo_cup(k+1) - zo_cup(k)   
         total_dz   = total_dz    + dz
         vert_int_w = vert_int_w  + dz * w(k) ! vert veloc m/s
         if( zo_cup(k+1) >= vert_height_threshold) exit loopK 
         
      end do loopK
      vert_int_w =  vert_int_w / (1.e-12+total_dz)
      
      v_ratio    = (vert_int_w / vert_int_w_threshold)**4.
 
   end subroutine fix_scale_dep

!------------------------------------------------------------------------------------
   subroutine get_zu_zd_pdf_orig(draft,ierr,kb,kt,zs,zuf,ztop,zu,kts,kte,ktf)

      implicit none
      integer, intent(in) ::kb,kt,kts,kte,ktf
      real, intent(in) :: Zs,Zuf,Ztop
      real, intent(inout) :: zu(kts:kte)
      integer, intent(inout) :: ierr
      character*(*), intent(in) ::draft

      !- local var
      integer :: add,i,nrec=0,k,kb_adj
      real ::zumax,ztop_adj
      real ::beta, alpha,kratio,tunning

      !- kb cannot be at 1st level
      kb_adj=max(kb,2)

      !-- fill zu with zeros
      zu=0.0

      if(draft == "UP" .or. draft == "up" ) then
         if(kt<=kb_adj) then
            !stop "ktop must be larger than kbcon"
            ierr=99
            return
         endif
         !beta=4.  !=> must larger than 1
                   !=> higher makes the profile sharper
                   !=> around the maximum zu
         add=0     !=> additional levels above kbcon, where
                   !=> the maximum zu will resides
         kb_adj=kb_adj+add
         kb_adj=max(10,kb_adj)

         !- this alpha constrains the location of the maximun ZU to be at
         !- "kb_adj" vertical level
         !alpha=1. + (beta-1.0)*(float(kb_adj)/float(kt+1))/(1.0-(float(kb_adj)/float(kt+1)))

         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         tunning = 0.6
         beta    = 2.0/tunning
         alpha   = tunning*beta

          !- Beta PDF
         do k=kts,min(kte,kt+1)
            kratio= float(k)/float(kt+1)

            zu(k) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      elseif(draft == "DOWN" .or. draft == "DOWNM") then
         add=0    !=> additional levels above kbcon, where
                  !=> the maximum zu will resides
         beta=4.  !=> must larger than 1
                  !=> higher makes the profile sharper
                  !=> around the maximum zu
         alpha= 0.25*beta

         !- for downdrafts kb = jmin(i)-levadj
         kb_adj=kb_adj+add

         !- 2nd approach for beta and alpha parameters
         !- the tunning parameter must be between 0.5 (low  level max zu)
         !-                                   and 1.5 (high level max zu)
         tunning = 1.
         beta    = 2.0/tunning
         alpha   = tunning*beta

         !- Beta PDF
         do k=kts,min(kte,kt)
            kratio= float(k)/float(kt)

            zu(k+1) = kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      elseif(draft == "shallow" .or. draft == "SHALLOW") then

         alpha= 3.
         beta = 2.*alpha
         kb_adj=1 ! level where mass flux starts

         !- Beta PDF
         do k=kts+kb_adj-1,min(kte,kt+1)
            kratio=float(k+1-kb_adj)/float(kt+1)  !-kb_adj+1)

            zu(k)=kratio**(alpha-1.0) * (1.0-kratio)**(beta-1.0)
         enddo

      else
         print*, "unknown type of flow" ,draft
         stop "routine get_zu_zd"

      endif

      !- normalize ZU
      zu(kts:min(kte,kt+1))= zu(kts:min(kte,kt+1))/ maxval(zu(kts:min(kte,kt+1)))

      !--- Sanity checks
      if(beta <= 1) stop "beta must be larger than 1"

      if(minval(zu(:)) < 0.0 ) then
         print*," zu has negative values for ", draft
         stop   " zu < zero"
      endif
      if(maxval(zu(:)) > 1.0 ) then
         print*," zu has values greater than 1 for ", draft
         stop   " zu  >  one"
      endif

      return

   !OPEN(19,FILE= 'zu.gra', FORM='unformatted',ACCESS='direct'&
   !       ,STATUS='unknown',RECL=4)
   ! DO k = kts,kte
   !    nrec=nrec+1
   !     WRITE(19,REC=nrec) zu(k)
   ! END DO
   !close (19)

   end subroutine get_zu_zd_pdf_orig

  !------------------------------------------------------------------------------------

   real function satvap(temp2)
      implicit none
      real :: temp2, temp, toot, toto, eilog, tsot,  &
         ewlog, ewlog2, ewlog3, ewlog4
      temp = temp2-273.155
      if (temp.lt.-20.) then   !!!! ice saturation
         toot = 273.16 / temp2
         toto = 1 / toot
         eilog = -9.09718 * (toot - 1) - 3.56654 * (log(toot) / &
            log(10.)) + .876793 * (1 - toto) + (log(6.1071) / log(10.))
         satvap = 10 ** eilog
      else
         tsot = 373.16 / temp2
         ewlog = -7.90298 * (tsot - 1) + 5.02808 * &
            (log(tsot) / log(10.))
         ewlog2 = ewlog - 1.3816e-07 * &
            (10 ** (11.344 * (1 - (1 / tsot))) - 1)
         ewlog3 = ewlog2 + .0081328 * &
            (10 ** (-3.49149 * (tsot - 1)) - 1)
         ewlog4 = ewlog3 + (log(1013.246) / log(10.))
         satvap = 10 ** ewlog4
      endif

   end function

!---------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------
end  module lib_module_cu_gf_monan
