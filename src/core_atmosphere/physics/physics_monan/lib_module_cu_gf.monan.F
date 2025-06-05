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
   !------------------------------------------------------------------------------------
   ! !
   ! !--- Implements Becker et al (2021) closure, part 1
   ! !
   ! if (DICYCLE == 2 .and. trim(cumulus) == 'deep') then
   !       call BeckerEtAlClosure(cumulus,ave_layer,its, itf, ite, kts, ktf, kte, k22, start_level, ktop, klcl, kbcon &
   !                            , t, tn, tn_adv, q, qo, qo_adv, psur, tsur,us, vs, zqexec, ztexec                     &
   !                            , x_add_buoy, xland, tpert, zuo, up_massdetro , up_massentro, zuo, p_liq_ice          &
   !                            , qrco, ierr, zo, po, z1,aa1_radpbl, aa1_adv,dicycle)
   ! end if
   ! !--- Adds Becker et al (2021) closure, part 2
   ! !
   ! if( dicycle==2 ) then
   !          call get_Qadv(cumulus,itf,ktf,its,ite,kts,kte,ierr,dtime,q,qo,qo_adv,po,po_cup &
   !                       ,qeso, q_adv,col_sat_adv,alpha_adv,tau_bl,zo_cup,kbcon,ktop)
   ! endif
   ! subroutine BeckerEtAlClosure(cumulus,ave_layer,its, itf, ite, kts, ktf, kte, k22, start_level, ktop, klcl, kbcon &
   !                            , t_in, tn, tn_adv , q_in, qo, qo_adv, psur, tsur,us, vs, zqexec, ztexec, x_add_buoy  &
   !                            , xland, tpert, zu, up_massdetro, up_massentro, zuo, p_liq_ice, qrco, ierr, zo, po, z1 &
   !                            , aa1_radpbl, aa1_adv,DICYCLE)
   !    !! ## Becker et Al closure
   !    !!
   !    !! Author: Saulo R Freitas [SRF]
   !    !!
   !    !! E-mail: <mailto:saulo.freitas@inpe.br>
   !    !!
   !    !! Date: 22Fevereiro2023 13:46
   !    !!
   !    !! #####Version: 0.1.0
   !    !!
   !    !! ---
   !    !! **Full description**:
   !    !!
   !    !! Becker et Al closure
   !    !!
   !    !! ** History**:
   !    !!
   !    !! --- 
   !    !! ** Licence **: Under the terms of the GNU General Public version 3
   !    !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !    !!

   !    implicit none
   !    ! Parameters:
   !    character(len=*), parameter :: p_procedure_name = 'beckerEtAlClosure' 
   !    !! subroutine name
   
   !    ! Variables (input, output, inout)
   !    character *(*)   ,intent (in) :: cumulus
   !    real   , intent(in) :: ave_layer
   !    integer, intent(in) :: its, itf, ite, kts, ktf, kte, DICYCLE
   !    integer, intent(in) :: k22(:)
   !    integer, intent(in) :: start_level(:)
   !    integer, intent(in) :: ktop(:)
   !    integer, intent(in) :: klcl(:)
   !    integer, intent(in) :: kbcon(:)
   !    real, intent(in) :: t_in(:,:)
   !    real, intent(in) :: tn(:,:)
   !    real, intent(in) :: tn_adv(:,:)
   !    real, intent(in) :: q_in(:,:)
   !    real, intent(in) :: qo(:,:)
   !    real, intent(in) :: qo_adv(:,:)
   !    real, intent(in) :: psur(:)      
   !    real, intent(in) :: tsur(:)
   !    real, intent(in) :: us(:, :)
   !    real, intent(in) :: vs(:, :)
   !    real, intent(in) :: zqexec(:)
   !    real, intent(in) :: ztexec(:)
   !    real, intent(in) :: x_add_buoy(:)
   !    real, intent(in) :: xland(:)
   !    real, intent(in) :: tpert(:,:)
   !    real, intent(in) :: zu(:,:)
   !    real, intent(in) :: up_massdetro(:,:)
   !    real, intent(in) :: up_massentro(:,:)
   !    real, intent(in) :: zuo(:,:)
   !    real, intent(in) :: p_liq_ice(:,:)
   !    real, intent(in) :: qrco(:,:)
   !    integer, intent(inout) :: ierr(:)
      
   !    real, intent(in ) :: zo(:,:)
   !    real, intent(in ) :: po(:, :)
   !    real, intent(in ) :: z1(:)
   !    real, intent(out) :: aa1_radpbl(:)
   !    real, intent(out) :: aa1_adv(:)
      
   !    ! Local variables:
   !    integer ::ierr_dummy(its:ite) 
   !    real :: hkbo_x      (its:ite)
   !    real :: qeso_x      (kts:kte,its:ite)
   !    real :: heo_x       (kts:kte,its:ite)
   !    real :: heso_x      (kts:kte,its:ite)
   !    real :: qeso_cup_x  (kts:kte,its:ite)
   !    real :: qo_cup_x    (kts:kte,its:ite)
   !    real :: heo_cup_x   (kts:kte,its:ite)
   !    real :: u_cup_x     (kts:kte,its:ite)
   !    real :: v_cup_x     (kts:kte,its:ite)
   !    real :: heso_cup_x  (kts:kte,its:ite)
   !    real :: zo_cup_x    (kts:kte,its:ite)
   !    real :: po_cup_x    (kts:kte,its:ite)
   !    real :: gammao_cup_x(kts:kte,its:ite)
   !    real :: tn_cup_x    (kts:kte,its:ite)
   !    real :: hco_x       (kts:kte,its:ite)
   !    real :: dbyo_x      (kts:kte,its:ite)
   !    real :: tn_x        (kts:kte,its:ite)
   !    real :: qo_x        (kts:kte,its:ite)
   !    integer :: ki, i, k
   !    real :: x_add, denom
      
   !    do ki = 1,2
   !       if(DICYCLE==2 .and. ki==2) cycle 

   !       if(ki==1) then
   !        !-- get the cloud work function for updrafts associated only with RAD + PBL
   !        tn_x = t_in + tn - tn_adv
   !        qo_x = q_in + qo - qo_adv
   !        !-- to check => aa1_radpbl=aa1
   !        !! tn_x = tn
   !        !! qo_x = qo
   !       endif

   !       if(ki==2) then
   !          !-- get the cloud work function for updrafts associated only with Qv-advection
   !         !tn_x = tn_adv  ! orig 
   !         tn_x = t_in        ! v2
   !         qo_x = qo_adv
   !       endif

   !       ierr_dummy=ierr

   !       call cup_env     (zo,qeso_x,heo_x,heso_x,tn_x,qo_x,po,z1,psur,ierr_dummy,-1,itf,ktf, its,ite, kts,kte)
   !       call cup_env_clev(tn_x,qeso_x,qo_x,heo_x,heso_x,zo,po,qeso_cup_x,qo_cup_x,heo_cup_x,us,vs   &
   !                        ,u_cup_x,v_cup_x,heso_cup_x,zo_cup_x,po_cup_x,gammao_cup_x,tn_cup_x,psur,tsur  &
   !                        ,ierr_dummy,z1,itf,ktf,its,ite, kts,kte)

   !       !--- get MSE
   !       do i=its,itf
   !          if(ierr_dummy(i) /= 0) cycle
   !          x_add = (c_alvl*zqexec(i)+c_cp*ztexec(i)) +  x_add_buoy(i)
   !          call get_cloud_bc(cumulus,ave_layer,kts,kte,ktf,xland(i),po(kts:kte,i),heo_cup_x(kts:kte,i)&
   !                           ,hkbo_x(i),k22(i),x_add,Tpert(kts:kte,i))
   !          hco_x (kts:start_level(i),i) = hkbo_x(i)

   !          do k=start_level(i)+1 , ktop(i)+1  ! mass cons option
   !             denom =(zu(k-1,i)-.5*up_massdetro (k-1,i)+up_massentro (k-1,i)) + 1.e-12
               
   !             hco_x(k,i)=(hco_x(k-1,i)*zuo(k-1,i)-.5*up_massdetro(k-1,i)*hco_x(k-1,i)+ &
   !                                                    up_massentro(k-1,i)*heo_x(k-1,i))/ denom
   !             if(k==start_level(i)+1) then
   !                   !x_add = (xlv*zqexec(i)+cp*ztexec(i)) +  x_add_buoy(i)
   !                   hco_x(k,i)= hco_x(k,i) + x_add*up_massentro(k-1,i)/denom
   !             endif
   !             !
   !             !- includes glaciation effects on HCO_X
   !             hco_x(k,i)= hco_x(k,i)+(1.-p_liq_ice(k,i))*qrco(k,i)*c_xlf
   !          enddo
   !          hco_x (ktop(i)+2:ktf,i)=heso_cup_x(ktop(i)+2:ktf,i)
   !       enddo

   !       call get_buoyancy(itf,ktf, its,ite, kts,kte,ierr_dummy,klcl,kbcon,ktop &
   !                        ,hco_x,heo_cup_x,heso_cup_x,dbyo_x,zo_cup_x)

   !       if(ki==1) &  ! RAD+PBL only
   !          call cup_up_aa0(aa1_radpbl,zo_cup_x,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x   &
   !                         ,k22,klcl,kbcon,ktop,ierr_dummy,itf,ktf,its,ite, kts,kte)
   !       !-- get AA1_ADV
   !       !aa1_adv = aa1 + aa0 - aa1_radpbl
         

   !       if(ki==2) & ! ADV of Qv only
   !       call cup_up_aa0(aa1_adv,zo_cup_x,zuo,dbyo_x,GAMMAo_CUP_x,tn_cup_x   &
   !                      ,k22,klcl,kbcon,ktop,ierr_dummy,itf,ktf,its,ite, kts,kte)
   !       !Observe that : 
   !       !aa1 ~ aa0 + (aa1_radpbl-aa0) + (aa1_adv-aa0)

   !    enddo ! ki   
   ! end subroutine BeckerEtAlClosure

   !---------------------------------------------------------------------------------------------!
   ! !
   ! !--- determine the entrainment dependent on environmental moist (here relative humidity)
   ! !--- also the controls of RH on the diurnal cycle (see Tian et al 2022 GRL)
   ! ! if(cumulus == 'deep') &
   ! !   call rh_controls(whoami_all,itf,ktf,its,ite,kts,kte,ierr,tn,po,qo,qeso,po_cup,cumulus,rh_entr_factor, &
   ! !                   rh_dicycle_fct,entr_rate_input, entr_rate ,xlons,dtime)         
   ! subroutine rh_controls(mynum,itf,ktf,its,ite,kts,kte,ierr,t,po,qo,qeso,po_cup,cumulus  &
   !                     ,rh_entr_factor,rh_dicycle_fct,entr_rate_input, entr_rate,xlons,dt)
           
   !    implicit none
   !    character *(*), intent (in)             :: cumulus
   !    integer  ,intent (in )                  :: mynum,itf,ktf, its,ite, kts,kte
   !    integer  ,intent (in )  ,dimension(:)   :: ierr
   !    real     ,intent (in )                  :: entr_rate_input,dt
   !    real     ,intent (in )  ,dimension(:)   :: xlons
   !    real     ,intent (in )  ,dimension(:,:) :: t,po,qo,po_cup,qeso
   !    real     ,intent (inout),dimension(:)   :: entr_rate 
   !    real     ,intent (inout),dimension(:)   :: rh_entr_factor,rh_dicycle_fct

   !    !--locals
   !    integer :: i,k,vtp_index
   !    real(8)  :: y,x
   !    real    :: dpg, trash, dayhr, p_start = 1000.
   !    real    ,dimension(its:ite) :: frh,dayhrr
   !    real    ,parameter :: ref_local_time = 8., ftun3=0.25
   !    logical ,parameter :: free_troposphere = .true. 
      
   ! if(moist_trigger > 1) then 
      
   !    !-- ave rh from 1000 -> 450 hPa, following Tian et al 2022 GRL.
   !    ! OR
   !    !-- ave rh from 800 -> 450 hPa accounts only for the ´free troposphere'
   !    if(free_troposphere) p_start = 800. 

   !    do vtp_index = get_num_elements(vec_ok),1,-1
   !       i = get_data_value(vec_ok,vtp_index)
         
   !       frh(i) = 0.
   !       trash  = 0.
   !       loopN:    do k=kts,ktf
   !          if( po(k,i) .gt. p_start .and. po(k,i) .lt. 450.) cycle loopN
   !          dpg=100.*(po_cup(k,i)-po_cup(k+1,i))/c_grav
   !          trash=trash+dpg
   !          frh(i)= frh(i) + (qo(k,i)/qeso(k,i))*dpg 
   !       enddo loopN
         
   !       !--average relative humidity
   !       frh(i) =   100.*frh(i)/(1.e-8+trash) ! no unit
   !       frh(i) = max(1., min(100., frh(i)))
   !       !
   !       !--- this is for the moist_trigger = 2
   !       x = dble(frh(i))          
   !       y = 9.192833D0 - 0.2529055D0*x + 0.002344832d0*x**2 &
   !         - 0.000007230408D0*x**3
         
   !       !--- if RH increases, rh_entr_factor decreases => entr decreases
   !       if(moist_trigger == 2) rh_entr_factor(i) = max(0.5,min(1.0,real(y,4)))
         
   !       !--- if RH increases, rh_entr_factor increases => entr increases
   !       if(moist_trigger == 3) rh_entr_factor(i) = max(0.1,min(1.0,frh(i)*0.01))

   !       entr_rate(i) = entr_rate(i) * rh_entr_factor(i)    
   !    enddo
   ! endif  
   ! if(rh_dicycle == 1) then 
   !    do vtp_index = get_num_elements(vec_ok),1,-1
   !       i = get_data_value(vec_ok,vtp_index)
   !        !--- local time
   !        dayhr  = (time_in / 3600. + float(itime1_in/100)+float(mod(itime1_in,100))/60.) 
   !        dayhrr(i) = mod(dayhr+xlons(i)/15.+24., 24.)

   !       if(abs ( dayhrr(i) - ref_local_time) < 1. .or. time_in < dt+1. ) &  
   !         !--- ftun3 controls the domain of the dicycle closure 
   !         !    ftun3 => 1 the closure is insensitive to the mean tropospheric RH
   !         !    ftun3 => 0 it depends on the RH, with mininum = ftun3
   !         rh_dicycle_fct(i) = ftun3 +(1. - ftun3)*&
   !                          (1.-(atan((frh(i)-55.)/10.)+atan(55.))/3.1016)
   !    enddo
   ! endif
   !    !-- code to test the atan function
   !    !  do i = 1 ,100 !relative humidity
   !    !      y = 0.25 +0.75*(1.-(atan((float(i)-60.)/10.)+atan(50.))/3.1016)/0.9523154
   !    !      print*,i,y
   !    !  enddo
      
   !    !print*,"FRH",maxval(frh),minval(frh),maxval(rh_dicycle_fct),minval(rh_dicycle_fct)
   !    !call flush(6)
   ! end subroutine rh_controls
 
!---------------------------------------------------------------------------------------------!
!    !---------------------------------------------------------------------------------------------------
!    !--- brams driver 
!    !---------------------------------------------------------------------------------------------------
!    subroutine convParGFDriver(mxp,myp,mzp,mtp,nmp, time, itime1 &
!                            ,ims,ime, jms,jme, kms,kme               &
!                            ,its,ite, jts,jte, kts,kte               &
!                            ,flip                                    &
!                            ,FSCAV                                   &
!                            ,mynum                 &
!                            ,dt                    &
!                            ,dx2d                  &
!                            ,stochastic_sig        &
!                            ,zm                    &
!                            ,zt                    &
!                            ,dm                    &

!                            ,lons                  &
!                            ,lats                  &
!                            ,aot500                &
!                            ,temp2m                &
!                            ,sflux_r               &
!                            ,sflux_t               &
!                            ,qexcp                 &
!                            ,hexcp                 &
!                            ,wlpool                &      
!                            ,topt                  &
!                            ,xland                 &
!                            ,sfc_press             &
!                            ,kpbl                  &
!                            ,tke_pbl               &
!                            ,turb_len_scale        & !NEW 26FEB2024

!                            ,col_sat               &
!                            ,u                     &
!                            ,v                     &
!                            ,w                     &
!                            ,temp                  &
!                            ,press                 &
!                            ,rvap                  &
!                            ,mp_ice                &
!                            ,mp_liq                &
!                            ,mp_cf                 &
!                            ,curr_rvap             &
!                            ,TRACER                &!-note: uses GEOS-5 data structure
!                            ,cnvcf                 & ! conv cloud fraction
!                            !---- forcings---
!                            ,buoy_exc              &
!                            ,rthften               &! gsf_t
!                            ,rqvften               &! gsf_q
!                            ,rth_advten            &!advf_t
!                            ,rthblten              &!sgsf_t
!                            ,rqvblten              &!sgsf_q
!                            !---- output ----
!                            ,conprr                &
!                            ,lightn_dens           &
!                            ,rh_dicycle_fct        &
!                            ,rthcuten              &
!                            ,rqvcuten              &
!                            ,rqccuten              &
!                            ,rnlcuten              &
!                            ,rnicuten              &
!                            ,rucuten               &
!                            ,rvcuten               &
!                            ,sub_mpqi              &
!                            ,sub_mpql              &
!                            ,sub_mpcf              &
!                            ,rbuoycuten            &
!                            ,rchemcuten            &
!                            ,revsu_gf              &
!                            ,prfil_gf              &
!                            !
!                            ,do_this_column        &
!                            ,ierr4d                &
!                            ,jmin4d                &
!                            ,klcl4d                &
!                            ,k224d                 &
!                            ,kbcon4d               &
!                            ,ktop4d                &
!                            ,kstabi4d              &
!                            ,kstabm4d              &
!                            ,cprr4d                &
!                            ,xmb4d                 &
!                            ,edt4d                 &
!                            ,pwav4d                &
!                            ,sigma4d               &
!                            ,pcup5d                &
!                            ,up_massentr5d         &
!                            ,up_massdetr5d         &
!                            ,dd_massentr5d         &
!                            ,dd_massdetr5d         &
!                            ,zup5d                 &
!                            ,zdn5d                 &
!                            ,prup5d                &
!                            ,prdn5d                &
!                            ,clwup5d               &
!                            ,tup5d                 &
!                            ,conv_cld_fr5d         &
!                            !-- for debug/diagnostic
!                            ,AA0,AA1,AA1_ADV,AA1_RADPBL,AA1_BL,AA2,AA3,AA1_CIN,TAU_BL,TAU_EC &
!                            ,VAR2d,VAR3d_aGF,VAR3d_bGF,VAR3d_cGF,VAR3d_dGF)

!       implicit none
!       !include "mpif.h"
!       !------------------------------------------------------------------------
!       integer, intent(in) :: ims,ime, jms,jme, kms,kme,    &
!          its,ite, jts,jte, kts,kte,    &
!          mynum,mzp,mxp,myp,mtp,nmp, itime1

!       real,    intent(in) :: DT, time

!       integer, intent(in), dimension(mzp) :: flip

!       real:: FSCAV(mtp)

!       real,    dimension(kts:kte,its:ite,jts:jte), intent(in)  ::       &
!          zm,        &
!          zt,        &
!          u,         &
!          v,         &
!          w,         &
!          rvap,      &
!          temp,      &
!          press,     &
!          dm,        &
!          curr_rvap, &
!          buoy_exc,  &
!          cnvcf,     &
!          qexcp,     &
!          hexcp 

!       integer, dimension(its:ite,jts:jte), intent(in) :: kpbl

!       !-- intent (in)
!       real,    dimension(its:ite,jts:jte) :: topt ,aot500 ,temp2m ,sfc_press &
!                                             ,sflux_r ,sflux_t,xland,lons,lats &
!                                             ,dx2d,col_sat,stochastic_sig,tke_pbl,wlpool

!       real,    dimension(kts:kte,its:ite,jts:jte), intent(in) :: rthften    &
!          ,rqvften    &
!          ,rth_advten &
!          ,rthblten   &
!          ,rqvblten   &
!          ,turb_len_scale

!       real,    dimension(its:ite,jts:jte),         intent(out  ) ::   CONPRR,LIGHTN_DENS
!       real,    dimension(its:ite,jts:jte),         intent(inout) ::   rh_dicycle_fct

!       real,    dimension(kts:kte,its:ite,jts:jte), intent(out) ::     rthcuten   &
!          ,rqvcuten   &
!          ,rqccuten   &
!          ,rnlcuten   &
!          ,rnicuten   &
!          ,rucuten    &
!          ,rvcuten    &
!          ,rbuoycuten &
!          ,revsu_gf   &
!          ,prfil_gf   &
!          ,var3d_agf  &
!          ,var3d_bgf  &
!          ,var3d_cgf  &
!          ,var3d_dgf


!       real,    dimension(nmp,kts:kte,its:ite,jts:jte), intent(in)  :: mp_ice  &
!          ,mp_liq  &
!          ,mp_cf

!       real,    dimension(nmp,kts:kte,its:ite,jts:jte), intent(out) :: sub_mpqi   &
!          ,sub_mpql   &
!          ,sub_mpcf

!       !-***** TRACER has different data structure   (i,j,k,ispc) *********
!       real,    dimension(its:ite,jts:jte,kts:kte,mtp), intent(in )  :: TRACER
!       !-***** rchemcuten uses the GF data structure (ispc,k,i,j) *********
!       real,    dimension(mtp,kts:kte,its:ite,jts:jte), intent(out)  :: rchemcuten

!       integer, dimension(its:ite,jts:jte), intent(inout) :: do_this_column

!       !- for convective transport and cloud/radiation (OUT)
!       integer,dimension(mxp,myp,maxiens)::  &
!          !  integer, dimension(its:ite,jts:jte,maxiens) , INTENT(OUT) ::    &
!           ierr4d                    &
!          ,jmin4d                    &
!          ,klcl4d                    &
!          ,k224d                     &
!          ,kbcon4d                   &
!          ,ktop4d                    &
!          ,kstabi4d                  &
!          ,kstabm4d

!       real,dimension(mxp,myp,maxiens)::     &
!          !   real,dimension(its:ite,jts:jte,maxiens)    , INTENT(OUT) ::    &
!           cprr4d                    &
!          ,xmb4d                     &
!          ,edt4d                     &
!          ,pwav4d                    &
!          ,sigma4d
!       real,dimension(mzp,mxp,myp,maxiens):: &
!          !   real,dimension(its:ite,jts:jte,kts:kte,maxiens), INTENT(OUT) ::    &
!           pcup5d                    &
!          ,up_massentr5d             &
!          ,up_massdetr5d             &
!          ,dd_massentr5d             &
!          ,dd_massdetr5d             &
!          ,zup5d                     &
!          ,zdn5d                     &
!          ,prup5d                    &
!          ,prdn5d                    &
!          ,clwup5d                   &
!          ,tup5d                     &
!          ,conv_cld_fr5d
!       !--for debug
!       real   ,dimension(mxp,myp)  ,intent(inout)  :: aa0,aa1,aa1_adv,aa1_radpbl    &
!                                                     ,aa2,aa3,aa1_bl,aa1_cin,tau_bl &
!                                                     ,tau_ec,var2d

!       !----------------------------------------------------------------------
!       ! LOCAL VARS

!       ! basic environmental input includes
!       ! omega (omeg), windspeed (us,vs), and a flag (aaeq) to turn off
!       ! convection for this call only and at that particular gridpoint
!       ! 
!       real,   dimension (its:ite,jts:jte) ::  rtgt

!       real,   dimension (kts:kte,its:ite) ::  &
!                                               zo,temp_old,qv_old,po,us,vs,rhoi,phil            &
!                                              ,temp_new_dp,qv_new_dp,temp_new_sh,qv_new_sh,z2d  &
!                                              ,tkeg,rcpg,dhdt,temp_new_md,qv_new_md             &
!                                              ,temp_new_bl,qv_new_bl,dm2d,temp_tendqv,qv_curr   &
!                                              ,buoy_exc2d,revsu_gf_2d,prfil_gf_2d,var3d_agf_2d  &
!                                              ,var3d_bgf_2d,temp_new,qv_new                     &
!                                              ,temp_new_adv,qv_new_adv,cnvcf2d,turb_len_scale2d &
!                                              ,piexner 

!       real,   dimension (kts:kte,its:ite,maxiens) ::  outt,outq,outqc,outu,outv,outbuoy &
!                                                      ,outnliq,outnice
!       real,   dimension (kts:kte,its:ite,maxiens) ::  subten_Q,subten_T,subten_U,subten_v

!       real,   dimension (mtp,kts:kte,its:ite)         :: se_chem
!       real,   dimension (mtp,kts:kte,its:ite,maxiens) :: out_chem

!       real,   dimension (nmp,kts:kte,its:ite)         :: mpqi,mpql,mpcf
!       real,   dimension (nmp,kts:kte,its:ite,maxiens) :: outmpqi,outmpql,outmpcf

!       real,   dimension (its:ite)   :: ter11, xlandi,pbl,zws,ccn,psur             &
!                                       ,ztexec,zqexec,h_sfc_flux,le_sfc_flux,tsur  &
!                                       ,xlons,xlats,fixout_qv,cum_ztexec,cum_zqexec&
!                                       ,zlcl_sfc,plcl_sfc,tlcl_sfc

!       real,   dimension (kts:kte,its:ite,1:ens4)      ::  omeg

!       real,   dimension (kts:kte) :: min_tend,distance
!       integer,dimension (its:ite) :: kpbli

!       integer :: i,j,k,kr,n,itf,jtf,ktf,ispc,zmax,status

!       real :: dp,dq, dtdt,pten,pqen,paph,zrho,pahfs,pqhfl,zkhvfl,pgeoh
!       real :: fixouts,dt_inv

!       real,   dimension(mxp,myp,-1:5) :: dummy_precip
!       integer :: imemory,irun,jlx,kk,kss,plume,ii_plume

!       !----------------------------------------------------------------------
!       !-do not change this
!       itf=ite
!       ktf=kte-1
!       jtf=jte
!       int_time   = int_time + dt
!       WHOAMI_ALL = mynum
!       time_in    = time
!       itime1_in  = itime1
!       !----------------------------------------------------------------------

!       !-- big loop over j dimension
!       do j = jts,jtf
!          JCOL = J
!           conprr     (:,j) = 0.0
!           lightn_dens(:,j) = 0.0
!           var2d      (:,j) = 0.0

!          !-- initialization
!          rtgt(:,:) = 1.0
!          ztexec   (:) = 0.0
!          zqexec   (:) = 0.0
!          fixout_qv(:) = 1.0
!             !
!             !--- (k,i)
!          revsu_gf_2d (:,:) = 0.0
!          prfil_gf_2d (:,:) = 0.0
!          var3d_agf_2d(:,:) = 0.0
!          var3d_bgf_2d(:,:) = 0.0
!          temp_tendqv (:,:) = 0.0
!             !
!          omeg   (:,:,:) = 0.0
!          !- tendencies (w/ maxiens)
!          outt   (:,:,:) = 0.0
!          outu   (:,:,:) = 0.0
!          outv   (:,:,:) = 0.0
!          outq   (:,:,:) = 0.0
!          outqc  (:,:,:) = 0.0
!          outnice(:,:,:) = 0.0
!          outnliq(:,:,:) = 0.0
!          outbuoy(:,:,:) = 0.0

!          if(APPLY_SUB_MP == 1) then
!             !- tendencies (w/ nmp and maxiens)
!                outmpqi(:,:,:,:) = 0.0
!                outmpql(:,:,:,:) = 0.0
!                outmpcf(:,:,:,:) = 0.0
!          endif

!          if(USE_TRACER_TRANSP==1) then
!             out_chem(:,:,:,:) = 0.0
!          endif
!          !
!          if(autoconv == 2) then
!             do i= its,itf
!                ccn(i) = max( 100., ( 370.37*(0.01+MAX(0.,aot500(i,j))))**1.555 )
!             enddo
!          else
!             do i= its,itf
!                ccn(i) = 100.
!             enddo
!          endif

!          do i=its,itf

!             xlandi(i) = xland(i,j)!flag < 1 para land
!                                   !flag  =1 para water
!             psur  (i) = sfc_press(i,j)*1.e-2 ! mbar
!             tsur  (i) = temp2m(i,j)
!             ter11 (i) = max(0.,topt(i,j))
!             kpbli (i) = kpbl(i,j)
!             xlons (i) = lons(i,j)*180./3.14159
!             xlats (i) = lats(i,j)*180./3.14159
!          enddo

!          do i=its,itf
!             do k=kts,ktf
!                kr=k   !+1   !<<<< only kr=k (the input was already converted to the BRAMS vertical grid,
!                             !                see cup_grell3.f90 routine)

!                !- heigths, current pressure, temp and water vapor mix ratio
!                zo      (k,i)  = zt(kr,i,j)*rtgt(i,j)+topt(i,j)
!                po      (k,i)  = press(kr,i,j)*1.e-2 !mbar
!                temp_old(k,i)  = temp(kr,i,j)

!                piexner (k,i)  = 1. ! for brams, change this later for the real exner function
!                qv_old  (k,i)  = rvap     (kr,i,j) ! @ begin of the timestep
!                qv_curr (k,i)  = curr_rvap(kr,i,j) ! current (after dynamics + physical processes called before GF)

!                !- air density, TKE and cloud liq water mixing ratio
!                rhoi    (k,i)  = 1.e2*po (k,i)/( 287.04*temp_old(k,i)*(1.+0.608*qv_old(k,i)))
!                tkeg    (k,i)  = 1.e-5
!                rcpg    (k,i)  = 0.
              
!                !- cloud fraction 
!                cnvcf2d  (k,i)   =  cnvcf(kr,i,j)

!                !- wind velocities
!                us      (k,i)  =  u (kr,i,j)
!                vs      (k,i)  =  v (kr,i,j)
!                dm2d    (k,i)  =  dm(kr,i,j)
!                omeg    (k,i,:)= -c_grav*rhoi(k,i)*w(kr,i,j)
!                !-buoyancy excess
!                buoy_exc2d(k,i)= buoy_exc(kr,i,j)
!                !- temp/water vapor modified only by advection
!                temp_new_ADV(k,i)= temp_old(k,i)  +  (rth_advten(kr,i,j) )*dt
!                qv_new_ADV  (k,i)=   qv_old(k,i)  +  (rqvften   (kr,i,j) )*dt
!                !-- turb length scale
!                turb_len_scale2d (k,i)  = turb_len_scale (kr,i,j)
!             enddo
!          enddo

!          if(APPLY_SUB_MP == 1) then
!             do i=its,itf
!                do k=kts,ktf
!                   kr=k   !+1   !<<<< only kr=k
!                   !- microphysics ice and liq mixing ratio, and cloud fraction of the host model
!                   !- (only subsidence is applied)
!                   mpqi   (:,k,i) = mp_ice  (:,kr,i,j) ! kg/kg
!                   mpql   (:,k,i) = mp_liq  (:,kr,i,j) ! kg/kg
!                   mpcf   (:,k,i) = mp_cf   (:,kr,i,j) ! 1
!                enddo
!             enddo
!          endif
!          if(USE_TRACER_TRANSP==1) then
!             do i=its,itf
!                do k=kts,kte
!                   kr=k !+1
!                   !- atmos composition
!                   do ispc=1,mtp
!                      se_chem(ispc,k,i) = max(p_mintracer, TRACER(i,j,flip(kr),ispc))
!                   enddo
!                enddo
!             enddo
!          endif
!          !- pbl  (i) = depth of pbl layer (m)
!          !- kpbli(i) = index of zo(k,i)
!          do i=its,itf
!             pbl  (i)  = zo(kpbli(i),i) - topt(i,j)
!          enddo

!          !- begin: for GATE soundings-------------------------------------------
!          !- this section is intended for model developments only and must
!          !- not be used for normal runs.
!          if(P_USE_GATE) then
!             if(CLEV_GRID == 0) stop "use_gate requires CLEV_GRID 1 or 2"
!             if(USE_TRACER_TRANSP==1) then
!                ispc_CO=1
!                if( .not. allocated(Hcts)) allocate(Hcts(mtp))
!                CHEM_NAME_MASK (:) = 1
!                !--- dummy initization FSCAV
!                do i=1,mtp
!                   !FSCAV(i) = 0.1  !km^-1

!                   FSCAV(i) = 1.e-5  !km^-1
!                   Hcts(i)%hstar  = 0.0 !8.300e+4! 2.4E+3 !59.
!                   Hcts(i)%dhr    = 0.0 !7400.   !5000.  !4200.
!                   Hcts(i)%ak0    = 0.0
!                   Hcts(i)%dak    = 0.0
!                    ! H2O2      0.00000      8.300e+4    7400.00000       0.00000       0.00000
!                    ! HNO3      0.00000      2.100e+5    8700.00000       0.00000       0.00000
!                    ! NH3       0.00000      59.00000    4200.00000       0.00000       0.00000
!                    ! SO2       0.00000      2.400e+3    5000.00000       0.00000       0.00000
!                enddo
!                do i=its,itf
!                   se_chem(1:mtp,kts:kpbli(i)-1,i) = 1.+1.e-6
!                   do k=kpbli(i),kte
!                      se_chem(1:mtp,k,i) = 1.*exp(-(max(0.,0.9*float(k-kpbli(i)))/float(kpbli(i))))+1.e-6
!                   enddo
!                   do k=kts+1,kte-1
!                      se_chem(1:mtp,k,i) = 1./3. *( se_chem(1:mtp,k,i) + se_chem(1:mtp,k-1,i) + se_chem(1:mtp,k+1,i))
!                   enddo
!                enddo
!             endif

!             !--- only for GATE soundingg
!             if(trim(RUNDATA) == "GATE.dat") then
!                jlx= jl
!                !jlx= 10 ! to run with only one soundings
!                !buoy_exc2d(:,:) = float(jl)*30.
!                !print*,"GATE",jl,jlx,buoy_exc2d(1,1)
!                !jlx= 42 ! to run with only one soundings

!                do i=its,itf
!                   do k=kts,kte
!                      po       (k,i) = 0.5*(ppres(jlx,k)+ppres(jlx,min(kte,k+1)))
!                      temp_old (k,i) = ptemp(jlx,k)+273.15
!                      qv_old   (k,i) = pq(jlx,k)/1000.
!                      us       (k,i) = pu(jlx,k)
!                      vs       (k,i) = pv(jlx,k)
!                      omeg     (k,i,:)=pvervel(jlx,k)
!                      phil     (k,i) = pgeo(jlx,k)*c_grav   !geo
!                      rhoi     (k,i) = 1.e2*po(k,i)/(c_rgas*temp_old(k,i))
!                   enddo

!                   do k=kts,kte
!                      mpql     (:,k,i) = 0.
!                      mpql     (:,k,i) = 0.
!                      mpcf     (:,k,i) = 0.
!                      if(po(k,i) > 900. .or. po(k,i)<300.) cycle
!                      pqen  =  exp((-3.e-5*(po(k,i)-550.)**2))
!                      pten  =  min(1., (max(0.,(temp_old(k,i)-c_Tice))/(c_T00-c_Tice))**2)
!                      mpql  (:,k,i) =3.*pqen* pten
!                      mpqi  (:,k,i) =3.*pqen*(1.- pten)
!                      mpcf  (:,k,i) = (mpqi  (:,k,i)+mpql  (:,k,i))*100.
!                   enddo

!                   do k=kts,kte
!                      zo       (k,i) = 0.5*(phil(k,i)+phil(min(kte,k+1),i))/c_grav    !meters
!                   enddo
!                   ter11(i)  = phil(1,i)/c_grav  ! phil is given in g*h.
!                   psur (i)  = ppres(jlx,1)
!                   tsur (i)  = temp2m(i,j) !temp_old(i,1)
!                   kpbli(i)  = 5
!                   pbl  (i)  = zo(kpbli(i),i)
!                   zws  (i)  = 1.0 ! wstar
!                   do k=kts,ktf
!                      temp_new(k,i) = temp_old(k,i) + dt *(zadvt(jlx,k)+zqr(jlx,k))/86400.
!                      qv_new  (k,i) = qv_old  (k,i) + dt * zadvq(jlx,k)

!                      temp_new_dp (k,i) = temp_old(k,i) + dt *(zadvt(jlx,k)+zqr(jlx,k))/86400.
!                      qv_new_dp   (k,i) = qv_old  (k,i) + dt * zadvq(jlx,k)

!                      temp_new_md (k,i) = temp_new_dp(k,i)
!                      qv_new_md   (k,i) = qv_new_dp  (k,i)
!                      temp_new_bl (k,i) = temp_new_dp(k,i)
!                      qv_new_bl   (k,i) = qv_new_dp  (k,i)
!                      temp_new_adv(k,i) = temp_old   (k,i) + dt * zadvt(jlx,k)/86400.
!                      qv_new_adv  (k,i) = qv_old     (k,i) + dt * zadvq(jlx,k)
!                   enddo
!                enddo
!             endif
!          endif !- end:   for GATE soundings-------------------------------------------
!          !
!          !- get excess T and Q for source air parcels
!          do i=its,itf
!             pten = temp_old(1,i)
!             pqen = qv_old  (1,i)
!             paph = 100.*psur(i)
!             zrho = paph/(287.04*(temp_old(1,i)*(1.+0.608*qv_old(1,i))))
!             !- sensible and latent sfc fluxes for the heat-engine closure
!             h_sfc_flux (i)=zrho*c_cp  *sflux_t(i,j)!W/m^2
!             le_sfc_flux(i)=zrho*c_alvl*sflux_r(i,j)!W/m^2
!             !
!             !- local le and h fluxes for W*
!             pahfs=-sflux_t(i,j) *zrho*1004.64  !W/m^2
!             pqhfl=-sflux_r(i,j)                !kg/m^2/s
!             !- buoyancy flux (h+le)
!             zkhvfl= (pahfs/1004.64+0.608*pten*pqhfl)/zrho ! K m s-1
!             !- depth of 1st model layer
!             !- (zo(1)-top is ~ 1/2 of the depth of 1st model layer, => mult by 2)
!             pgeoh =  2.*( zo(1,i)-topt(i,j) )*c_grav ! m+2 s-2
!             !-convective-scale velocity w*
!             !- in the future, change 0.001 by ustar^3
!             zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten) ! m+3 s-3
!             !
!             !-- get LCL properties for parcels from surface level
!             call calc_lcl(pten,paph,pqen,tlcl_sfc(i),plcl_sfc(i),zlcl_sfc(i))
!             zlcl_sfc(i) = max(zlcl_sfc(i), 0.)
!             !print*,'lcl',minval(zlcl_sfc(its:itf)),maxval(zlcl_sfc(its:itf))
!             !
!             if(zws(i) > tiny(pgeoh)) then
!                !-convective-scale velocity w*
!                zws(i) = 1.2*zws(i)**.3333
!                !- temperature excess
!                ztexec(i)     = max(0.,-1.5*pahfs/(zrho*zws(i)*1004.64)) ! K
!                !print*,"exce1=",pahfs,zrho,ztexec(i),zws(i),pgeoh,zo(1,i),topt(i,j)
!                !call flush(6)
!                !- moisture  excess
!                zqexec(i)     = max(0.,-1.5*pqhfl/(zrho*zws(i)))        !kg kg-1
!             endif   ! zws > 0
!             !
!             !- zws for shallow convection closure (Grant 2001)
!             !- depth of the pbl
!             pgeoh = pbl(i)*c_grav
!             !-convective-scale velocity W* (m/s)
!             zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten)
!             zws(i) = 1.2*zws(i)**.3333
!          enddo
         
!          !
!          !------ CALL CUMULUS PARAMETERIZATION
!          !

!          do ii_plume = 1, maxiens

!             if(ii_plume == 1) then
!                plume = shal
!                c0 = c0_shal
!             endif
!             if(ii_plume == 2) then
!                plume = deep
!                c0 = c0_deep
!             endif
!             if(ii_plume == 3) then
!                plume = mid
!                c0 = c0_mid
!             endif

!             if(icumulus_gf(plume) == OFF ) cycle
            
!             use_smooth_tend=  cum_use_smooth_tend(plume)
!             hei_down_land  =  cum_hei_down_land  (plume)
!             hei_down_ocean =  cum_hei_down_ocean (plume)
!             hei_updf_land  =  cum_hei_updf_land  (plume)
!             hei_updf_ocean =  cum_hei_updf_ocean (plume)
!             max_edt_land   =  cum_max_edt_land   (plume)
!             max_edt_ocean  =  cum_max_edt_ocean  (plume)
!             fadj_massflx   =  cum_fadj_massflx   (plume)
!             use_excess     =  cum_use_excess     (plume)
!             ave_layer      =  cum_ave_layer      (plume)
!             T_star         =  cum_t_star         (plume)
!             fr_min_entr    =  cum_fr_min_entr    (plume)
!             !print*,"plume=",plume,shal,mid,deep

!             !-- set minimum/max for excess of T and Q
!             if(use_excess == 0) then
!                cum_ztexec(:)= 0.
!                cum_zqexec(:)= 0.
!             elseif (use_excess == 1) then
!                cum_ztexec(:)= ztexec(:)
!                cum_zqexec(:)= zqexec(:)
!             elseif (use_excess == 2) then
!                do i=its,itf
!                   cum_zqexec(i)=min(5.e-4, max(1.e-4,zqexec(i)))! kg kg^-1
!                   cum_ztexec(i)=min(0.5,   max(0.2  ,ztexec(i)))! Kelvin
!                enddo
!             else
!                do i=its,itf
!                   if(xlandi(i) > 0.98) then ! ocean
!                      cum_zqexec(i)=min(8.e-4, max(5.e-4,zqexec(i)))! kg kg^-1
!                      cum_ztexec(i)=min(1.,    max(0.5  ,ztexec(i)))! Kelvin
!                   else                      ! land
!                      cum_ztexec(i)= ztexec(i)
!                      cum_zqexec(i)= zqexec(i)
!                   endif
!                enddo
!             endif
!             !
!             !--- replace 'q' and 't' excess in case of use of the cold pool scheme
!             !
!             if(convection_tracer == 1 .and. plume == deep) then           
!                  if(use_gustiness == 1 ) then
!                     cum_ztexec(:)= ztexec(:)
!                     cum_zqexec(:)= zqexec(:)
!                   endif

!             !   if(use_gustiness == 1 ) then           
!             !      k = 2 ! surface in brams
!             !      do i=its,itf
!             !       cum_ztexec(i)= (hexcp(k,i,j) - qexcp(k,i,j)*c_alvl)/cp
!             !       cum_zqexec(i)= qexcp(k,i,j)
!             !      enddo
!             !   else 
!             !       cum_ztexec(:)= 0.
!             !       cum_zqexec(:)= 0.
!             !   endif
!             endif
!             !
!             !
!             !-- shallow convection
!             !
!             if(plume == shal) then
!                do i=its,itf
!                   do k=kts,ktf
!                      kr=k!+1 <<<<
!                      if(p_use_gate) then
!                         dhdt    (k,i)= c_cp*(temp_new_dp(k,i)-temp_old(k,i))+c_alvl*(qv_new_dp(k,i)-qv_old(k,i))
!                         temp_new(k,i)= temp_new_dp(k,i)
!                         qv_new  (k,i)= qv_new_dp  (k,i)
!                      else

!                         temp_new(k,i)=temp_old(k,i) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
!                         qv_new  (k,i)=  qv_old(k,i) + (rqvblten(kr,i,j)+rqvften(kr,i,j))*dt
!                         qv_new  (k,i)= max(p_smaller_qv,qv_new  (k,i))

!                         !- only pbl forcing changes moist static energy
!                         dhdt(k,i)= c_cp   *(rthblten(kr,i,j)) +  c_alvl  *(rqvblten(kr,i,j))

!                         !- all forcings change moist static energy
!                         dhdt(k,i)=dhdt(k,i) + c_cp*rthften(kr,i,j) + c_alvl*rqvften(kr,i,j)

!                      endif
!                   enddo
!                enddo
!             endif
!             !
!             !--- deep convection
!             if(plume == deep) then

!                if(p_use_gate) then
!                   do k=kts,ktf
!                      do i=its,itf
!                         temp_new(k,i) = temp_new_dp(k,i)
!                         qv_new  (k,i) = qv_new_dp  (k,i)
!                      enddo
!                   enddo
!                else
!                   do k=kts,ktf
!                      do i=its,itf
!                         kr=k!+1 <<<<
!                         temp_new    (k,i)= temp_old(k,i)  +  (rthblten(kr,i,j) + rthften (kr,i,j))*dt
!                         qv_new      (k,i)=   qv_old(k,i)  +  (rqvblten(kr,i,j) + rqvften (kr,i,j))*dt

!                         temp_new_BL (k,i)= temp_old(k,i)  +  (rthblten(kr,i,j) )*dt
!                         qv_new_BL   (k,i)=   qv_old(k,i)  +  (rqvblten(kr,i,j) )*dt
!                      enddo
!                   enddo
!                endif
!             endif
!             !
!             !--- mid/congestus type convection
!             if(plume == mid) then

!                if(p_use_gate) then
!                   do k=kts,ktf
!                      do i=its,itf
!                         temp_new(k,i) = temp_new_dp(k,i)
!                         qv_new  (k,i) = qv_new_dp  (k,i)
!                      enddo
!                   enddo
!                else
!                   do i=its,itf
!                      do k=kts,ktf
!                         kr=k!+1 <<<<

!                         temp_new(k,i)=temp_old(k,i) + (rthblten(kr,i,j)+rthften(kr,i,j))*dt
!                         qv_new  (k,i)=  qv_old(k,i) + (rqvblten(kr,i,j)+rqvften(kr,i,j))*dt
!                         qv_new  (k,i)= max(p_smaller_qv,qv_new  (k,i))

!                         !- only pbl forcing changes moist static energy
!                         dhdt(k,i)= c_cp   *(rthblten(kr,i,j)) +  c_alvl  *(rqvblten(kr,i,j))

!                         !- all forcings change moist static energy
!                         dhdt(k,i)=dhdt(k,i) + c_cp*rthften(kr,i,j) + c_alvl*rqvften(kr,i,j)

!                         !- temp/water vapor modified only by bl processes
!                         temp_new_BL(k,i)= temp_old(k,i)  +  (rthblten(kr,i,j) )*dt
!                         qv_new_BL  (k,i)= qv_old  (k,i)  +  (rqvblten(kr,i,j) )*dt

!                      enddo
!                   enddo
!                endif
!             endif
!             !

!             call CUP_GF(its,ite,kts,kte, itf,ktf, mtp, nmp, FSCAV  &
!                         ,cumulus_type  (plume)            &
!                         ,closure_choice(plume)            &
!                         ,cum_entr_rate (plume)            &
!                         ,cum_use_excess(plume)            &
!                         !- input data
!                         ,dx2d          (:,j)              &
!                         ,stochastic_sig(:,j)              &
!                         ,col_sat       (:,j)              &
!                         ,tke_pbl       (:,j)              &
!                         ,rh_dicycle_fct(:,j)              &
!                         ,wlpool        (:,j)              &
!                         ,dt                               &
!                         ,kpbli                            &
!                         ,cum_ztexec                       &
!                         ,cum_zqexec                       &
!                         ,ccn                              &
!                         ,rhoi                             &
!                         ,omeg                             &
!                         ,temp_old                         &
!                         ,qv_old                           &
!                         ,ter11                            &
!                         , h_sfc_flux                      &
!                         ,le_sfc_flux                      &
!                         ,zlcl_sfc                         &
!                         ,xlons                            &
!                         ,xlats                            &
!                         ,xlandi                           &
!                         ,temp_new                         &
!                         ,qv_new                           &
!                         ,temp_new_BL                      &
!                         ,qv_new_BL                        &
!                         ,temp_new_ADV                     &
!                         ,qv_new_ADV                       &
!                         ,zo                               &
!                         ,po                               &
!                         ,piexner                          &
!                         ,tsur                             &
!                         ,psur                             &
!                         ,us                               &
!                         ,vs                               &
!                         ,dm2d                             &
!                         ,se_chem                          &
!                         ,zws                              &
!                         ,dhdt                             &
!                         ,buoy_exc2d                       &
!                         ,cnvcf2d                          &
!                         ,turb_len_scale2d                 &
!                         ,mpqi                             &
!                         ,mpql                             &
!                         ,mpcf                             &
!                         !output data
!                         ,outt                 (:,:,plume) &
!                         ,outq                 (:,:,plume) &
!                         ,outqc                (:,:,plume) &
!                         ,outu                 (:,:,plume) &
!                         ,outv                 (:,:,plume) &
!                         ,subten_Q             (:,:,plume) &
!                         ,subten_T             (:,:,plume) &
!                         ,subten_U             (:,:,plume) &
!                         ,subten_V             (:,:,plume) &
!                         ,outnliq              (:,:,plume) &
!                         ,outnice              (:,:,plume) &
!                         ,outbuoy              (:,:,plume) &
!                         ,outmpqi            (:,:,:,plume) &
!                         ,outmpql            (:,:,:,plume) &
!                         ,outmpcf            (:,:,:,plume) &
!                         ,out_chem           (:,:,:,plume) &
!                         !- for convective transport
!                         ,ierr4d               (:,j,plume) &
!                         ,jmin4d               (:,j,plume) &
!                         ,klcl4d               (:,j,plume) &
!                         ,k224d                (:,j,plume) &
!                         ,kbcon4d              (:,j,plume) &
!                         ,ktop4d               (:,j,plume) &
!                         ,kstabi4d             (:,j,plume) &
!                         ,kstabm4d             (:,j,plume) &
!                         ,cprr4d               (:,j,plume) &
!                         ,xmb4d                (:,j,plume) &
!                         ,edt4d                (:,j,plume) &
!                         ,pwav4d               (:,j,plume) &
!                         ,sigma4d              (:,j,plume) &
!                         ,pcup5d             (:,:,j,plume) &
!                         ,up_massentr5d      (:,:,j,plume) &
!                         ,up_massdetr5d      (:,:,j,plume) &
!                         ,dd_massentr5d      (:,:,j,plume) &
!                         ,dd_massdetr5d      (:,:,j,plume) &
!                         ,zup5d              (:,:,j,plume) &
!                         ,zdn5d              (:,:,j,plume) &
!                         ,prup5d             (:,:,j,plume) &
!                         ,prdn5d             (:,:,j,plume) &
!                         ,clwup5d            (:,:,j,plume) &
!                         ,tup5d              (:,:,j,plume) &
!                         ,conv_cld_fr5d      (:,:,j,plume) &
!                         !-- for debug/diag
!                         ,AA0(:,j),AA1(:,j),AA1_ADV(:,j),AA1_RADPBL(:,j),AA2(:,j),AA3(:,j) &
!                         ,AA1_BL(:,j),AA1_CIN(:,j),TAU_BL(:,j),TAU_EC(:,j) &
!                         !-- for diag
!                         ,lightn_dens  (:,j)               &
!                         ,var2d        (:,j)               &
!                         ,revsu_gf_2d                      &
!                         ,prfil_gf_2d                      &
!                         ,var3d_agf_2d                     &
!                         ,var3d_bgf_2d                     &
!                         )
!          enddo !- plume

!          !--- reset ierr4d to value different of zero in case the correspondent
!          !--- plume (shalllow, congestus, deep) was not actually used
!          do n=1,maxiens
!             if(icumulus_gf(n) == OFF ) ierr4d (:,j,n) = -99
!          enddo

!          do i=its,itf
!             do_this_column(i,j) = 0
!             loop1:  do n=1,maxiens
!                if(ierr4d (i,j,n) == 0 ) then
!                   do_this_column(i,j) = 1
!                   exit loop1
!                endif
!             enddo loop1
!          enddo
!          !----------- check for negative water vapor mix ratio
!          do i=its,itf
!             if(do_this_column(i,j) == 0) cycle
!             do k = kts,ktf
!                temp_tendqv(k,i)= outq (k,i,shal) + outq (k,i,deep) + outq (k,i,mid )
!             enddo

!             do k = kts,ktf
!                distance(k)= qv_curr(k,i) + temp_tendqv(k,i) * dt
!             enddo

!             if(minval(distance(kts:ktf)) < 0.) then
!                zmax   =  MINLOC(distance(kts:ktf),1)

!                if( abs(temp_tendqv(zmax,i) * dt) <  p_mintracer) then
!                   fixout_qv(i)= 0.999999
!                 !fixout_qv(i)= 0.
!                else
!                   fixout_qv(i)= ( (p_smaller_qv - qv_curr(zmax,i))) / (temp_tendqv(zmax,i) *dt)
!                endif
!                fixout_qv(i)=max(0.,min(fixout_qv(i),1.))
!             endif
!          enddo
!          !------------ feedback
!          !-- deep convection
!          do i=its,itf
!             if(do_this_column(i,j) == 0) cycle
!             cprr4d(i,j,deep) =  cprr4d(i,j,deep)* fixout_qv(i)
!             cprr4d(i,j,mid)  =  cprr4d(i,j,mid) * fixout_qv(i)
!             cprr4d(i,j,shal) =  cprr4d(i,j,shal)* fixout_qv(i)
!             conprr(i,j)= (cprr4d(i,j,deep) + cprr4d(i,j,mid) + cprr4d(i,j,shal))
!             conprr(i,j)= max(0.,conprr(i,j))
!          enddo

!          !-- deep + shallow + mid convection
!          do i = its,itf
!             if(do_this_column(i,j) == 0) cycle
!             do k = kts,kte
!                kr=k!+1
!                !- feedback the tendencies from convection
!                RTHCUTEN (kr,i,j)= (outt (k,i,shal) + outt (k,i,deep) + outt (k,i,mid)) *fixout_qv(i)

!                RQVCUTEN (kr,i,j)= (outq (k,i,shal) + outq (k,i,deep) + outq (k,i,mid)) *fixout_qv(i)

!                RQCCUTEN (kr,i,j)= (outqc(k,i,shal) + outqc(k,i,deep) + outqc(k,i,mid)) *fixout_qv(i)

!                REVSU_GF (kr,i,j)= revsu_gf_2d(k,i)*fixout_qv(i) !-- already contains deep and mid amounts.

!                !---these arrays are only for the deep plume mode
!                PRFIL_GF (kr,i,j)= prfil_gf_2d (k,i)*fixout_qv(i) !-- ice/liq prec flux of the deep plume
!                !VAR3d_aGF(kr,i,j)= var3d_gf_2d(k,i)              !-- vertical velocity of the deep plume
!                VAR3d_aGF(kr,i,j)= outt (k,i,mid)*fixout_qv(i)   !--
!                VAR3d_bGF(kr,i,j)= outq (k,i,mid)*fixout_qv(i)   !--

!                if(icumulus_gf(shal) == OFF) then
!                   VAR3d_cGF(kr,i,j)= outqc (k,i,deep)*fixout_qv(i)  !--
!                   VAR3d_dGF(kr,i,j)= outqc (k,i,mid) *fixout_qv(i)  !--
!                else
!                   VAR3d_cGF(kr,i,j)= outt (k,i,shal)*fixout_qv(i)   !--
!                   VAR3d_dGF(kr,i,j)= outq (k,i,shal)*fixout_qv(i)   !--
!                endif

!             enddo
!          enddo
!          if(USE_MOMENTUM_TRANSP > 0) then
!             do i = its,itf
!                if(do_this_column(i,j) == 0) cycle
!                do k = kts,kte
!                   kr=k!+1
!                   RUCUTEN (kr,i,j)= (outU(k,i,deep)+outU(k,i,mid)+outU(k,i,shal)) *fixout_qv(i)
!                   RVCUTEN (kr,i,j)= (outV(k,i,deep)+outV(k,i,mid)+outV(k,i,shal)) *fixout_qv(i)
!                enddo
!             enddo
!          endif

!          if(APPLY_SUB_MP == 1) then
!             do i = its,itf
!                if(do_this_column(i,j) == 0) cycle
!                do k = kts,kte
!                   kr=k!+1
!                   SUB_MPQL (:,kr,i,j)= (outmpql(:,k,i,deep)+outmpql(:,k,i,mid)+outmpql(:,k,i,shal)) *fixout_qv(i)
!                   SUB_MPQI (:,kr,i,j)= (outmpqi(:,k,i,deep)+outmpqi(:,k,i,mid)+outmpqi(:,k,i,shal)) *fixout_qv(i)
!                   SUB_MPCF (:,kr,i,j)= (outmpcf(:,k,i,deep)+outmpcf(:,k,i,mid)+outmpcf(:,k,i,shal)) *fixout_qv(i)
!                enddo
!             enddo
!          endif

!          if(LIQ_ICE_NUMBER_CONC == 1) then
!             do i = its,itf
!                if(do_this_column(i,j) == 0) cycle
!                do k = kts,kte
!                   kr=k!+1
!                   RNICUTEN (kr,i,j)= (outnice(k,i,shal) + outnice(k,i,deep) + outnice(k,i,mid)) *fixout_qv(i)
!                   RNLCUTEN (kr,i,j)= (outnliq(k,i,shal) + outnliq(k,i,deep) + outnliq(k,i,mid)) *fixout_qv(i)
!                enddo
!             enddo
!          endif

!          if(USE_TRACER_TRANSP==1) then
!             do i = its,itf
!                if(do_this_column(i,j) == 0) cycle
!                do k = kts,kte
!                   kr=k!+1
!                   RCHEMCUTEN (:,kr,i,j)= (out_CHEM(:,k,i,deep) +out_CHEM(:,k,i,mid)+out_CHEM(:,k,i,shal)) *fixout_qv(i)
!                enddo
!             enddo

!             !- constrain positivity for tracers
!             do i = its,itf
!                if(do_this_column(i,j) == 0) cycle

!                do ispc=1,mtp
!                   if(CHEM_NAME_MASK (ispc) == 0 ) cycle

!                   do k=kts,ktf
!                      distance(k)= se_chem(ispc,k,i) + RCHEMCUTEN(ispc,k,i,j)* dt
!                   enddo

!                   !-- fixer for mass of tracer
!                   if(minval(distance(kts:ktf)) < 0.) then
!                      zmax   =  MINLOC(distance(kts:ktf),1)

!                      if( abs(RCHEMCUTEN(ispc,zmax,i,j)*dt) <  p_mintracer) then
!                         fixouts= 0.999999
!                       !fixouts= 0.
!                      else
!                         fixouts=  ( (p_mintracer - se_chem(ispc,i,zmax))) / (RCHEMCUTEN(ispc,zmax,i,j)*dt)
!                      endif
!                      if(fixouts > 1. .or. fixouts <0.)fixouts=0.

!                      RCHEMCUTEN(ispc,kts:ktf,i,j)=fixouts*RCHEMCUTEN(ispc,kts:ktf,i,j)
!                   endif
!                enddo
!             enddo
!          endif

!          if(CONVECTION_TRACER==1) then
!             do i = its,itf
!                if(do_this_column(i,j) == 0) cycle
!                do k = kts,kte
!                   kr=k!+1
!                   RBUOYCUTEN (kr,i,j)= (outbuoy(k,i,deep)+outbuoy(k,i,mid)+outbuoy(k,i,shal)) *fixout_qv(i)
!                enddo
!             enddo
!          endif

!   !do i = its,itf
!   ! AA1_ADV(i,j)=zt(ktop4d(i,j,deep),i,j)*rtgt(i,j)
!   ! !AA0    (i,j)=zlcl_sfc(i)
!   !enddo

!       enddo
!    end subroutine convParGFDriver

!---------------------------------------------------------------------------------------------------
end  module lib_module_cu_gf_monan
