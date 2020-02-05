module opacity_table_module

  integer         , save :: nromax, ntmax, nyemax, ngmax
  double precision, save :: romin, romax, tmin, tmaxi
  double precision, save :: yemin, yemax
  double precision, allocatable, save :: enuk(:)
  real, allocatable, save :: absopac(:,:,:,:)
  real, allocatable, save :: emis(:,:,:,:)
  real, allocatable, save :: scaopac(:,:,:,:)
  real, allocatable, save :: sdelta(:,:,:,:)
  integer    , parameter :: phi_nmom=2,phi_nfreq=12,phi_neta=30,phi_nt=30
  integer    , parameter :: nnt=30,nneta=30,nie=12,nio=12
  double precision, parameter :: phie_eta1=0.d0, phie_eta2=50.d0
  double precision, parameter :: phie_t1=-1.d0,phie_t2=1.5d0,phie_t12=-1.d0,phie_t22=1.5d0
  integer, save :: ngroups, ng0, ng1, ng2
  double precision, allocatable, save :: egroupMeV(:)
  double precision, allocatable, save :: phieetable(:,:,:,:,:)
  double precision, allocatable, save :: phiaetable(:,:,:,:,:)
  double precision, allocatable, save :: phimetable(:,:,:,:,:)
  double precision, allocatable, save :: phientable(:,:,:,:,:)
  double precision, allocatable, save :: phiantable(:,:,:,:,:)
  double precision, allocatable, save :: phimntable(:,:,:,:,:)
  double precision, allocatable, save :: phieptable(:,:,:,:,:)
  double precision, allocatable, save :: phiaptable(:,:,:,:,:)
  double precision, allocatable, save :: phimptable(:,:,:,:,:)
!  double precision, save :: freqe(phi_nfreq), freqa(phi_nfreq), freqm(phi_nfreq)

contains

  ! dummy subroutine

!     ============================================
!
!     opacity for a given temperature, density, and ye computed
!     by an interpolation of the precalculated opacity table
!
!     this is a simplified routine with all interpolations linear
!
!     input:   er - neutrino energy (mev)
!              temp   - temperature (mev)
!              yein  - electron fraction
!              rho (in state vector) - density     (g cm^-3)
!     output:  ab  - absorptive opacity  (1/cm)
!              sc  - scattering opacity  (1/cm)
!              delta - scattering anisotropy factor
!              eta   - emissivity
!
  subroutine get_opacity_emissivity( &
       kappa, sigma, delta, eta,           &
       dkappa, dsigma, ddelta, deta, &
       rho, yein, temp, er, ngr, ntot, &
       get_k, get_s, get_d, get_j, get_dk, &
       get_ds, get_dd, get_dj)

    implicit none

    integer, intent(in) :: ntot, get_k, get_s, get_d, get_j
    integer, intent(in) :: get_dk, get_ds, get_dd, get_dj
    integer, intent(in) :: ngr(3)
    double precision, intent(in) :: rho, yein, temp, er(ntot)
    double precision, intent(out) :: kappa(ntot), sigma(ntot), delta(ntot), eta(ntot)
    double precision, intent(out) :: dkappa(ntot), dsigma(ntot), ddelta(ntot), deta(ntot)
    double precision :: coff(8)
    
    integer jye, jr, jt, jf, g, grp, nu, ng

    real*8 ye

    real*8 tl, rl !, el
    !real*8 frmin, frmax
    real*8 deltaye, deltar, deltat, deltaf
    real*8 r1i, r2i, dri
    real*8 t1i, t2i, dti
    real*8 ye1i, ye2i, dyei
    !real*8 f1i, f2i, dfi
    real*8 dfr
	real*8 yesave
    

    
    ! LHH hack to go beyond high end of table
    real*8 yetrue

    yesave = yein

    if( (rho.gt.romax) .or. (rho.lt.romin) ) then
!       print*,'rho out of range ',romin,rho,romax,iflag
!       stop'done'
    endif
    
    if( (temp.gt.tmaxi) .or. (temp.lt.tmin) ) then
!       print*,'temp out of range ',tmin,temp,tmaxi,iflag
!       stop'done'
    endif

!           if( (yein.gt.yemax) .or. (yein.lt.yemin) ) then
!              print*,'ye out of range ',yemin,yein,yemax,iflag
!              stop'done'
!           endif

    if (yein .gt. 1.d0) then
       print *, 'Ye grossly out of range', yein
       stop
    endif

    tl = log(temp)
    rl = log(rho)
    
    ye = yein
    
    ! hack to hold ye within range of table, see also getgderivsYe below
    ye = max(min(ye, yemax), yemin)
    tl = max(min(tl, log(tmaxi)),log(tmin))
    rl = max(min(rl, log(romax)),log(romin))
    
!    frmin = log(enuk(1))
!    frmax = log(enuk(ngmax))
!    dfr = log(enuk(ngmax)/enuk(1))/(ngmax-1.d0)
    deltaye = (ye-yemin)/(yemax-yemin)*float(nyemax-1)
    deltar=(rl-dlog(romin))/(dlog(romax/romin))*float(nromax-1)
    deltat=(tl-dlog(tmin))/(dlog(tmaxi/tmin))*float(ntmax-1)
    jye = 1 + idint(deltaye)
    jr  = 1 + idint(deltar)
    jt  = 1 + idint(deltat)
    if(jr.lt.1) jr = 1
    if(jr.gt.(nromax-1)) jr = nromax-1
    if(jt.lt.1) jt = 1
    if(jt.gt.(ntmax-1)) jt = ntmax-1
    if(jye.le.1) jye = 1
    if(jye.gt.(nyemax-1)) jye = nyemax-1
    
    r1i=dlog(romin) + (dlog(romax)-dlog(romin))*dble(jr-1)/ &
         dble(nromax-1)
    r2i=dlog(romin) + (dlog(romax)-dlog(romin))*dble(jr)/ &
         dble(nromax-1)
    dri=(rl-r1i)/(r2i-r1i)
    if(dri .lt. 0.d0) dri = 0.d0
    
    t1i=dlog(tmin) + (dlog(tmaxi)-dlog(tmin))*dble(jt-1)/ &
         dble(ntmax-1)
    t2i=dlog(tmin) + (dlog(tmaxi)-dlog(tmin))*dble(jt)/ &
         dble(ntmax-1)
    dti=(tl-t1i)/(t2i-t1i)
    if(dti .lt. 0.d0) dti = 0.d0

    ye1i=yemin + (yemax-yemin)*dble(jye-1)/dble(nyemax-1)
    ye2i=yemin + (yemax-yemin)*dble(jye)/dble(nyemax-1)
    dyei=(yein-ye1i)/(ye2i-ye1i)

    !if(dyei .lt. 0.d0) dyei = 0.d0
    
    ! LHH hack to go beyond high end of table
    !dyei=(ye-ye1i)/(ye2i-ye1i)

    !write(*,*) rho, exp(r1i), exp(r2i), dri, jr
    !write(*,*) temp, exp(t1i), exp(t2i), dti, jt
    !write(*,*) ye, ye1i, ye2i, dyei, jye

    !if(dri .gt. 1.d0 .OR. dri .lt. 0.d0) then
    !   write(*,*) "dri, dti, dyei: ", dri, dti, dyei
    !endif
    !if(dti .gt. 1.d0 .OR. dti .lt. 0.d0) then
    !   write(*,*) "dri, dti, dyei: ", dri, dti, dyei
    !endif
    !if(dyei .gt. 1.d0 .OR. dyei .lt. 0.d0) then
    !   write(*,*) "dri, dti, dyei: ", dri, dti, dyei
    !endif
    
    !
    ! absorption
    !
    coff(1) = (1.d0-dri)*(1.d0-dti)*(1.d0-dyei)
    coff(2) = dri*(1.d0-dti)*(1.d0-dyei)
    coff(3) = (1.d0-dri)*dti*(1.d0-dyei)
    coff(4) = dri*dti*(1.d0-dyei)
    coff(5) = (1.d0-dri)*(1.d0-dti)*dyei
    coff(6) = dri*(1.d0-dti)*dyei
    coff(7) = (1.d0-dri)*dti*dyei
    coff(8) = dri*dti*dyei


    if(get_k .eq. 1) then
       kappa(:) = coff(1)*absopac(:,jr,jt,jye)
       kappa(:) = kappa(:) + coff(2)*absopac(:,jr+1,jt,jye)
       kappa(:) = kappa(:) + coff(3)*absopac(:,jr,jt+1,jye)
       kappa(:) = kappa(:) + coff(4)*absopac(:,jr+1,jt+1,jye)
       kappa(:) = kappa(:) + coff(5)*absopac(:,jr,jt,jye+1)
       kappa(:) = kappa(:) + coff(6)*absopac(:,jr+1,jt,jye+1)
       kappa(:) = kappa(:) + coff(7)*absopac(:,jr,jt+1,jye+1)
       kappa(:) = kappa(:) + coff(8)*absopac(:,jr+1,jt+1,jye+1)
       kappa(:) = exp(kappa(:))*rho
    endif

    if(get_s .eq. 1) then
       sigma(:) = coff(1)*scaopac(:,jr,jt,jye)
       sigma(:) = sigma(:) + coff(2)*scaopac(:,jr+1,jt,jye)
       sigma(:) = sigma(:) + coff(3)*scaopac(:,jr,jt+1,jye)
       sigma(:) = sigma(:) + coff(4)*scaopac(:,jr+1,jt+1,jye)
       sigma(:) = sigma(:) + coff(5)*scaopac(:,jr,jt,jye+1)
       sigma(:) = sigma(:) + coff(6)*scaopac(:,jr+1,jt,jye+1)
       sigma(:) = sigma(:) + coff(7)*scaopac(:,jr,jt+1,jye+1)
       sigma(:) = sigma(:) + coff(8)*scaopac(:,jr+1,jt+1,jye+1)
       sigma(:) = exp(sigma(:))*rho
    endif

    if(get_d .eq. 1) then
       delta(:) = coff(1)*sdelta(:,jr,jt,jye)
       delta(:) = delta(:) + coff(2)*sdelta(:,jr+1,jt,jye)
       delta(:) = delta(:) + coff(3)*sdelta(:,jr,jt+1,jye)
       delta(:) = delta(:) + coff(4)*sdelta(:,jr+1,jt+1,jye)
       delta(:) = delta(:) + coff(5)*sdelta(:,jr,jt,jye+1)
       delta(:) = delta(:) + coff(6)*sdelta(:,jr+1,jt,jye+1)
       delta(:) = delta(:) + coff(7)*sdelta(:,jr,jt+1,jye+1)
       delta(:) = delta(:) + coff(8)*sdelta(:,jr+1,jt+1,jye+1)
    endif

    if(get_j .eq. 1) then
! .and. jye .gt. 1) then
       eta(:) = coff(1)*emis(:,jr,jt,jye)
       eta(:) = eta(:) + coff(2)*emis(:,jr+1,jt,jye)
       eta(:) = eta(:) + coff(3)*emis(:,jr,jt+1,jye)
       eta(:) = eta(:) + coff(4)*emis(:,jr+1,jt+1,jye)
       eta(:) = eta(:) + coff(5)*emis(:,jr,jt,jye+1)
       eta(:) = eta(:) + coff(6)*emis(:,jr+1,jt,jye+1)
       eta(:) = eta(:) + coff(7)*emis(:,jr,jt+1,jye+1)
       eta(:) = eta(:) + coff(8)*emis(:,jr+1,jt+1,jye+1)
       eta(:) = exp(eta(:))*rho
    endif

!    if(get_j .eq. 1 .and. jye .le. 1) then
!       sfunc(:,nu) = coff(1)*(emis(:,jr,jt,jye,nu)-absopac(:,jr,jt,jye,nu))
!       sfunc(:,nu) = sfunc(:,nu) + coff(2)*(emis(:,jr+1,jt,jye,nu)-absopac(:,jr+1,jt,jye,nu))
!       sfunc(:,nu) = sfunc(:,nu) + coff(3)*(emis(:,jr,jt+1,jye,nu)-absopac(:,jr,jt+1,jye,nu))
!       sfunc(:,nu) = sfunc(:,nu) + coff(4)*(emis(:,jr+1,jt+1,jye,nu)-absopac(:,jr+1,jt+1,jye,nu))
!       sfunc(:,nu) = sfunc(:,nu) + coff(5)*(emis(:,jr,jt,jye+1,nu)-absopac(:,jr,jt,jye+1,nu))
!       sfunc(:,nu) = sfunc(:,nu) + coff(6)*(emis(:,jr+1,jt,jye+1,nu)-absopac(:,jr+1,jt,jye+1,nu))
!       sfunc(:,nu) = sfunc(:,nu) + coff(7)*(emis(:,jr,jt+1,jye+1,nu)-absopac(:,jr,jt+1,jye+1,nu))
!       sfunc(:,nu) = sfunc(:,nu) + coff(8)*(emis(:,jr+1,jt+1,jye+1,nu)-absopac(:,jr+1,jt+1,jye+1,nu))
       !emi(:,nu) = sfunc(:,nu) + opac(:,nu)
!    endif

    !do ng=1,ngmax
    !   write(*,*) enuk(ng), emi(ng,1), opac(ng,1)
    !enddo
    !stop

  end subroutine get_opacity_emissivity


  subroutine get_inelastic_opacity_emissivity(je,he,t,r,y,etael, &
                                              xxn,xxp,source,sink,scatt)

    implicit none


    double precision, intent(in) :: t, r, y, etael, xxn, xxp
    double precision, intent(in) :: je(ngroups), he(ngroups*NDIM)
    double precision, intent(out) :: source(ngroups), sink(ngroups), scatt(ngroups)

    double precision :: xxj, xxe, xxeb, cut
    double precision :: sourcee, sinke, scatte
    double precision :: sourcen, sinkn, scattn
    integer :: nfreq,n0,n1,g,nh0,nh1
    integer :: inu 
    character*2 setscat
    character*4 nutype

    !setscat = 'bt'
    setscat = 'es'

    do inu=1,ngroups

        if(inu .le. ng0)then
          nutype = 'enu '
          n0 = 1
          n1 = ng0
          nh0 = 1
          nh1 = nh0 + ng0*NDIM - 1
          nfreq = ng0
        else if(inu .le. ng0 + ng1)then
          nutype = 'aenu'
          n0 = ng0+1
          n1 = n0+ng1-1
          nh0 = ng0*NDIM+1
          nh1 = nh0 + ng1*NDIM - 1
          nfreq = ng1
        else
          nutype = 'unu '
          n0 = ng0+ng1+1
          n1 = n0+ng2-1
          nh0 = (ng0+ng1)*NDIM+1
          nh1 = nh0 + ng2*NDIM - 1
          nfreq = ng2
        end if

        xxj = je(inu)
        xxe = egroupMeV(inu)

        xxeb = 50.d0

        sourcee = 0.d0
        sinke = 0.d0
        scatte = 0.d0
        sourcen = 0.d0
        sinkn = 0.d0
        scattn = 0.d0
        !source = 0.d0
        !sink = 0.d0
        !scatt = 0.d0

        if(r.gt.1.d8.and.setscat.ne.'no'.and.xxj.gt.0.d0)then
           if(setscat.eq.'es' .or. setscat.eq.'bt')then
              call resourcesink(nutype,nfreq,inu-n0+1,je(n0:n1),he(nh0:nh1), &
                              egroupMeV(n0:n1),t,r,y,etael,sourcee,sinke,scatte)
           endif
           if(setscat.eq.'ns' .or. setscat.eq.'bt')then
             if(xxe.lt.xxeb)then
               call xxnucresourcesink(nutype,nfreq,inu-n0+1,je(n0:n1), &
                              he(nh0:nh1),egroupMeV(n0:n1),t,r,y,xxn,xxp, &
                              sourcen,sinkn,scattn)
             else
               sourcen = 0.d0
               sinkn   = 0.d0
               scattn  = 0.d0
             endif
           endif
!           if(setscat.eq.'bt')then
!             call resourcesink(nutype,nfreq,inu-n0+1,je(n0:n1),he(nh0:nh1), &
!                              egroupMeV(n0:n1),t,r,y,etael,sourcee,sinke,scatte)
!             if(xxe.lt.xxeb)then
!               call xxnucresourcesink(nutype,nfreq,inu-n0+1,je(n0:n1), &
!                              he(nh0:nh1),egroupMeV(n0:n1),t,r,y,xxn,xxp, &
!                              sourcen,sinkn,scattn) 
!             else
!               sourcen = 0.d0
!               sinkn   = 0.d0
!               scattn  = 0.d0
!             endif
!          endif
        else
          source(inu) = 0.d0
          sink(inu)   = 0.d0
          scatt(inu)  = 0.d0
        endif

        source(inu) = sourcee + sourcen
        sink(inu)   = sinke   + sinkn
        scatt(inu)  = scatte  !  + scattn

        cut=1.d0
        !if(r.gt.1.d13)then
        !  cut=r/1.d13
        !endif   
    
        source(inu) = source(inu) / cut
        sink(inu)   = sink(inu)   / cut
        scatt(inu)  = scatt(inu)  / cut

        source(inu) = max(source(inu), 0.d0)
        sink(inu) = max(sink(inu), 0.d0)

        !source = 0.d0
        !sink = 0.d0
        !scatt = 0.d0

        !write(6,*) g, egroupMev(g), je(g), he(g), source, sink
    enddo
  
  end subroutine get_inelastic_opacity_emissivity

  subroutine resourcesink(nutype,nfreq,nf,je,he,freqe,t,r,y,etael,srce,sinke,scatt)

    implicit none

    character*4, intent(in) :: nutype
    integer, intent(in) :: nfreq, nf
    double precision, intent(in) :: t,r,y,etael
    double precision, intent(in) :: freqe(nfreq),je(nfreq),he(nfreq*NDIM)
    double precision, intent(out) :: srce, sinke, scatt

    double precision, parameter :: hbar=6.582122d-22,bigG=3.937d-17,clt=2.99792458d10
    double precision, parameter :: fourpi=12.56637061d0,pi=3.1415926535898d0
    double precision, parameter :: factf=(2.d0*pi*hbar*clt)**3/clt/1.60217733d-6
    double precision, parameter :: hc2pi=(2.d0*pi*hbar*clt)
    double precision, parameter :: constin=fourpi*bigG**2/hc2pi**6*1.60217733d-6
    double precision, parameter :: constout=fourpi*bigG**2/hc2pi**3/clt
    
    integer :: id,nfp,nfpp, nfpm
    double precision :: cut, czero, xje, xhe(NDIM), sumin, sumout, ssum
    double precision :: dnue, xjpe, xhpe(NDIM), omegae, expe, phi0, phi1
    double precision :: term, enu, fdotf, phi0ee, phi1ee
  
    cut   = 1.d0
    czero = 1.d0
    if(nutype.eq.'unu '.or.nutype.eq.'aunu') then
      cut = 4.d0
      czero = 0.d0
    endif
    if(nutype.ne.'unu '.and.nutype.ne.'aunu'.and.nutype.ne.'aenu'.and.nutype.ne.'enu ') then
       print*,'Wrong type in Resourcesink'
       stop
    endif

    xje    = max(0.d0,min(je(nf)*factf/freqe(nf)**3/cut,1.d0))
    !xhe    = max(0.d0,min(he(nf)*factf/freqe(nf)**3/cut,1.d0))
    do id=1,NDIM
      xhe(id) = he((nf-1)*NDIM + id)*factf/freqe(nf)**3/cut
      if (xhe(id) .ne. xhe(id)) then
         xhe(id) = 0.d0
      end if
    end do


    sumin  = 0.d0
    sumout = 0.d0
    ssum   = 0.d0

    do nfp = 1, nfreq
       nfpp   = min(nfp + 1, nfreq)
       nfpm   = max(nfp - 1, 1    )
       dnue   = (freqe(nfpp)-freqe(nfpm))/2.d0
       xjpe   = max(0.d0,min(je(nfp)*factf/freqe(nfp)**3/cut,1.d0))
       !xhpe   = max(0.d0,min(he(nfp)*factf/freqe(nfp)**3/cut,1.d0))
       fdotf  = 0.d0
       do id=1,NDIM
         xhpe(id) = he((nfp-1)*NDIM + id)*factf/freqe(nfp)**3/cut
         if (xhpe(id) .ne. xhpe(id)) then
            xhpe(id) = 0.d0
         end if
         fdotf = fdotf + xhe(id)*xhpe(id)
       end do
       omegae = freqe(nf)-freqe(nfp)
       expe   = max(min(exp(-omegae/t),1.d90),1.d-90)
       call phifind_interp(nutype,nf,nfp,t,etael,phi0ee,phi1ee)
       phi0   = phi0ee
       phi1   = phi1ee
       term   = freqe(nfp)**2*dnue
       sumin  = sumin  + term*expe*(0.5d0*phi0*xjpe*(1.d0-xje)- &
                         czero*3.d0*0.5d0*phi1*fdotf)
       sumout = sumout + term*(0.5d0*phi0*(1.d0-xjpe)- &
                         czero*3.d0*0.5d0*phi1*fdotf/xje)
       ssum   = ssum   + term*(0.5d0*phi0*(1.d0-xjpe+expe*xjpe))
    end do 
    enu=freqe(nf)
    sumin=sumin*cut

    srce  = constin  * enu**3 * sumin 
    sinke = constout * sumout  
    scatt = constout * ssum

    return
  end subroutine resourcesink

  subroutine phifind_interp(nutype,nf,nfp,xtemp,etalep,phi0e,phi1e)

    implicit none

    character*4, intent(in) :: nutype
    integer, intent(in) :: nf,nfp
    double precision, intent(in) :: xtemp, etalep
    double precision, intent(out) :: phi0e, phi1e

    integer :: jeta, jq
    double precision :: etal,tl,alpha,beta,ql,delta,sp,sq
    double precision :: t1,t2,t12,t22,eta1,eta2

    t1 = phie_t1
    t2 = phie_t2
    t12 = phie_t12
    t22 = phie_t22
    eta1 = phie_eta1
    eta2 = phie_eta2

    etal  = etalep
    tl    = dlog10(xtemp)
    alpha = t1+(etal-eta1)/(eta2-eta1)*(t12-t1)
    beta  = t2-t1+((t22-t12)-(t2-t1))*(etal-eta1)/(eta2-eta1)
    ql    = (tl - alpha)/beta
    delta = (etal-eta1)/(eta2-eta1)*dble(phi_neta)
    jeta  = 1 + idint(delta)
    jq    = 1 + idint(dble(phi_nt)*ql)
    if( jeta .lt. 2    ) jeta = 2
    if( jq   .lt. 2    ) jq   = 2
    if( jeta .ge. phi_neta ) jeta = phi_neta-1
    if( jq   .ge. phi_nt   ) jq   = phi_nt-1
    sp = delta - (jeta-1)
    sq = dble(phi_nt)*ql - (jq-1)
    if(nutype.eq.'enu ')then
       phi0e = 0.5D0*sQ*(sQ-1.D0)*phieetable(1,nf,nfp,jeta,jq-1) &
             + 0.5D0*sP*(sP-1.D0)*phieetable(1,nf,nfp,jeta-1,jq) &
             + (1.D0+sP*sQ-sP*sP-sQ*sQ)*phieetable(1,nf,nfp,jeta,jq) &
             + 0.5D0*sP*(sP-2.D0*sQ+1.D0)*phieetable(1,nf,nfp,jeta+1,jq) &
             + 0.5D0*sQ*(sQ-2.D0*sP+1.D0)*phieetable(1,nf,nfp,jeta,jq+1) &
             + sP*sQ*phieetable(1,nf,nfp,jeta+1,jq+1)
       phi1e = 0.5D0*sQ*(sQ-1.D0)*phieetable(2,nf,nfp,jeta,jq-1) &
             + 0.5D0*sP*(sP-1.D0)*phieetable(2,nf,nfp,jeta-1,jq) &
             + (1.D0+sP*sQ-sP*sP-sQ*sQ)*phieetable(2,nf,nfp,jeta,jq) &
             + 0.5D0*sP*(sP-2.D0*sQ+1.D0)*phieetable(2,nf,nfp,jeta+1,jq) &
             + 0.5D0*sQ*(sQ-2.D0*sP+1.D0)*phieetable(2,nf,nfp,jeta,jq+1) &
             + sP*sQ*phieetable(2,nf,nfp,jeta+1,jq+1)
    elseif(nutype.eq.'aenu')then
       phi0e = 0.5D0*sQ*(sQ-1.D0)*phiaetable(1,nf,nfp,jeta,jq-1) &
             + 0.5D0*sP*(sP-1.D0)*phiaetable(1,nf,nfp,jeta-1,jq) &
             + (1.D0+sP*sQ-sP*sP-sQ*sQ)*phiaetable(1,nf,nfp,jeta,jq) &
             + 0.5D0*sP*(sP-2.D0*sQ+1.D0)*phiaetable(1,nf,nfp,jeta+1,jq) &
             + 0.5D0*sQ*(sQ-2.D0*sP+1.D0)*phiaetable(1,nf,nfp,jeta,jq+1) &
             + sP*sQ*phiaetable(1,nf,nfp,jeta+1,jq+1)
       phi1e = 0.5D0*sQ*(sQ-1.D0)*phiaetable(2,nf,nfp,jeta,jq-1) &
             + 0.5D0*sP*(sP-1.D0)*phiaetable(2,nf,nfp,jeta-1,jq) &
             + (1.D0+sP*sQ-sP*sP-sQ*sQ)*phiaetable(2,nf,nfp,jeta,jq) &
             + 0.5D0*sP*(sP-2.D0*sQ+1.D0)*phiaetable(2,nf,nfp,jeta+1,jq) &
             + 0.5D0*sQ*(sQ-2.D0*sP+1.D0)*phiaetable(2,nf,nfp,jeta,jq+1) &
             + sP*sQ*phiaetable(2,nf,nfp,jeta+1,jq+1)
    elseif(nutype.eq.'unu '.or.nutype.eq.'aunu')then
       phi0e = 0.5D0*sQ*(sQ-1.D0)*phimetable(1,nf,nfp,jeta,jq-1) &
             + 0.5D0*sP*(sP-1.D0)*phimetable(1,nf,nfp,jeta-1,jq) &
             + (1.D0+sP*sQ-sP*sP-sQ*sQ)*phimetable(1,nf,nfp,jeta,jq) &
             + 0.5D0*sP*(sP-2.D0*sQ+1.D0)*phimetable(1,nf,nfp,jeta+1,jq) &
             + 0.5D0*sQ*(sQ-2.D0*sP+1.D0)*phimetable(1,nf,nfp,jeta,jq+1) &
             + sP*sQ*phimetable(1,nf,nfp,jeta+1,jq+1)
       phi1e = 0.5D0*sQ*(sQ-1.D0)*phimetable(2,nf,nfp,jeta,jq-1) &
             + 0.5D0*sP*(sP-1.D0)*phimetable(2,nf,nfp,jeta-1,jq) &
             + (1.D0+sP*sQ-sP*sP-sQ*sQ)*phimetable(2,nf,nfp,jeta,jq) & 
             + 0.5D0*sP*(sP-2.D0*sQ+1.D0)*phimetable(2,nf,nfp,jeta+1,jq) &
             + 0.5D0*sQ*(sQ-2.D0*sP+1.D0)*phimetable(2,nf,nfp,jeta,jq+1) &
             + sP*sQ*phimetable(2,nf,nfp,jeta+1,jq+1)
    endif 
    return
  end subroutine phifind_interp

  subroutine xxnucresourcesink(nutype,nfreq,nf,je,he,freqe,t,r,y,xxn,xxp,sourcen,sinkn,scattn) 

    implicit none

    character*4, intent(in) :: nutype
    integer, intent(in) :: nfreq, nf
    double precision, intent(in) :: freqe(nfreq),je(nfreq),he(nfreq*NDIM)
    double precision, intent(in) :: t, r, y, xxn, xxp
    double precision, intent(out) :: sourcen, sinkn, scattn

    double precision, parameter :: hbar=6.582122d-22,bigG=3.937d-17,clt=2.99792458d10
    double precision, parameter :: fourpi=12.56637061d0,pi=3.1415926535898d0
    double precision, parameter :: factf=(2.d0*pi*hbar*clt)**3/clt/1.60217733d-6
    double precision, parameter :: hc2pi=(2.d0*pi*hbar*clt)
    double precision, parameter :: constin=fourpi*bigG**2/hc2pi**6*1.60217733d-6
    double precision, parameter :: constout=fourpi*bigG**2/hc2pi**3/clt

    character*1 :: tgt
    integer :: iln, ihn, ilp, ihp, jetan, jqn, jetap, jqp
    double precision :: etaxn, etaxp, sumin, sumout, ssum
    double precision :: suminn, suminp, sumoutn, sumoutp
    double precision :: ssumn, ssump, enu
    double precision :: sinn1, soutn1, ssn1
    double precision :: sinn2, soutn2, ssn2
    double precision :: sinn3, soutn3, ssn3
    double precision :: sinn4, soutn4, ssn4
    double precision :: sinn5, soutn5, ssn5
    double precision :: sinn6, soutn6, ssn6
    double precision :: sinp1, soutp1, ssp1
    double precision :: sinp2, soutp2, ssp2
    double precision :: sinp3, soutp3, ssp3
    double precision :: sinp4, soutp4, ssp4
    double precision :: sinp5, soutp5, ssp5
    double precision :: sinp6, soutp6, ssp6


! .. calculate chemical potentials

    etaxn  = max(dlog((2.d0*pi/939.57d0/t)**1.5d0/2.d0*xxn* &
                r/1.675d-24*(hbar*clt)**3.),-20.d0)
    etaxp  = max(dlog((2.d0*pi/938.27d0/t)**1.5d0/2.d0*xxp* &
                r/1.675d-24*(hbar*clt)**3.),-20.d0)

! .. find ilows and ihighs for given etas and temps

    call xifind(t,etaxn,etaxp,jetan,jqn,jetap,jqp)

    iln=max(nf-5,1)
    ihn=min(nf+5,nfreq)

    ilp=max(nf-5,1)
    ihp=min(nf+5,nfreq)

! .. integrate source, sink, and scattering opacity

    tgt='n'
    call xiint(nutype,tgt,nf,nfreq,je,he,freqe,iln,ihn,jetan  ,jqn  ,sinn1,soutn1,ssn1)
    call xiint(nutype,tgt,nf,nfreq,je,he,freqe,iln,ihn,jetan+1,jqn  ,sinn2,soutn2,ssn2)
    call xiint(nutype,tgt,nf,nfreq,je,he,freqe,iln,ihn,jetan  ,jqn+1,sinn3,soutn3,ssn3)
    call xiint(nutype,tgt,nf,nfreq,je,he,freqe,iln,ihn,jetan+1,jqn+1,sinn6,soutn6,ssn6)

    tgt='p'
    call xiint(nutype,tgt,nf,nfreq,je,he,freqe,ilp,ihp,jetap  ,jqp  ,sinp1,soutp1,ssp1)
    call xiint(nutype,tgt,nf,nfreq,je,he,freqe,ilp,ihp,jetap+1,jqp  ,sinp2,soutp2,ssp2)
    call xiint(nutype,tgt,nf,nfreq,je,he,freqe,ilp,ihp,jetap  ,jqp+1,sinp3,soutp3,ssp3)
    call xiint(nutype,tgt,nf,nfreq,je,he,freqe,ilp,ihp,jetap+1,jqp+1,sinp6,soutp6,ssp6)

!    if (je(nf) .gt. 0.d0) then
!      write(6,931) "n ", nutype, nf, freqe(nf), je(nf), he(nf), sinn1, sinn2, sinn3, sinn6, t, etaxn, jetan, jqn, iln, ihn
!      write(6,931) "p ", nutype, nf, freqe(nf), je(nf), he(nf), sinp1, sinp2, sinp3, sinp6, t, etaxp, jetap, jqp, ilp, ihp
!    end if
!931 format (A2, A4, 2X, I3, 9(2X, E12.5), 4(2X, I3))

! .. interpolate each source, sink, and scatt integral

    call xiinterp(t,etaxn,sinn1,sinn2,sinn3, &
          sinn4,sinn5,sinn6,suminn)
    call xiinterp(t,etaxp,sinp1,sinp2,sinp3, &
          sinp4,sinp5,sinp6,suminp)

    call xiinterp(t,etaxn,soutn1,soutn2,soutn3, &
          soutn4,soutn5,soutn6,sumoutn)
    call xiinterp(t,etaxp,soutp1,soutp2,soutp3, &
          soutp4,soutp5,soutp6,sumoutp)

    call xiinterp(t,etaxn,ssn1,ssn2,ssn3, &
          ssn4,ssn5,ssn6,ssumn)
    call xiinterp(t,etaxp,ssp1,ssp2,ssp3, &
          ssp4,ssp5,ssp6,ssump)

! .. sum neutron and proton contributions

    sumin  = suminn  + suminp
    sumout = sumoutn + sumoutp
    ssum   = ssumn   + ssump

! .. construct total source, sink, and scattering opacity

!      if(type.eq.'enu '                  ) enu = freqe(nf)
!      if(type.eq.'aenu'                  ) enu = freqe(nf)
!      if(type.eq.'unu '.or.type.eq.'aunu') 
    enu = freqe(nf)

    sumin  = constin  * enu**3 * sumin     ! erg/MeV/s/cm^3/Str
    sumout = constout * sumout             ! cm^-1
    ssum   = constout * ssum               ! cm^-1

    sourcen  = sumin
    sinkn   = sumout
    scattn   = ssum

    return
  end subroutine xxnucresourcesink

  subroutine xiint(nutype,tgt,nf,nfreq,je,he,freqe,ilow,ihigh,jeta,jq, &
                   sumin,sumout,ssum)

    implicit none

    character*4, intent(in) :: nutype
    character*1, intent(in) :: tgt
    integer, intent(in) :: nf, ilow, ihigh, jeta, jq,nfreq
    double precision, intent(in) :: je(nfreq), he(nfreq), freqe(nfreq*NDIM)
    double precision, intent(out) :: sumin, sumout, ssum
    double precision, parameter :: hbar=6.582122d-22,bigG=3.937d-17,clt=2.99792458d10
    double precision, parameter :: fourpi=12.56637061d0,pi=3.1415926535898d0
    double precision, parameter :: factf=(2.d0*pi*hbar*clt)**3/clt/1.60217733d-6
    double precision, parameter :: hc2pi=(2.d0*pi*hbar*clt)
    double precision, parameter :: constin=fourpi*bigG**2/hc2pi**6*1.60217733d-6
    double precision, parameter :: constout=fourpi*bigG**2/hc2pi**3/clt

    integer :: nn,id,nfp, nfpp
    double precision :: xjp1,xjp2,xhp1(NDIM),xhp2(NDIM),aj,bj,ah,bh
    double precision :: xje,xhe(NDIM),sum20,sum30,sum2e0,sum3e0
    double precision :: aaa,bbb,ccc

    sumin  = 0.d0
    sumout = 0.d0
    ssum   = 0.d0

    if(nutype.eq.'enu ')then
       !xje    = max(0.d0,min(je(nf)*factf/freqe(nf)**3,1.d0))
       !xhe    = max(0.d0,min(he(nf)*factf/freqe(nf)**3,1.d0))
       !do id=1,flux_dim
       !   xhe(id) = he((nf-1)*flux_dim + id)
       !end do
       do nn=1,(ihigh-ilow)
          nfp    = ilow+nn-1
          nfpp   = nfp+1
          xjp1   = min(je(nfp)*factf/freqe(nfp)**3,1.d0)
          xjp2   = min(je(nfpp)*factf/freqe(nfpp)**3,1.d0)
          !do id=1, flux_dim
          !   xhp1(id) = he((nfp-1)*flux_dim + id)*factf/freqe(nfp)**3
          !   xhp2(id) = he((nfpp-1)*flux_dim + id)*factf/freqe(nfpp)**3
          !end do
          !xhp1   = min(he(nfp)*factf/freqe(nfp)**3,1.d0)
          !xhp2   = min(he(nfpp)*factf/freqe(nfpp)**3,1.d0)
          aj     = (xjp1-xjp2)/(freqe(nfp)-freqe(nfpp))
          bj     = xjp1-aj*freqe(nfp)
          !ah     = (xhp1-xhp2)/(freqe(nfp)-freqe(nfpp))
          !bh     = xhp1-ah*freqe(nfp)
          if(tgt.eq.'n')then
             sum20  = phientable(3 ,nf,nn,jeta,jq)
             sum30  = phientable(4 ,nf,nn,jeta,jq)
             sum2e0 = phientable(5 ,nf,nn,jeta,jq)
             sum3e0 = phientable(6 ,nf,nn,jeta,jq)
          elseif(tgt.eq.'p')then
             sum20  = phieptable(3 ,nf,nn,jeta,jq)
             sum30  = phieptable(4 ,nf,nn,jeta,jq)
             sum2e0 = phieptable(5 ,nf,nn,jeta,jq)
             sum3e0 = phieptable(6 ,nf,nn,jeta,jq)
          else
             print*,'wrong tgt'
             stop
          endif
          sumin  = sumin  + (aj*sum3e0+bj*sum2e0)!(1.d0-xje)*(aj*sum3e0+bj*sum2e0)!-
!     .                         xhe*(ah*sum3e1+bh*sum2e1)
          sumout = sumout + sum20-(aj*sum30+bj*sum20)+(aj*sum3e0+bj*sum2e0)!-
!     .                         (xhe/xje)*(ah*sum31+bh*sum21)
          ssum   = ssum   + sum20!-(aj*sum30+bj*sum20)!+
!     .                         (aj*sum3e0+bj*sum2e0)
       enddo

    elseif(nutype.eq.'aenu')then
       !xje    = max(0.d0,min(je(nf)*factf/freqe(nf)**3,1.d0))
       !xhe    = max(0.d0,min(he(nf)*factf/freqe(nf)**3,1.d0))
       do nn=1,(ihigh-ilow)
          nfp    = ilow+nn-1
          nfpp   = nfp+1
          xjp1   = min(je(nfp)*factf/freqe(nfp)**3,1.d0)
          xjp2   = min(je(nfpp)*factf/freqe(nfpp)**3,1.d0)
          !xhp1   = min(he(nfp)*factf/freqe(nfp)**3,1.d0)
          !xhp2   = min(he(nfpp)*factf/freqe(nfpp)**3,1.d0)
          aj     = (xjp1-xjp2)/(freqe(nfp)-freqe(nfpp))
          bj     = xjp1-aj*freqe(nfp)
          !ah     = (xhp1-xhp2)/(freqe(nfp)-freqe(nfpp))
          !bh     = xhp1-ah*freqe(nfp)
          if(tgt.eq.'n')then
             sum20  = phiantable(3 ,nf,nn,jeta,jq)
             sum30  = phiantable(4 ,nf,nn,jeta,jq)
             sum2e0 = phiantable(5 ,nf,nn,jeta,jq)
             sum3e0 = phiantable(6 ,nf,nn,jeta,jq)
          elseif(tgt.eq.'p')then
             sum20  = phiaptable(3 ,nf,nn,jeta,jq)
             sum30  = phiaptable(4 ,nf,nn,jeta,jq)
             sum2e0 = phiaptable(5 ,nf,nn,jeta,jq)
             sum3e0 = phiaptable(6 ,nf,nn,jeta,jq)
          else
             print*,'wrong tgt'
             stop
          endif
          sumin  = sumin  + (aj*sum3e0+bj*sum2e0)!(1.d0-xje)*(aj*sum3e0+bj*sum2e0)!-
!     .                         xhe*(ah*sum3e1+bh*sum2e1)
          sumout = sumout + sum20-(aj*sum30+bj*sum20)+(aj*sum3e0+bj*sum2e0)!!-
!     .                         (xhe/xje)*(ah*sum31+bh*sum21)
          ssum   = ssum   + sum20!-(aj*sum30+bj*sum20)!+
!     .                         (aj*sum3e0+bj*sum2e0)
       enddo

    elseif(nutype.eq.'unu '.or.nutype.eq.'aunu')then
       xje    = min(max(0.d0,je(nf)*factf/freqe(nf)**3/4.d0),1.d0)!/4.d0)
       !xhe    = min(max(0.d0,he(nf)*factf/freqe(nf)**3/4.d0),1.d0)!/4.d0)
       do nn=1,(ihigh-ilow)
          nfp    = ilow+nn-1
          nfpp   = nfp+1
          xjp1   = min(max(0.d0,je(nfp)*factf/freqe(nfp)**3/4.d0),1.d0)
          xjp2   = min(max(0.d0,je(nfpp)*factf/freqe(nfpp)**3/4.d0),1.d0)
          !xhp1   = min(max(0.d0,he(nfp)*factf/freqe(nfp)**3/4.d0),1.d0)
          !xhp2   = min(max(0.d0,he(nfpp)*factf/freqe(nfpp)**3/4.d0),1.d0)
          aj     = (xjp1-xjp2)/(freqe(nfp)-freqe(nfpp))
          bj     = xjp1-aj*freqe(nfp)
          !ah     = (xhp1-xhp2)/(freqe(nfp)-freqe(nfpp))
          !bh     = xhp1-ah*freqe(nfp)

          if(tgt.eq.'n')then
             sum20  = phimntable(3 ,nf,nn,jeta,jq)
             sum30  = phimntable(4 ,nf,nn,jeta,jq)
             sum2e0 = phimntable(5 ,nf,nn,jeta,jq)
             sum3e0 = phimntable(6 ,nf,nn,jeta,jq)
          elseif(tgt.eq.'p')then
             sum20  = phimptable(3 ,nf,nn,jeta,jq)
             sum30  = phimptable(4 ,nf,nn,jeta,jq)
             sum2e0 = phimptable(5 ,nf,nn,jeta,jq)
             sum3e0 = phimptable(6 ,nf,nn,jeta,jq)
          else
             print*,'wrong tgt'
             stop
          endif

          if(nfp.le.25)then
             aaa=(1.d0-xje)*sum2e0*xjp2
             bbb=sum20-sum20*xjp2
             ccc=sum20
          else
             aaa=0.d0
             bbb=0.d0
             ccc=0.d0
          endif

          sumin  = sumin  + aaa!(aj*sum3e0+bj*sum2e0)!-
!     .                         xhe*(ah*sum3e1+bh*sum2e1)
          sumout = sumout + bbb!(aj*sum30+bj*sum20)!+4.d0*(aj*sum3e0+bj*sum2e0)!!-
!     .                         (xhe/xje)*(ah*sum31+bh*sum21)
          ssum   = ssum   + ccc!-(aj*sum30+bj*sum20)!+
!     .                         (aj*sum3e0+bj*sum2e0)

       enddo
       sumin  = 4.d0*sumin
       sumout = sumout
    else
       print*,'wrong type in xiint'
       stop
    endif

    return
  end subroutine xiint


  subroutine xifind(t,etaxn,etaxp,jetan,jqn,jetap,jqp)

    implicit none

    integer, intent(out) :: jetan, jqn, jetap, jqp
    double precision, intent(in) :: t, etaxn, etaxp

    integer :: nt, neta, jeta, jq
    double precision :: tl, eta1, eta2, t1, t2, t12, t22
    double precision :: etal, alpha, beta, ql, delta, sp, sq

    nt=nnt
    neta=nneta

    tl     =  dlog10(t)
! for the nucleons, not the electrons
    eta1   = -15.d0  
    eta2   =  5.d0   
    t1     =  -1.d0
    t2     =  1.4d0
    t12    =  -1.d0
    t22    =  1.4d0

! .. neutron scattering
    etal   =  etaxn
    alpha  = t1+(etal-eta1)/(eta2-eta1)*(t12-t1)
    beta   = t2-t1+((t22-t12)-(t2-t1))*(etal-eta1)/(eta2-eta1)
    ql     = (tl - alpha)/beta
    delta  = (etal-eta1)/(eta2-eta1)*dble(neta)
    jeta   = 1 + idint(delta)
    jq     = 1 + idint(dble(nt)*ql)
    sp     = delta - (jeta-1)
    sq     = dble(nt)*ql - (jq-1)

    jetan = jeta
    jqn   = jq

! .. proton scattering

    etal   =  etaxp
    alpha  = t1+(etal-eta1)/(eta2-eta1)*(t12-t1)
    beta   = t2-t1+((t22-t12)-(t2-t1))*(etal-eta1)/(eta2-eta1)
    ql     = (tl - alpha)/beta
    delta  = (etal-eta1)/(eta2-eta1)*dble(neta)
    jeta   = 1 + idint(delta)
    jq     = 1 + idint(dble(nt)*ql)
    sp     = delta - (jeta-1)
    sq     = dble(nt)*ql - (jq-1)

    jetap = jeta
    jqp   = jq

    if( jetan .le. 1     ) jetan = 1+1
    if( jetan .ge. nneta ) jetan = nneta-1
    if( jetap .le. 1     ) jetap = 1+1
    if( jetap .ge. nneta ) jetap = nneta-1
    if( jqn   .le. 1     ) jqn   = 1+1
    if( jqn   .ge. nnt   ) jqn   = nnt-1
    if( jqp   .le. 1     ) jqp   = 1+1
    if( jqp   .ge. nnt   ) jqp   = nnt-1

    return
  end subroutine xifind
 
  subroutine xiinterp(t,etan,x1,x2,x3,x4,x5,x6,xpt)

    implicit none

    double precision, intent(in) :: t, etan, x1, x2, x3, x4, x5, x6
    double precision, intent(out) :: xpt

    integer :: nt, neta, jeta, jq
    double precision :: tl, etal, eta1, eta2, t1, t2, t12, t22
    double precision :: alpha, beta, ql, delta, sp, sq

    nt=nnt
    neta=nneta
    tl    =  dlog10(t)
    etal  =  etan
    eta1  = -15.d0
    eta2  =  5.d0
    t1    =  -1.d0
    t2    =  1.4d0
    t12   =  -1.d0
    t22   =  1.4d0
    alpha = t1+(etal-eta1)/(eta2-eta1)*(t12-t1)
    beta  = t2-t1+((t22-t12)-(t2-t1))*(etal-eta1)/(eta2-eta1)
    ql    = (tl - alpha)/beta
    delta = (etal-eta1)/(eta2-eta1)*dble(neta)
    jeta  = 1 + idint(delta)
    jq    = 1 + idint(dble(nt)*ql)
    sp = delta - (jeta-1)
    sq = dble(nt)*ql - (jq-1)
    xpt = (1.D0-sP)*(1.d0-sQ)*x1+sP*(1.D0-sQ)*x2+ &
          sQ*(1.D0-sP)*x3+sP*sQ*x6
    return
  end subroutine xiinterp

end module opacity_table_module

  subroutine init_opacity_table(egrp,ngrp,ntot,opac_param_name,opac_name)

    use opacity_table_module

    implicit none
    
    integer, intent(in) :: ntot, ngrp(3)
    double precision, intent(in) :: egrp(ntot)
    integer i, j, k, ng, jf
    integer irecl
    integer nu, nu_max, igrp
    real*8 f1i, f2i, dfi, el
    real*8 dfr, frmin, frmax

    real, allocatable :: temp(:,:,:,:,:)
    
    !character*80 rstfil
    character(*), intent(in) :: opac_param_name, opac_name
    integer :: ncpu,mype,info,ierr

#if(USE_MPI==TRUE)
    include 'mpif.h'
    ! Init MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,mype,info)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,info)
#else
    ncpu = 1
    mype = 0
#endif

    if (mype == 0) then
        open(3, file=opac_param_name, &
             status='old')
        read(3,*) ngmax
    endif
#if(USE_MPI==TRUE)
    call MPI_BCAST(ngmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#endif
    allocate(enuk(ngmax), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_opacity_table]:  proc',mype,'failed to allocate enuk!'
    if(mype == 0) then
        read(3,*) romin, romax, nromax, tmin, tmaxi, ntmax, &
             yemin, yemax, nyemax, ngmax, (enuk(i), i=1,ngmax)
        close(3)
    endif
#if(USE_MPI==TRUE)
    call MPI_BCAST(romin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(romax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(nromax,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(tmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(tmaxi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(ntmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(yemin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(yemax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(nyemax,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(ngmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
    call MPI_BCAST(enuk,ngmax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
#endif

    frmin = log(enuk(1))
    frmax = log(enuk(ngmax))
    dfr = log(enuk(ngmax)/enuk(1))/(ngmax-1.d0)

    allocate(absopac(ntot,nromax,ntmax,nyemax), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_opacity_table]:  proc',mype,'failed to allocate absopac!'
    allocate(scaopac(ntot,nromax,ntmax,nyemax), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_opacity_table]:  proc',mype,'failed to allocate scaopac!'
    allocate(emis(ntot,nromax,ntmax,nyemax), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_opacity_table]:  proc',mype,'failed to allocate emis!'
    allocate(sdelta(ntot,nromax,ntmax,nyemax), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_opacity_table]:  proc',mype,'failed to allocate sdelta!'


    if(mype == 0) then
        allocate(temp(nromax, ntmax, nyemax,ngmax,3), stat=ierr)
        if (ierr /= 0) write (*,*) '[init_opacity_table]:  proc',mype,'failed to allocate temp!'
        irecl = 4*nromax*ntmax*nyemax*ngmax + 10
        open(2,file=opac_name,form='unformatted',access='direct', &
             recl=irecl,status='old')
        nu_max = 3

        ! read in the absorption, interp onto runs energy grid
        do nu = 1, nu_max
           read(2, rec = nu) ((((temp(i,j,k,ng,nu), i = 1,nromax), &
                j=1,ntmax),k=1,nyemax), ng=1,ngmax)
        enddo
        do k=1,nyemax
            do j=1,ntmax
                do i=1,nromax
                    do ng=1,ntot
                        igrp = 1
                        if(ng > ngrp(1)) then
                            if(ng > ngrp(1)+ngrp(2)) then
                                igrp = 3
                            else
                                igrp = 2
                            end if
                        end if
                        el = log(egrp(ng))
                        jf = idint((el-frmin)/dfr) + 1
                        if(jf .lt. 1) then
                            jf = 1
                            dfi = 0.d0
                        else if(jf .gt. ngmax-1) then
                            jf = ngmax-1
                            dfi = 1.d0
                        else
                            f1i = frmin + (jf-1)*dfr
                            f2i = frmin + jf*dfr
                            dfi = (el - f1i)/dfr
                        end if
                        absopac(ng,i,j,k) = (1.d0 - dfi)*temp(i,j,k,jf,igrp) + dfi*temp(i,j,k,jf+1,igrp)
                    enddo
                enddo
            enddo
        enddo


        do nu = 1, nu_max
            read(2, rec =  3 + nu) ((((temp(i,j,k,ng,nu), &
                i = 1,nromax), j=1,ntmax), k=1,nyemax), ng=1,ngmax)
        enddo
        do k=1,nyemax
            do j=1,ntmax
                do i=1,nromax
                    do ng=1,ntot
                        igrp = 1
                        if(ng > ngrp(1)) then
                            if(ng > ngrp(1)+ngrp(2)) then
                                igrp = 3
                            else
                                igrp = 2
                            end if
                        end if
                        el = log(egrp(ng))
                        jf = idint((el-frmin)/dfr) + 1
                        if(jf .lt. 1) then
                            jf = 1
                            dfi = 0.d0
                        else if(jf .gt. ngmax-1) then
                            jf = ngmax-1
                            dfi = 1.d0
                        else
                            f1i = frmin + (jf-1)*dfr
                            f2i = frmin + jf*dfr
                            dfi = (el - f1i)/dfr
                        end if
                        scaopac(ng,i,j,k) = (1.d0 - dfi)*temp(i,j,k,jf,igrp) + dfi*temp(i,j,k,jf+1,igrp)
                    enddo
                enddo
            enddo
        enddo

        do nu = 1, nu_max
            read(2, rec = 6 + nu) ((((temp(i,j,k,ng,nu), i = 1,nromax), &
                j=1,ntmax), k=1,nyemax), ng=1,ngmax)
        enddo
        do k=1,nyemax
            do j=1,ntmax
                do i=1,nromax
                    do ng=1,ntot
                        igrp = 1
                        if(ng > ngrp(1)) then
                            if(ng > ngrp(1)+ngrp(2)) then
                                igrp = 3
                            else
                                igrp = 2
                            end if
                        end if
                        el = log(egrp(ng))
                        jf = idint((el-frmin)/dfr) + 1
                        if(jf .lt. 1) then
                            jf = 1
                            dfi = 0.d0
                        else if(jf .gt. ngmax-1) then
                            jf = ngmax-1
                            dfi = 1.d0
                        else
                            f1i = frmin + (jf-1)*dfr
                            f2i = frmin + jf*dfr
                            dfi = (el - f1i)/dfr
                        end if
                        emis(ng,i,j,k) = (1.d0 - dfi)*temp(i,j,k,jf,igrp) + dfi*temp(i,j,k,jf+1,igrp)
                    enddo
                enddo
            enddo
        enddo


        do nu = 1, nu_max
            read(2, rec = 9 + nu) ((((temp(i,j,k,ng,nu), &
                i = 1,nromax), j=1,ntmax), k=1,nyemax), ng=1,ngmax)
        enddo
        do k=1,nyemax
            do j=1,ntmax
                do i=1,nromax
                    do ng=1,ntot
                        igrp = 1
                        if(ng > ngrp(1)) then
                            if(ng > ngrp(1)+ngrp(2)) then
                                igrp = 3
                            else
                                igrp = 2
                            end if
                        end if
                        el = log(egrp(ng))
                        jf = idint((el-frmin)/dfr) + 1
                        if(jf .lt. 1) then
                            jf = 1
                            dfi = 0.d0
                        else if(jf .gt. ngmax-1) then
                            jf = ngmax-1
                            dfi = 1.d0
                        else
                            f1i = frmin + (jf-1)*dfr
                            f2i = frmin + jf*dfr
                            dfi = (el - f1i)/dfr
                        end if
                        sdelta(ng,i,j,k) = (1.d0 - dfi)*temp(i,j,k,jf,igrp) + dfi*temp(i,j,k,jf+1,igrp)
                    enddo
                enddo
            enddo
        enddo

        close(2)
        deallocate(temp)
    endif

#if(USE_MPI==TRUE)  
    call MPI_BCAST(absopac, size(absopac), MPI_REAL, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(scaopac, size(scaopac), MPI_REAL, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(emis   , size(emis   ), MPI_REAL, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(sdelta , size(sdelta ), MPI_REAL, 0, MPI_COMM_WORLD, info)
#endif
    !do ng=1,ngmax
    !  write(*,*) enuk(ng), dexp(emis(ng,36,20,15,1))*1.389495494373d12, &
    !                dexp(absopac(ng,36,20,15,1))*1.389495494373d12
    !enddo

    !do k=1,nyemax
    !   write(*,"(f6.4)",advance='no') yemin + (k-1)*(yemax-yemin)/(nyemax-1)
    !   do ng=1,ngmax
    !      write(*,"(E14.7)",advance='no') exp(emis(ng,44,37,k,1) - absopac(ng,44,37,k,1))
    !   enddo
    !   write(*,*)
    !enddo

    !stop
    
  end subroutine init_opacity_table

  subroutine init_inelastic_opacity_table(egrp, ngrp, ntot, inelastic_root)

    use opacity_table_module

    implicit none

    integer, intent(in) :: ngrp(3),ntot
    double precision, intent(in) :: egrp(ntot)
    character(*), intent(in) :: inelastic_root
    integer :: irecl,irecln,irec,m,ke,le,jeta,jt
    integer :: ncpu,mype,info,ierr
    character*3 :: nfe, nfa, nfm
    double precision, allocatable :: fe(:), fa(:), fm(:)

#if(USE_MPI==TRUE)
    include 'mpif.h'
    ! Init MPI
    call MPI_COMM_RANK(MPI_COMM_WORLD,mype,info)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,info)
#else
    ncpu = 1
    mype = 0
#endif

    ng0 = ngrp(1)
    ng1 = ngrp(2)
    ng2 = ngrp(3)
    ngroups = ntot

    allocate(egroupMeV(ngroups), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate egroupMeV!'

    do m=1, ngroups
        egroupMeV(m) = egrp(m)
    enddo


    write(nfe,'(I3.3)') ng0
    write(nfa,'(I3.3)') ng1
    write(nfm,'(I3.3)') ng2

    allocate(phieetable(2,ng0,ng0,phi_neta,phi_nt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phieetable!'
    allocate(phiaetable(2,ng1,ng1,phi_neta,phi_nt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phiaetable!'
    allocate(phimetable(2,ng2,ng2,phi_neta,phi_nt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phimetable!'
    allocate(phientable(6,ng0,nio,nneta,nnt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phientable!'    
    allocate(phiantable(6,ng1,nio,nneta,nnt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phiantable!'
    allocate(phimntable(6,ng2,nio,nneta,nnt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phimntable!'
    allocate(phieptable(6,ng0,nie,nneta,nnt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phieptable!'    
    allocate(phiaptable(6,ng1,nio,nneta,nnt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phiaptable!'
    allocate(phimptable(6,ng2,nio,nneta,nnt), stat=ierr)
    if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate phimptable!'

    if (mype.eq.0) then

        allocate(fe(ng0), stat=ierr)
        if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate fe!'
        allocate(fa(ng1), stat=ierr)
        if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate fa!'
        allocate(fm(ng2), stat=ierr)
        if (ierr /= 0) write (*,*) '[init_inelastic_opacity_table]:  proc',mype,'failed to allocate fm!'

        irecl = 2*8*ng0*ng0*phi_neta*phi_nt
        irecln = 6*8*ng0*nio*nneta*nnt

        open(20,file=inelastic_root//'phi_ee'//nfe//'_030x030_a.dat', &
             form='unformatted', access='direct',recl=irecl)
        irec = 1
        read(20,rec=irec)(fe(ke),ke = 1,ng0)
        irec = 2
        read(20,rec=irec)(((((phieetable(m,ke,le,jeta,jt),jt = 1,phi_nt), &
             jeta = 1,phi_neta),le = 1,ng0),ke = 1,ng0),m = 1,2)
        close(20)

        open(20,file=inelastic_root//'phi_en'//nfe//'_030x030_a.dat', &
             status='unknown')
        irec = 1
        read(20,*)(fe(ke),ke = 1,ng0)
        irec = 2
        read(20,*)(((((phientable(m,ke,le,jeta,jt),jt = 1,nnt), &
             jeta = 1,nneta),le = 1,nio),ke = 1,ng0),m = 1,6)
        close(20)

        open(20,file=inelastic_root//'phi_ep'//nfe//'_030x030_a.dat', &
             status='unknown')
        irec = 1
        read(20,*)(fe(ke),ke = 1,ng0)
        irec = 2
        read(20,*)(((((phieptable(m,ke,le,jeta,jt),jt = 1,nnt), &
             jeta = 1,nneta),le = 1,nie),ke = 1,ng0),m = 1,6)
        close(20)

        irecl = 2*8*ng1*ng1*phi_neta*phi_nt
        irecln = 6*8*ng1*nio*nneta*nnt

        open(20,file=inelastic_root//'phi_ae'//nfa//'_030x030_a.dat', &
             form='unformatted', access='direct',recl=irecl)
        irec = 1
        read(20,rec=irec)(fa(ke),ke = 1,ng1)
        irec = 2
        read(20,rec=irec)(((((phiaetable(m,ke,le,jeta,jt),jt = 1,phi_nt), &
             jeta = 1,phi_neta),le = 1,ng1),ke = 1,ng1),m = 1,2)
        close(20)

        open(20,file=inelastic_root//'phi_an'//nfa//'_030x030_a.dat', &
             status='unknown')
        irec = 1
        read(20,*)(fa(ke),ke = 1,ng1)
        irec = 2
        read(20,*)(((((phiantable(m,ke,le,jeta,jt),jt = 1,nnt), &
             jeta = 1,nneta),le = 1,nio),ke = 1,ng1),m = 1,6)
        close(20)

        open(20,file=inelastic_root//'phi_ap'//nfa//'_030x030_a.dat', &
              status='unknown')
        irec = 1
        read(20,*)(fa(ke),ke = 1,ng1)
        irec = 2
        read(20,*)(((((phiaptable(m,ke,le,jeta,jt),jt = 1,nnt), &
             jeta = 1,nneta),le = 1,nio),ke = 1,ng1),m = 1,6)
        close(20)

        irecl = 2*8*ng2*ng2*phi_neta*phi_nt
        irecln = 6*8*ng2*nio*nneta*nnt

        open(20,file=inelastic_root//'phi_me'//nfm//'_030x030_a.dat', &
             form='unformatted', access='direct',recl=irecl)
        irec = 1
        read(20,rec=irec)(fm(ke),ke = 1,ng2)
        irec = 2
        read(20,rec=irec)(((((phimetable(m,ke,le,jeta,jt),jt = 1,phi_nt), &
             jeta = 1,phi_neta),le = 1,ng2),ke = 1,ng2),m = 1,2)
        close(20)

        open(20,file=inelastic_root//'phi_mn'//nfm//'_030x030_a.dat', &
              status='unknown')
        irec = 1
        read(20,*)(fm(ke),ke = 1,ng2)
        irec = 2
        read(20,*)(((((phimntable(m,ke,le,jeta,jt),jt = 1,nnt), &
             jeta = 1,nneta),le = 1,nio),ke = 1,ng2),m = 1,6)
        close(20)

        open(20,file=inelastic_root//'phi_mp'//nfm//'_030x030_a.dat', &
              status='unknown')
        irec = 1
        read(20,*)(fm(ke),ke = 1,ng2)
        irec = 2
        read(20,*)(((((phimptable(m,ke,le,jeta,jt),jt = 1,nnt), &
             jeta = 1,nneta),le = 1,nio),ke = 1,ng2),m = 1,6)
        close(20)

        deallocate(fe)
        deallocate(fa)
        deallocate(fm)

    endif

#if(USE_MPI==TRUE)
    call MPI_BCAST(phieetable, size(phieetable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(phiaetable, size(phiaetable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(phimetable, size(phimetable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(phientable, size(phientable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(phiantable, size(phiantable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(phimntable, size(phimntable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(phieptable, size(phieptable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(phiaptable, size(phiaptable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(phimptable, size(phimptable), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
#endif

  end subroutine init_inelastic_opacity_table

