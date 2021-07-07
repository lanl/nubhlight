! ======================================================================
! copyright 2020. Triad National Security, LLC. All rights
! reserved. This program was produced under U.S. Government contract
! 89233218CNA000001 for Los Alamos National Laboratory (LANL), which
! is operated by Triad National Security, LLC for the U.S. Department
! of Energy/National Nuclear Security Administration. All rights in
! the program are reserved by Triad National Security, LLC, and the
! U.S. Department of Energy/National Nuclear Security
! Administration. The Government is granted for itself and others
! acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display
! publicly, and to permit others to do so.
! ======================================================================

! Originally written by Josh Dolence
! Table reader for opacity/emissivity tables by Adam Burrows
! First published in Burrows, Reddy, Thompson.
! DOI: 10.1016/j.nuclphysa.2004.06.012

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
