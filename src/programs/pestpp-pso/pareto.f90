subroutine psopareto(basnam)
!========================================================================================
!==== This subroutine executes basic pso in estimation mode.                         ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================
  use psodat
#ifdef __INTEL_COMPILER
  use ifport
#endif  

  implicit none

! specifications:
!----------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------- 
! interfaces
!----------------------------------------------------------------------------------------
  interface 
  subroutine readrst(basnam)
  integer::ipart,iparm,irep,iobs,igp,ipto
  integer,dimension(:),allocatable::measgp
  
  character(len=100), optional, intent(in)::basnam
  character(len=100)::scrc,fnam
  end subroutine readrst
  end interface
!---------------------------------------------------------------------------------------- 
! external routines for run management via PANTHER
!----------------------------------------------------------------------------------------
  external rmif_run
  integer rmif_run
  external rmif_get_num_failed_runs
  integer rmif_get_num_failed_runs
  external rmif_get_failed_run_ids
  integer rmif_get_failed_run_ids
  external rmif_get_num_total_runs
  integer rmif_get_num_total_runs
  external rmif_delete
  integer rmif_delete
  
!----------------------------------------------------------------------------------------
! local PSO-main variables
!----------------------------------------------------------------------------------------
  integer::err,irun,iiter,ipart,iparm,iparm1,iparm2,repred,iptout,iobs,ipto,irep,fail,&
           inpar
  
  double precision::alpha
  
  character(len=100),intent(in)::basnam
  integer :: n_seed, ierr
  integer, dimension(:), allocatable :: a_seed
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

  err      = 0
  repred   = 0
  iptout   = 1
  fail     = 0
  alpha    = 0.0d+00
  !
! remove all repository output files from previous runs
#ifdef __windows__  
!  ierr = system('del *.rep')
#else
!  ierr = system('rm ./*.rep')
#endif  
  !
! restart from previous run if requested
  if (rstpso == 1) then
    !
    call readrst(basnam)
    !
    if (nrepact > 1) then
      !
!     manage repository and calculate fitness
      if (repmode == 1) then
        !
!       use grid-based method described by Coello et al, 2004
        call repgrid(repred)
        !
      else if (repmode == 2) then
        !
!       use loneliness method 
        call loneliness(repred,alpha)
        !
      end if
      !
    else if (nrepact == 1) then
      !
      do irep=1,npop+nrep
        if (repindx(irep) == 1) fitness(irep) = 1.00d+00
      end do
      !
    end if
    !
!   write initial conditions to record file
    call listipt(basnam)
    !
    call listpto(-1,iptout,basnam,alpha)
    !
  else
    !
!   write initial conditions to record file
    call listipt(basnam)
    !
!-- generate random parameter sets and velocities, set pbest, and add corresponding runs to the queue  
    call random_seed (size = n_seed)
    allocate(a_seed(1:n_seed))
    a_seed = iseed
    call random_seed(put = a_seed)
    deallocate(a_seed)
    !
    inpar = 0
    !
    do ipart=1,npop
      !
      if (initp == 1 .and. ipart == 1) then
        !
        do iparm=1,npar
          partval(ipart,iparm) = parval1(iparm)
        end do
        !
        call unisamp(ipart,0)
        !
      else if (initp == 2) then
        !
        inpar = inpar + 1
        !
        do iparm1=1,npar
          !
          do iparm2=1,npar
            !
            if (trim(parnme(iparm1)) == trim(iprnme(iparm2))) then
              partval(ipart,iparm1) = (iprval(iparm2,inpar) - offset(iparm1))/scale(iparm1)
              exit
            end if
            !
            if (iparm2 == npar) then
              !
              write(*,'(A,A/A)')'Parameter names in Pareto parameter file do not match',&
                                ' in control file','-- stopping execution --'
              stop
              !
            end if
            !
          end do
          !
        end do
        !
        call unisamp(ipart,0)
        !
        if (inpar == nitp) initp = 0
        !
      else
        !
        call unisamp(ipart,1)
        !
      end if
      !
      do iparm=1,npar
        parval(iparm) = scale(iparm)*partval(ipart,iparm) + offset(iparm)
      end do
      !
!     add model runs to the queue
      call modelrm(1,irun,fail)
      !
    end do
    !
!-- execute model runs
    err = rmif_run()
    !
!-- evaluate each initial objective value for each particle
    do ipart=1,npop
      !
      modfail(ipart) = 0
      !
!     get model run results
      irun = ipart-1
      call modelrm(0,irun,fail)
      !
      modfail(ipart) = fail
      !
      if (fail == 0) then
        !
        call ptoobjeval(ipart)
        !
!       save observations for recording purposes
        do iobs=1,nobs
          mobssav(ipart,iobs) = mobsval(iobs)
        end do
        !
      else
        !
        do ipto=1,nptogp
          ptogp(ipart,ipto) = 1.00d+30
        end do
        !
        if (nptocon > 0) violate(ipart) = 1.00d+30
        !
      end if
      !
      if (nptocon > 0) then
        !
!       since this is the initial condition, vioopt set directly to initial violate
        vioopt(ipart) = violate(ipart)
        !
      end if
      !
!     since this is the initial condition, ptogpopts set directly to their objs
      do ipto=1,nptogp
        ptogpopt(ipart,ipto) = ptogp(ipart,ipto)
      end do
      !
    end do
    !
!   check if maximum allowable model failures has been exceeded
    fail = 0
    !
    do ipart=1,npop
      fail = fail + modfail(ipart)
    end do
    !
    if (fail > nforg) then
      write(*,'(I0,A,I0)')fail,' failed model runs greater than limit set by user, ',nforg
      write(*,'(A)')'-- stopping execution --'
      stop
    end if
    !
!   update pareto repository
    call repcont(0,repred)
    !
    if (nrepact > 1) then
      !
!     manage repository and calculate fitness
      if (repmode == 1) then
        !
!       use grid-based method described by Coello et al, 2004
        call repgrid(repred)
        !
      else if (repmode == 2) then
        !
!       use loneliness method 
        call loneliness(repred,alpha)
        !
      end if
      !
    else if (nrepact == 1) then
      !
      do irep=1,npop+nrep
        if (repindx(irep) == 1) fitness(irep) = 1.00d+00
      end do
      !
    end if
    !
!   list output for pareto
    call listpto(0,iptout,basnam,alpha)
    !
!-- write restart data if requested
    if (trim(rstfle) == 'restart') then
      call wrparerst(basnam)
    end if
    !
!   stop if noptmax = 0
    if (noptmax == 0) then
      !
      err = rmif_delete()
      !
      call wrfinpareto(basnam)
      !
      stop
      !
    end if
    !
  end if
 

! begin the PSO iterative procedure
! ---------------------------------
  do iiter=1,noptmax
    !
!-- write message to terminal
    call witmess(iiter)
    !
!-- calculate velocities and update particle positions  
    call pertpareto(iiter)
    !
!-- reinitialize run manager and make another set of runs
    call initialrm(1)
    !
!-- add model runs to the queue  
    do ipart=1,npop
      !
      do iparm=1,npar
        parval(iparm) = scale(iparm)*partval(ipart,iparm) + offset(iparm)
      end do
      !
      call modelrm(1,irun,fail)
      !
    end do
    !
!-- execute model runs
    err = rmif_run()
    !
!-- evaluate objective value for each particle, find pbest, gbest and gindex
    do ipart=1,npop
      !
      modfail(ipart) = 0
      !
!---- get model run results
      irun = ipart-1
      call modelrm(0,irun,fail)
      !
      modfail(ipart) = fail
      !
      if (fail == 0) then
        !
        call ptoobjeval(ipart)
        !
        do iobs=1,nobs
          mobssav(ipart,iobs) = mobsval(iobs)
        end do
        !
      else
        !
        do ipto=1,nptogp
          ptogp(ipart,ipto) = 1.00d+30
        end do
        !
        if (nptocon > 0) violate = 1.00d+30
        !
      end if
      !
    end do
    !
!   check if maximum allowable model failures has been exceeded
    fail = 0
    !
    do ipart=1,npop
      fail = fail + modfail(ipart)
    end do
    !
    if (fail > nforg) then
      write(*,'(I0,A,I0)')fail,' failed model runs greater than limit set by user, ',nforg
      write(*,'(A)')'-- stopping execution --'
      stop
    end if
    !
!-- update pareto repository
    call repcont(iiter,repred)
    !
    if (nrepact > 1) then
      !
!     manage repository and calculate fitness
      if (repmode == 1) then
        !
!       use grid-based method described by Coello et al, 2004
        call repgrid(repred)
        !
      else if (repmode == 2) then
        !
!       use loneliness method 
        call loneliness(repred,alpha)
        !
      end if
      !
    else if (nrepact == 1) then
      !
      do irep=1,npop+nrep
        if (repindx(irep) == 1) fitness(irep) = 1.00d+00
      end do
      !
    end if
    !
!-- set pbest
    call pbestpareto(iiter)
    !
!-- list output for pareto
    call listpto(iiter,iptout,basnam,alpha)
    !
!-- write restart data if requested
    if (trim(rstfle) == 'restart') then
      call wrparerst(basnam)
    end if
    !
  end do
  !
! end main loop of PSO iterative procedure
!-----------------------------------------
  !
! write out final repository data
  call wrfinpareto(basnam)  
  
end subroutine psopareto
















