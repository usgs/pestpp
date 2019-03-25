subroutine psoest(basnam)
!========================================================================================
!==== This subroutine executes basic pso in estimation mode.                         ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================
  use psodat

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
  integer::err,irun,iiter,ipart,iparm,iparm1,iparm2,inpar,gindex,igp,fail,iphistp
  
  double precision::gbest,gmbest,gpbest,objmin
  
  character(len=100),intent(in)::basnam
  integer :: n_seed
  integer, dimension(:), allocatable :: a_seed
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

  err     = 0
  fail    = 0
  gindex  = 0
  iphistp = 0
  gbest   = 0.0d+00
  gmbest  = 0.0d+00
  gpbest  = 0.0d+00
  
! restart from previous run if requested
  if (rstpso == 1) then
    !
    call readrst()
    !
!-- set gbest
    call getgbest(gindex,gbest,gmbest,gpbest)
    !
    objmin = gbest
    !
!   write initial conditions to record file
    call listini(basnam,gindex)
    !
  else
    !
!   write initial conditions to record file
    call listini(basnam,gindex)
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
        call estobjeval(ipart)
        !
!       since this is the initial condition, objopts set direct to their objs
        objmopt(ipart) = objm(ipart)
        objpopt(ipart) = objp(ipart)
        objopt(ipart)  = obj(ipart)
        do igp=1,nobsgp
          objgpopt(ipart,igp) = objgp(ipart,igp)
        end do
        !
      else
        !
!       since this is the initial condition, objopts set to a massive value
        objmopt(ipart) = 1.00d+30
        objpopt(ipart) = 1.00d+30
        objopt(ipart)  = 1.00d+30
        do igp=1,nobsgp
          objgpopt(ipart,igp) = 1.00d+30
        end do
        !
      end if
      !
    end do
    !
!   set gbest
    call getgbest(gindex,gbest,gmbest,gpbest)
    !
    objmin = gbest
    !  
!   list output for initial conditions
    call listout(0,gbest,gmbest,gpbest,gindex,objmin,iphistp)
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
      err = rmif_delete()
      stop
    end if
    !
!   stop if noptmax = 0
    if (noptmax == 0) then
      err = rmif_delete()
      stop
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
    call pertpart(iiter,gindex)
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
        call estobjeval(ipart)
        !
      else
        !
        obj(ipart)  = 1.00d+30
        objm(ipart) = 1.00d+30
        objp(ipart) = 1.00d+30
        do igp=1,nobsgp
          objgp(ipart,igp) = 1.00d+30
        end do
        !
      end if
      !
    end do
    !
!-- set pbest
    call getpbest()
    !
!-- set gbest
    call getgbest(gindex,gbest,gmbest,gpbest)
    !
!-- write restart data if requested
    if (trim(rstfle) == 'restart') then
      call writerst(basnam)
    end if
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
!   check relative phi reduction for termination criteria
    if ((1.0d+00 - gbest/objmin) <= phiredstp) then
      iphistp = iphistp + 1
    else
      iphistp = 0
    end if
    !
!-- list iteration output
    call listout(iiter,gbest,gmbest,gpbest,gindex,objmin,iphistp)
    !
    if (iphistp == nphistp) then
      !
!     swarm is no longer reducing phi significantly, terminate
      write(*,'(A)')'termination due to phiredstp criteria'
      write(*,'(A)')'-- stopping execution --'
      !
      exit
      !
    end if
    !
    objmin = gbest
    !
  end do
! end main loop of PSO iterative procedure
!-----------------------------------------
  
  
! running model one last time with best parameters
! reinitialize run manager
  call initialrm(1)
  !  
! add model run to the queue  
  do iparm=1,npar
    parval(iparm) = scale(iparm)*pbest(gindex,iparm) + offset(iparm)
  end do
  !
  call modelrm(1,irun,fail)
  !    
! execute model run
  err = rmif_run()
  !
! get model run results
  irun = 0
  call modelrm(0,irun,fail)
  !
! write pbest and gbest to file
  call writebest(gindex,basnam)
  
  
end subroutine psoest