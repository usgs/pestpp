subroutine parunc(basnam)
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
! external routines for run management via YAMR
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
  integer::err,irun,iiter,ipart,iparm,iparm1,iparm2,inpar,gindex,igp,fail,iphistp,nsav,&
           istop,iobs
  
  double precision::gbest,gmbest,gpbest
  
  character(len=100),intent(in)::basnam
  character(len=100)::rdwr
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------

  err     = 0
  fail    = 0
  gindex  = 0
  iphistp = 0
  nsav    = 0
  istop   = 0
  gbest   = 0.0d+00
  gmbest  = 0.0d+00
  gpbest  = 0.0d+00
  !
! read Jacobian file if using SVD
  if (sdim > 0) then
    !
    call readjac()
    !
    call svdjac()
    !
  end if
  !
! restart from previous run if requested
  if (rstpso == 1) then
    !
    call readrst()
    !
!-- set gbest
    call getgbest(gindex,gbest,gmbest,gpbest)
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
    call random_seed()
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
        if (sdim > 0) then
!           !
!           open(71,file='test_pso-randpar1.dat')
!           write(71,*)'single point'
!           do iparm=1,npar
!             write(71,*)parnme(iparm),partval(ipart,iparm),1.00d+00,0.0d+00
!           end do
!           close(71)
!           !
          call transpar(ipart,1)
          call transbnd(1)
          !
!-------- project random parameter perturbation onto null space
          call nullproj(ipart)
          !
          call transpar(ipart,0)
          call transbnd(0)
          !
!-------- retie projected child parameters to their parents
          do iparm=1,npar
            if (trim(partrans(iparm)) == 'tied') then
              partval(ipart,iparm) = tiedrat(iparm)*partval(ipart,partied(iparm))
            end if
          end do
          !
!-------- reset pbest to projected starting particle position
          do iparm=1,npar
            pbest(ipart,iparm) = partval(ipart,iparm)
          end do
!           !
!           open(71,file='test-pso-prandpar1.dat')
!           do iparm=1,npar
!             write(71,*)partval(ipart,iparm)
!           end do
!           close(71)
!           !
!           stop
        end if
        !
      end if
      !
      do iparm=1,npar
        parval(iparm) = scale(iparm)*partval(ipart,iparm) + offset(iparm)
      end do
      !
!---- add model runs to the queue
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
      call modelrm(0,ipart-1,fail)
      !
      modfail(ipart) = fail
      !
      if (fail == 0) then
        !
        call estobjeval(ipart)
        !
!       save observations for recording purposes
        do iobs=1,nobs
          mobssav(ipart,iobs) = mobsval(iobs)
        end do
        !
!       since this is the initial condition, objopts set directly to their objs
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
!-- set gbest
    call getgbest(gindex,gbest,gmbest,gpbest)
    !  
!-- list general iteration output
    call listout(0,gbest,gmbest,gpbest,gindex,0.0d+00,iphistp)
    !
!-- check if maximum allowable model failures has been exceeded
    fail = 0
    !
    do ipart=1,npop
      fail = fail + modfail(ipart)
    end do
    !
    if (fail > nforg) then
      write(rdwr,'(I0,A,I0)')fail,' failed model runs greater than limit set by user, ',nforg
      call killrm(trim(rdwr))
    end if
    !
!-- save any calibrated particles
    do ipart=1,npop
      !
      if (objopt(ipart) <= thresh) then
        !
        nsav = nsav + 1
        !
        do iparm=1,npar
          parsav(iparm,nsav) = pbest(ipart,iparm)
        end do
        !
        do igp=1,nobsgp
          objgpsav(igp,nsav) = objgpopt(ipart,igp)
        end do
        !
        do iobs=1,nobs
          punobs(iobs,nsav) = mobssav(ipart,iobs)
        end do
        !
        if (nsav == nreal) then
          !
          istop = 1
          !
          exit
          !
        end if
        !
      end if
      !
    end do
    !
!-- write details to record file
    call listpun(nsav)
    !
!-- write parsav and objgpsav to file
    call wparunc(basnam,nsav)
    !
!   stop if noptmax = 0
    if (noptmax == 0) call killrm('Intitial swarm complete (NOPTMAX = 0)')
    !
  end if
  !
  !
  !
  !
! Run the no-calibration NSMC process (i.e., no PSO is used)
  if (nopso == 1) then
    !
    do iiter=1,noptmax
      !
      !
      !
      if (istop == 1) then
        !
        write(rdwr,'(I0,A)')nreal,' realizations obtained.'
        call killrm(trim(rdwr))
        !
      end if
      !
      do ipart=1,npop
        !
        call unisamp(ipart,1)
        !
        if (sdim > 0) then
          !
          call transpar(ipart,1)
          call transbnd(1)
          !
!-------- project random parameter perturbation onto null space
          call nullproj(ipart)
          !
          call transpar(ipart,0)
          call transbnd(0)
          !
!-------- retie projected child parameters to their parents
          do iparm=1,npar
            if (trim(partrans(iparm)) == 'tied') then
              partval(ipart,iparm) = tiedrat(iparm)*partval(ipart,partied(iparm))
            end if
          end do
          !
        end if
        !
        do iparm=1,npar
          parval(iparm) = scale(iparm)*partval(ipart,iparm) + offset(iparm)
        end do
        !
!------ add model runs to the queue
        call modelrm(1,irun,fail)
        !
      end do
      !
!---- execute model runs
      err = rmif_run()
      !
!---- evaluate each initial objective value for each particle
      do ipart=1,npop
        !
        modfail(ipart) = 0
        !
!       get model run results
        call modelrm(0,ipart-1,fail)
        !
        modfail(ipart) = fail
        !
        if (fail == 0) then
          !
          call estobjeval(ipart)
          !
!         save observations for recording purposes
          do iobs=1,nobs
            mobssav(ipart,iobs) = mobsval(iobs)
          end do
          !
        end if
        !
      end do
      !
!---- check if maximum allowable model failures has been exceeded
      fail = 0
      !
      do ipart=1,npop
        fail = fail + modfail(ipart)
      end do
      !
      if (fail > nforg) then
        write(rdwr,'(I0,A,I0)')fail,' failed model runs greater than limit set by user, ',nforg
        call killrm(trim(rdwr))
      end if
      !
!---- save any calibrated particles
      do ipart=1,npop
        !
        if (obj(ipart) <= thresh) then
          !
          nsav = nsav + 1
          !
          do iparm=1,npar
            parsav(iparm,nsav) = partval(ipart,iparm)
          end do
          !
          do igp=1,nobsgp
            objgpsav(igp,nsav) = objgp(ipart,igp)
          end do
          !
          do iobs=1,nobs
            punobs(iobs,nsav) = mobssav(ipart,iobs)
          end do
          !
          if (nsav == nreal) then
            !
            istop = 1
            !
            exit
            !
          end if
          !
        end if
        !
      end do
      !
!---- write details to record file
      call listpun(nsav)
      !
!---- write parsav and objgpsav to file
      call wparunc(basnam,nsav)
      !
      !
      !
    end do
    !
    write(rdwr,'(I0,A)')noptmax,' NSMC iterations reached'
    call killrm(trim(rdwr))
    !
  end if
  




! begin the PSO iterative procedure
! ---------------------------------
  do iiter=1,noptmax
    !
!-- write message to terminal
    call witmess(iiter)
    !
    if (istop == 1) then
      !
      write(rdwr,'(I0,A)')nreal,' realizations obtained.'
      call killrm(trim(rdwr))
      !
    end if
    !
    do ipart=1,npop
      !
      if (objopt(ipart) <= thresh .or. nopso == 1) then
        !
        call unisamp(ipart,1)
        !
        if (sdim > 0) then
          !
          call transpar(ipart,1)
          call transbnd(1)
          !
!-------- project random parameter perturbation onto null space
          call nullproj(ipart)
          !
          call transbnd(0)
          call transpar(ipart,0)
          !
!-------- retie projected child parameters to their parents
          do iparm=1,npar
            if (trim(partrans(iparm)) == 'tied') then
              partval(ipart,iparm) = tiedrat(iparm)*partval(ipart,partied(iparm))
            end if
          end do
          !
!-------- reset pbest to projected starting particle position
          do iparm=1,npar
            pbest(ipart,iparm) = partval(ipart,iparm)
          end do
          !
        end if
        !
        objopt(ipart) = 1.00d+15
        !
      end if
      !
    end do
    !
!-- reset gbest after calibrated particles have been replaced with random NSMC realizations
    call getgbest(gindex,gbest,gmbest,gpbest)
    !
    if (neibr == 1) then
      call wrtneib()
    else
      write(unit(1),'(A,I0)')'g-best position is now ',gindex
    end if
    !
!-- calculate velocities and update particle positions  
    call pertunc(iiter,gindex)
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
      call modelrm(0,ipart-1,fail)
      !
      modfail(ipart) = fail
      !
      if (fail == 0) then
        !
        call estobjeval(ipart)
        !
!       save observations for recording purposes
        do iobs=1,nobs
          mobssav(ipart,iobs) = mobsval(iobs)
        end do
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
      write(rdwr,'(I0,A,I0)')fail,' failed model runs greater than limit set by user, ',nforg
      call killrm(trim(rdwr))
    end if
    !
!-- get gbest for listing output
    call getgbest(gindex,gbest,gmbest,gpbest)
    !
!-- list general iteration output
    call listout(iiter,gbest,gmbest,gpbest,gindex,0.0d+00,iphistp)
    !
!-- save any calibrated particles
    do ipart=1,npop
      !
      if (objopt(ipart) <= thresh) then
        !
        nsav = nsav + 1
        !
        do iparm=1,npar
          parsav(iparm,nsav) = pbest(ipart,iparm)
        end do
        !
        do igp=1,nobsgp
          objgpsav(igp,nsav) = objgpopt(ipart,igp)
        end do
        !
        do iobs=1,nobs
          punobs(iobs,nsav) = mobssav(ipart,iobs)
        end do
        !
        if (nsav == nreal) then
          !
          istop = 1
          !
          exit
          !
        end if
        !
      end if
      !
    end do
    !
!-- write details to record file
    call listpun(nsav)
    !
!-- write parsav and objgpsav to file
    call wparunc(basnam,nsav)
    !
  end do
! end main loop of PSO iterative procedure
!-----------------------------------------
  
  
end subroutine parunc
