subroutine unisamp(ipart,velonly)
!========================================================================================
!==== This subroutine develops initial random particle positions and velocities      ====
!==== using basic uniform random number sampling.                                    ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::ipart,velonly
  integer::iparm
  
  double precision::tpar,tlbnd,tubnd,r
!----------------------------------------------------------------------------------------

   do iparm=1,npar
     !
     if (trim(partrans(iparm)) == 'fixed') then
       !
       partval(ipart,iparm) = parval1(iparm)
       !
     else
       !
       if (velonly == 1) then
         !
         if (trim(partrans(iparm)) == 'log') then
           !
           call random_number(r)
           tpar = (log10(parubnd(iparm)) - log10(parlbnd(iparm)))*r + &
                  log10(parlbnd(iparm))
           !
           partval(ipart,iparm) = 1.0d+01**(tpar)
           !
         else if (trim(partrans(iparm)) == 'eqlog') then
           !
           tlbnd = (log10(parlbnd(iparm)))/(log10(base(iparm)))
           tubnd = (log10(parubnd(iparm)))/(log10(base(iparm)))
           !
           call random_number(r)
           tpar = (tubnd - tlbnd)*r + tlbnd
           !
           partval(ipart,iparm) = base(iparm)**(tpar)
           !
         else if (trim(partrans(iparm)) == 'none') then
           !
           call random_number(r)
           partval(ipart,iparm) = (parubnd(iparm) - parlbnd(iparm))*r + &
                                   parlbnd(iparm)
           !
         end if
         !
       end if
       !
       if (trim(partrans(iparm)) /= 'tied') then
         !
!        partvel values and hence, vmax, are based on transformed parameter values
         !
!        if none  is used, vmax is specified based on base parameter values' ranges
!        if log10 is used, vmax should be specified based on orders of magnitude of 
!                             parameter ranges
!        if eqlog is used, vmax is specified as a fraction of maxrange (which was calc'd 
!                             previously in the readpst subroutine)
         !
         call random_number(r)
         partvel(ipart,iparm) = 2.0d+00*vmax(iparm)*r - vmax(iparm)
         !
       end if
       !
     end if
     !
   end do
   !
!  tied parameters
   do iparm=1,npar
     !
     if (trim(partrans(iparm)) == 'tied' .and. velonly == 1) then
       partval(ipart,iparm) = tiedrat(iparm)*partval(ipart,partied(iparm))
     end if
     !
   end do
   !
!  since this is an initialized particle, it's pbest position is set to the initial value
   do iparm=1,npar
     pbest(ipart,iparm) = partval(ipart,iparm)
   end do
  
end subroutine unisamp



subroutine svdjac()
!========================================================================================
!==== This subroutine conducts SVD on Jacobian.                                      ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!---------------------------------------------------------------------------------------- 
  integer::iparm,iobs,lwork,lwmax,info,inull,isol
  
  double precision,dimension(:),allocatable::work,sp
  double precision,dimension(:,:),allocatable::up,vpt,v2t,v1t,jacsvd
!----------------------------------------------------------------------------------------

  lwmax = 100000
  !
  allocate(jacsvd(nobs,nadj),sp(nadj),up(nobs,nadj),vpt(nadj,nadj),v2t(nadj-sdim,nadj),&
           work(lwmax),v1t(sdim,nadj),nproj(nadj,nadj),sproj(nadj,nadj))
  !
  v2t   = 0.0d+00
  v1t   = 0.0d+00
  nproj = 0.0d+00
  sproj = 0.0d+00
  !
  do iparm=1,nadj
    do iobs=1,nobs
      jacsvd(iobs,iparm) = jac(iobs,iparm)
    end do
  end do
!   !
!   open(71,file='test-jac.dat')
!   do iobs=1,nobs
!     write(71,*)(jacsvd(iobs,iparm),iparm=1,100)
!   end do
!   close(71)
  !
  lwork = -1   !<-- query for optimal performance
  info = 0
  work = 0.0d+00
  sp   = 0.0d+00
  up   = 0.0d+00
  vpt  = 0.0d+00
  !
! undertake SVD
  call dgesvd('N','A',nobs,nadj,jacsvd,nobs,sp,up,1,vpt,nadj,&
                  work,lwork,info)
  !
  if (info /= 0) then
    !
    write(*,*)'SVD of Jacobian failed'
    write(*,*)'-- stopping execution --'
    stop
    !
  end if
  !
  lwork = min(lwmax,int(work(1)))
  !
  if (lwmax < int(work(1))) then
    !
    write(*,*)'Need to increase memory for SVD work array'
    write(*,*)'work(1) = ',work(1)
    write(*,*)'-- stopping execution --'
    stop
    !
  end if
  !
  info = 0
  work = 0.0d+00
  sp   = 0.0d+00
  up   = 0.0d+00
  vpt  = 0.0d+00
  !
  call dgesvd('N','A',nobs,nadj,jacsvd,nobs,sp,up,1,vpt,nadj,&
                  work,lwork,info)
  !
  if (info /= 0) then
    !
    write(*,*)'SVD of Jacobian failed'
    write(*,*)'-- stopping execution --'
    stop
    !
  end if
  !
! basis functions spanning the null space
  do inull=1,nadj-sdim
    do iparm=1,nadj
      v2t(inull,iparm) = vpt(inull+sdim,iparm)
    end do
  end do
  !
! basis functions spanning the solution space
  do isol=1,sdim
    do iparm=1,nadj
      v1t(isol,iparm) = vpt(isol,iparm)
    end do
  end do
!   !
!   open(71,file='test-vt.dat')
!   do inull=1,nadj
!     write(71,*)(vpt(iparm,inull),iparm=301,414)
!   end do
!   close(71)
!   !
!   open(71,file='test-s.dat')
!   do iparm=1,nadj
!     write(71,*)sp(iparm)**2
!   end do
!   close(71)
!   stop
  !
! projection matrix for null space
  call dsyrk('U','T',nadj,nadj-sdim,1.0d+00,v2t,nadj-sdim,0.00d+00,nproj,nadj)
  !
! projection matrix for solution space
  call dsyrk('U','T',nadj,sdim,1.0d+00,v1t,sdim,0.00d+00,sproj,nadj)
!   !
!   open(71,file='test-v2v2t.dat')
!   do inull=1,nadj
!     write(71,*)(nproj(inull,iparm),iparm=1,100)
!   end do
!   close(71)
  !
  deallocate(work,sp,up,vpt,v2t,v1t,jacsvd)

end subroutine svdjac



subroutine nullproj(ipart)
!========================================================================================
!==== This subroutine projects randomly generated parameter pertubrations onto       ====
!====   the null space of the inverse problem.                                       ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::ipart
  integer::iparm,icnt,iobs
  
  double precision,dimension(:),allocatable::dp,pdp
!----------------------------------------------------------------------------------------

  allocate(dp(nadj),pdp(nadj))
  !
  dp   = 0.0d+00
  pdp  = 0.0d+00
  icnt = 0
  !
!   open(71,file='test-dp.dat')
! calculate the perturbation in the adjustable parameters
  do iparm=1,npar
    !
    if (trim(partrans(iparm)) == 'log' .or. trim(partrans(iparm)) == 'eqlog' .or. &
        trim(partrans(iparm)) == 'none') then
      !
      icnt = icnt + 1
      !
      dp(icnt) = parval1(iparm) - partval(ipart,iparm)
!       write(71,*)dp(icnt)
      !
    end if
    !
  end do
!   close(71)
  !
! project parameter perturbation (from calibrated parameter set) onto null space
  call dsymv('U',nadj,1.00d+00,nproj,nadj,dp,1,0.0d+00,pdp,1)
  !
  icnt = 0
  !
!   open(71,file='test-pdp.dat')
! calculate the parameter set that results from the projected perturbation
  do iparm=1,npar
    !
    if (trim(partrans(iparm)) == 'log' .or. trim(partrans(iparm)) == 'eqlog' .or. &
        trim(partrans(iparm)) == 'none') then
      !
      icnt = icnt + 1
      !
      partval(ipart,iparm) = parval1(iparm) - pdp(icnt)
!       write(71,*)pdp(icnt)
      !
    end if
    !
  end do
!   close(71)
  !
! for now we deal with parameters going beyond bounds by setting them to their bounds
  do iparm=1,npar
    !
    if (partval(ipart,iparm) > parubnd(iparm)) partval(ipart,iparm) = parubnd(iparm)
    if (partval(ipart,iparm) < parlbnd(iparm)) partval(ipart,iparm) = parlbnd(iparm)
    !
  end do
  !
  deallocate(dp,pdp)

end subroutine nullproj



subroutine solproj(ipart)
!========================================================================================
!==== This subroutine projects PSO generated parameter pertubrations onto the        ====
!====   the solution space of the inverse problem.                                   ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::ipart
  integer::iparm,icnt
  
  double precision,dimension(:),allocatable::dp,pdp
!----------------------------------------------------------------------------------------

  allocate(dp(nadj),pdp(nadj))
  !
  dp   = 0.0d+00
  pdp  = 0.0d+00
  icnt = 0
  !
!   open(71,file='test-dp.dat')
! calculate the perturbation in the adjustable parameters
  do iparm=1,npar
    !
    if (trim(partrans(iparm)) == 'log' .or. trim(partrans(iparm)) == 'eqlog' .or. &
        trim(partrans(iparm)) == 'none') then
      !
      icnt = icnt + 1
      !
      dp(icnt) = parval1(iparm) - partval(ipart,iparm)
!       write(71,*)dp(icnt)
      !
    end if
    !
  end do
!   close(71)
  !
! project parameter perturbation (from calibrated parameter set) onto null space
  call dsymv('U',nadj,1.00d+00,sproj,nadj,dp,1,0.0d+00,pdp,1)
  !
  icnt = 0
  !
!   open(71,file='test-pdp.dat')
! calculate the parameter set that results from the projected perturbation
  do iparm=1,npar
    !
    if (trim(partrans(iparm)) == 'log' .or. trim(partrans(iparm)) == 'eqlog' .or. &
        trim(partrans(iparm)) == 'none') then
      !
      icnt = icnt + 1
      !
      partval(ipart,iparm) = parval1(iparm) - pdp(icnt)
!       write(71,*)pdp(icnt)
      !
    end if
    !
  end do
!   close(71)
  !
! for now we deal with parameters going beyond bounds by setting them to their bounds
  do iparm=1,npar
    !
    if (partval(ipart,iparm) > parubnd(iparm)) partval(ipart,iparm) = parubnd(iparm)
    if (partval(ipart,iparm) < parlbnd(iparm)) partval(ipart,iparm) = parlbnd(iparm)
    !
  end do
  !
  deallocate(dp,pdp)

end subroutine solproj



subroutine transpar(ipart,trans)
!========================================================================================
!==== This subroutine transforms parameters.                                         ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::ipart,trans
  integer::irep,iparm
!----------------------------------------------------------------------------------------

  do iparm=1,npar
    !
    if (trim(partrans(iparm)) == 'log' .and. trans == 1) then
      !
      partval(ipart,iparm) = log10(partval(ipart,iparm))
      pbest(ipart,iparm)   = log10(pbest(ipart,iparm))
      !
    else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 1) then
      !
      partval(ipart,iparm) = (log10(partval(ipart,iparm)))/(log10(base(iparm)))
      pbest(ipart,iparm) = (log10(pbest(ipart,iparm)))/(log10(base(iparm)))
      !
    end if
    !
    if (trim(partrans(iparm)) == 'log' .and. trans == 0) then
      !
      partval(ipart,iparm) = 1.0d+01**(partval(ipart,iparm))
      pbest(ipart,iparm)   = 1.0d+01**(pbest(ipart,iparm))
      !
    else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 0) then
      !
      partval(ipart,iparm) = base(iparm)**(partval(ipart,iparm))
      pbest(ipart,iparm)   = base(iparm)**(pbest(ipart,iparm))
      !
    end if
    !
  end do

end subroutine transpar



subroutine transrep(trans)
!========================================================================================
!==== This subroutine transforms parameters.                                         ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::trans
  integer::irep,iparm
!----------------------------------------------------------------------------------------

  do iparm=1,npar
    !
    if (trim(partrans(iparm)) == 'log' .and. trans == 1) then
      !
      do irep=1,nrep+npop
        if (repindx(irep) == 1) then
          reposit(irep,iparm) = log10(reposit(irep,iparm))
        end if
      end do
      !
    else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 1) then
      !
      do irep=1,nrep+npop
        if (repindx(irep) == 1) then
          reposit(irep,iparm) = (log10(reposit(irep,iparm)))/(log10(base(iparm)))
        end if
      end do
      !
    end if
    !
    if (trim(partrans(iparm)) == 'log' .and. trans == 0) then
      !
      do irep=1,nrep+npop
        if (repindx(irep) == 1) then
          reposit(irep,iparm) = 1.0d+01**(reposit(irep,iparm))
        end if
      end do
      !
    else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 0) then
      !
      do irep=1,nrep+npop
        if (repindx(irep) == 1) then
          reposit(irep,iparm) = base(iparm)**(reposit(irep,iparm))
        end if
      end do
      !
    end if
    !
  end do

end subroutine transrep



subroutine transbnd(trans)
!========================================================================================
!==== This subroutine transforms parameters.                                         ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::trans
  integer::iparm
!----------------------------------------------------------------------------------------

  do iparm=1,npar
    !
    if (trim(partrans(iparm)) == 'log' .and. trans == 1) then
      !
      parlbnd(iparm)       = log10(parlbnd(iparm))
      parubnd(iparm)       = log10(parubnd(iparm))
      !
      parval1(iparm)       = log10(parval1(iparm))
      !
    else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 1) then
      !
      parlbnd(iparm) = (log10(parlbnd(iparm)))/(log10(base(iparm)))
      parubnd(iparm) = (log10(parubnd(iparm)))/(log10(base(iparm)))
      !
      parval1(iparm) = (log10(parval1(iparm)))/(log10(base(iparm)))
      !
    end if
    !
    if (trim(partrans(iparm)) == 'log' .and. trans == 0) then
      !
      parlbnd(iparm)       = 1.0d+01**(parlbnd(iparm))
      parubnd(iparm)       = 1.0d+01**(parubnd(iparm))
      !
      parval1(iparm)       = 1.0d+01**(parval1(iparm))
      !
    else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 0) then
      !
      parlbnd(iparm)       = base(iparm)**(parlbnd(iparm))
      parubnd(iparm)       = base(iparm)**(parubnd(iparm))
      !
      parval1(iparm)       = base(iparm)**(parval1(iparm))
      !
    end if
    !
  end do
  
end subroutine transbnd


subroutine transparunc(nsav,trans)
!========================================================================================
!==== This subroutine transforms parameters.                                         ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::trans,nsav
  integer::iparm,isav
!----------------------------------------------------------------------------------------

  do iparm=1,npar
    !
    if (trim(partrans(iparm)) == 'log' .and. trans == 1) then
      !
      do isav=1,nsav
        parsav(iparm,isav) = log10(parsav(iparm,isav))
      end do
      !
    else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 1) then
      !
      do isav=1,nsav
        parsav(iparm,isav) = (log10(parsav(iparm,isav)))/(log10(base(iparm)))
      end do
      !
    end if
    !
    if (trim(partrans(iparm)) == 'log' .and. trans == 0) then
      !
      do isav=1,nsav
        parsav(iparm,isav) = 1.0d+01**(parsav(iparm,isav))
      end do
      !
    else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 0) then
      !
      do isav=1,nsav
        parsav(iparm,isav) = base(iparm)**(parsav(iparm,isav))
      end do
      !
    end if
    !
  end do
  
end subroutine transparunc







