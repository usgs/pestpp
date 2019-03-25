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
  
  double precision::tpar,tlbnd,tubnd, r1
!----------------------------------------------------------------------------------------

   do iparm=1,npar
     !
     if (trim(partrans(iparm)) == 'fixed') then
       !
       partval(ipart,iparm) = parval1(iparm)
       pbest(ipart,iparm) = parval1(iparm)       
       !
     else
       !
       if (velonly == 1) then
         !
         if (trim(partrans(iparm)) == 'log') then
           !
           call random_number(r1)
           tpar = (log10(parubnd(iparm)) - log10(parlbnd(iparm)))*r1 + &
                  log10(parlbnd(iparm))
           !
           partval(ipart,iparm) = 1.0d+01**(tpar)
           !
         else if (trim(partrans(iparm)) == 'eqlog') then
           !
           tlbnd = (log10(parlbnd(iparm)))/(log10(base(iparm)))
           tubnd = (log10(parubnd(iparm)))/(log10(base(iparm)))
           !
           call random_number(r1)
           tpar = (tubnd - tlbnd)*r1 + tlbnd
           !
           partval(ipart,iparm) = base(iparm)**(tpar)
           !
         else if (trim(partrans(iparm)) == 'none') then
           !
           call random_number(r1)
           partval(ipart,iparm) = (parubnd(iparm) - parlbnd(iparm))*r1 + &
                                   parlbnd(iparm)
           !
         end if
         !
       end if
       !
       if (trim(partrans(iparm)) /= 'tied') then
         !
!        since this is the initial condition, set pbest = partval
         pbest(ipart,iparm) = partval(ipart,iparm)
         !
!        partvel values and hence, vmax, are based on transformed parameter values
         !
!        if none  is used, vmax is specified based on base parameter values' ranges
!        if log10 is used, vmax should be specified based on orders of magnitude of 
!                             parameter ranges
!        if eqlog is used, vmax is specified as a fraction of maxrange (which was calc'd 
!                             previously in the readpst subroutine)
         !
         call random_number(r1)
         partvel(ipart,iparm) = 2.0d+00*vmax(iparm)*r1 - vmax(iparm)
         !
       end if
       !
     end if
     !
   end do
   
   do iparm=1,npar
     !
     if (trim(partrans(iparm)) == 'tied' .and. velonly == 1) then
       !
       partval(ipart,iparm) = tiedrat(iparm)*partval(ipart,partied(iparm))
       pbest(ipart,iparm) = partval(ipart,iparm)
       !
     end if
     !
   end do
  
end subroutine unisamp



subroutine transpar(iparm,trans)
!========================================================================================
!==== This subroutine transforms parameters.                                         ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::iparm,trans
  integer::ipart,irep
!----------------------------------------------------------------------------------------

  if (trim(partrans(iparm)) == 'log' .and. trans == 1) then
    !
    do ipart=1,npop
      partval(ipart,iparm) = log10(partval(ipart,iparm))
      pbest(ipart,iparm)   = log10(pbest(ipart,iparm))
    end do
    !
    parlbnd(iparm)         = log10(parlbnd(iparm))
    parubnd(iparm)         = log10(parubnd(iparm))
    !
  else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 1) then
    !
    do ipart=1,npop
      partval(ipart,iparm) = (log10(partval(ipart,iparm)))/(log10(base(iparm)))
      pbest(ipart,iparm) = (log10(pbest(ipart,iparm)))/(log10(base(iparm)))
    end do
    !
    parlbnd(iparm) = (log10(parlbnd(iparm)))/(log10(base(iparm)))
    parubnd(iparm) = (log10(parubnd(iparm)))/(log10(base(iparm)))
    !
  end if
  
  if (trim(partrans(iparm)) == 'log' .and. trans == 0) then
    !
    do ipart=1,npop
      partval(ipart,iparm) = 1.0d+01**(partval(ipart,iparm))
      pbest(ipart,iparm)   = 1.0d+01**(pbest(ipart,iparm))
    end do
    !
    parlbnd(iparm)         = 1.0d+01**(parlbnd(iparm))
    parubnd(iparm)         = 1.0d+01**(parubnd(iparm))
    !
  else if (trim(partrans(iparm)) == 'eqlog' .and. trans == 0) then
    !
    do ipart=1,npop
      partval(ipart,iparm) = base(iparm)**(partval(ipart,iparm))
      pbest(ipart,iparm)   = base(iparm)**(pbest(ipart,iparm))
    end do
    !
    parlbnd(iparm)         = base(iparm)**(parlbnd(iparm))
    parubnd(iparm)         = base(iparm)**(parubnd(iparm))
    !
  end if
  
  if (trim(pestmode) == 'pareto') then
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
  end if
 
  
end subroutine transpar







