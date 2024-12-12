subroutine pertpart(iiter,gindex)
!========================================================================================
!==== This subroutine perturbs particles randomly based on PSO algorithm.            ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------
  integer,intent(in)::gindex,iiter
  integer::ipart,iparm,spin
  
  double precision::local,global,r
!----------------------------------------------------------------------------------------    

! calculate current value for inertia
  if (iiter < inerti) then
    inertia = iinert + (finert - iinert)*(dble(iiter)/dble(inerti))
  else
    inertia = finert
  end if
  !
  !
! transform all parameters for all particles
  call transbnd(1)
  !
  do ipart=1,npop
    call transpar(ipart,1)
  end do
  !
  !
! perturb particles in their transformed rhealm
  do ipart=1,npop
    !
    do iparm=1,npar
      !
      if (trim(partrans(iparm)) == 'none' .or. trim(partrans(iparm)) == 'log' .or. &
          trim(partrans(iparm)) == 'eqlog') then
        !
        if (neibr == 0) then
          !
          call random_number(r)
          local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r)
          global = c2*r*(pbest(gindex,iparm) - partval(ipart,iparm))
          !
        else
          !
          call random_number(r)
          local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r)
          global = c2*r*(pbest(gneibr(ipart),iparm) - partval(ipart,iparm))
          !
        end if
        !
!       handle precision issues for when particle stops moving
        if (dabs(local) < 1.00d-32) then
          local = 0.0d+00
        end if
        if (dabs(global) < 1.00d-32) then
          global = 0.0d+00
        end if
        if (dabs(partvel(ipart,iparm)) < 1.0d-32) then
          partvel(ipart,iparm) = 0.0d+00
        end if
        !
        !
!       calculate particle velocity
        partvel(ipart,iparm) = inertia*partvel(ipart,iparm) + local + global 
        !
        if (dabs(partvel(ipart,iparm)) > vmax(iparm)) then
          if (partvel(ipart,iparm) > 0.0d+00) partvel(ipart,iparm) = vmax(iparm)
          if (partvel(ipart,iparm) < 0.0d+00) partvel(ipart,iparm) = -vmax(iparm)
        end if
        !
        !
!       initial particle perturbation
        partval(ipart,iparm) = partval(ipart,iparm) + partvel(ipart,iparm)
        !
        !
!       keep perturbing particles if they violate parameter bounds
        if (partval(ipart,iparm) > parubnd(iparm) .or. &
          partval(ipart,iparm) < parlbnd(iparm)) then
          !
          spin = 0
          !
          do
            ! 
            spin = spin + 1
            !
            if (neibr == 0) then
              !
              call random_number(r)
              local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
              call random_number(r)
              global = c2*r*(pbest(gindex,iparm) - partval(ipart,iparm))
              !
            else
              !
              call random_number(r)
              local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
              call random_number(r)
              global = c2*r*(pbest(gneibr(ipart),iparm) - partval(ipart,iparm))
              !
            end if
            !
            partvel(ipart,iparm) = inertia*partvel(ipart,iparm) + local + global
            !
            if (dabs(partvel(ipart,iparm)) > vmax(iparm)) then
              if (partvel(ipart,iparm) > 0.0d+00) partvel(ipart,iparm) = vmax(iparm)
              if (partvel(ipart,iparm) < 0.0d+00) partvel(ipart,iparm) = -vmax(iparm)
            end if
            !
            partval(ipart,iparm) = partval(ipart,iparm) + partvel(ipart,iparm)
            !
            if (partval(ipart,iparm) < parubnd(iparm) .and. &
                partval(ipart,iparm) > parlbnd(iparm)) exit
            !
            if (spin > 1000) then
              !
              write(*,*)"I'm spinning my wheels trying to find a feasible parameter value:"
              write(*,*)'-- stopping execution --'
              write(*,*)'ipart ',ipart,'  iparm ',iparm
              write(*,*)'partval ',partval(ipart,iparm)
              write(*,*)'partvel ',partvel(ipart,iparm)
              stop
              !
            end if
            !
          end do
          !
        end if
        !
      end if
      !
      !
    end do
    !
  end do
  !
! transform all parameters for all particles back to original form
  do ipart=1,npop
    call transpar(ipart,0)
  end do
  !
  call transbnd(0)
  !
  !
! tie child parameters to their parents
  do ipart=1,npop
    do iparm=1,npar
      if (trim(partrans(iparm)) == 'tied') then
        !
        partval(ipart,iparm) = tiedrat(iparm)*partval(ipart,partied(iparm))
        !
      end if
    end do
  end do          
    
end subroutine pertpart




subroutine pertpareto(iiter)
!========================================================================================
!==== This subroutine perturbs particles randomly based on MOPSO algorithm.          ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------
  integer,intent(in)::iiter
  integer::i,ipart,iparm,irep,repcnt,spin,k1,k2,gindex
  
  double precision::local,global,r1,r2,viomin,r
!----------------------------------------------------------------------------------------    

! calculate current value for inertia
  if (iiter < inerti) then
    inertia = iinert + (finert - iinert)*(dble(iiter)/dble(inerti))
  else
    inertia = finert
  end if
  !  
! transform parameters
  call transbnd(1)
  call transrep(1)
  !
  do ipart=1,npop
    call transpar(ipart,1)
  end do
  !
  if (nrepact == 1) then
    !
!   there is only a single particle in the repository, set gbest to this particle
    do irep=1,nrep+npop
      if (repindx(irep) == 1) gindex = irep
    end do
    !
  end if
  !
  if (nrepact == 0) then
    !
    if (nptocon > 0) then
      !
!     there are no feasible particles, so the repository is empty. set gbest as the 
!     -- pbest with the smallest vioopt
      i = 0
      viomin = 0.0d+00
      !
      do ipart=1,npop
        !
        i = i + 1
        !
        if (vioopt(ipart) < viomin .or. i == 1) then
          !
          viomin = vioopt(ipart)
          gindex = ipart
          !
        end if
        !
      end do
      !
    else
      !
      write(*,*)'All particle positions violate objective limits (PTOUB)'
      write(*,*)'-- stopping execution --'
      stop
      !
    end if
    !
  end if
  !
! perturb particles in their transformed rhealm
  do ipart=1,npop
    !
    if (nrepact > 1) then
      !
!     select gbest particle randomly from the repository using roulette-wheel based on fitness
      do
        !
        repcnt = 0
        !
!       select random repository position
        call random_number(r)
        r1 = r
        k1 = int(r1*dble(nrepact)) + 1
        !
        do irep=1,nrep+npop
          !
          if (repindx(irep) == 1) repcnt = repcnt + 1
          !
          if (repcnt == k1) then
            k2 = irep
            exit
          end if
          !
        end do
        !
!       accept or reject repository position based on random number and fitness
        call random_number(r)
        r2 = r
        !
        if (r2 < fitness(k2)) then
          !
          gindex = k2
          exit
          !
        end if
        !
      end do
      !
    end if
    !
    do iparm=1,npar
      !
      if (trim(partrans(iparm)) == 'none' .or. trim(partrans(iparm)) == 'log' .or. &
          trim(partrans(iparm)) == 'eqlog') then
        !
        if (nrepact == 0) then
          !
          call random_number(r)
          local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r)
          global = c2*r*(pbest(gindex,iparm) - partval(ipart,iparm))
          !
        else
          !
          call random_number(r)
          local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r)
          global = c2*r*(reposit(gindex,iparm) - partval(ipart,iparm))
          !
        end if
        !
        partvel(ipart,iparm) = inertia*partvel(ipart,iparm) + local + global 
        !
        if (dabs(partvel(ipart,iparm)) > vmax(iparm)) then
          if (partvel(ipart,iparm) > 0.0d+00) partvel(ipart,iparm) = vmax(iparm)
          if (partvel(ipart,iparm) < 0.0d+00) partvel(ipart,iparm) = -vmax(iparm)
        end if
        !
        partval(ipart,iparm) = partval(ipart,iparm) + partvel(ipart,iparm)
        !
!       keep perturbing particles if they violate parameter bounds
        if (partval(ipart,iparm) > parubnd(iparm) .or. &
          partval(ipart,iparm) < parlbnd(iparm)) then
          !
          spin = 0
          !
          do
            ! 
            spin = spin + 1
            !
            if (nrepact == 0) then
              !
              call random_number(r)
              local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
              call random_number(r)
              global = c2*r*(pbest(gindex,iparm) - partval(ipart,iparm))
              !
            else
              !
              call random_number(r)
              local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
              call random_number(r)
              global = c2*r*(reposit(gindex,iparm) - partval(ipart,iparm))
              !
            end if
            !
            partvel(ipart,iparm) = inertia*partvel(ipart,iparm) + local + global
            !
            if (dabs(partvel(ipart,iparm)) > vmax(iparm)) then
              if (partvel(ipart,iparm) > 0.0d+00) partvel(ipart,iparm) = vmax(iparm)
              if (partvel(ipart,iparm) < 0.0d+00) partvel(ipart,iparm) = -vmax(iparm)
            end if
            !
            partval(ipart,iparm) = partval(ipart,iparm) + partvel(ipart,iparm)
            !
            if (partval(ipart,iparm) < parubnd(iparm) .and. &
                partval(ipart,iparm) > parlbnd(iparm)) exit
            !
            if (spin > 1000) then
              !
              write(*,*)"I'm spinning my wheels trying to find a feasible parameter value:"
              write(*,*)'-- stopping execution --'
              write(*,*)'ipart ',ipart,'  iparm ',iparm
              write(*,*)'partval ',partval(ipart,iparm)
              write(*,*)'partvel ',partvel(ipart,iparm)
              stop
              !
            end if
            !
          end do
          !
        end if
        !
      end if
      !
    end do
    !
  end do
  !
! transform parameters back to native values
  do ipart=1,npop
    call transpar(ipart,0)
  end do
  !
  call transbnd(0)
  call transrep(0)
  !
! tie child parameters to their parents
  do ipart=1,npop
    do iparm=1,npar
      if (trim(partrans(iparm)) == 'tied') then
        !
        partval(ipart,iparm) = tiedrat(iparm)*partval(ipart,partied(iparm))
        !
      end if
    end do
  end do          
    
end subroutine pertpareto




subroutine pertunc(iiter,gindex)
!========================================================================================
!==== This subroutine perturbs particles randomly based on PSO algorithm.            ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------
  integer,intent(in)::gindex,iiter
  integer::ipart,iparm,spin,gsav
  
  double precision::local,global,r
!----------------------------------------------------------------------------------------    

! calculate current value for inertia
!   if (iiter < inerti) then
!     inertia = iinert + (finert - iinert)*(dble(iiter)/dble(inerti))
!   else
!     inertia = finert
!   end if
  !
  call transbnd(1)
  !
  do ipart=1,npop
    call transpar(ipart,1)
  end do
  !
! perturb particles in their transformed rhealm
  do ipart=1,npop
    !
    do iparm=1,npar
      !
      if (trim(partrans(iparm)) == 'none' .or. trim(partrans(iparm)) == 'log' .or. &
          trim(partrans(iparm)) == 'eqlog') then
        !
        if (neibr == 0) then
          !
          call random_number(r)
          local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r)
          global = c2*r*(pbest(gindex,iparm) - partval(ipart,iparm))
          !
        else
          !
          call random_number(r)
          local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r)
          global = c2*r*(pbest(gneibr(ipart),iparm) - partval(ipart,iparm))
          !
        end if
        !
        !
!       handle precision issues for when particle stops moving
        if (dabs(local) < 1.00d-32) then
          local = 0.0d+00
        end if
        if (dabs(global) < 1.00d-32) then
          global = 0.0d+00
        end if
        if (dabs(partvel(ipart,iparm)) < 1.0d-32) then
          partvel(ipart,iparm) = 0.0d+00
        end if
        !
        !
!       reduce inertia if particle is close to calibration
        if (objopt(ipart) <= inerthr) then
          inertia = finert
        else
          inertia = iinert
        end if
        !
        !
!       calculate particle velocity
        partvel(ipart,iparm) = inertia*partvel(ipart,iparm) + local + global 
        !
        if (dabs(partvel(ipart,iparm)) > vmax(iparm)) then
          if (partvel(ipart,iparm) > 0.0d+00) partvel(ipart,iparm) = vmax(iparm)
          if (partvel(ipart,iparm) < 0.0d+00) partvel(ipart,iparm) = -vmax(iparm)
        end if
        !
        !
!       initial particle perturbation
        partval(ipart,iparm) = partval(ipart,iparm) + partvel(ipart,iparm)
        !
        !
!       keep perturbing particles if they violate parameter bounds (for basic parunc - not nsmc)
        if (sdim == 0) then
          !
          if (partval(ipart,iparm) > parubnd(iparm) .or. &
            partval(ipart,iparm) < parlbnd(iparm)) then
            !
            spin = 0
            !
            do
              ! 
              spin = spin + 1
              !
              if (neibr == 0) then
                !
                call random_number(r)
                local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
                call random_number(r)
                global = c2*r*(pbest(gindex,iparm) - partval(ipart,iparm))
                !
              else
                !
                call random_number(r)
                local  = c1*r*(pbest(ipart,iparm) - partval(ipart,iparm))
                call random_number(r)
                global = c2*r*(pbest(gneibr(ipart),iparm) - partval(ipart,iparm))
                !
              end if
              !
              partvel(ipart,iparm) = inertia*partvel(ipart,iparm) + local + global
              !
              if (dabs(partvel(ipart,iparm)) > vmax(iparm)) then
                if (partvel(ipart,iparm) > 0.0d+00) partvel(ipart,iparm) = vmax(iparm)
                if (partvel(ipart,iparm) < 0.0d+00) partvel(ipart,iparm) = -vmax(iparm)
              end if
              !
              partval(ipart,iparm) = partval(ipart,iparm) + partvel(ipart,iparm)
              !
              if (partval(ipart,iparm) < parubnd(iparm) .and. &
                  partval(ipart,iparm) > parlbnd(iparm)) exit
                !
              if (spin > 1000) then
                !
                write(*,*)"I'm spinning my wheels trying to find a feasible parameter value:"
                write(*,*)'-- stopping execution --'
                write(*,*)'ipart ',ipart,'  iparm ',iparm
                write(*,*)'partval ',partval(ipart,iparm)
                write(*,*)'partvel ',partvel(ipart,iparm)
                stop
                !
              end if
              !
            end do
            !
          end if
          !
        end if
        !
        !
      end if
      !
    end do
    !
    !
    if (sdim > 0) then
      !
!     project perturbation onto null space
      call nullproj(ipart)
      !
    end if
    !
    !
  end do
  !
  !
! transform parameters back to native values
  do ipart=1,npop
    call transpar(ipart,0)
  end do
  !
  call transbnd(0)
  !
! tie child parameters to their parents
  do ipart=1,npop
    do iparm=1,npar
      if (trim(partrans(iparm)) == 'tied') then
        partval(ipart,iparm) = tiedrat(iparm)*partval(ipart,partied(iparm))
      end if
    end do
  end do          
    
end subroutine pertunc







































