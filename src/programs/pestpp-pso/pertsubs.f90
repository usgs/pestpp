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
  
  double precision::local,global,r1
!----------------------------------------------------------------------------------------    

! calculate current value for inertia
  if (iiter < inerti) then
    inertia = iinert + (finert - iinert)*(dble(iiter)/dble(inerti))
  else
    inertia = finert
  end if
  !  
! transform all parameters for all particles
  do iparm=1,npar
    !
    call transpar(iparm,1)
    !
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
        if (neibr == 1) then
          call random_number(r1)
          local  = c1*r1*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r1)
          global = c2*r1*(pbest(gneibr(ipart),iparm) - partval(ipart,iparm))
        else
          call random_number(r1)
          local  = c1*r1*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r1)
          global = c2*r1*(pbest(gindex,iparm) - partval(ipart,iparm))
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
            if (neibr == 1) then
              call random_number(r1)
              local  = c1*r1*(pbest(ipart,iparm) - partval(ipart,iparm))
              call random_number(r1)
              global = c2*r1*(pbest(gneibr(ipart),iparm) - partval(ipart,iparm))
            else
              call random_number(r1)
              local  = c1*r1*(pbest(ipart,iparm) - partval(ipart,iparm))
              call random_number(r1)
              global = c2*r1*(pbest(gindex,iparm) - partval(ipart,iparm))
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
! transform all parameters for all particles back to original form
  do iparm=1,npar
    !
    call transpar(iparm,0)
    !
  end do
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
  
  double precision::local,global,r1,r2,r3, viomin
!----------------------------------------------------------------------------------------    

! calculate current value for inertia
  if (iiter < inerti) then
    inertia = iinert + (finert - iinert)*(dble(iiter)/dble(inerti))
  else
    inertia = finert
  end if
  !  
! transform all parameters for all particles
  do iparm=1,npar
    !
    call transpar(iparm,1)
    !
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
        call random_number(r1)
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
        call random_number(r2)
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
        if (nrepact == 1) then
          !
          call random_number(r3)
          local  = c1*r3*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r3)
          global = c2*r3*(pbest(gindex,iparm) - partval(ipart,iparm))
          !
        else
          !
          call random_number(r3)  
          local  = c1*r3*(pbest(ipart,iparm) - partval(ipart,iparm))
          call random_number(r3)
          global = c2*r3*(reposit(gindex,iparm) - partval(ipart,iparm))
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
              call random_number(r3)
              local  = c1*r3*(pbest(ipart,iparm) - partval(ipart,iparm))
              call random_number(r3)
              global = c2*r3*(pbest(gindex,iparm) - partval(ipart,iparm))
              !
            else
              !
              call random_number(r3)
              local  = c1*r3*(pbest(ipart,iparm) - partval(ipart,iparm))
              call random_number(r3)
              global = c2*r3*(reposit(gindex,iparm) - partval(ipart,iparm))
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
! transform all parameters for all particles back to original form
  do iparm=1,npar
    !
    call transpar(iparm,0)
    !
  end do
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








































