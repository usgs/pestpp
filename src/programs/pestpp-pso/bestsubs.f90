subroutine getgbest(gindex,gbest,gmbest,gpbest)
!========================================================================================
!==== This subroutine determines the global best particle within the neighborhood of ====
!==== each particle.
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(out)::gindex
  integer::i,j,ipart,left,right,ineibr,iineibr
  
  double precision,intent(out)::gbest,gmbest,gpbest
  double precision::gnbest
!----------------------------------------------------------------------------------------

  if (neibr == 1) then
    !
!   find gbest within each neighborhood
    gnbest = 0.0d+00
    !
    if (mod(nneibr,2) > 0) then
      left  = (nneibr - 1)/2
      right = left
    else
      left  = nneibr/2
      right = nneibr/2 - 1
    end if
    !
    do ipart=1,npop
      !
      i = 0
      j = 0
      !
      do ineibr=1,nneibr
        !
        iineibr = ipart - (left - ineibr + 1)
        !
!       wrap index around if range goes beyond maximum or minimum index
        if (iineibr < 1)    iineibr = ipart +   (ineibr - 1)    + (npop - left)
        if (iineibr > npop) iineibr = ipart + (ineibr - nneibr) + (right - npop)
        !
        i = i + 1
        !
        if (objopt(iineibr) < gnbest .or. i == 1) then
          !
          gnbest = objopt(iineibr)
          gneibr(ipart) = iineibr
          !
        end if  
        !
      end do
      !
    end do
    !
  end if
  !
! find gbest of the swarm
  !
  i = 0
  !
  do ipart=1,npop
    !
    i = i + 1
    !
    if (objopt(ipart) < gbest .or. i == 1) then
      gbest  = objopt(ipart)
      gmbest = objmopt(ipart)
      gpbest = objpopt(ipart)
      gindex = ipart
    end if
    !
  end do
  
  
end subroutine getgbest  



subroutine getpbest()
!========================================================================================
!==== This subroutine determines the personal best particle position.                ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer::ipart,iparm,igp
!----------------------------------------------------------------------------------------

! if running standard PSO, set pbest based on composite objective function
  do ipart=1,npop
    !
    if (obj(ipart) < objopt(ipart)) then
      !
      objopt(ipart)  = obj(ipart)
      objmopt(ipart) = objm(ipart)
      objpopt(ipart) = objp(ipart)
      !
      do igp=1,nobsgp
        objgpopt(ipart,igp) = objgp(ipart,igp)
      end do
      !
      do iparm=1,npar
        pbest(ipart,iparm) = partval(ipart,iparm)
      end do
      !
    end if
    !
  end do


end subroutine getpbest



subroutine pbestpareto(iiter)
!========================================================================================
!==== This subroutine determines the personal best particle positions for either     ====
!==== constrained or unconstrained MOPSO.                                            ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::iiter
  integer::ipart,iparm,ipto
!----------------------------------------------------------------------------------------
  
  if (nptocon > 0) then
    !
    do ipart=1,npop
      !
      if (feas(ipart) == 0) then
        !
!       particle has never been feasible, set pbest based on violate
        if (violate(ipart) < vioopt(ipart)) then
          !
          vioopt(ipart) = violate(ipart)
          !
          do iparm=1,npar
            pbest(ipart,iparm) = partval(ipart,iparm)
          end do
          !
          do ipto=1,nptogp
            ptogpopt(ipart,ipto) = ptogp(ipart,ipto)
          end do
          !
        end if
        !
      else
        !
!       particle has been feasible at some stage
        if (violate(ipart) > 0.0d+00) then
          !
!         particle was feasible but is not feasible now, do not set pbest
          cycle
          !
        else
          !
          if (vioopt(ipart) > 0.0d+00) then
            !
!           particle is feasible now for the first time, update pbest and zero vioopt
            vioopt(ipart) = 0.0d+00
            !
            do iparm=1,npar
              pbest(ipart,iparm) = partval(ipart,iparm)
            end do
            !
            do ipto=1,nptogp
              ptogpopt(ipart,ipto) = ptogp(ipart,ipto)
            end do
            !
          else
            !
!           particle is feasible now and has been before, update pbest based on MOPSO procedure
            call pdominance(ipart,iiter)
            !
          end if
          !
        end if
        !
      end if
      !
    end do
    !
  else
    !
    do ipart=1,npop
      call pdominance(ipart,iiter)
    end do
    !
  end if

end subroutine pbestpareto



subroutine pdominance(ipart,iiter)
!========================================================================================
!==== This subroutine assigns the personal best position to a particle in the swarm  ====
!==== for basic MOPSO.                                                               ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------
  integer,intent(in)::iiter,ipart
  integer::iparm,ipto,dom,irep,repposs,repposp,coin
  
  double precision::r1,ptoubvr,ptoubvs,dist,mindisp,mindiss,repfits,repfitp
  double precision,dimension(:),allocatable::minrobj,maxrobj
!----------------------------------------------------------------------------------------

! check if pbest is dominated by current particle position
  dom = 0
  !
  do ipto=1,nptogp
    !
    if (ptogpopt(ipart,ipto) > ptogp(ipart,ipto)) then
      !
      dom = 1
      !
    else
      !
      dom = 0
      exit
      !
    end if
    !
  end do
  !
  if (dom == 1) then
    !
!   we know current particle position is dominant --> update pbest
    do iparm=1,npar
      pbest(ipart,iparm) = partval(ipart,iparm)
    end do
    !
    do ipto=1,nptogp
      ptogpopt(ipart,ipto) = ptogp(ipart,ipto)
    end do
    !
  end if
  !
  if (dom == 0) then
    !
!   check if pbest is dominant over current particle position
    dom = 0
    !
    do ipto=1,nptogp
      !
      if (ptogp(ipart,ipto) > ptogpopt(ipart,ipto)) then
        !
        dom = 1
        !
      else
        !
        dom = 0
        exit
        !
      end if
      !
    end do
    !
    if (dom == 1) then
      !
!     we know that pbest is dominant --> do nothing
      return
      !
    else
      !
!     we know that neither particle dominates, choose one carefully based on PTOUB
      ptoubvs = 0.0d+00
      ptoubvr = 0.0d+00
      !
      do ipto=1,nptogp
        !
        if (ptogp(ipart,ipto) > ptoub(ipto)) then
          ptoubvs = ptoubvs + ((ptogp(ipart,ipto) - ptoub(ipto))/ptoub(ipto))**2
        end if
        !
        if (ptogpopt(ipart,ipto) > ptoub(ipto)) then
          ptoubvr = ptoubvr + ((ptogpopt(ipart,ipto) - ptoub(ipto))/ptoub(ipto))**2
        end if
        !
      end do
      !
      ptoubvs = sqrt(ptoubvs)
      ptoubvr = sqrt(ptoubvr)
      !
      if (ptoubvs > ptoubvr) then
        !
!       swarm position is closer to ptogp upper bound - update pbest
        do iparm=1,npar
          pbest(ipart,iparm) = partval(ipart,iparm)
        end do
        !
        do ipto=1,nptogp
          ptogpopt(ipart,ipto) = ptogp(ipart,ipto)
        end do
        !
      end if
      !
      if (ptoubvr - ptoubvs < 1.0d-16) then
        !
!       either both positions are within ptogp upper bounds, or current position
!       is nearly as good as the pbest position - choose one randomly 
        call random_number(r1)
        !
        if (r1 < 5.0d-01) then
          !
          do iparm=1,npar
            pbest(ipart,iparm) = partval(ipart,iparm)
          end do
          !
          do ipto=1,nptogp
            ptogpopt(ipart,ipto) = ptogp(ipart,ipto)
          end do
          !
        end if
!         !
!         !
! !       either both positions are within ptogp upper bounds, or current position
! !       is nearly as good as the pbest position - choose one based on repository fitness 
!         allocate(minrobj(nptogp),maxrobj(nptogp))
!         repfits = 0.0d+00
!         repfitp = 0.0d+00
!         maxrobj = -1.0d+30
!         minrobj = 1.0d+30
!         mindisp = 1.0d+30
!         mindiss = 1.0d+30
!         !
! !       find the max and min of the objs in the repository for normalizing distances
!         do irep=1,npop+nrep
!           !
!           if (repindx(irep) == 1) then
!             !
!             do ipto=1,nptogp
!               if (repoobj(irep,ipto) < minrobj(ipto)) minrobj(ipto) = repoobj(irep,ipto)
!               if (repoobj(irep,ipto) > maxrobj(ipto)) maxrobj(ipto) = repoobj(irep,ipto)
!             end do
!             !
!           end if
!           !
!         end do
!         !
! !       find the fitness associated with the repository position closest to both the 
! !       pbest (repfitp) and current swarm (repfits) position
!         do irep=1,npop+nrep
!           !
!           if (repindx(irep) == 1) then
!             !
! !           current position
!             dist = 0.0d+00
!             do ipto=1,nptogp
!               dist = dist + ((ptogp(ipart,ipto) - repoobj(irep,ipto))/&
!                   (maxrobj(ipto) - minrobj(ipto)))**2
!             end do
!             dist = sqrt(dist)
!             if (iiter == 4 .and. ipart == 1) write(*,*)dist,mindiss
!             !
!             if (dist < mindiss) then
!               mindiss = dist
!               repfits = fitness(irep)
!               repposs = irep
!             end if
!             !
! !           pbest position
!             dist = 0.0d+00
!             do ipto=1,nptogp
!               dist = dist + ((ptogpopt(ipart,ipto) - repoobj(irep,ipto))/&
!                   (maxrobj(ipto) - minrobj(ipto)))**2
!             end do
!             dist = sqrt(dist)
!             if (iiter == 4 .and. ipart == 1) write(*,*)dist,mindisp
!             !
!             if (dist < mindisp) then
!               mindisp = dist
!               repfitp = fitness(irep)
!               repposp = irep
!             end if
!             !
!           end if
!           !
!         end do
!         !
! !       update pbest if the closest repository position to the current swarm position has
! !       a greater fitness than that of the pbest position
! !         if (repposs == repposp) then
! !           !
! ! !         both pbest and current position are closest to the same repository position - 
! ! !         choose one randomly
! !           coin = int(1 + 2*rand(0))
! !           !
! !           if (coin == 1) then
! !             !
! !             do iparm=1,npar
! !               pbest(ipart,iparm) = partval(ipart,iparm)
! !             end do
! !             !
! !             do ipto=1,nptogp
! !               ptogpopt(ipart,ipto) = ptogp(ipart,ipto)
! !             end do
! !             !
! !           end if
! !           !
! !           return
! !           !
! !         end if
!         !
!         if (repfits > repfitp) then
!           !
!           do iparm=1,npar
!             pbest(ipart,iparm) = partval(ipart,iparm)
!           end do
!           !
!           do ipto=1,nptogp
!             ptogpopt(ipart,ipto) = ptogp(ipart,ipto)
!           end do
!           !
!         end if
!         !
!         deallocate(minrobj,maxrobj)
        !
      end if
      !
    end if
    !
  end if
  

end subroutine pdominance













