subroutine repcont(iiter,repred)
!========================================================================================
!==== This subroutine manages dominance in the pareto repository.                    ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::iiter
  integer,intent(inout)::repred
  integer::ipart,ipart1,ipart2,iparm,ipto,dom,irep,irep1,irep2,repcnt,iobs,nptoubv
  integer,dimension(:),allocatable::nondom,ptoubv
!----------------------------------------------------------------------------------------

  allocate(nondom(npop),ptoubv(npop))
  !  
  repred = 0
  repcnt = 0
  ptoubv = 0  
  !
! check repository and make sure it is not too big
  do irep=1,nrep+npop
    if (repindx(irep) == 1) repcnt = repcnt + 1
  end do
  !
  if (repcnt > nrep) then
    !
    write(*,*)'Repository is too big - there is a bug somewhere'
    write(*,*)'--stopping execution--'
    stop
    !
  end if
  !
! reset the dominance recorder for the current swarm
  do ipart=1,npop
    nondom(ipart) = 1
  end do
  !
! check if particle violates limits on ptogp's
  nptoubv = 0
  !
  do ipart=1,npop
    do ipto=1,nptogp
      !
      if (ptogp(ipart,ipto) > ptoub(ipto)) then
        ptoubv(ipart) = 1
        nondom(ipart) = 0
      end if
      !
    end do
    !
    if (nondom(ipart) == 0) nptoubv = nptoubv + 1  
    !
  end do
  !
  if (nptoubv == npop .and. nrepact == 0) then
    !
    write(*,*)'All particles are violating upper bounds on objectives'
    write(*,*)'--stopping execution--'
    stop
    !
  end if
  !
! if particles currently violate the constraints, do not allow them to contribute
! -- to the dominance calculation or enter the repository
  if (nptocon > 0) then
    !
!   setting nondom=0 for infeasible particles will prevent them from entering the repository
    do ipart=1,npop
      if (violate(ipart) > 0.0d+00) nondom(ipart) = 0
    end do
    !
!   determine dominance relationships in the current swarm only among feasible particles
    do ipart1=1,npop
      !
      if (violate(ipart1) > 0.0d+00 .or. ptoubv(ipart1) == 1) then
        !
        cycle
        !
      else
        !
        do ipart2=1,npop
          !
          if (violate(ipart2) > 0.0d+00 .or. ptoubv(ipart2) == 1) then
            !
            cycle
            !
          else
            !
            if (ipart1 /= ipart2) then
              !
              dom = 0
              !
              do ipto=1,nptogp
                !
                if (ptogp(ipart1,ipto) >= ptogp(ipart2,ipto)) then
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
              if (dom == 1) nondom(ipart1) = 0
              !
            end if
            !
          end if
          !
        end do
        !
      end if
      !
    end do
    !
  else
    !
!   determine dominance relationships in the entire current swarm
    do ipart1=1,npop
      !
      if (ptoubv(ipart1) == 1) then
        !
        cycle
        !
      else
        !
        do ipart2=1,npop
          !
          if (ipart1 /= ipart2 .and. ptoubv(ipart2) == 0) then
            !
            dom = 0
            !
            do ipto=1,nptogp
              !
              if (ptogp(ipart1,ipto) >= ptogp(ipart2,ipto)) then
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
            if (dom == 1) nondom(ipart1) = 0
            !
          end if
          !
        end do
        !
      end if
      !
    end do
    !
  end if
  !  
! update repository
  if (iiter == 0) then
    !
!   This is the initial condition - initiate the repository based on starting particle positions
    irep = 0
    !
    do ipart=1,npop
      !
      if (nondom(ipart) == 1) then
        !
        irep = irep + 1
        repindx(irep) = 1
        !
        do iparm=1,npar
          reposit(irep,iparm) = partval(ipart,iparm)
        end do
        !
        do ipto=1,nptogp
          repoobj(irep,ipto) = ptogp(ipart,ipto)
        end do
        !
        do ipto=1,nptocon
          repcon(irep,ipto) = ptocon(ipart,ipto)
        end do
        !
        do iobs=1,nobs
          repoobs(irep,iobs) = mobssav(ipart,iobs)
        end do
        !
      end if
      !
    end do
    !
  else
    !
!   Determine if there are particles that should enter the repository
    do ipart=1,npop
      !
      if (nondom(ipart) == 1) then
        !
        do irep=1,nrep
          !
          if (repindx(irep) == 1) then
            !
            dom = 0
            !
            do ipto=1,nptogp
              !
              if (ptogp(ipart,ipto) >= repoobj(irep,ipto)) then
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
              nondom(ipart) = 0
              exit
              !
            end if
            !
          end if
          !
        end do
        !
      end if
      !
    end do
    !
!   add particles to the repository as needed
    do ipart=1,npop
      !
      if (nondom(ipart) == 1) then
        !
        do irep=1,nrep+npop
          !
          if (repindx(irep) == 0) then
            !
            do iparm=1,npar
              reposit(irep,iparm) = partval(ipart,iparm)
            end do
            !
            do ipto=1,nptogp
              repoobj(irep,ipto) = ptogp(ipart,ipto)
            end do
            !
            do ipto=1,nptocon
              repcon(irep,ipto) = ptocon(ipart,ipto)
            end do
            !
            do iobs=1,nobs
              repoobs(irep,iobs) = mobssav(ipart,iobs)
            end do
            !
            repindx(irep) = 1
            !
            exit
            !
          end if
          !
          if (irep == nrep+npop) then
            !
            write(*,*)'There is an issue with the repository management routine'
            write(*,*)'--stopping execution--'
            stop
            !
          end if
          !
        end do
        !
      end if
      !
    end do    
    !
!   check repository for dominance relationships and remove particles dominated by incoming 
!   -- particles from the swarm
    do irep1=1,nrep+npop
      !
      if (repindx(irep1) == 1) then
        !
        do irep2=1,nrep+npop
          !
          if (repindx(irep2) == 1 .and. irep1 /= irep2) then
            !
            dom = 0
            !
            do ipto=1,nptogp
              !
              if (repoobj(irep1,ipto) >= repoobj(irep2,ipto)) then
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
              repindx(irep1) = 0
              exit
              !
            end if
            !
          end if
          !
        end do
        !
      end if
      !
    end do
    !
  end if
  !
! check how many particles are in the repository
  nrepact  = 0
  do irep=1,nrep+npop
    if (repindx(irep) == 1) nrepact = nrepact + 1
  end do
  !
! check if repository is now too big
  if (nrepact > nrep) repred = 1
  !
  deallocate(nondom,ptoubv)
  
end subroutine repcont



subroutine repgrid(repred)
!========================================================================================
!==== This subroutine maintains the size of the pareto repository using a grid-based ====
!==== based on Coello et al, 2004                                                    ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(inout)::repred
  integer::i,irep,irep1,irep2,ipto,igrd,hyp,nhypmax,nhypmin,k1,repcnt,rmcnt,icnt,imin
  integer,dimension(:),allocatable::reprm,nhyp,repmin
  
  double precision::lgrd,rgrd,r1,normhyp
  double precision,dimension(:),allocatable::objmax,objmin,delobj
!----------------------------------------------------------------------------------------

  allocate(reprm(nrep+npop),nhyp(nrep+npop))
  allocate(objmax(nptogp),objmin(nptogp),delobj(nptogp))
  allocate(repmin(nptogp))
  objmax = 0.0d+00
  objmin = 0.0d+00
  grd    = 0
  nhyp   = 0
  repmin = 0

! get particles associated with max and min objective
  do ipto=1,nptogp
    !
    i = 0
    !
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1) then
        !
        i = i + 1
        !
        if (repoobj(irep,ipto) > objmax(ipto) .or. i == 1) then
          !
          objmax(ipto) = repoobj(irep,ipto)
          !
        end if
        !
        if (repoobj(irep,ipto) < objmin(ipto) .or. i == 1) then
          !
          objmin(ipto) = repoobj(irep,ipto)
          repmin(ipto) = irep
          !
        end if
        !
      end if
      !
    end do
    !
  end do
  
! set up adaptive grid
  do ipto=1,nptogp
    !
    delobj(ipto) = (objmax(ipto) - objmin(ipto))/dble(ngrid-1)
    !
  end do
  
! assign hypercube to each particle
  do irep=1,nrep+npop
    !
    do ipto=1,nptogp
      !
      lgrd = objmin(ipto) - 5.0d-01*delobj(ipto)
!      if (lgrd < 0.0d+00) lgrd = 0.0d+00
      rgrd = objmin(ipto) + 5.0d-01*delobj(ipto)
!      if (rgrd < 0.0d+00) rgrd = 0.0d+00
      !
      do igrd=1,ngrid
        !
        if (repoobj(irep,ipto) > lgrd .and. repoobj(irep,ipto) <= rgrd) then
          !
          grd(irep,ipto) = igrd
          !
        end if
        !
        lgrd = rgrd
        rgrd = rgrd + delobj(ipto)
        !
      end do
      !
    end do
    !
  end do
  
  
  if (repred == 1) then
    !
!   determine all of the particles that lie in the hypercube(s) with the maximum 
!   number of particles in it(them)
    nhypmax = 0
    nhyp    = 0
    reprm   = 0
    icnt    = 0
    !
    do irep1=1,nrep+npop
      !
      imin = 0
      !
      if (repindx(irep1) == 1) then
        !
!       check if particle is a minimum for any objective
!       if so, do not allow it to be susceptible to removal
        do ipto=1,nptogp
          !
          if (irep1 == repmin(ipto)) then
            imin = 1
            exit
          end if
          !
        end do
        !
        if (imin == 1) cycle
        !
        icnt = icnt + 1
        !
        nhyp(irep1) = 1
        !
        do irep2=1,nrep+npop
          !
          hyp = 0
          !
          if (repindx(irep2) == 1 .and. irep1 /= irep2) then
            !
!           check if particle 1 is in same hypercube as particle 2
            do ipto=1,nptogp
              !
              if (grd(irep1,ipto) == grd(irep2,ipto)) then
                hyp = 1
              else
                hyp = 0
                exit
              end if
              !
            end do
            !
          end if
          !
          nhyp(irep1) = nhyp(irep1) + hyp
          !
        end do
        !
!       see if new hypercube is discovered with the most number of particles in it
!       thus far
        if (nhyp(irep1) > nhypmax .or. icnt == 1) then
          !
          do irep=1,nrep+npop
            reprm(irep) = 0
          end do
          !
          nhypmax = nhyp(irep1)
          !
          reprm(irep1) = 1
          !
        end if
        !
!       if the current particle location is tied with the maximum, record it
        !
        if (nhyp(irep1) == nhypmax) then
          !
          reprm(irep1) = 1
          !
        end if
        !
      end if
      !
    end do
    ! 
!   remove selected particles at random, based on density, until repository size is nrep
    icnt = 0
    do
      !
      icnt = icnt + 1
      !
      repcnt = 0
      rmcnt  = 0
      !
!     count the number of particles in the repository
!     also count (among those particles in the repository) the ones susceptible to removal
      do irep=1,nrep+npop
        !
        if (repindx(irep) == 1) then
          !
          repcnt = repcnt + 1
          !
          if (reprm(irep) == 1) rmcnt = rmcnt + 1
          !
        end if
        !
      end do
      !
      if (rmcnt == 0) then
        !
!       if there are no more particles set for removal, make all particles susceptible to removal.
!       this can happen if the hypercube with the maximum number of particles contains fewer particles
!       --  than need to be removed in order to get repcnt down to nrep; or, if there are the same number
!       --  of particles in all hypercubes containing particles, which most likely happens when there 
!       --  are a small number of hypercubes containing particles and each hypercube only contains a 
!       --  single particle.
        !
        do irep=1,nrep+npop
          reprm(irep) = 1
        end do
        !
!       set rmcnt equal to the number of particles in the repository minus the minimum objective particles
        rmcnt = repcnt - nptogp
        !
!       prevent minimum objective particles from being susceptible to removal
        do ipto=1,nptogp
          reprm(repmin(ipto)) = 0
        end do
        !
      end if
      !
      if (repcnt > nrep) then 
        !
!       there are too many particles in the repository still, so let's remove one      
!       select a particle susceptible to removal randomly
        call random_number(r1)
        k1 = int(r1*dble(rmcnt)) + 1
        !
!       count our way up until we arrive at the randomly selected particle
        rmcnt = 0
        !
        do irep=1,nrep+npop
          !
          if (reprm(irep) == 1 .and. repindx(irep) == 1) rmcnt = rmcnt + 1
          !
          if (rmcnt == k1) then
            repindx(irep) = 0
            exit
          end if
          !
        end do
        !
      else if (repcnt == nrep) then
        !
!       the repository is down to nrep, so exit
        exit
        !
      else
        !
        write(*,*)'something is wrong here'
        !
      end if
      !
      if (icnt > npop) then
        write(*,*)'something is wrong'
        write(*,*)k1,repcnt
        stop
      end if
      !
    end do
    !
    repred = 0
    !
  end if
  !
! run a final check on the repository size
  repcnt = 0
  do irep=1,nrep+npop
    if (repindx(irep) == 1) repcnt = repcnt + 1
  end do
  if (repcnt > nrep) then
    write(*,*)'something is wrong here'
    write(*,*)'repcnt = ',repcnt
    stop
  end if
  !
! recalculate grid-based fitness
! re-evaluate the number of particles sharing a hypercube with each particle
  nhyp    = 0
  nrepact = 0
  !
  do irep1=1,nrep+npop
    !
    if (repindx(irep1) == 1) then
      !
      nhyp(irep1) = 1
      nrepact = nrepact + 1
      !
      do irep2=1,nrep+npop
        !
        hyp = 0
        !
        if (repindx(irep2) == 1 .and. irep1 /= irep2) then
          !
!         check if particle 1 is in same hypercube as particle 2
          do ipto=1,nptogp
            !
            if (grd(irep1,ipto) == grd(irep2,ipto)) then
              hyp = 1
            else
              hyp = 0
              exit
            end if
            !
          end do
          !
        end if
        !
        nhyp(irep1) = nhyp(irep1) + hyp
        !
      end do
      !
    end if
    !
  end do
  !
! develop fitness array for particle movement
  fitness = 0.0d+00
  normhyp = 0.0d+00
  !
! find the number of particles in the hypercube with the least amount of particles in it
  icnt = 0
  !
  do irep=1,nrep+npop
    !
    if (repindx(irep) == 1) then
      !
      icnt = icnt + 1
      !
      if (nhyp(irep) < nhypmin .or. icnt == 1) then
        nhypmin = nhyp(irep)
      end if
      !
    end if
    !
  end do
  !
! calculate grid-based fitness
  do irep=1,nrep+npop
    !
    if (repindx(irep) == 1) then
      grdfit(irep) = (dble(nhypmin)/dble(nhyp(irep)))**rfit
      fitness(irep) = grdfit(irep)
    end if
    !
  end do
  !
  deallocate(objmax,objmin,delobj,reprm,nhyp,repmin)

end subroutine repgrid




subroutine loneliness(repred,alpha)
!========================================================================================
!==== This subroutine adjusts the fitness of the grid-based approach based on        ====
!==== clustering.                                                                    ==== 
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!---------------------------------------------------------------------------------------- 
  integer,intent(inout)::repred
  integer::irep,irep1,irep2,ipto,ipto1,ipto2,icnt,reprm
  integer,dimension(:),allocatable::objmin
  integer,dimension(:,:),allocatable::mindom
  
  double precision,intent(out)::alpha
  double precision::sum,maxloner,minloner,robjmin,rise,run,perfull,a
  double precision,dimension(:),allocatable::edist,minobj,maxobj
  double precision,dimension(2)::slope
!----------------------------------------------------------------------------------------

  allocate(edist(nrep+npop))
  allocate(mindom(nrep+npop,nptogp))
  allocate(objmin(nrep+npop))
  allocate(minobj(nptogp),maxobj(nptogp))
!   write(*,*)noweak,nptogp,'poop'
  !
  edist    = 0.0d+00
  lonely   = 0.0d+00
  fitness  = 0.0d+00
  minobj   = 0.0d+00
  maxobj   = 0.0d+00
  robjmin  = 0.0d+00
  maxloner = 0.0d+00
  minloner = 0.0d+00
  reprm    = 0
  objmin   = 0
  mindom   = 0
  !
  if (nrepact > 2) then
    !
!   calculate rfit
    if (nrepact < nrep .and. nrepact > 3) then
      !
      perfull = (dble(nrepact) - 3.00d+00)/(dble(nrep) - 3.0d+00)
      a = 1.00d+00/(exp(rramp) - 1.0d+00)
      alpha = 1.00d+00 + (a*(exp(rramp*perfull) - 1.00d+00))*(rfit - 1.00d+00)
      !
    else if (nrepact == 3) then
      !
      alpha = 1.00d+00
      !
    else if (nrepact >= nrep) then
      !
      alpha = rfit
      !
    end if
    !
!   for each position in the repository, find the minimum dominated position for each objective
    do irep1=1,nrep+npop
      !
      if (repindx(irep1) == 1) then
        !
        do ipto=1,nptogp
          !
          icnt = 0
          !
          do irep2=1,nrep+npop
            !
            if (repindx(irep2) == 1 .and. irep1 /= irep2) then              
              !
              if (repoobj(irep1,ipto) < repoobj(irep2,ipto)) then
                !
                icnt = icnt + 1
                !
                if (repoobj(irep2,ipto) < robjmin .or. icnt == 1) then
                  !
                  robjmin = repoobj(irep2,ipto)
                  mindom(irep1,ipto) = irep2
                  !
                end if
                !
              end if
              !
            end if
            !
          end do
          !
!         if icnt = 0, we are the tail of the front, i.e., at the minimum value of the current objective
!         since the tails (or extremes) are important, the fitness assigned to them will be 1.00
          if (icnt == 0) then
            objmin(irep1) = 1
          end if
          !
        end do
        !
      end if
      !
    end do
    !
!   find extreme values on the current pareto front
    do ipto=1,nptogp
      !
      icnt = 0
      !
      do irep=1,nrep+npop
        !
        if (repindx(irep) == 1) then
          !
          icnt = icnt + 1
          !
          if (repoobj(irep,ipto) < minobj(ipto) .or. icnt == 1) then
            minobj(ipto) = repoobj(irep,ipto)
          end if
          !
          if (repoobj(irep,ipto) > maxobj(ipto) .or. icnt == 1) then
            maxobj(ipto) = repoobj(irep,ipto)
          end if
          !
        end if
        !
      end do
      !
    end do
    !
!   calculate euclidean distance among each repository position's nearest neighbors
    icnt = 0
    !
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1 .and. objmin(irep) == 0) then
        !
        do ipto1=1,nptogp
          !
          sum = 0.0d+00
          !
          if (lonorm == 1) then
            !
            do ipto2=1,nptogp
              sum = sum + ((repoobj(irep,ipto2) - repoobj(mindom(irep,ipto1),ipto2))/&
                  (maxobj(ipto2) - minobj(ipto2)))**2 
            end do
            !
          else
            !
            do ipto2=1,nptogp
              sum = sum + (repoobj(irep,ipto2) - repoobj(mindom(irep,ipto1),ipto2))**2 
            end do
            !
          end if
          !
          edist(irep) = edist(irep) + sqrt(sum)
          !
        end do
        !
!       find loneliest repository position
        icnt = icnt + 1
        !
        if (edist(irep) > maxloner .or. icnt == 1) then
          maxloner = edist(irep)
        end if
        !
      end if
      !
    end do
    !
!   calculate fitness based on loneliness
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1) then
        !
        if (objmin(irep) == 1) then
          !
          if (repred == 1) then
            lonely(irep) = 1.0d+00
          else
            lonely(irep) = 1.0d+00
          end if
          !
          fitness(irep) = lonely(irep)
          !
        else
          !
          lonely(irep) = (edist(irep)/maxloner)**alpha
          fitness(irep) = lonely(irep)
          !
        end if
        !
      end if
      !
    end do
    !
  else if (nrepact == 2) then
    !
    do irep=1,nrep+npop
      if (repindx(irep) == 1) fitness(irep) = 1.00d+00
    end do
    !
  end if
  !
! remove repository positions until repository size is nrep
  if (repred == 1) then
    !
    do 
      !
!     find particle with minimum loneliness
      icnt = 0
      !
      do irep=1,nrep+npop
        !
        if (repindx(irep) == 1 .and. objmin(irep) == 0) then
          !
          icnt = icnt + 1
          !
          if (lonely(irep) < minloner .or. icnt == 1) then
            !
            minloner = lonely(irep)
            reprm = irep
            !
          end if
          !
        end if
        !
      end do
      !
      repindx(reprm) = 0
      !
      icnt = 0
      !
      do irep=1,nrep+npop
        if (repindx(irep) == 1) icnt = icnt + 1
      end do
      !
      if (icnt == nrep) exit
      !
    end do
    !
!   recalculate fitness based on reduced repository
!   for each position in the repository, find the minimum dominated position for each objective
    edist    = 0.0d+00
    lonely   = 0.0d+00
    robjmin  = 0.0d+00
    maxloner = 0.0d+00
    mindom   = 0
    !
    do irep1=1,nrep+npop
      !
      if (repindx(irep1) == 1) then
        !
        do ipto=1,nptogp
          !
          icnt = 0
          !
          do irep2=1,nrep+npop
            !
            if (repindx(irep2) == 1 .and. irep1 /= irep2) then              
              !
              if (repoobj(irep1,ipto) < repoobj(irep2,ipto)) then
                !
                icnt = icnt + 1
                !
                if (repoobj(irep2,ipto) < robjmin .or. icnt == 1) then
                  !
                  robjmin = repoobj(irep2,ipto)
                  mindom(irep1,ipto) = irep2
                  !
                end if
                !
              end if
              !
            end if
            !
          end do
          !
        end do
        !
      end if
      !
    end do
    !
!   calculate euclidean distance among each repository position's nearest neighbors
    icnt = 0
    !
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1 .and. objmin(irep) == 0) then
        !
        do ipto1=1,nptogp
          !
          sum = 0.0d+00
          !
          if (lonorm == 1) then
            !
            do ipto2=1,nptogp
              sum = sum + ((repoobj(irep,ipto2) - repoobj(mindom(irep,ipto1),ipto2))/&
                  (maxobj(ipto2) - minobj(ipto2)))**2 
            end do
            !
          else
            !
            do ipto2=1,nptogp
              sum = sum + (repoobj(irep,ipto2) - repoobj(mindom(irep,ipto1),ipto2))**2 
            end do
            !
          end if
          !
          edist(irep) = edist(irep) + sqrt(sum)
          !
        end do
        !
! !       calculate percent slope with neighbors for two-d problems
!         if (noweak == 1 .and. nptogp == 2) then
!           !
!           rise = ((repoobj(irep,2) - repoobj(mindom(irep,1),2))/&
!               repoobj(irep,2))
!           run  = ((repoobj(irep,1) - repoobj(mindom(irep,1),1))/&
!               repoobj(irep,1))
!           slope(1) = rise/run
!           !
!           rise = ((repoobj(irep,1) - repoobj(mindom(irep,2),1))/&
!               repoobj(irep,1))
!           run  = ((repoobj(irep,2) - repoobj(mindom(irep,2),2))/&
!               repoobj(irep,2))
!           slope(2) = rise/run
!           !
!         end if
! !       this is only working if repred = 1 - need to update
!         write(*,*)'slope',slope(1),slope(2)
!         write(*,*)'objs ',repoobj(irep,1),repoobj(irep,2)
        !
!       find loneliest repository position
        icnt = icnt + 1
        !
        if (edist(irep) > maxloner .or. icnt == 1) then
          maxloner = edist(irep)
        end if
        !
      end if
      !
    end do
    !
!   calculate fitness based on loneliness
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1) then
        !
        if (objmin(irep) == 1) then
          !
          lonely(irep) = 1.0d+00
          fitness(irep) = lonely(irep)
          !
        else
          !
          lonely(irep) = (edist(irep)/maxloner)**alpha
          fitness(irep) = lonely(irep)
          !
        end if
        !
      end if
      !
    end do
    !
  end if
  !
! final check of repository size
  nrepact = 0
  do irep=1,nrep+npop
    if (repindx(irep) == 1) nrepact = nrepact + 1
  end do
  if (nrepact > nrep) then
    write(*,*)'Something is wrong with lonliness-based repository management'
    write(*,*)'--stopping execution--'
    stop
  end if
  !
  repred = 0
  !
  deallocate(mindom,minobj,maxobj,objmin,edist)
  
    
end subroutine loneliness




subroutine combgrlo()
!========================================================================================
!==== This subroutine combines the fitness of the grid-based and loneliness approach ====
!==== into one fitness criteria.                                                     ==== 
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!---------------------------------------------------------------------------------------- 
  integer::irep
!----------------------------------------------------------------------------------------

  fitness = 0.0d+00
  !
  do irep=1,nrep+npop
    !
    if (repindx(irep) == 1) then
      !
      fitness(irep) = (grdfit(irep) + lonely(irep))/2.0d+00
      !
    end if
    !
  end do
  
end subroutine combgrlo





