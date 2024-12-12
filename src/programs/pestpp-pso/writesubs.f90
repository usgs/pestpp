subroutine witmess(iiter)
!========================================================================================
!==== This subroutine writes iteration output to the terminal.                       ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::iiter
!----------------------------------------------------------------------------------------

  write(*,'(A,I0)')'OPTIMISATION ITERATION: ',iiter
  
end subroutine witmess


subroutine listini(basnam,gindex)
!========================================================================================
!==== This subroutine writes list output at the very beginning of the simulation.    ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::gindex
  integer::iparm,ipart,igp
  
  character(len=100),intent(in)::basnam
  character(len=9)::iniparm
!----------------------------------------------------------------------------------------
  
  write(unit(1),'(A)')      '   Particle Swarm Optimization for PEST -- Version 0.0.0'
  write(unit(1),'(/,A,2/)') '                      by Adam Siade'     
  !
  write(unit(1),'(A)')   '                   PEST++ Version 3.0.0'                
  write(unit(1),'(/,A)') '                      by Dave Welter'
  write(unit(1),'(A,2/)')'          Computational Water Resource Engineering'
  !     
  write(unit(1),'(A,A,/)')'PEST-PSO RUN RECORD: CASE ',trim(basnam)
  write(unit(1),'(A,/)')'PEST-PSO run mode:-'
  write(unit(1),'(3X,A,/)')trim(pestmode)
  if (rstpso == 1) write(unit(1),'(A,/)')'PEST-PSO restarted with pest-pso.rst file'
  write(unit(1),'(A,/)')'Case dimensions:-'
  write(unit(1),'(3X,A,T40,I6)')'Number of parameters',npar
  write(unit(1),'(3X,A,T40,I6)')'Number of adjustable parameters',nadj
  write(unit(1),'(3X,A,T40,I6)')'Number of observations',nobs
  write(unit(1),'(3X,A,T40,I6,/)')'Number of prior estimates',nprior
  !
  write(unit(1),'(A,/)')'Parameter definitions:-'
  !
  if (ntied > 0) then
    !
    write(unit(1),'(A)')'Tied parameters were specified as follows,'
    do iparm=1,npar
      if (partied(iparm) > 0) then
        write(unit(1),'(3X,3(A))')trim(parnme(iparm)),' <--> ',trim(parnme(partied(iparm)))
      end if
    end do
    !
    write(unit(1),'(A,/)')'Note that upper and lower bounds are ignored for tied parameters'
    !
  end if
  !
  if (rstpso == 1) then
    write(unit(1),'(A/)')'PEST-PSO restarted with best particle value as follows:'
  end if
  write(unit(1),'(3X,A,T19,A,T34,A,T49,A,T64,A,T79,A,T94,A,T109,A)')&
      'Name','Trans-','Initial','Lower','Upper','Vmax','Log','Range'
  write(unit(1),'(T19,A,T34,A,T49,A,T64,A,T94,A,T94,A)')&
      'formation','value','bound','bound','base'
  !
  if (rstpso == 1) then
    !
    do iparm=1,npar
      if (trim(partrans(iparm)) == 'log') then
        !
        write(unit(1),'(3X,A,T19,A,T29,6(1X,1PE13.6,1X))')parnme(iparm),partrans(iparm),&
          parval1(iparm),parlbnd(iparm),parubnd(iparm),vmax(iparm),1.0d+01,&
          dabs(log10(parubnd(iparm)) - log10(parlbnd(iparm)))
        !
      else if (trim(partrans(iparm)) == 'eqlog') then
        !
        write(unit(1),'(3X,A,T19,A,T29,6(1X,1PE13.6,1X))')parnme(iparm),partrans(iparm),&
          parval1(iparm),parlbnd(iparm),parubnd(iparm),vmax(iparm),base(iparm),maxrange
        !
      else
        !
        write(unit(1),'(3X,A,T19,A,T29,4(1X,1PE13.6,1X),T94,A,T109,A)')parnme(iparm),&
          partrans(iparm),parval1(iparm),parlbnd(iparm),parubnd(iparm),vmax(iparm),&
          '--','--'
        !
      end if
    end do
    !
  else
    !
    do iparm=1,npar
      if (trim(partrans(iparm)) == 'log') then
        !
        write(unit(1),'(3X,A,T19,A,T29,6(1X,1PE13.6,1X))')parnme(iparm),partrans(iparm),&
          parval1(iparm),parlbnd(iparm),parubnd(iparm),vmax(iparm),1.0d+01,&
          dabs(log10(parubnd(iparm)) - log10(parlbnd(iparm)))
        !
      else if (trim(partrans(iparm)) == 'eqlog') then
        !
        write(unit(1),'(3X,A,T19,A,T29,6(1X,1PE13.6,1X))')parnme(iparm),partrans(iparm),&
          parval1(iparm),parlbnd(iparm),parubnd(iparm),vmax(iparm),base(iparm),maxrange
        !
      else
        !
        write(unit(1),'(3X,A,T19,A,T29,4(1X,1PE13.6,1X),T94,A,T109,A)')parnme(iparm),&
          partrans(iparm),parval1(iparm),parlbnd(iparm),parubnd(iparm),vmax(iparm),&
          '--','--'
        !
      end if
    end do
    !
  end if
  !
  write(unit(1),'(/,A,/)')'Control settings:-'
  write(unit(1),'(3X,A,T63,I9)')'Maximum number of optimisation iterations',noptmax    
  write(unit(1),'(3X,A,T63,I9)')'Number of particles in swarm',npop
  write(unit(1),'(3X,A,T63,1PE9.2)')'Value for cognitive parameter, c1',c1
  write(unit(1),'(3X,A,T63,1PE9.2)')'Value for social parameter, c2',c2
  iniparm = 'no'
  if (initp == 1) iniparm = 'yes'
  write(unit(1),'(3X,A,T63,A)')'Initial parameter values used?',adjustr(iniparm)
  write(unit(1),'(3X,A,T63,1PE9.2)')'Initial inertia',iinert
  write(unit(1),'(3X,A,T63,1PE9.2)')'Final inertia',finert
  !
  if (neibr == 1) then
    write(unit(1),'(/,3X,A)')'Neighborhoods are being implemented'
    write(unit(1),'(3X,A,T63,I9)')'Neghborhood size',nneibr
  end if
  !
  if (nprior > 0) then
    write(unit(1),'(/,3X,A)')'Prior information is included in objective function'
  end if
  !
  if (trim(pestmode) == 'parunc') then
    write(unit(1),'(/,3X,A)')'Parameter uncertainty mode is active:-'
    write(unit(1),'(3X,A,T63,I9)')'Number of calibrated particles desired',nreal
    write(unit(1),'(3X,A,T63,1PE9.2)')'Obj threshold for calibration',thresh
  end if
  !
  write(unit(1),'(3/)')
  !
  if (rstpso == 1) then
    !
    write(unit(1),'(A,A)')           'OPTIMISATION ITERATION NO.        : ','restart'
    write(unit(1),'(A,I0)')          '   Model calls so far             : ',0
    write(unit(1),'(A,1PE10.3)')     '   Best particle value            : ',objopt(gindex)
    if (nprior > 0) then
      write(unit(1),'(A,1PE10.3)')   '   Best meas obj value            : ',objmopt(gindex)
      write(unit(1),'(A,1PE10.3)')   '   Best pri inf obj value         : ',objpopt(gindex)
    end if
    write(unit(1),'(A,I0)')          '   Best particle no.              : ',gindex
    do igp=1,nobsgp
      write(unit(1),'(3(A),1PE10.3)')'    Meas/pri obj. ',adjustl(obgnme(igp)),' : ',&
        objgpopt(gindex,igp)
    end do
    write(unit(1),'(/,A)')           '   Particle positions upon restart:'
    !
    if (suppart > 0) then
    !
!     write particle identifiers  
      do ipart=1,npop
        if (ipart == 1) write(unit(1),'(25X,A)',advance='no')' '
        if (ipart < npop) then
          write(unit(1),'(A,I3.3,3X)',advance='no')'particle-',ipart
        else
          write(unit(1),'(A,I3.3,3X)')'particle-',ipart
        end if
      end do
      !
!     write best composite objective function for each particle
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best comp obj   -->'
        !
        if (ipart < npop) then
          write(unit(1),'(5X,1PE10.3)',advance='no')objopt(ipart)
        else
          write(unit(1),'(5X,1PE10.3)')objopt(ipart)
        end if
        !
      end do
      !
!     write current composite objective function for each particle
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'current comp obj-->'
        !
        if (ipart < npop) then
         write(unit(1),'(5X,1PE10.3)',advance='no')obj(ipart)
        else
          write(unit(1),'(5X,1PE10.3)')obj(ipart)
        end if
        !
      end do
      !
      if (nprior > 0) then
        !
!       write best measurement objective function for each particle
        do ipart=1,npop
          !
          if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best meas obj   -->'
          !
          if (ipart < npop) then
            write(unit(1),'(5X,1PE10.3)',advance='no')objmopt(ipart)
          else
            write(unit(1),'(5X,1PE10.3)')objmopt(ipart)
          end if
          !
        end do
        !
!       write current measurement objective function for each particle
        do ipart=1,npop
          !
          if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'current meas obj-->'
          !
          if (ipart < npop) then
            write(unit(1),'(5X,1PE10.3)',advance='no')objm(ipart)
          else
            write(unit(1),'(5X,1PE10.3)')objm(ipart)
          end if
          !
        end do
        !
!       write best prior information objective function for each particle
        do ipart=1,npop
          !
          if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best pri inf obj-->'
          !
          if (ipart < npop) then
            write(unit(1),'(5X,1PE10.3)',advance='no')objpopt(ipart)
          else
            write(unit(1),'(5X,1PE10.3)')objpopt(ipart)
          end if
          !
        end do
        !
!       write current prior information objective function for each particle
        do ipart=1,npop
          !
          if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'current pri obj -->'
          !
          if (ipart < npop) then
            write(unit(1),'(5X,1PE10.3)',advance='no')objp(ipart)
          else
            write(unit(1),'(5X,1PE10.3)')objp(ipart)
          end if
          !
       end do
        ! 
      end if
      !
!     write best neighbor particle identifier if neighborhoods are used
      if (neibr > 0) then
        !
        do ipart=1,npop
          !
          if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best neighbor   -->'
          !
          if (ipart < npop) then
            write(unit(1),'(8X,I4.4,3X)',advance='no')gneibr(ipart)
          else
            write(unit(1),'(8X,I4.4,3X)')gneibr(ipart)
          end if
          !
        end do
        !
      end if
      !
      if (suppart > 1) then
        !
!       write initial parameter values for each particle
        write(unit(1),'()')
        !
        do iparm=1,npar
          do ipart=1,npop
            if (ipart == 1) write(unit(1),'(3X,A19)',advance='no')adjustl(parnme(iparm))
            if (ipart < npop) then
              write(unit(1),'(5X,1PE10.3)',advance='no')scale(iparm)*partval(ipart,iparm)&
                + offset(iparm)
            else
              write(unit(1),'(5X,1PE10.3)')scale(iparm)*partval(ipart,iparm)&
                + offset(iparm)
            end if
          end do
        end do
        !
      end if
      !
    end if
    !
  end if
  !
  write(unit(1),'()')

end subroutine listini



subroutine listipt(basnam)
!========================================================================================
!==== This subroutine writes list output at the very beginning of the simulation.    ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer::iparm
  
  character(len=100),intent(in)::basnam
  character(len=9)::iniparm,minnorm
!----------------------------------------------------------------------------------------

  write(unit(1),'(A,A,/)')'PEST-PSO RUN RECORD: CASE ',trim(basnam)
  write(unit(1),'(A,/)')'PEST-PSO run mode:-'
  write(unit(1),'(3X,A,/)')trim(pestmode)
  if (rstpso == 1) write(unit(1),'(A,/)')'PEST-PSO restarted with ',trim(basnam) // '.rst file'
  write(unit(1),'(A,/)')'Case dimensions:-'
  write(unit(1),'(3X,A,T40,I6)')'Number of parameters',npar
  write(unit(1),'(3X,A,T40,I6)')'Number of adjustable parameters',nadj
  write(unit(1),'(3X,A,T40,I6)')'Number of observations',nobs
  write(unit(1),'(3X,A,T40,I6,/)')'Number of prior estimates',nprior
  !
  write(unit(1),'(A,/)')'Parameter definitions:-'
  !
  if (ntied > 0) then
    !
    write(unit(1),'(A)')'Tied parameters were specified as follows,'
    do iparm=1,npar
      if (partied(iparm) > 0) then
        write(unit(1),'(3X,3(A))')trim(parnme(iparm)),' <--> ',trim(parnme(partied(iparm)))
      end if
    end do
    !
    write(unit(1),'(A,/)')'Note that upper and lower bounds are ignored for tied parameters'
    !
  end if
  !
  if (rstpso == 1) then
    write(unit(1),'(A/)')'PEST-PSO restarted with parameter definitions as follows:'
  end if
  write(unit(1),'(3X,A,T19,A,T34,A,T49,A,T64,A,T79,A,T94,A)')&
      'Name','Trans-','Lower','Upper','Vmax','Log','Range'
  write(unit(1),'(T19,A,T34,A,T49,A,T79,A)')&
      'formation','bound','bound','base'
  !
  if (rstpso == 1) then
    !
    do iparm=1,npar
      if (trim(partrans(iparm)) == 'log') then
        !
        write(unit(1),'(3X,A,T19,A,T29,5(1X,1PE13.6,1X))')parnme(iparm),partrans(iparm),&
          parlbnd(iparm),parubnd(iparm),vmax(iparm),1.0d+01,&
          dabs(log10(parubnd(iparm)) - log10(parlbnd(iparm)))
        !
      else if (trim(partrans(iparm)) == 'eqlog') then
        !
        write(unit(1),'(3X,A,T19,A,T29,5(1X,1PE13.6,1X))')parnme(iparm),partrans(iparm),&
          parlbnd(iparm),parubnd(iparm),vmax(iparm),base(iparm),maxrange
        !
      else
        !
        write(unit(1),'(3X,A,T19,A,T29,3(1X,1PE13.6,1X),T79,A,T94,A)')parnme(iparm),&
          partrans(iparm),parlbnd(iparm),parubnd(iparm),vmax(iparm),&
          '--','--'
        !
      end if
    end do
    !
  else
    !
    do iparm=1,npar
      if (trim(partrans(iparm)) == 'log') then
        !
        write(unit(1),'(3X,A,T19,A,T29,5(1X,1PE13.6,1X))')parnme(iparm),partrans(iparm),&
          parlbnd(iparm),parubnd(iparm),vmax(iparm),1.0d+01,&
          dabs(log10(parubnd(iparm)) - log10(parlbnd(iparm)))
        !
      else if (trim(partrans(iparm)) == 'eqlog') then
        !
        write(unit(1),'(3X,A,T19,A,T29,5(1X,1PE13.6,1X))')parnme(iparm),partrans(iparm),&
          parlbnd(iparm),parubnd(iparm),vmax(iparm),base(iparm),maxrange
        !
      else
        !
        write(unit(1),'(3X,A,T19,A,T29,3(1X,1PE13.6,1X),T79,A,T94,A)')parnme(iparm),&
          partrans(iparm),parlbnd(iparm),parubnd(iparm),vmax(iparm),&
          '--','--'
        !
      end if
    end do
    !
  end if
  !
  write(unit(1),'(/,A,/)')'Control settings:-'
  write(unit(1),'(3X,A,T63,I9)')'Maximum number of optimisation iterations',noptmax    
  write(unit(1),'(3X,A,T63,I9)')'Number of particles in swarm',npop
  write(unit(1),'(3X,A,T63,1PE9.2)')'Value for cognitive parameter, c1',c1
  write(unit(1),'(3X,A,T63,1PE9.2)')'Value for social parameter, c2',c2
  if (rstpso == 0) then
    iniparm = 'no'
    if (initp == 1) iniparm = 'yes'
    write(unit(1),'(3X,A,T63,A)')'Initial parameter values used?',adjustr(iniparm)
  end if
  write(unit(1),'(3X,A,T63,1PE9.2)')'Initial inertia',iinert
  write(unit(1),'(3X,A,T63,1PE9.2)')'Final inertia',finert
  write(unit(1),'(3X,A,T63,I9)')'Maximum repository size',nrep
  write(unit(1),'(3X,A)',advance='no')'Repository management method'
  if (repmode == 1) then
    write(unit(1),'(T31,A)')'      grid'
    write(unit(1),'(3X,A,T63,I9)')'Number of grid divisions',ngrid
  else if (repmode == 2) then
    write(unit(1),'(T31,A)')'loneliness'
    if (lonorm == 1) then
      minnorm = 'yes'
    else
      minnorm = 'no'
    end if
    write(unit(1),'(3X,A,T63,A)')'Min-Normalized objective?',adjustr(minnorm)
  end if
  write(unit(1),'(3X,A,T63,1PE9.2)')'Alpha',rfit
  !
  write(unit(1),'(2/)')

end subroutine listipt



subroutine listout(iiter,gbest,gmbest,gpbest,gindex,objmin,iphistp)
!========================================================================================
!==== This subroutine writes list output at the end of each iteration.               ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::iiter,gindex,iphistp
  integer::ipart,iparm,igp
  
  double precision,intent(in)::gbest,gmbest,gpbest,objmin
  double precision::gsum,gave
!----------------------------------------------------------------------------------------

! calculate basic stats on objectives
  gsum = 0.0d+00
  !
  do ipart=1,npop
    gsum = gsum + objopt(ipart)
  end do
  !
  gave = gsum/dble(npop)
  !
! write listing output
  write(unit(1),'(A,I10)')                      'OPTIMISATION ITERATION NO.        : ',iiter
  write(unit(1),'(A,I10)')                      '   Model calls so far             : ',modeval
  if (iiter > 0) then
    write(unit(1),'(A,1PE10.3)')                '   Current inertia value          : ',inertia  
  end if
  write(unit(1),'(A,1PE10.3)')                  '   Best particle value            : ',gbest
  write(unit(1),'(A,1PE10.3)')                  '   Average among pbest positions  : ',gave
  if (iiter > 0 .and. trim(pestmode) == 'estimation') then
    write(unit(1),'(A,F10.3)')                  '   Percent of previous iter.      : ',&
        (gbest/objmin)*1.0d+02
    write(unit(1),'(A,I10)')                    '   Iters. since last sig. reduct. : ',iphistp
  end if
  if (nprior > 0) then
    write(unit(1),'(A,1PE10.3)')                '   Best measurement obj.          : ',gmbest
    write(unit(1),'(A,1PE10.3)')                '   Best prior information obj.    : ',gpbest
  end if
  write(unit(1),'(A,I10)')                      '   Best particle no.              : ',gindex
  do igp=1,nobsgp
    write(unit(1),'(3(A),1PE10.3)')             '    Meas/pri obj. ',adjustl(obgnme(igp)),' : ',&
        objgpopt(gindex,igp)
  end do
  if (suppart > 0) then
    write(unit(1),'(/,A)')                      '   Particle positions evaluated at current iteration:'
  end if
  !
  if (suppart > 0) then
    !
!   write particle identifiers  
    do ipart=1,npop
      if (ipart == 1) write(unit(1),'(25X,A)',advance='no')' '
      if (ipart < npop) then
        write(unit(1),'(A,I3.3,3X)',advance='no')'particle-',ipart
      else
        write(unit(1),'(A,I3.3,3X)')'particle-',ipart
      end if
    end do
    !
!   write best composite objective function for each particle
    do ipart=1,npop
      !
      if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best comp obj   -->'
      !
      if (ipart < npop) then
        write(unit(1),'(5X,1PE10.3)',advance='no')objopt(ipart)
      else
        write(unit(1),'(5X,1PE10.3)')objopt(ipart)
      end if
      !
    end do
    !
!   write current composite objective function for each particle 
    do ipart=1,npop
      !
      if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'current comp obj-->'
      !
      if (ipart < npop) then
        !
        if (modfail(ipart) == 0) then
          write(unit(1),'(5X,1PE10.3)',advance='no')obj(ipart)
        else
          write(unit(1),'(5X,A)',advance='no')'  failed  '
        end if
        !
      else
        !
        if (modfail(ipart) == 0) then
          write(unit(1),'(5X,1PE10.3)')obj(ipart)
        else
          write(unit(1),'(5X,A)')'  failed  '
        end if
        !
      end if
      !
    end do
    !
    if (nprior > 0) then
      !
!     write best measurement objective function for each particle
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best meas obj   -->'
        !
        if (ipart < npop) then
          write(unit(1),'(5X,1PE10.3)',advance='no')objmopt(ipart)
        else
          write(unit(1),'(5X,1PE10.3)')objmopt(ipart)
        end if
        !
      end do
      !
!     write current measurement objective function for each particle (if not reassigned)  
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'current meas obj-->'
        !
        if (ipart < npop) then
          !
          if (modfail(ipart) == 0) then
            write(unit(1),'(5X,1PE10.3)',advance='no')objm(ipart)
          else
            write(unit(1),'(15X)',advance='no')
          end if
          !
        else
          !
          if (modfail(ipart) == 0) then
            write(unit(1),'(5X,1PE10.3)')objm(ipart)
          else
            write(unit(1),'(15X)')
          end if
          !
        end if
        !
      end do
      !
!     write best prior info objective function for each particle
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best pri inf obj-->'
        !
        if (ipart < npop) then
         write(unit(1),'(5X,1PE10.3)',advance='no')objpopt(ipart)
        else
          write(unit(1),'(5X,1PE10.3)')objpopt(ipart)
        end if
        !
      end do
      !  
!     write current prior info objective function for each particle 
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'current pri obj -->'
        !
        if (ipart < npop) then
          write(unit(1),'(5X,1PE10.3)',advance='no')objp(ipart)  
        else
         write(unit(1),'(5X,1PE10.3)')objp(ipart)
         end if
        !
      end do
      !
    end if
    !
!   write best neighbor particle identifier if neighborhoods are used
    if (neibr == 1 .and. trim(pestmode) == 'estimation') then
      !
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best neighbor   -->'
        !
        if (ipart < npop) then
          write(unit(1),'(8X,I4.4,3X)',advance='no')gneibr(ipart)
        else
          write(unit(1),'(8X,I4.4,3X)')gneibr(ipart)
        end if
        !
      end do 
      !
    end if
    !
!   write extra blank line
    write(unit(1),'(A)')' '
    !
    if (suppart > 1) then
      !
!     write current particle positions
      do iparm=1,npar
        do ipart=1,npop
          if (ipart == 1) write(unit(1),'(3X,A19)',advance='no')adjustl(parnme(iparm))
          if (ipart < npop) then
            write(unit(1),'(5X,1PE10.3)',advance='no')scale(iparm)*partval(ipart,iparm)&
                + offset(iparm)
          else
            write(unit(1),'(5X,1PE10.3)')scale(iparm)*partval(ipart,iparm)&
                + offset(iparm)
          end if
        end do
      end do
      !
    end if
    !
  end if
  !
  write(unit(1),'(A,/)')' '
  
end subroutine listout



subroutine writebest(gindex,basnam)
!========================================================================================
!==== This subroutine writes pbest and gbest to file at the end of the simulation.   ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::gindex
  integer::ipart,iparm,iobs,igp
  
  double precision::sum
  
  character(len=100),intent(in)::basnam
  character(len=100)::fnam
!----------------------------------------------------------------------------------------

  fnam = trim(basnam) // '.gbs'
  open(unit(2),file=trim(fnam))
  fnam = trim(basnam) // '.pbs'
  open(unit(3),file=trim(fnam))
  fnam = trim(basnam) // '.obs'
  open(unit(4),file=trim(fnam))
  
! write swarm size
  write(unit(3),'(I4.4)',advance='no')npop
  !
! write particle identifiers  
  do ipart=1,npop
    if (ipart == 1) write(unit(3),'(19X,A)',advance='no')' '
    if (ipart < npop) then
      write(unit(3),'(A,I4.4,3X)',advance='no')'particle-',ipart
    else
      write(unit(3),'(A,I4.4,3X)')'particle-',ipart
    end if
  end do
  
! write pbest
  do iparm=1,npar
    do ipart=1,npop
      if (ipart == 1) write(unit(3),'(3X,A19)',advance='no')adjustl(parnme(iparm))
      if (ipart < npop) then
        write(unit(3),'(1X,1PE15.8)',advance='no')scale(iparm)*pbest(ipart,iparm)&
            + offset(iparm)
      else
        write(unit(3),'(1X,1PE15.8)')scale(iparm)*pbest(ipart,iparm)&
            + offset(iparm)
      end if
    end do
  end do
  !
  write(unit(3),'(A)')' '
  !
! write objective values for pbest
  do igp=1,nobsgp
    do ipart=1,npop
      if (ipart == 1) write(unit(3),'(3X,A19)',advance='no')adjustl(obgnme(igp))
      if (ipart < npop) then
        write(unit(3),'(1X,1PE15.8)',advance='no')objgpopt(ipart,igp)
      else
        write(unit(3),'(1X,1PE15.8)')objgpopt(ipart,igp)
      end if
    end do
  end do
  !
  do ipart=1,npop
    if (ipart == 1) write(unit(3),'(A20,2X)',advance='no')adjustl('composite-obj')
    if (ipart < npop) then
      write(unit(3),'(1X,1PE15.8)',advance='no')objopt(ipart)
    else
      write(unit(3),'(1X,1PE15.8)')objopt(ipart)
    end if
  end do
  
! write gbest
  write(unit(2),'(I0)')1
  do iparm=1,npar
    write(unit(2),'(A19,3X,1PE15.8)')adjustl(parnme(iparm)),scale(iparm)*pbest(gindex,iparm)&
        + offset(iparm)
  end do
  write(unit(2),'(/)')
  do igp=1,nobsgp
    write(unit(2),'(A19,3X,1PE15.8)')adjustl(obgnme(igp)),objgpopt(gindex,igp)
  end do
  
! write observations
  sum = 0.0d+00
  do iobs=1,nobs
    if (iobs == 1) then
      write(unit(4),'(17X,A,4(7X,A))')'simulated ',' observed ',' residual ',' wght res ',&
          'wght res^2'
    end if
    write(unit(4),'(A15)',advance='no')adjustl(obsnme(iobs))
    write(unit(4),'(1PE14.7,3X)',advance='no')mobsval(iobs)
    write(unit(4),'(1PE14.7,3X)',advance='no')obsval(iobs)
    write(unit(4),'(1PE14.7,3X)',advance='no')mobsval(iobs) - obsval(iobs)
    write(unit(4),'(1PE14.7,3X)',advance='no')weight(iobs)*(mobsval(iobs) - obsval(iobs))
    write(unit(4),'(1PE14.7,3X)')(weight(iobs)*(mobsval(iobs) - obsval(iobs)))**2
    sum = sum + (weight(iobs)*(mobsval(iobs) - obsval(iobs)))**2
  end do
  write(unit(4),'(49X,A,1PE14.7)')'weighted sum of squared residuals:',sum
  
  close(unit(2))
  close(unit(3))
  close(unit(4))
  
end subroutine writebest



subroutine writerst(basnam)
!========================================================================================
!==== This subroutine writes all data necessary for a full restart.                  ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer::ipart,iparm,igp
  
  character(len=100),intent(in)::basnam
  character(len=100)::fnam
!----------------------------------------------------------------------------------------

  fnam = trim(basnam) // '.rst'
  open(unit(6),file=trim(fnam))
  !
! write partval
  do iparm=1,npar
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')partval(ipart,iparm)
      else
        write(unit(6),'(1X,1PE15.8)')partval(ipart,iparm)
      end if
    end do
  end do
  !
! write partvel
  do iparm=1,npar
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')partvel(ipart,iparm)
      else
        write(unit(6),'(1X,1PE15.8)')partvel(ipart,iparm)
      end if
    end do
  end do
  !
! write pbest
  do iparm=1,npar
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')pbest(ipart,iparm)
      else
        write(unit(6),'(1X,1PE15.8)')pbest(ipart,iparm)
      end if
    end do
  end do
  !
! write current obj's
  do igp=1,nobsgp
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')objgp(ipart,igp)
      else
        write(unit(6),'(1X,1PE15.8)')objgp(ipart,igp)
      end if
    end do
  end do
  !
! write best obj's
  do igp=1,nobsgp
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')objgpopt(ipart,igp)
      else
        write(unit(6),'(1X,1PE15.8)')objgpopt(ipart,igp)
      end if
    end do
  end do
  !
  close(unit(6))
  
end subroutine writerst



subroutine listpto(iiter,iptout,basnam,alpha)
!========================================================================================
!==== This subroutine writes list output at the end of each iteration for pareto.    ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::iiter
  integer,intent(inout)::iptout
  integer::ipart,ipto,irep,iparm,repcnt
  
  double precision,intent(in)::alpha
  
  character(len=100),intent(in)::basnam
  character(len=100)::fnam,rdwr
!----------------------------------------------------------------------------------------

  if (iiter == -1) then
    !
    write(unit(1),'(A,A)')       'OPTIMISATION ITERATION NO.        : ','restart'
    write(unit(1),'(A,I0)')      '   Model calls so far             : ',modeval
    write(unit(1),'(A,I0)')      '   Repository size                : ',nrepact
    write(unit(1),'(A,1PE9.3,/)')'   Alpha value                    : ',alpha
    !
  else
    !
    write(unit(1),'(A,I0)')      'OPTIMISATION ITERATION NO.        : ',iiter
    write(unit(1),'(A,I0)')      '   Model calls so far             : ',modeval
    write(unit(1),'(A,I0)')      '   Repository size                : ',nrepact
    write(unit(1),'(A,1PE9.3,/)')'   Alpha value                    : ',alpha
    !
  end if
  !
  if (suppart > 0) then
    !
    write(unit(1),'(A,/)')'   Particle positions:'
    !
!   write particle identifiers  
    do ipart=1,npop
      if (ipart == 1) write(unit(1),'(T26,A)',advance='no')' '
      if (ipart < npop) then
        write(unit(1),'(A,I3.3,3X)',advance='no')'particle-',ipart
      else
        write(unit(1),'(A,I3.3,3X,/)')'particle-',ipart
      end if
    end do
    !
  end if
  !
! write parameter values for each particle
  if (suppart > 1) then
    !
    do iparm=1,npar
      do ipart=1,npop
        if (ipart == 1) write(unit(1),'(T8,A)',advance='no')adjustl(parnme(iparm))
        if (ipart < npop) then
          write(unit(1),'(5X,1PE10.3)',advance='no')scale(iparm)*partval(ipart,iparm)&
            + offset(iparm)
        else
          write(unit(1),'(5X,1PE10.3)')scale(iparm)*partval(ipart,iparm)&
            + offset(iparm)
        end if
      end do
    end do
    !
  end if
  !
  if (suppart > 0) then
    !
    write(unit(1),'(/,A,/)')'   Particle objectives:'
    !
!   write each objective function for each particle
    do ipto=1,nptogp
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(T8,A)',advance='no')adjustl(ptonme(ipto))
        !
        if (ipart < npop) then
          !
          if (modfail(ipart) == 0) then
            write(unit(1),'(5X,1PE10.3)',advance='no')ptogp(ipart,ipto)
          else
            write(unit(1),'(5X,A)',advance='no')'  failed  '
          end if
          !
        else
          !
          if (modfail(ipart) == 0) then
            write(unit(1),'(5X,1PE10.3)')ptogp(ipart,ipto)
          else
            write(unit(1),'(5X,A)')'  failed  '
          end if
          !
        end if
        !
      end do
    end do
    !
    if (nptocon > 0) then
      !
      write(unit(1),'(/,A,/)')'   Particle constraints:'
      !
!     write particle constraints
      do ipto=1,nptocon
        do ipart=1,npop
          !
          if (ipart == 1) write(unit(1),'(T8,A)',advance='no')adjustl(ptcnme(ipto))
          !
          if (ipart < npop) then
            write(unit(1),'(5X,1PE10.3)',advance='no')ptocon(ipart,ipto)
          else
            write(unit(1),'(5X,1PE10.3)')ptocon(ipart,ipto)
          end if
          !
        end do
      end do
      !
    end if
    !
    write(unit(1),'(/,A,/)')'   Particle personal best objectives:'
    !
!   write each pbest objective function for each particle
    do ipto=1,nptogp
      do ipart=1,npop
        !
        if (ipart == 1) write(unit(1),'(T8,A)',advance='no')adjustl(ptonme(ipto))
        !
        if (ipart < npop) then
          write(unit(1),'(5X,1PE10.3)',advance='no')ptogpopt(ipart,ipto)
        else
          write(unit(1),'(5X,1PE10.3)')ptogpopt(ipart,ipto)
        end if
        !
      end do
    end do
    !
  end if
  !
  if (suprep > 0) then
    !
    write(unit(1),'(/,A,/)')'   Repository positions:'
    !
!   write particle identifiers  
    repcnt = 0
    !
    do irep=1,nrep+npop
      !
      if (irep == 1) write(unit(1),'(T27,A)',advance='no')' '
      !
      if (repindx(irep) == 1) then
        !
        repcnt = repcnt + 1
        !
        if (repcnt > nrepact) then
          write(*,*)'Error: nrepact not calculated correctly'
          write(*,*)'--stopping execution--'
         stop
        end if
        !
        if (repcnt < nrepact) then
          write(unit(1),'(A,I3.3,4X)',advance='no')'reposit-',repcnt
        else if (repcnt == nrepact) then
          write(unit(1),'(A,I3.3,4X)')'reposit-',repcnt
        end if
        !
      end if
      !
    end do
    !
    write(unit(1),'()')
    !
!   write repository positions
    if (suprep > 1) then
      !
      do iparm=1,npar
        !
        repcnt = 0
        !
        do irep=1,nrep+npop
          !
          if (irep == 1) write(unit(1),'(T8,A)',advance='no')adjustl(parnme(iparm))
          !
          if (repindx(irep) == 1) then
            !
            repcnt = repcnt + 1
            !
            if (repcnt > nrepact) then
              write(*,*)'Error: nrepact not calculated correctly'
              write(*,*)'--stopping execution--'
              stop
            end if
            !
            if (repcnt < nrepact) then
              write(unit(1),'(5X,1PE10.3)',advance='no')scale(iparm)*reposit(irep,iparm)&
                + offset(iparm)
            else if (repcnt == nrepact) then
              write(unit(1),'(5X,1PE10.3)')scale(iparm)*reposit(irep,iparm)&
                + offset(iparm)
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
    write(unit(1),'(/,A,/)')'   Repository objectives:'
    !
!   write repository objectives
    do ipto=1,nptogp
      !
      repcnt = 0
      !
      do irep=1,nrep+npop
        !
        if (irep == 1) write(unit(1),'(T8,A)',advance='no')adjustl(ptonme(ipto))
        !
        if (repindx(irep) == 1) then
          !
          repcnt = repcnt + 1
          !
          if (repcnt > nrepact) then
            write(*,*)'Error: nrepact not calculated correctly'
            write(*,*)'--stopping execution--'
            stop
          end if
          !
          if (repcnt < nrepact) then
            write(unit(1),'(5X,1PE10.3)',advance='no')repoobj(irep,ipto) 
          else if (repcnt == nrepact) then
            write(unit(1),'(5X,1PE10.3)')repoobj(irep,ipto) 
          end if
          !
        end if
        !
      end do
      !
    end do
    !
!   write repository constraints
    if (nptocon > 0) then
      !
      write(unit(1),'(/,A,/)')'   Repository constraints'
      !
      do ipto=1,nptocon
        !
        repcnt = 0
        !
        do irep=1,nrep+npop
          !
          if (irep == 1) write(unit(1),'(T8,A)',advance='no')adjustl(ptcnme(ipto))
          !
          if (repindx(irep) == 1) then
            !
            repcnt = repcnt + 1
            !
            if (repcnt > nrepact) then
              write(*,*)'Error: nrepact not calculated correctly'
              write(*,*)'--stopping execution--'
              stop
            end if
            !
            if (repcnt < nrepact) then
              write(unit(1),'(5X,1PE10.3)',advance='no')repcon(irep,ipto) 
            else if (repcnt == nrepact) then
              write(unit(1),'(5X,1PE10.3)')repcon(irep,ipto) 
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
    write(unit(1),'(/,A,/)')'   Repository fitness:'
    !
!   write repository fitness
    repcnt = 0
    ! 
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1) then
        !
        repcnt = repcnt + 1
        !
        if (repcnt > nrepact) then
          write(*,*)'Error: nrepact not calculated correctly'
          write(*,*)'--stopping execution--'
          stop
        end if
        !
        if (repcnt == 1) write(unit(1),'(22X)',advance='no')
        !
        if (repcnt < nrepact) then
          write(unit(1),'(5X,1PE10.3)',advance='no')fitness(irep) 
        else if (repcnt == nrepact) then
          write(unit(1),'(5X,1PE10.3)')fitness(irep) 
        end if
        !
      end if
      !
    end do
    !
!   write repository grid locations
    if (repmode == 1) then
      !
      write(unit(1),'(/,A,/)')'   Hypercube locations:'
      !
      do ipto=1,nptogp
        !
        repcnt = 0
        !
        do irep=1,nrep+npop
          !
          if (irep == 1) write(unit(1),'(T8,A)',advance='no')adjustl(ptonme(ipto))
          !
          if (repindx(irep) == 1) then
            !
            repcnt = repcnt + 1
            !
            if (repcnt > nrepact) then
              write(*,*)'Error: nrepact not calculated correctly'
              write(*,*)'--stopping execution--'
              stop
            end if
            !
            if (repcnt < nrepact) then
              write(unit(1),'(5X,I10)',advance='no')grd(irep,ipto)
            else if (repcnt == nrepact) then
              write(unit(1),'(5X,I10)')grd(irep,ipto)
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
    write(unit(1),'(4/)')
    !
  end if
  !
! write repository output for specified iterations
  if (iiter > 0) then
    !
    if (niterout < 0 .or. iiter == iterout(iptout)) then
      !
      if (niterout > 0) iptout = iptout + 1
      !
      write(rdwr,'(A,A,I4.4,A)')trim(basnam),'_',iiter,'.rep'
      read(rdwr,'(A)')fnam
      open(unit(5),file=trim(fnam))
      !
      do irep=1,nrep+npop
        !
        if (repindx(irep) == 1) then
          !
          do ipto=1,nptogp
            !
            if (ipto < nptogp) then
              write(unit(5),'(1PE15.8,1X)',advance='no')repoobj(irep,ipto)
            else
              write(unit(5),'(1PE15.8)')repoobj(irep,ipto)
            end if
            !
          end do
          !
        end if
        !
      end do
      !
      close(unit(5))
      !
    end if
    !
  else if (iiter == 0) then
    !
    write(rdwr,'(A,A,I4.4,A)')trim(basnam),'_',iiter,'.rep'
    read(rdwr,'(A)')fnam
    open(unit(5),file=trim(fnam))
    !
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1) then
        !
        do ipto=1,nptogp
          !
          if (ipto < nptogp) then
            write(unit(5),'(1PE15.8,1X)',advance='no')repoobj(irep,ipto)
          else
            write(unit(5),'(1PE15.8)')repoobj(irep,ipto)
          end if
          !
        end do
        !
      end if
      !
    end do
    !
    close(unit(5))
    !
  end if
  
end subroutine listpto
  
  
  
  
subroutine wrfinpareto(basnam)
!========================================================================================
!==== This subroutine writes out final repository data.                              ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer::irep,irep1,irep2,ipto,iobs,reprec,iparm
  integer,dimension(:),allocatable::reporder
  
  double precision::minobj
  
  character(len=100),intent(in)::basnam
  character(len=100)::fnam
!----------------------------------------------------------------------------------------
   
  fnam = trim(basnam) // '.rep'
  open(unit(5),file=trim(fnam))
  !
  allocate(reporder(nrep))
  !
  reporder = 0
  minobj   = 0.0d+00
  !
! organize repository index by magnitude of first objective function
  do irep1=1,nrepact
    !
    reprec = 0
    !
    do irep2=1,nrep+npop
      !
      if (repindx(irep2) == 1) then
        !
        if (repoobj(irep2,1) < minobj .or. reprec == 0) then
          !
          minobj = repoobj(irep2,1)
          reprec = irep2
          !
        end if
        !
      end if
      !
    end do
    !
    reporder(irep1) = reprec
    repindx(reprec) = 0
    !
  end do    
  !
  do irep1=1,nrepact
    !
    do ipto=1,nptogp
      !
      if (ipto < nptogp) then
        write(unit(5),'(1PE15.8,1X)',advance='no')repoobj(reporder(irep1),ipto)
      else
        write(unit(5),'(1PE15.8)')repoobj(reporder(irep1),ipto)
      end if
      !
    end do
    !
  end do
  !
  close(unit(5))
  !
  fnam = trim(basnam) // '.obs'
  open(unit(5),file=trim(fnam))
  !
  do ipto=1,nptogp 
    !
    do irep=1,nrepact
      !
      if (irep == 1) write(unit(5),'(T8,A)',advance='no')adjustl(ptonme(ipto))
      !
      if (irep < nrepact) then
        write(unit(5),'(1PE15.8,1X)',advance='no')repoobj(reporder(irep),ipto)
      else
        write(unit(5),'(1PE15.8)')repoobj(reporder(irep),ipto)
      end if
      !
    end do
    !
  end do
  !
  do iobs=1,nobs
    !
    do irep=1,nrepact
      !
      if (irep == 1) write(unit(5),'(T8,A)',advance='no')adjustl(obsnme(iobs))
      !
      if (irep > nrepact) then
        write(*,*)'Error: nrepact not calculated correctly'
        write(*,*)'--stopping execution--'
        stop
      end if
      !
      if (irep < nrepact) then
        write(unit(5),'(1PE15.8,1X)',advance='no')repoobs(reporder(irep),iobs) 
      else if (irep == nrepact) then
        write(unit(5),'(1PE15.8)')repoobs(reporder(irep),iobs) 
      end if
      !
    end do
    !
  end do
  !
  close(unit(5))
  !
  fnam = trim(basnam) // '.par'
  open(unit(5),file=trim(fnam))
  !
  write(unit(5),'(I0)')nrepact
  !
  do iparm=1,npar
    !
    do irep=1,nrepact
      !
      if (irep == 1) write(unit(5),'(T8,A)',advance='no')adjustl(parnme(iparm))
      !
      if (irep < nrepact) then
        write(unit(5),'(1PE15.8,1X)',advance='no')scale(iparm)*reposit(reporder(irep),iparm) &
            + offset(iparm)
      else if (irep == nrepact) then
        write(unit(5),'(1PE15.8)')scale(iparm)*reposit(reporder(irep),iparm) + offset(iparm)
      end if
      !
    end do
    !
  end do
  !
  close(unit(5))
  
  deallocate(reporder)
  
end subroutine wrfinpareto



subroutine wrparerst(basnam)
!========================================================================================
!==== This subroutine writes all MOPSO data necessary for a full restart.            ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer::ipart,iparm,irep,ipto,repcnt,iobs
  
  character(len=100),intent(in)::basnam
  character(len=100)::fnam
!----------------------------------------------------------------------------------------

  fnam = trim(basnam) // '.rst'
  open(unit(6),file=trim(fnam))
  !
! write partval
  write(unit(6),'(A)')'* partval'
  do iparm=1,npar
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')partval(ipart,iparm)
      else
        write(unit(6),'(1X,1PE15.8)')partval(ipart,iparm)
      end if
    end do
  end do
  !
! write partvel
  write(unit(6),'(A)')'* partvel'
  do iparm=1,npar
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')partvel(ipart,iparm)
      else
        write(unit(6),'(1X,1PE15.8)')partvel(ipart,iparm)
      end if
    end do
  end do
  !
! write obj's
  write(unit(6),'(A)')'* ptogp'
  do ipto=1,nptogp
    !
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')ptogp(ipart,ipto)
      else
        write(unit(6),'(1X,1PE15.8)')ptogp(ipart,ipto)
      end if
    end do
    !
  end do
  !
! write pbest
  write(unit(6),'(A)')'* pbest'
  do iparm=1,npar
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')pbest(ipart,iparm)
      else
        write(unit(6),'(1X,1PE15.8)')pbest(ipart,iparm)
      end if
    end do
  end do
  !
! write ptogpopt's
  write(unit(6),'(A)')'* ptogpopt'
  do ipto=1,nptogp
    !
    do ipart=1,npop
      if (ipart < npop) then
        write(unit(6),'(1X,1PE15.8)',advance='no')ptogpopt(ipart,ipto)
      else
        write(unit(6),'(1X,1PE15.8)')ptogpopt(ipart,ipto)
      end if
    end do
    !
  end do
  !
! write repository positions
  write(unit(6),'(A)')'* reposit'
  write(unit(6),'(I4)')nrepact
  do iparm=1,npar
    !
    repcnt = 0
    !
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1) then
        !
        repcnt = repcnt + 1
        !
        if (repcnt < nrepact) then
          write(unit(6),'(1X,1PE15.8)',advance='no')reposit(irep,iparm)
        else
          write(unit(6),'(1X,1PE15.8)')reposit(irep,iparm)
        end if
        !
      end if
      !
    end do
    !
  end do
  !
! write repository objectives
  write(unit(6),'(A)')'* repoobj'
  do ipto=1,nptogp
    !
    repcnt = 0
    !
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1) then
        !
        repcnt = repcnt + 1
        !
        if (repcnt < nrepact) then
          write(unit(6),'(1X,1PE15.8)',advance='no')repoobj(irep,ipto)
        else
          write(unit(6),'(1X,1PE15.8)')repoobj(irep,ipto)
        end if
        !
      end if
      !
    end do
    !
  end do
  !
! write repository model observations
  write(unit(6),'(A)')'* repoobs'
  do iobs=1,nobs
    !
    repcnt = 0
    !
    do irep=1,nrep+npop
      !
      if (repindx(irep) == 1) then
        !
        repcnt = repcnt + 1
        !
        if (repcnt < nrepact) then
          write(unit(6),'(1X,1PE15.8)',advance='no')repoobs(irep,iobs)
        else
          write(unit(6),'(1X,1PE15.8)')repoobs(irep,iobs)
        end if
        !
      end if
      !
    end do
    !
  end do
  !
  close(unit(6))
  
end subroutine wrparerst



subroutine wparunc(basnam,nsav)
!========================================================================================
!==== This subroutine writes all calibrated particles and their obj's.               ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::nsav
  integer::ipart,iparm,isav,igp,iobs
  
  character(len=100),intent(in)::basnam
  character(len=100)::fnam
  
  double precision::objsav
!----------------------------------------------------------------------------------------

  fnam = trim(basnam) // '.pun'
  open(unit(9),file=trim(fnam))
  !
! write parsav
  write(unit(9),'(I0)')nsav
  !
  do iparm=1,npar
    !
    do isav=1,nsav
      !
      if (isav == 1) write(unit(9),'(T8,A)',advance='no')adjustl(parnme(iparm))
      !
      if (isav < nsav) then
        write(unit(9),'(1PE15.8,1X)',advance='no')scale(iparm)*parsav(iparm,isav) &
            + offset(iparm)
      else if (isav == nsav) then
        write(unit(9),'(1PE15.8)')scale(iparm)*parsav(iparm,isav) + offset(iparm)
      end if
      !
    end do
    !
  end do
  !
! write objgpsav
  write(unit(9),'(/)')
  !
  do igp=1,nobsgp
    !
    do isav=1,nsav
      !
      if (isav == 1) write(unit(9),'(T8,A)',advance='no')adjustl(obgnme(igp))
      !
      if (isav < nsav) then
        write(unit(9),'(1PE15.8,1X)',advance='no')objgpsav(igp,isav)
      else if (isav == nsav) then
        write(unit(9),'(1PE15.8)')objgpsav(igp,isav)
      end if
      !
    end do
    !
  end do
  !
! write total objective for each parsav
!   write(unit(9),'(/)')
  !
  do isav=1,nsav
    !
    objsav = 0.0d+00
    !
    if (isav == 1) write(unit(9),'(T8,A)',advance='no')'total obj      '
    !
    do igp=1,nobsgp
      objsav = objsav + objgpsav(igp,isav)
    end do
    !
    if (isav < nsav) then
      write(unit(9),'(1PE15.8,1X)',advance='no')objsav
    else if (isav == nsav) then
      write(unit(9),'(1PE15.8)')objsav
    end if
    !
  end do
  !
! write observations
  write(unit(9),'(/)')
  !
  do iobs=1,nobs
    !
    do isav=1,nsav
      !
      if (isav == 1) write(unit(9),'(T8,A)',advance='no')adjustl(obsnme(iobs))
      !
      if (isav < nsav) then
        write(unit(9),'(1PE15.8,1X)',advance='no')punobs(iobs,isav)
      else if (isav == nsav) then
        write(unit(9),'(1PE15.8)')punobs(iobs,isav)
      end if
      !
    end do
    !
  end do
  !
  close(unit(9))
  
end subroutine wparunc
  


subroutine listpun(nsav)
!========================================================================================
!==== This subroutine writes all calibrated particles and their obj's.               ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::nsav
!----------------------------------------------------------------------------------------  

  backspace(unit(1))
  !
  write(unit(1),'(3X,A,I0,A)')'PARUNC: ',nsav,' calibrated particles have been saved'
  !
  if (neibr == 0) write(unit(1),'(A,/)')' '
  
end subroutine listpun




subroutine wrtneib()
!========================================================================================
!==== This subroutine writes best neihbors for parunc mode.                          ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer::ipart
!----------------------------------------------------------------------------------------  

  do ipart=1,npop
    !
    if (ipart == 1) write(unit(1),'(3X,A)',advance='no')'best neighbor   -->'
    !
    if (ipart < npop) then
      write(unit(1),'(8X,I4.4,3X)',advance='no')gneibr(ipart)
    else
      write(unit(1),'(8X,I4.4,3X)')gneibr(ipart)
    end if
    !
  end do 
  !
! write extra blank line
  write(unit(1),'(A)')' '
  
end subroutine wrtneib













