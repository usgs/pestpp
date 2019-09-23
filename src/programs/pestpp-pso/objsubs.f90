subroutine estobjeval(ipart)
!========================================================================================
!==== This subroutine calculates the objective function for a given particle         ====
!==== position in estimation mode.                                                   ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::ipart
  integer::iobs,igp,gpnum,icon
  
  double precision::priobj
!----------------------------------------------------------------------------------------

! initialization
  obj(ipart) = 0.0d+00
  do igp=1,nobsgp
    objgp(ipart,igp) = 0.0d+00
  end do
  if (nptocon > 0) then
    do icon=1,nptocon
      ptocon(ipart,icon) = 0.0d+00
    end do
    violate(ipart) = 0.0d+00
  end if
  !
! calculate measurement objective for each observation group
  do iobs=1,nobs
    do igp=1,nobsgp
      !
      if (trim(obsgp(iobs)) == trim(obgnme(igp))) then
        !
        gpnum = igp
        !
!       toggle objective
        if (objmeth(igp) == 1) then
          objgp(ipart,gpnum) = objgp(ipart,gpnum) + (weight(iobs)*(obsval(iobs) - mobsval(iobs)))**2
        else if (objmeth(igp) == 2) then
          objgp(ipart,gpnum) = objgp(ipart,gpnum) + weight(iobs)*mobsval(iobs)
        end if
        !
        exit
        !
      end if
      !
    end do
  end do
  !
! calculate the objective functions associated with prior information
  if (nprior > 0) then
    call priobjeval(ipart,priobj)
  end if
  !
! calculate overall objective function
  do igp=1,nobsgp
    obj(ipart) = obj(ipart) + objgp(ipart,igp)
  end do
  !
! evaluate constraints
  if (nptocon > 0) then
    !
    do iobs=1,nobs
      do icon=1,nptocon
        !
        if (trim(obsgp(iobs)) == trim(ptcnme(icon))) then
          !
          gpnum = icon
          !
!         toggle objective
          if (conmeth(gpnum) == 1) then
            ptocon(ipart,gpnum) = ptocon(ipart,gpnum) + (weight(iobs)*(obsval(iobs) - mobsval(iobs)))**2
          else if (conmeth(gpnum) == 2) then
            ptocon(ipart,gpnum) = ptocon(ipart,gpnum) + weight(iobs)*mobsval(iobs)
          end if
          !
          exit
          !
        end if
        !
      end do
    end do
    !
!   check if constraints are violated
    do icon=1,nptocon
      !
      if (ptocon(ipart,icon) > ptoeps(icon)) then
        violate(ipart) = violate(ipart) + (ptocon(ipart,icon) - ptoeps(icon))
      end if
      !
    end do
    !
!   if no constraints are violated, particle is feasible
    if (violate(ipart) < 1.00d-16) then
      feas(ipart) = 1
    end if
    !
  end if
  !
  !
  !
end subroutine estobjeval




subroutine ptoobjeval(ipart)
!========================================================================================
!==== This subroutine calculates the objective function for a given particle         ====
!==== position in pareto mode.                                                       ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::ipart
  integer::iobs,igp,gpnum,ipto,icon
  
  double precision::priobj
!----------------------------------------------------------------------------------------

  do igp=1,nobsgp
    objgp(ipart,igp) = 0.0d+00
  end do
  do ipto=1,nptogp
    ptogp(ipart,ipto) = 0.0d+00
  end do
  if (nptocon > 0) then
    do ipto=1,nptocon
      ptocon(ipart,ipto) = 0.0d+00
    end do
    violate(ipart) = 0.0d+00
  end if
  !
! calculate measurement objective for each observation group
  do iobs=1,nobs
    do igp=1,nobsgp
      !
      if (trim(obsgp(iobs)) == trim(obgnme(igp))) then
        !
        gpnum = igp
        !
!       toggle objective
        if (objmeth(igp) == 1) then
          objgp(ipart,gpnum) = objgp(ipart,gpnum) + (weight(iobs)*(obsval(iobs) - mobsval(iobs)))**2
        else if (objmeth(igp) == 2) then
          objgp(ipart,gpnum) = objgp(ipart,gpnum) + weight(iobs)*mobsval(iobs)
        end if
        !
        exit
        !
      end if
      !
    end do
  end do
  !
! ! calculate measurement objective for each observation group
!   do iobs=1,nobs
!     !
!     do igp=1,nobsgp
!       if (trim(obsgp(iobs)) == trim(obgnme(igp))) gpnum = igp
!     end do
!     !
! !   toggle objective
!     if (objmeth(iobs) == 1) then
!       objgp(ipart,gpnum) = objgp(ipart,gpnum) + (weight(iobs)*(obsval(iobs) - mobsval(iobs)))**2
!     else if (objmeth(iobs) == 2) then
!       objgp(ipart,gpnum) = objgp(ipart,gpnum) + weight(iobs)*mobsval(iobs)
!     end if
!     !
!   end do
  !
! calculate the objective functions associated with prior information
  if (nprior > 0) then
    call priobjeval(ipart,priobj)
  end if
  !
! calculate composite objectives for pareto analysis
  do ipto=1,nptogp
    do igp=1,nobsgp
      !
      if (trim(ptonme(ipto)) == trim(ptgpnme(igp))) then
        ptogp(ipart,ipto) = ptogp(ipart,ipto) + ptow(igp)*objgp(ipart,igp)
      end if
      !
    end do
  end do
  !
! evaluate constraints
  if (nptocon > 0) then
    !
    do iobs=1,nobs
      do icon=1,nptocon
        !
        if (trim(obsgp(iobs)) == trim(ptcnme(icon))) then
          !
          gpnum = icon
          !
!         toggle objective
          if (conmeth(gpnum) == 1) then
            ptocon(ipart,gpnum) = ptocon(ipart,gpnum) + (weight(iobs)*(obsval(iobs) - mobsval(iobs)))**2
          else if (conmeth(gpnum) == 2) then
            ptocon(ipart,gpnum) = ptocon(ipart,gpnum) + weight(iobs)*mobsval(iobs)
          end if
          !
          exit
          !
        end if
        !
      end do
    end do
    !
!   check if constraints are violated
    do icon=1,nptocon
      !
      if (ptocon(ipart,icon) > ptoeps(icon)) then
        violate(ipart) = violate(ipart) + (ptocon(ipart,icon) - ptoeps(icon))
      end if
      !
    end do
    !
!   if no constraints are violated, particle is feasible
    if (violate(ipart) < 1.00d-16) then
      feas(ipart) = 1
    end if
    !
  end if
  !
! ! evaluate epsilon-constraints
!   if (nptocon > 0) then
!     !
!     do ipto=1,nptocon
!       !
!       do igp=1,nobsgp
!         !
!         if (trim(ptcnme(ipto)) == trim(ptgpnme(igp))) then
!           ptocon(ipart,ipto) = ptocon(ipart,ipto) + ptow(igp)*objgp(ipart,igp)
!         end if
!         !
!       end do
!       !
! !     check for constraint violation
!       if (ptocon(ipart,ipto) > ptoeps(ipto)) then
!         violate(ipart) = violate(ipart) + (ptocon(ipart,ipto) - ptoeps(ipto))
!       end if
!       !
!     end do
!     !
!     if (violate(ipart) < 1.00d-16) then
!       feas(ipart) = 1
!     end if
!     !
!   end if
    
end subroutine ptoobjeval



subroutine priobjeval(ipart,priobj)
!========================================================================================
!==== This subroutine calculates the objective function associated with prior info   ====
!==== for a given particle position.                                                 ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  integer,intent(in)::ipart
  integer::ipi,iparm,igp
  
  double precision,intent(out)::priobj
  double precision::objpi,w,targ
  
  character(len=100)::rdwr,wrd,prigp
  character(len=1)::pm
!----------------------------------------------------------------------------------------

  priobj = 0.0d+00
  !
  do ipi=1,nprior
    !
    objpi = 0.0d+00
    !
    call getwrd(pri(ipi),wrd,1)
    !
    pilbl(ipi) = trim(wrd)
    !
    call getwrd(pri(ipi),wrd,2)
    write(rdwr,'(A)')wrd
    read(rdwr,*)w
    !
    call getwrd(pri(ipi),wrd,4)
    !
    if (wrd(1:3) == 'log') then
      !
      do iparm=1,npar
        if (trim(parnme(iparm)) == wrd(5:len(trim(wrd))-1)) then
          objpi = w*log10(partval(ipart,iparm))
        end if
      end do
      !
    else
      !
      do iparm=1,npar
        if (trim(parnme(iparm)) == trim(wrd)) then
          objpi = w*partval(ipart,iparm)
        end if
      end do
      !
    end if
    !
    call getwrd(pri(ipi),wrd,5)
    !
    if (trim(wrd) == '+' .or. trim(wrd) == '-') then
      !
      pm = trim(wrd)
      !
      call getwrd(pri(ipi),wrd,6)
      write(rdwr,'(A)')wrd
      read(rdwr,*)w
      !
      call getwrd(pri(ipi),wrd,8)
      !
      if (wrd(1:3) == 'log') then
        !
        do iparm=1,npar
          !
          if (trim(parnme(iparm)) == wrd(5:len(trim(wrd))-1)) then
            !
            if (pm == '-') objpi = objpi - w*log10(partval(ipart,iparm))
            if (pm == '+') objpi = objpi + w*log10(partval(ipart,iparm))
            !
          end if
          !
        end do
        !
      else
        !
        do iparm=1,npar
          !
          if (trim(parnme(iparm)) == trim(wrd)) then
            !
            if (pm == '-') objpi = objpi - w*partval(ipart,iparm)
            if (pm == '+') objpi = objpi + w*partval(ipart,iparm)
            !
          end if
          !
        end do
        !
      end if
      !
      call getwrd(pri(ipi),wrd,10)
      !
      write(rdwr,'(A)')wrd
      read(rdwr,*)targ
      !
      call getwrd(pri(ipi),wrd,11)
      !
      write(rdwr,'(A)')wrd
      read(rdwr,*)w
      !
      call getwrd(pri(ipi),wrd,12)
      prigp = wrd
      !
    else
      !
      call getwrd(pri(ipi),wrd,6)
      !
      write(rdwr,'(A)')wrd
      read(rdwr,*)targ
      !
      call getwrd(pri(ipi),wrd,7)
      !
      write(rdwr,'(A)')wrd
      read(rdwr,*)w
      !
      call getwrd(pri(ipi),wrd,8)
      prigp = wrd
      !
    end if
    !
    do igp=1,nobsgp
      if (trim(prigp) == trim(obgnme(igp))) then
        objgp(ipart,igp) = objgp(ipart,igp) + (w*(objpi - targ))**2
      end if
    end do
    !
    priobj = priobj + (w*(objpi - targ))**2
    !
  end do

    
end subroutine priobjeval





















