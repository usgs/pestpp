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
  integer::iobs,igp,gpnum
  
  double precision::priobj
!----------------------------------------------------------------------------------------

  obj(ipart)   = 0.0d+00
  objm(ipart)  = 0.0d+00
  objp(ipart)  = 0.0d+00
  do igp=1,nobsgp
    objgp(ipart,igp) = 0.0d+00
  end do
  !
! calculate overall measurement objective
  do iobs=1,nobs
    !
    if (objmeth(iobs) == 1) then
      objm(ipart) = objm(ipart) + (weight(iobs)*(obsval(iobs) - mobsval(iobs)))**2
    else if (objmeth(iobs) == 2) then
      objm(ipart) = objm(ipart) + (weight(iobs)*mobsval(iobs))
    end if
    !
  end do
  !
! calculate measurement objective for each observation group
  do iobs=1,nobs
    !
    do igp=1,nobsgp
      if (trim(obsgp(iobs)) == trim(obgnme(igp))) gpnum = igp
    end do
    !
    if (objmeth(iobs) == 1) then
      objgp(ipart,gpnum) = objgp(ipart,gpnum) + (weight(iobs)*(obsval(iobs) - mobsval(iobs)))**2
    else if (objmeth(iobs) == 2) then
      objgp(ipart,gpnum) = objgp(ipart,gpnum) + mobsval(iobs)
    end if
    !
  end do
  !
! calculate the objective functions associated with prior information
  if (nprior > 0) then
    !
    call priobjeval(ipart,priobj)
    !
    objp(ipart) = priobj
    !
  end if
  !
! calculate composite objective function used by basic PSO
  obj(ipart) = objm(ipart) + objp(ipart)
    
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
  integer::iobs,igp,gpnum,ipto
  
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
    !
    do igp=1,nobsgp
      if (trim(obsgp(iobs)) == trim(obgnme(igp))) gpnum = igp
    end do
    !
!   toggle objective
    if (objmeth(iobs) == 1) then
      objgp(ipart,gpnum) = objgp(ipart,gpnum) + (weight(iobs)*(obsval(iobs) - mobsval(iobs)))**2
    else if (objmeth(iobs) == 2) then
      objgp(ipart,gpnum) = objgp(ipart,gpnum) + weight(iobs)*mobsval(iobs)
    end if
    !
  end do
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
! evaluate epsilon-constraints
  if (nptocon > 0) then
    !
    do ipto=1,nptocon
      !
      do igp=1,nobsgp
        !
        if (trim(ptcnme(ipto)) == trim(ptgpnme(igp))) then
          ptocon(ipart,ipto) = ptocon(ipart,ipto) + ptow(igp)*objgp(ipart,igp)
        end if
        !
      end do
      !
!     check for constraint violation
      if (ptocon(ipart,ipto) > ptoeps(ipto)) then
        violate(ipart) = violate(ipart) + (ptocon(ipart,ipto) - ptoeps(ipto))
      end if
      !
    end do
    !
    if (violate(ipart) < 1.00d-16) then
      feas(ipart) = 1
    end if
    !
  end if
    
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





















