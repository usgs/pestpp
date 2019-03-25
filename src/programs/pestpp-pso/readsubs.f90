module psodat

! base data for all PSO operations
!----------------------------------------------------------------------------------------
  integer::npar,nobs,npargp,nprior,nobsgp,ntplfle,ninsfle,noptmax,npop,iseed,&
      modeval,initp,rstpso,ntied,nadj,nforg,inerti,nphistp,suppart,suprep,nitp
  integer,dimension(:),allocatable::partied,unit,modfail,objmeth
  !
  double precision::c1,c2,iinert,finert,inertia,maxrange,phiredstp
  double precision,dimension(:),allocatable::parval1,parlbnd,parubnd,obsval,weight,&
      parval,mobsval,vmax,scale,offset,base,tiedrat
  double precision,dimension(:,:),allocatable::pbest,partvel,partval,predval,objgp,&
      objgpopt,iprval
      !
  character(len=15)::precis,dpoint,rstfle,pestmode
  character(len=50),dimension(:),allocatable::tempfle,infle,insfle,outfle
  character(len=15),dimension(:),allocatable::parnme,partrans,obgnme,obsnme,obsgp,&
      pilbl,iprnme
  character(len=20),dimension(:),allocatable::comline
  character(len=100),dimension(:),allocatable::pri
!----------------------------------------------------------------------------------------
  
  
! data for estimation mode
!----------------------------------------------------------------------------------------
  integer::neibr,nneibr
  integer,dimension(:),allocatable::gneibr
  !
  double precision,dimension(:),allocatable::objmopt,objpopt,objopt,obj,objm,objp
!----------------------------------------------------------------------------------------
  
  
! data for pareto mode
!----------------------------------------------------------------------------------------
  integer::nrep,repmode,nrepact,niterout,nptogp,nptocon,noweak
  integer,dimension(:),allocatable::repindx,iterout,repmem,repenter,feas
  !
  double precision::rfit,rramp
  double precision,dimension(:),allocatable::fitness,ptow,ptoeps,violate,vioopt,ptoub
  double precision,dimension(:,:),allocatable::reposit,repoobj,repoobs,mobssav,&
      ptogp,ptogpopt,ptocon,repcon
      !
  character(len=15),dimension(:),allocatable::ptonme,ptgpnme,ptcnme
  !
! data for repmode=1
  integer::ngrid
  integer,dimension(:,:),allocatable::grd
  !
  double precision,dimension(:),allocatable::grdfit
  !
! data for repmode=2
  integer::lonorm
  !
  double precision,dimension(:),allocatable::lonely
!----------------------------------------------------------------------------------------
  
end module psodat



subroutine readpst(pstnam)
!========================================================================================
!==== This subroutine opens PEST control file and reads in the details of the        ====
!==== parameter estimation problem. Some of the wheel has be re-invented here, so    ====
!==== this subroutine will likely be deleted upon inclusion within PEST suite.       ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================
  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------
  integer::i,j,k,idat,eof,iparm,ipart
  
  double precision::rangelog,scl,off
  
  logical::fex
  
  character(len=100),intent(in)::pstnam
  character(len=100)::scrc,line,wrd,tied1,tied2
!----------------------------------------------------------------------------------------
  
  allocate(unit(10))
  unit = 0
  
! open PST file
  open(19,file=trim(pstnam))
  
  eof = 0
  
! read basic parameter estimation control data
  do i=1,2
    read(19,*)scrc
  end do
  read(19,*)rstfle,pestmode
  if (trim(pestmode) == 'pareto') then
    read(19,*)npar,nobs,npargp,nprior,nobsgp,nptogp,nptocon
  else
    read(19,*)npar,nobs,npargp,nprior,nobsgp
  end if
  read(19,*)ntplfle,ninsfle,precis,dpoint
  do i=1,3
    read(19,*)scrc
  end do
  read(19,*)noptmax,phiredstp,nphistp

! allocate and initialize
  allocate(comline(1))
  allocate(parval1(npar),parval(npar),parlbnd(npar),parubnd(npar),scale(npar),offset(npar))
  allocate(base(npar),partied(npar),tiedrat(npar))
  allocate(obsval(nobs),weight(nobs),mobsval(nobs))
  allocate(parnme(npar),partrans(npar))
  allocate(obsnme(nobs),obsgp(nobs),objmeth(nobs))
  allocate(obgnme(nobsgp))
  allocate(pri(nprior),pilbl(nprior))
  allocate(tempfle(ntplfle),infle(ntplfle))
  allocate(insfle(ninsfle),outfle(ninsfle))
  allocate(iterout(noptmax))
  wrd      = ' '
  parnme   = ' '
  partrans = ' '
  obsnme   = ' '
  obgnme   = ' '
  pilbl    = ' '
  obsgp    = ' '
  tempfle  = ' '
  infle    = ' '
  insfle   = ' '
  outfle   = ' '
  pri      = ' '
  idat     = 0
  ntied    = 0
  nadj     = 0
  modeval  = 0
  partied  = 0
  iterout  = 0
  lonorm   = 0
  objmeth  = 0
  tiedrat  = 0.0d+00
  parval1  = 0.0d+00
  parval   = 0.0d+00
  parlbnd  = 0.0d+00
  parubnd  = 0.0d+00
  obsval   = 0.0d+00
  mobsval  = 0.0d+00
  weight   = 0.0d+00
  base     = 0.0d+00
  
  do

    read(19,'(A)',iostat=eof)line
    if (eof < 0) exit
    
!   read PSO control data
    if (trim(line) == '* pso   ') then
      !
      read(19,*)rstpso,nforg,suppart
      read(19,*)npop,c1,c2,iseed
      read(19,*)initp,iinert,finert,inerti
      !
      if (initp == 2) then
        !
        read(19,'(A)')line
        inquire(file=trim(line),exist=fex)
        if (fex) then
          open(unit(8),file=trim(line))
        else
          write(*,'(3(A)/A)')'Pareto parameter file, ',trim(line),', does not exist',&
                             '-- stopping execution --'
          stop
        end if
        !
        read(unit(8),*)nitp
        !
        allocate(iprnme(npar),iprval(npar,nitp))
        iprnme = ' '
        iprval = 0.0d+00
        !
        do iparm=1,npar
          read(unit(8),*)iprnme(iparm),(iprval(iparm,ipart),ipart=1,nitp)
        end do
        !
        close(unit(8))
        !
      end if
      !
      if (trim(pestmode) == 'estimation') then
        !
        read(19,*)neibr,nneibr
        !
      else if (trim(pestmode) == 'pareto') then
        !
        read(19,*)nrep,repmode,rfit,rramp
        noweak = 1
        !
        if (repmode == 1) then
          read(19,*)ngrid
        else if (repmode == 2) then
          read(19,*)lonorm
        end if
        !
        read(19,*)suprep,niterout
        !
        if (niterout > 0) then
          do i=1,niterout
            read(19,*)iterout(i)
            if (iterout(i) > noptmax) then
              write(*,*)'Error: Check iteration output in pso section'
              write(*,*)'--stopping execution--'
              stop
            end if
          end do
        end if
        !
      end if
      !
      idat = 1
      !
!     allocate and initialize pso data
      allocate(pbest(npop,npar),partvel(npop,npar),partval(npop,npar))
      allocate(vmax(npar))
      allocate(objgp(npop,nobsgp))
      allocate(modfail(npop))
      !
      modfail = 0
      inertia = 0.0d+00
      pbest   = 0.0d+00
      partvel = 0.0d+00
      partval = 0.0d+00
      objgp   = 0.0d+00   
      !
      if (trim(pestmode) == 'estimation') then
        !
        allocate(objmopt(npop),objpopt(npop),objopt(npop),obj(npop),objm(npop),objp(npop))
        allocate(gneibr(npop))
        allocate(objgpopt(npop,nobsgp))
        !
        gneibr   = 0
        objmopt  = 0.0d+00
        objpopt  = 0.0d+00
        objopt   = 0.0d+00
        obj      = 0.0d+00
        objm     = 0.0d+00
        objp     = 0.0d+00
        objgpopt = 0.0d+00        
        !
      else if (trim(pestmode) == 'pareto') then
        !
        allocate(reposit(nrep+npop,npar))
        allocate(repoobj(nrep+npop,nptogp))
        allocate(repindx(nrep+npop),repmem(nrep+npop),repenter(nrep+npop))
        allocate(repoobs(nrep+npop,nobs))
        allocate(mobssav(npop,nobs))
        allocate(fitness(nrep+npop))
        allocate(ptow(nobsgp),ptgpnme(nobsgp))
        allocate(ptonme(nptogp))
        allocate(ptogp(npop,nptogp),ptogpopt(npop,nptogp))
        allocate(ptoub(nptogp))
        !
        if (repmode == 1) then
          allocate(grdfit(nrep+npop))
          allocate(grd(nrep+npop,nptogp))
          grd    = 0
          grdfit = 0.0d+00
        end if
        !
        if (repmode == 2) then
          allocate(lonely(nrep+npop))
          lonely = 0.0d+00
        end if
        !
        if (nptocon > 0) then
          allocate(ptcnme(nptocon),ptoeps(nptocon))
          allocate(violate(npop),vioopt(npop),feas(npop))
          allocate(ptocon(npop,nptocon))
          allocate(repcon(nrep+npop,nptocon))
          feas    = 0
          ptoeps  = 0.0d+00
          ptocon  = 0.0d+00
          repcon  = 0.0d+00
          violate = 0.0d+00
          vioopt  = 0.0d+00
          ptcnme  = ' '
        end if
        !
        nrepact  = 0
        repindx  = 0
        repmem   = 0
        repenter = 0
        reposit  = 0.0d+00
        repoobj  = 0.0d+00
        repoobs  = 0.0d+00
        mobssav  = 0.0d+00
        fitness  = 0.0d+00
        ptow     = 0.0d+00
        ptogp    = 0.0d+00
        ptogpopt = 0.0d+00
        ptoub    = 0.0d+00
        ptgpnme  = ' '
        ptonme   = ' '       
        !
      end if
      !
    end if
    !
    if (idat == 0 .and. line(1:18) == '* parameter groups') then
      write(*,'(A)')'The PSO section in the PEST control file does not exist.'
      write(*,'(A)')'-- stopping execution --'
      stop
    end if
    
!   read parameter data
    if (trim(line) == '* parameter groups') then
      do i=1,npargp
        read(19,*)scrc!--> read in parameter groups stuff - not worried about it atm.
                      !--> may want to make vmax and other pso control parameters group-specific.
      end do
    end if
    !
    if (trim(line) == '* parameter data') then
      !
      j = 0
      !
      do i=1,npar
        !
        read(19,*)parnme(i),partrans(i),scrc,parval1(i),parlbnd(i),parubnd(i),scrc,&
            scale(i),offset(i),vmax(i)
        !
        if (trim(partrans(i)) == 'none' .or. trim(partrans(i)) == 'log' .or. &
            trim(partrans(i)) == 'eqlog') nadj = nadj + 1
        !
        if (parval1(i) < parlbnd(i) .or. parval1(i) > parubnd(i)) then
          write(*,*)'Error: Parameter is beyond its bounds, check control file:'
          write(*,*)trim(parnme(i))
          write(*,*)'-- stopping execution --'
          stop
        end if
        !
        if (trim(partrans(i)) == 'log' .or. trim(partrans(i)) == 'eqlog') then
          !
          if (parlbnd(i) <= 0.0d+00 .or. parubnd(i) <= 0.0d+00) then
            write(*,*)'Error: Cannot have zero or negative bound when using log transform'
            write(*,*)'--stopping execution--'
            stop
          end if
          !
        end if
        !
!       determine maximum log10 range if partrans = eqlog is used
        if (trim(partrans(i)) == 'eqlog') then
          !
          j = j + 1
          !
          rangelog = log10(parubnd(i)/parlbnd(i))
          !
          if (rangelog > maxrange .or. j == 1) then
            maxrange = rangelog
          end if
          !
        end if
        !
!       count the number of tied parameters
        if (trim(partrans(i)) == 'tied') then
          !
          ntied = ntied + 1
          !
        end if
        !
      end do
      !
!     get maxrange and calculate vmax if partrans = eqlog is being used
      if (j > 0) then
        !
        do i=1,npar
          !
          if (trim(partrans(i)) == 'eqlog') then
            !
            base(i) = (parubnd(i)/parlbnd(i))**(1.0d+00/maxrange)
            vmax(i) = maxrange*vmax(i)
            !
          end if
          !
        end do
        !
      end if
      !
!     determined tied parameter relationships
      if (ntied > 0) then
        !
        do i=1,ntied
          !
          read(19,*)tied1,tied2
          if (trim(tied1) == '* observation groups') then
            write(*,*)'tied parameters not specified correctly, check control file'
            write(*,*)'-- stopping execution --'
            stop
          end if
          !
          do j=1,npar
            if (trim(parnme(j)) == trim(tied1)) then
              do k=1,npar
                if (trim(parnme(k)) == trim(tied2)) then
                  !
                  partied(j) = k
                  tiedrat(j) = parval1(j)/parval1(k)
                  !
                end if
              end do
            end if
          end do
          !
        end do
        !
      end if
      !
    end if
    
    
!   read pareto group information
    if (trim(pestmode) == 'pareto') then
      !
      if (trim(line) == '* pareto groups') then
        do i=1,nptogp
          read(19,*)ptonme(i),ptoub(i)
        end do
      end if
      !
    end if
    
!   read pareto constraint information
    if (trim(pestmode) == 'pareto') then
      !
      if (trim(line) == '* pareto constraints') then
        do i=1,nptocon
          read(19,*)ptcnme(i),ptoeps(i)
        end do
      end if
      !
    end if
    
!   read observation data
    if (trim(line) == '* observation groups') then
      do i=1,nobsgp
        if (trim(pestmode) == 'pareto') then
          read(19,*)obgnme(i),ptgpnme(i),ptow(i)
        else
          read(19,*)obgnme(i)
        end if
      end do
    end if
    !
    if (trim(line) == '* observation data') then
      do i=1,nobs
        read(19,*)obsnme(i),obsval(i),weight(i),obsgp(i),objmeth(i)
      end do
    end if
    
!   read model command line
    if (trim(line) == '* model command line') then
      read(19,'(A)')comline(1)
    end if
    
!   read model input/output control
    if (trim(line) == '* model input/output') then
      !
      do i=1,ntplfle
        !
        read(19,'(A)')line
        !
        call getwrd(line,wrd,1)
        !
        tempfle(i) = trim(wrd)
        !
        call getwrd(line,wrd,2)
        !
        infle(i) = trim(wrd)
        !
      end do
      !
      do i=1,ninsfle
        !
        read(19,'(A)')line
        !
        call getwrd(line,wrd,1)
        !
        insfle(i) = trim(wrd)
        !
        call getwrd(line,wrd,2)
        !
        outfle(i) = trim(wrd)
      end do
      !
    end if
   
!   read prior information section
    if (trim(line) == '* prior information') then
      !
      do i=1,nprior
        !
        read(19,'(A)')line
        !
        pri(i) = line
        !
      end do
      !
    end if
    
  end do
  
  close(19)
  
end subroutine readpst



subroutine readrst(basnam)
!========================================================================================
!==== This subroutine reads restart data from "pest-pso.rst".                        ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================
  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------
  integer::ipart,iparm,irep,iobs,igp,ipto
  integer,dimension(:),allocatable::measgp
  
  character(len=100), optional, intent(in)::basnam
  character(len=100)::scrc,fnam
!----------------------------------------------------------------------------------------

  fnam = trim(basnam) // '.rst'
  open(unit(7),file=trim(fnam))
  !
  allocate(measgp(nobsgp))
  measgp = 0
  !
  if (trim(pestmode) == 'estimation') then
    !
    do iparm=1,npar
      read(unit(7),*)(partval(ipart,iparm),ipart=1,npop)
    end do
    !
    do iparm=1,npar
      read(unit(7),*)(partvel(ipart,iparm),ipart=1,npop)
    end do
    !
    do iparm=1,npar
      read(unit(7),*)(pbest(ipart,iparm),ipart=1,npop)
    end do
    ! 
    do igp=1,nobsgp
      read(unit(7),*)(objgp(ipart,igp),ipart=1,npop)
    end do
    !
    do igp=1,nobsgp
      read(unit(7),*)(objgpopt(ipart,igp),ipart=1,npop)
    end do
    !
!   determine which obj's are associated with measurements
    do igp=1,nobsgp
      do iobs=1,nobs
        if (trim(obgnme(igp)) == trim(obsgp(iobs))) then
          measgp(igp) = 1
          exit
        end if
      end do
    end do
    !
!   calculate obj's
    do ipart=1,npop
      !
      do igp=1,nobsgp
        if (measgp(igp) == 1) then
          objm(ipart)    = objm(ipart)    + objgp(ipart,igp)
          objmopt(ipart) = objmopt(ipart) + objgpopt(ipart,igp)
        else
          objp(ipart)    = objp(ipart)    + objgp(ipart,igp)
          objpopt(ipart) = objpopt(ipart) + objgpopt(ipart,igp) 
        end if
      end do
      !
      obj(ipart)    = objm(ipart) + objp(ipart)
      objopt(ipart) = objmopt(ipart) + objpopt(ipart)
      !
    end do
    !
    close(unit(7))
    !
  else if (trim(pestmode) == 'pareto') then
    !
    read(unit(7),*)scrc
    do iparm=1,npar
      read(unit(7),*)(partval(ipart,iparm),ipart=1,npop)
    end do
    !
    read(unit(7),*)scrc
    do iparm=1,npar
      read(unit(7),*)(partvel(ipart,iparm),ipart=1,npop)
    end do
    !
    read(unit(7),*)scrc
    do ipto=1,nptogp
      read(unit(7),*)(ptogp(ipart,ipto),ipart=1,npop)
    end do
    !
    read(unit(7),*)scrc
    do iparm=1,npar
      read(unit(7),*)(pbest(ipart,iparm),ipart=1,npop)
    end do
    !
    read(unit(7),*)scrc
    do ipto=1,nptogp
      read(unit(7),*)(ptogpopt(ipart,ipto),ipart=1,npop)
    end do
    !
    read(unit(7),*)scrc
    read(unit(7),*)nrepact
    !
    do irep=1,nrepact
      repindx(irep) = 1
    end do
    !
    do iparm=1,npar
      read(unit(7),*)(reposit(irep,iparm),irep=1,nrepact)
    end do
    !
    read(unit(7),*)scrc
    do ipto=1,nptogp
      read(unit(7),*)(repoobj(irep,ipto),irep=1,nrepact)
    end do
    !
    read(unit(7),*)scrc
    do iobs=1,nobs
      read(unit(7),*)(repoobs(irep,iobs),irep=1,nrepact)
    end do
    !
    close(unit(7))
    !
  end if
  !
  deallocate(measgp)
  
  
end subroutine readrst