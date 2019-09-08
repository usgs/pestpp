module psodat

! base data for all PSO operations
!----------------------------------------------------------------------------------------
  integer::ujac,uipar
  integer::npar,nobs,npargp,nprior,nobsgp,ntplfle,ninsfle,noptmax,npop,iseed,&
      modeval,initp,rstpso,ntied,nadj,nforg,inerti,nphistp,suppart,suprep,nitp,fint,&
      sdim,nopso,pnobsgp
  integer,dimension(:),allocatable::partied,unit,modfail,objmeth
  !
  double precision::c1,c2,iinert,finert,inertia,maxrange,phiredstp,constvmax
  double precision,dimension(:),allocatable::parval1,parlbnd,parubnd,obsval,weight,&
      parval,mobsval,scale,offset,base,tiedrat,vmax
  double precision,dimension(:,:),allocatable::pbest,partvel,partval,predval,objgp,&
      objgpopt,iprval,mobssav
      !
  character(len=15)::precis,dpoint,rstfle,pestmode
!   character(len=50),dimension(:),allocatable::tempfle,infle,insfle,outfle
  character(len=15),dimension(:),allocatable::parnme,partrans,obgnme,obsnme,obsgp,&
      pilbl,iprnme,pobsnme
!   character(len=20),dimension(:),allocatable::comline
  character(len=100),dimension(:),allocatable::pri
!----------------------------------------------------------------------------------------
  
  
! data for estimation mode
!----------------------------------------------------------------------------------------
  integer::neibr,nneibr
  integer,dimension(:),allocatable::gneibr
  !
  double precision,dimension(:),allocatable::objmopt,objpopt,objopt,obj,objm,objp
!----------------------------------------------------------------------------------------


! data for parunc mode
!----------------------------------------------------------------------------------------
  integer::nreal
  !
  double precision::thresh,inerthr
  double precision,dimension(:,:),allocatable::parsav,objgpsav,jac,nproj,sproj,punobs
  !
  character(len=100)::jacfle
!----------------------------------------------------------------------------------------
  
  
! data for pareto mode
!----------------------------------------------------------------------------------------
  integer::nrep,repmode,nrepact,niterout,nptogp,nptocon,noweak
  integer,dimension(:),allocatable::repindx,iterout,repmem,repenter,feas,conmeth
  !
  double precision::rfit,rramp
  double precision,dimension(:),allocatable::fitness,ptow,ptoeps,violate,vioopt,ptoub
  double precision,dimension(:,:),allocatable::reposit,repoobj,repoobs,&
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



subroutine readpso(pstnam)
!========================================================================================
!==== This subroutine searches through the PEST control file to find the PSO control ====
!==== file which is then opened and PSO-specific input data is read. The PSO control ====
!==== file should be input as: ++PSO("path-to-pso-control-file")                     ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================
  use psodat
  !
  implicit none
  !
! specifications
!----------------------------------------------------------------------------------------
  integer::i,j,k,l,eof,endpath,foundat,iparm,ipart
  !
  character(len=100),intent(in)::pstnam
  character(len=100)::line,line2,psonam
  !
  double precision::rangelog
  !
  logical::fex
!----------------------------------------------------------------------------------------
  !
  !
  !
! open PST file
  open(19,file=trim(pstnam))
  !
  eof = 0
  endpath = 0
  line2 = ' '
  !
  !
  !
! scan PEST control file for PSO control file name
! *************************************************************************************************
  do
    !
!-- read control file line by line
    read(19,'(A)',iostat=eof)line
    !
!-- end program if pso control file not found
    if (eof < 0) then
      write(*,*)'PSO control file not specified in PEST++ control file'
      write(*,*)'-- stopping execution --'
      stop
    end if
    !
!-- if line begins with ++ read PSO control file
    do i=1,99
      !
      if (line(i:i+1) == '++') then
        !
        do j=i+2,97
          !
          if (line(j:j+3) == 'PSO(') then
            !
!           deal with user variability on spaces
            !
            l = 0
            !
            do k=j+4,99
              !
              if (line(k:k) == ')') then
                endpath = 1
                exit
              end if
              !
              if (line(k:k) /= ' ') then
                l = l + 1
                line2(l:l) = line(k:k)
              end if
              !
            end do
            !
            if (endpath == 0) then
              write(*,*)'PSO control file path too long (> 92 characters)'
              write(*,*)'-- stopping execution --'
              stop
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
    if (l > 0) exit
    !
  end do
  !
  !
  !
! open pso control file
  inquire(file=trim(line2),exist=fex)
  if (fex) then
    open(99,file=trim(line2))
  else
    write(*,*)'PSO control file, ',trim(line2),', does not exist'
    write(*,*)'-- stopping execution --'
    stop
  end if
  !
  !
  !
! read pso control data
! *************************************************************************************************
  rewind(99)
  eof = 0
  foundat = 0
  !
  do
    !
    read (99,'(A)',iostat=eof)line
    if (eof < 0 .and. foundat == 0) then
      write(*,*)'PSO control data missing.'
      write(*,*)'-- stopping execution --'
      stop
    end if
    !
    !
    if (trim(line) == '* control data') then
      !
      foundat = 1
      !
!     control data common to all modes
      read(99,*)rstpso,nobsgp,nptocon,nforg,suppart
      read(99,*)npop,c1,c2,iseed
      read(99,*)initp,constvmax,iinert,finert,inerti
      !
      if (nobsgp + nptocon /= pnobsgp) then
         write(*,*)'Number of objectives plus constraints must sum to the number of ',&
           'observation groups in PEST control file'
         write(*,*)'-- stopping execution --'
         stop
      end if
      !
!     check if a parameter list file is used to initialize part of the swarm
      if (initp == 2) then
        !
        read(99,'(A)')line
        inquire(file=trim(line),exist=fex)
        if (fex) then
          open(uipar,file=trim(line))
        else
          write(*,'(3(A)/A)')'PSO parameter file, ',trim(line),', does not exist'
          write(*,*)'-- stopping execution --'
          stop
        end if
        !
        read(uipar,*)nitp
        !
        allocate(iprnme(npar),iprval(npar,nitp))
        iprnme = ' '
        iprval = 0.0d+00
        !
        do iparm=1,npar
          read(uipar,*)iprnme(iparm),(iprval(iparm,ipart),ipart=1,nitp)
        end do
        !
        close(uipar)
        !
      end if
      !
!     read control data specific to each mode
      if (trim(pestmode) == 'estimation' .or. trim(pestmode) == 'parunc') read(99,*)neibr,nneibr
      !
      if (nneibr > npop .and. neibr > 0) then
        !
        write(*,*)'Neighborhood is larger than swarm, switching neihborhoods off...'
        neibr = 0
        !
      end if
      !
      if (trim(pestmode) == 'parunc') then
        !
        read(99,*)nreal,thresh,inerthr,sdim,nopso
        !
        if (sdim > 0) read(99,'(A)')jacfle
        !
      end if
      !
      if (trim(pestmode) == 'pareto') then
        !
        read(99,*)nrep,repmode,rfit,rramp
        noweak = 1
        !
        if (repmode == 1) then
          read(99,*)ngrid
        else if (repmode == 2) then
          read(99,*)lonorm
        end if
        !
        read(99,*)suprep,niterout
        !
        if (niterout > 0) then
          do i=1,niterout
            read(99,*)iterout(i)
            if (iterout(i) > noptmax) then
              write(*,*)'Error: Check controls for iteration output in pso section'
              write(*,*)'--stopping execution--'
              stop
            end if
          end do
        end if
        !
      end if
      !
!     allocate and initialize some global pso data
      allocate(vmax(npar))
      allocate(objgp(npop,nobsgp))
      allocate(modfail(npop))
      allocate(mobssav(npop,nobs))
      allocate(objmeth(nobsgp))
      !
      modfail = 0
      objmeth = 0
      inertia = 0.0d+00
      objgp   = 0.0d+00
      vmax    = 0.0d+00
      mobssav = 0.0d+00
      !
!     allocate and initialize data specific to each mode
      if (trim(pestmode) == 'estimation' .or. trim(pestmode) == 'parunc') then
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
        if (trim(pestmode) == 'parunc') then
          !
          allocate(parsav(npar,nreal),objgpsav(nobsgp,nreal),punobs(nobs,nreal))
          !
          objgpsav = 0.0d+00
          parsav   = 0.0d+00
          punobs   = 0.0d+00
          !
        end if
        !
      else if (trim(pestmode) == 'pareto') then
        !
!       pareto dimension fixed at 2, but may expand this in the future if there's interest
        nptogp = 2
        !
        allocate(reposit(nrep+npop,npar))
        allocate(repoobj(nrep+npop,nptogp))
        allocate(repindx(nrep+npop),repmem(nrep+npop),repenter(nrep+npop))
        allocate(repoobs(nrep+npop,nobs))
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
        nrepact  = 0
        repindx  = 0
        repmem   = 0
        repenter = 0
        reposit  = 0.0d+00
        repoobj  = 0.0d+00
        repoobs  = 0.0d+00
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
      exit
      !
    end if
    !
  end do
  !
  !
  !
! set up parameter transformations - currently, set to eqlog as default, but may want to open this
! up later on...
! *************************************************************************************************
  do i=1,npar
    !
    partrans(i) = 'eqlog'
    !
!   determine maximum log10 range if eqlog is used
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
  end do
  !
  !
! get maxrange and calculate vmax if eqlog is being used
  if (j > 0) then
    !
    do i=1,npar
      !
      if (trim(partrans(i)) == 'eqlog') then
        !
        base(i) = (parubnd(i)/parlbnd(i))**(1.0d+00/maxrange)
        vmax(i) = maxrange*constvmax
        !
      end if
      !
    end do
    !
  end if
  !
! allocate some memory
  allocate(pbest(npop,npar),partvel(npop,npar),partval(npop,npar))
  !
  pbest   = 0.0d+00
  partvel = 0.0d+00
  partval = 0.0d+00
  !
  if (trim(pestmode) == 'parunc' .and. sdim > 0) then
    !
    allocate(jac(nobs,nadj))
    !
    jac = 0.0d+00
    !
  end if
  !
  !
  !
! read pso pareto groups
! *************************************************************************************************
  if (pestmode == 'pareto') then
    !
    rewind(99)
    eof = 0
    foundat = 0
    !
    do
      !
      read (99,'(A)',iostat=eof)line
      if (eof < 0 .and. foundat == 0) then
        write(*,*)'Pareto groups missing (PESTMODE set to pareto)'
        write(*,*)'-- stopping execution --'
        stop
      end if
      !
      !
      if (trim(line) == '* pareto groups') then
        !
        do i=1,nptogp
          read(99,*)ptonme(i),ptoub(i)
        end do
        !
        exit
        !
      end if
      !
    end do
    !
  end if
  !
  !
  !
! read pso objective data
! *************************************************************************************************
  rewind(99)
  eof = 0
  foundat = 0
  !
  do
    !
    read (99,'(A)',iostat=eof)line
    if (eof < 0 .and. foundat == 0) then
      write(*,*)'PSO objective data missing.'
      write(*,*)'-- stopping execution --'
      stop
    end if
    !
    !
    if (trim(line) == '* objective data') then
      !
      foundat = 1
      !
      do i=1,nobsgp
        if (pestmode == 'estimation') read(99,*)obgnme(i),objmeth(i)
        if (pestmode == 'pareto') read(99,*)obgnme(i),objmeth(i),ptgpnme(i),ptow(i)
      end do
      !
      exit
      !
    end if
    !
  end do
  !
  !
  !
! read pso constraint data
! *************************************************************************************************
  if (nptocon > 0) then
    !
    rewind(99)
    eof = 0
    foundat = 0
    !
    do
      !
      read (99,'(A)',iostat=eof)line
      if (eof < 0 .and. foundat == 0) then
        write(*,*)'PSO constraint data missing and NCONGP > 0.'
        write(*,*)'-- stopping execution --'
        stop
      end if
      !
      !
      if (trim(line) == '* constraint data') then
        !
        foundat = 1
        !
        allocate(ptcnme(nptocon),ptoeps(nptocon),conmeth(nptocon))
        allocate(violate(npop),vioopt(npop),feas(npop))
        allocate(ptocon(npop,nptocon))
        allocate(repcon(nrep+npop,nptocon)) !<- pareto
        feas    = 0
        conmeth = 0
        ptoeps  = 0.0d+00
        ptocon  = 0.0d+00
        repcon  = 0.0d+00 !<- pareto
        violate = 0.0d+00
        vioopt  = 0.0d+00
        ptcnme  = ' '
        !
        do i=1,nptocon
          read(99,*)ptcnme(i),conmeth(i),ptoeps(i)
        end do
        !
        exit
        !
      end if
      !
    end do
    !
  end if
  !
  !
  !
end subroutine readpso


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

  double precision::scl,off

  logical::fex

  character(len=100),intent(in)::pstnam
  character(len=100)::scrc,line,wrd,tied1,tied2
!----------------------------------------------------------------------------------------
  
! open PST file
  open(19,file=trim(pstnam))
  !
  eof = 0
  !
! read basic parameter estimation control data
  do i=1,2
    read(19,*)scrc
  end do
  !
  read(19,*)rstfle,pestmode
  read(19,*)npar,nobs,npargp,nprior,pnobsgp
  read(19,*)ntplfle,ninsfle,precis,dpoint
  !
  do i=1,3
    read(19,*)scrc
  end do
  read(19,*)noptmax,phiredstp,nphistp
  !
  !
  !
! allocate and initialize
  allocate(parval1(npar),parval(npar),parlbnd(npar),parubnd(npar),scale(npar),offset(npar))
  allocate(base(npar),partied(npar),tiedrat(npar))
  allocate(obsval(nobs),weight(nobs),mobsval(nobs))
  allocate(parnme(npar),partrans(npar))
  allocate(obsnme(nobs),obsgp(nobs))
  allocate(obgnme(nobsgp),pobsnme(nobsgp))
  allocate(pri(nprior),pilbl(nprior))
  allocate(iterout(noptmax))
  wrd      = ' '
  parnme   = ' '
  partrans = ' '
  obsnme   = ' '
  obgnme   = ' '
  pobsnme  = ' '
  pilbl    = ' '
  obsgp    = ' '
  pri      = ' '
  idat     = 0
  ntied    = 0
  nadj     = 0
  modeval  = 0
  partied  = 0
  iterout  = 0
  tiedrat  = 0.0d+00
  parval1  = 0.0d+00
  parval   = 0.0d+00
  parlbnd  = 0.0d+00
  parubnd  = 0.0d+00
  obsval   = 0.0d+00
  mobsval  = 0.0d+00
  weight   = 0.0d+00
  base     = 0.0d+00
  !
  !
  !
  do
    !
!-- read control file line by line
    read(19,'(A)',iostat=eof)line
    if (eof < 0) exit
    !
    if (trim(line) == '* parameter data') then
      !
      j = 0
      fint = 0
      !
      do i=1,npar
        !
        read(19,*)parnme(i),scrc,scrc,parval1(i),parlbnd(i),parubnd(i),scrc,&
            scale(i),offset(i)
            !
        if (parval1(i) < parlbnd(i) .or. parval1(i) > parubnd(i)) then
          write(*,*)'Error: Parameter is beyond its bounds, check control file:'
          write(*,*)trim(parnme(i))
          write(*,*)'-- stopping execution --'
          stop
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
    !
    if (trim(line) == '* observation data') then
      do i=1,nobs
        read(19,*)obsnme(i),obsval(i),weight(i),obsgp(i)
      end do
    end if
    !
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
    !
  end do
  !
  close(19)
  !
  !
  !
end subroutine readpst



subroutine readjac()
!========================================================================================
!==== This subroutine reads a PEST Jacobian file.                                    ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================
  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------
  integer::iparm,iparm1,iparm2,iobs,nterm,icnt,itemp1,itemp2,j,bn,ies,irow,iadd
  
  double precision::dtemp
  double precision,dimension(:,:),allocatable::jacsav
  
  logical::fex
  
  character(len=12)::aapar
  character(len=20)::aaobs
!----------------------------------------------------------------------------------------

  inquire(file=trim(jacfle),exist=fex)
  !
  if (fex) then
    open(ujac,file=trim(jacfle),form='unformatted',status='old',access='stream')
  else
    write(*,'(3(A)/A)')'Jacobain file, ',trim(jacfle),', does not exist',&
                       '-- stopping execution --'
    stop
  end if
  !
  allocate(jacsav(nobs,nadj))
  jacsav = 0.0d+00
  !
! read Jacobian matrix from file
  read(ujac)itemp1,itemp2
  read(ujac)nterm
  !
  bn = iabs(itemp2)
  !
  do icnt=1,nterm
    !
    read(ujac)j,dtemp
    !
    ies  = (j-1)/bn + 1
    irow = j - (ies-1)*bn
    jacsav(irow,ies) = dtemp
    !
  end do
  !
! remove columns from the Jacobian (parameters) with zero sensitivity 
! (for numerical stability of DGESVD)
  ies = 0
  !
  do iparm=1,nadj
    !
    read(ujac)aapar
    !
    call lowcas(aapar)
    !
    iadd = 0
    !
    do iobs=1,nobs
      !
      if (dabs(jacsav(iobs,iparm)) > 1.00d-16) then
        !
        iadd = 1
        exit
        !
      end if
      !
    end do
    !
    if (iadd == 1) then
      !
      ies = ies + 1
      !
      do iobs=1,nobs
        jac(iobs,ies) = weight(iobs)*jacsav(iobs,iparm)
      end do
      !
    else
      !
      write(*,'(3(A))')'parameter ',trim(aapar),&
                       ' has zero sensitivity, and has therefore become fixed.'
      !
      do iparm2=1,npar
        if (trim(aapar) == trim(parnme(iparm2))) partrans(iparm2) = 'fixed'
      end do
      !
    end if
    !
  end do
  !
  nadj = ies
  !
  deallocate(jacsav)

end subroutine readjac


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
  
  character(len=100),intent(in)::basnam
  character(len=100)::scrc,fnam
!----------------------------------------------------------------------------------------

  fnam = trim(basnam) // '.rst'
  open(unit(7),file=trim(fnam))
  !
  allocate(measgp(nobsgp))
  measgp = 0
  !
  if (trim(pestmode) == 'estimation' .or. trim(pestmode) == 'parunc') then
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
