      subroutine readins(nins,insfle,outfle,nobs,onam,oval,ifail)
      implicit none
      
      
      integer                   :: inst,immf,nins,nobs,n,nblbmx,asize,
     +                             numl,ninstr,nblc,j,isum,ins
      integer                   :: ifail,jfail,ierr
      integer,allocatable       :: lcins(:),obsn1(:),obsn2(:),ll(:),
     +                             iiobs(:)
      double precision          :: oval(nobs)
      character*1,allocatable   :: a(:),mrkdel(:)
      character*50              :: onam(nobs),insfle(nins),outfle(nins)
      character*2000            :: aline,buf
      logical                   :: lexist,eof
      
      ifail=0
      inst=35                   !unit number for instruction files
      immf=36                   !unit number for model output files
      
c     DEBUG - write list of ins files and observation names to file
c     open(62,file='OBS_DEBUG.dat')
c      write(62,'(a)') 'List of instruction file names ----->'
c      write(62,'(3a50)') ('/'//trim(insfle(n))//'\',n=1,nins)
c      write(62,'(a)') 'List of observation names ----->'
c      write(62,'(3a50)') ('/'//trim(onam(n))//'\',n=1,nobs)
c      close(62)
      
c     ensure all INS files exist
      do n=1,nins
        inquire(file=insfle(n),exist=lexist)
        if(.not.lexist) then
          write(*,10) trim(insfle(n))
          ifail=1
        end if
      end do
      if(ifail.ne.0) return
10    format(2/,3x,'Error INS file:',1x,a,1x,'does not exist.')

c     the following is from PEST routine ioctl which reads the instruction files
c     and allocates space for internal storage of instructions      
      allocate(mrkdel(nins))
      nblbmx=0
      asize=0
      numl=0
      ninstr=0
      do n=1,nins
        eof=.false.
        open(unit=inst,file=insfle(n),iostat=ierr)
        if(ierr.ne.0) then
          write(*,20) trim(insfle(n))
          ifail=2
          exit
        end if
        read(inst,'(a)',iostat=ierr) aline
        if(ierr.ne.0) then
          write(*,50) trim(insfle(n))
          ifail=3
          exit
        end if
        call tabrem(aline)
        call lowcas(aline(1:3))
        if((aline(1:3).ne.'pif').and.(aline(1:3).ne.'jif')) then
          write(*,30) trim(insfle(n))
          ifail=4
          exit
        end if
        mrkdel(n)=aline(5:5)
        if(mrkdel(n).eq.' ') then
          write(*,40) trim(insfle(n))
          ifail=5
          exit
        end if
        do
          read(inst,'(a)',iostat=ierr) aline
          if(ierr.ne.0) then
            eof=.true.
            exit
          end if
          call tabrem(aline)
          if(index(aline,mrkdel(n)).eq.0) call cmprss(aline)
          nblc=len_trim(aline)
          if(nblc.ne.0) then
            if(nblc.gt.nblbmx) nblbmx=nblc
            ninstr=ninstr+1
            do j=1,nblc
              if(aline(j:j).ne.' ') then
                if((aline(j:j).eq.'l').or.(aline(j:j).eq.'L')) then
                  numl=numl+1
                end if
                exit
              end if
            end do
            asize=asize+nblc
          end if
        end do
        close(inst)
      end do
      if(ifail.ne.0) return
      !get size of model output file names      
      nblbmx=nblbmx+1
      do n=1,nins
        asize=asize+2+len_trim(outfle(n))
      end do
      ninstr=ninstr+nins
20    format(2/,3x,'Error opening instruction file:',1x,a)
30    format(2/,3x,'Illegal header specified in instruction file:',1x,a,
     +        /,3x,'Acceptable header is PIF.')
40    format(2/,3x,'Marker delimiter not specified in file:',1x,a,
     +        /,3x,'It must be the 5th character on the first line.')
50    format(2/,3x,'Unexpected end of file:',1x,a)

c     get instructions for use with OUTRD routine      
      allocate(a(asize),lcins(ninstr))
      jfail=0
      ins=0
      isum=0
      do n=1,nins
c        eof=.false.
        open(unit=inst,file=insfle(n),iostat=ierr)
        read(inst,'(a)') aline
        ins=ins+1
        aline(1:1)=CHAR(2)
        aline(2:2)=' '
        aline(3:len(aline))=' '
        aline(3:len(aline))=outfle(n)
        lcins(ins)=isum+1
        nblc=len_trim(aline)
        do j=1,nblc
          A(j+isum)=aline(j:j)
        end do
        isum=isum+nblc
        do
          read(inst,'(a)',iostat=ierr) aline
          if(ierr.ne.0) then
c            eof=.true.
            exit
          end if
          call tabrem(aline)
          if(index(aline,mrkdel(n)).eq.0) call cmprss(aline)
          nblc=len_trim(aline)
          if(nblc.ne.0) then
            ins=ins+1
            lcins(ins)=isum+1
            do j=1,nblc
              A(j+isum)=aline(j:j)
            end do
            isum=isum+nblc
          end if
c          if(eof) exit
        end do
        close(unit=inst)
      end do

c     read model output files            
      allocate(obsn1(nobs),obsn2(nobs),ll(numl),iiobs(nobs))
      call outrd(jfail,ninstr,nins,asize,numl,nobs,nblbmx,lcins,ll,
     +           obsn1,obsn2,iiobs,oval,onam,a,mrkdel,aline,buf)
      deallocate(obsn1,obsn2,ll,iiobs)
      if (jfail.ne.0) ifail = jfail
      deallocate(mrkdel)
      deallocate(a, lcins)
      end subroutine readins