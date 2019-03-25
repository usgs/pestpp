c     ******************************************************************
      subroutine wrttpl(ntpl,tplfle,infle,npar,pnam,pval,ifail)
      implicit none
c     CTM Mar 2011
c     Routine to write model input files from template files
c     Modified from pINWRIT routine of PEST by J. Doherty

c     IFAIL  =  0  - Routine completed successfully
c     IFAIL  =  1  - TPL file does not exist
c     IFAIL  =  2  - Not all parameters listed in TPL file
c     IFAIL  =  3  - Illegal header specified in TPL file
c     IFAIL  =  4  - Error getting parameter name from template
c     IFAIL  =  5  - Error getting parameter name from template
c     IFAIL  = 10  - Error writing to model input file
c     ******************************************************************

      integer                   :: n,ntpl,npar,inst,imif,iline,
     +                             lc,j,j1,j2,ipar,precis,nopnt,
     +                             jfail,ifail,ierr,nw(npar)
      double precision          :: pval(npar),tval
      character*1               :: mark
      character*50              :: pword(npar),ftyp,tpar,pnam(npar)
      character*50              :: tplfle(ntpl),infle(ntpl)
      character*2000            :: aline
      logical                   :: lexist,eof
      
      inst=35                   !unit number for template files
      imif=40                   !unit number for model input files
      ifail=0                   !if routine fails, ifail is non-zero
      precis=0                  !0:single precision  1: double precision
      nopnt=0                   !0:point    1:no point
      j1=0
      
c     DEBUG - write list of tpl files and parameter names to file
c      open(10,file='PAR_DEBUG.dat')
c      write(10,'(a)') 'List of template file names ----->'
c      write(10,'(3a50)') ('/'//trim(tplfle(n))//'\',n=1,ntpl)
c      write(10,'(a)') 'List of parameter names ----->'
c      write(10,'(3a50)') ('/'//trim(pnam(n))//'\',n=1,npar)
c      close(10)
      
c     ensure all TPL files exist
      do n=1,ntpl
        inquire(file=tplfle(n),exist=lexist)
        if(.not.lexist) then
          write(*,10) trim(tplfle(n))
          ifail=1
        end if
      end do
      if(ifail.ne.0) return
10    format(2/,3x,'Error TPL file:',1x,a,1x,'does not exist.')

c     get the smallest field width for each parameter from the tpl files
c     to ensue each instance of a number is written exactly the same
      nw=1000
      do n=1,ntpl
        ifail=0
        ipar=1
        eof=.false.
        open(unit=inst,file=tplfle(n),iostat=ierr)
        iline=1
        read(inst,*,iostat=ierr) ftyp,mark
        if(ierr.ne.0) then
          write(*,20) trim(tplfle(n))
          ifail=3
          exit
        end if
        call lowcas(ftyp)
        if(ftyp.ne.'ptf'.and.ftyp.ne.'jtf') then
          write(*,20) trim(tplfle(n))
          ifail=3
          exit
        end if
        do
          read(inst,'(a)',iostat=ierr) aline
          if(ierr.ne.0) then
            eof=.true.
            exit
          end if
          lc=len_trim(aline)
          j2=0
          do
            if(j2.ge.lc.or.ifail.ne.0) exit
            j1=index(aline(j2+1:lc),mark)
            if(j1.ne.0) then
              j1=j1+j2
              j2=index(aline(j1+1:lc),mark)
              j2=j2+j1
              call parnam(jfail,j1,j2,tpar,aline)
              call lowcas(aline(j1:j2))
              if(jfail.ne.0) then
                write(*,30) aline(j1:j2),trim(tplfle(n))
                ifail=4
                eof=.true.
                exit
              end if
              call which1(jfail,npar,ipar,pnam,tpar)
              if(jfail.ne.0) then
                write(*,40) trim(tpar),trim(tplfle(n))
                ifail=5
                eof=.true.
                exit
              end if
              j=j2-j1+1
              nw(ipar)=min(nw(ipar),j)
            else
              exit
            end if
          end do !multiple templates on a line loop
          if(eof.or.ifail.ne.0) exit
        end do !file read loop
        close(inst)
        if(ifail.ne.0) exit
      end do !tpl file loop      
      if(ifail.ne.0) return
20    format(2/,3x,'Illegal header in TPL file:',1x,a)
30    format(2/,3x,'Error determining parameter name in template:',1x,a,
     +        /,3x,'Error in file:',1x,a)
40    format(2/,3x,'Error finding parameter name:',1x,a,
     +          1x,'in list of parameters.',
     +        /,3x,'Error in file:',1x,a)
     
c     ensure that each parameter is listed in the TPL files
c     and write the significant figures to the file
      do n=1,npar
        call lowcas(pnam(n))
        if(nw(n).eq.1000) then
          write(*,50) trim(pnam(n))
          ifail=2
        else
          call wrtsig(jfail,pval(n),pword(n),nw(n),precis,tval,nopnt)
          pval(n)=tval
        end if
      end do     
      if(ifail.ne.0) return
50    format(3x,'Error: Parameter:',1x,a,
     +       1x,'not listed in any TPL file.')

c     write the input files      
      do n=1,ntpl
        ifail=0
        ipar=1
        eof=.false.
c        write(*,*) tplfle(n)
        open(unit=inst,file=tplfle(n),iostat=ierr)
        open(unit=imif,file=infle(n),iostat=ierr)
        read(inst,*) ftyp,mark
        do
          read(inst,'(a)',iostat=ierr) aline
          if(ierr.ne.0) then
            eof=.true.
            exit
          end if
          lc=len_trim(aline)
          j2=0
          do
            if(j2.ge.lc.or.ifail.ne.0) exit
            j1=index(aline(j2+1:lc),mark)
            if(j1.ne.0) then
              j1=j1+j2
              j2=index(aline(j1+1:lc),mark)
              j2=j2+j1
              call parnam(jfail,j1,j2,tpar,aline)
              call lowcas(aline(j1:j2))
              call which1(jfail,npar,ipar,pnam,tpar)
c              do j=j1,j2
                aline(j1:j2)=' '
c              end do
              j=len_trim(pword(ipar))
              aline(j2-j+1:j2)=pword(ipar)(1:j)
            else
              exit
            end if
          end do
          write(imif,'(a)',iostat=ierr)
     +                      aline(:max(len_trim(aline),1))
          if(ierr.ne.0) then
            ifail=10
          end if
          if(ifail.ne.0) exit
        end do
        close(imif)
        close(inst)
      end do      

      end subroutine