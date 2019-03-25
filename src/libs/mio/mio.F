module model_input_output_interface
      use iso_c_binding
      external hash_add_table
      external hash_which1_lower
      external hash_free
      private
       
! -- Variables

! -- Size

      integer                            :: numinfile=0  ! Number of input files
      integer                            :: numoutfile=0 ! Number of output files
      integer                            :: name_len=200 !

! -- Counters and flags

      integer                            :: mcall                ! number of model calls
      integer                            :: template_status=0    ! status of template loading and reading
      integer                            :: instruction_status=0 ! status of instruction loading and reading

! -- Instruction-related variables

      integer                            :: asize=0    ! size of "a" array
      integer                            :: numl       ! size of "ll" array
      integer                            :: ninstr     ! size of "lcins" array
      integer, allocatable, dimension(:) :: ll         ! holds line advance data
      integer, allocatable, dimension(:) :: lcins      ! pointer to instructions
      integer, allocatable, dimension(:) :: obsn1(:)   ! stores part of observation instruction
      integer, allocatable, dimension(:) :: obsn2(:)   ! stores part of observation instruction
      integer, allocatable, dimension(:) :: iiobs(:)   ! stores part of observation instruction
      character (len=1), allocatable,   dimension(:)  :: a          ! holds compressed instruction set

! -- Template-related variables

      integer                                         :: precis     ! Precision protocol
      integer                                         :: nopnt      ! decimal point protocol
      integer, allocatable, dimension(:)              :: nw         ! minimum word length of a parameter
      character (len=1),   allocatable, dimension(:)  :: mrkdel     ! Marker delimiters
      character (len=1),   allocatable, dimension(:)  :: pardel     ! Parameter delimiters
      character (len=200), allocatable,  dimension(:)  :: pword      ! Word for each parameter

! -- Filenames

      character (len=180), allocatable, dimension(:)  :: tempfile   ! Template files
      character (len=180), allocatable, dimension(:)  :: modinfile  ! Model input files
      character (len=180), allocatable, dimension(:)  :: insfile    ! Instruction files
      character (len=180), allocatable, dimension(:)  :: modoutfile ! Model output files

      character (len=200)                :: afile        ! Temporary character storage
      character (len=80)                 :: errsub       ! Character string for error header
      character*500               :: amessage=' ' ! Error message string
      character (len=20000)               :: dline        ! Character string for text storage

! -- SUBROUTINES

! -- Visible subroutines

      public mio_initialise,                        &
             mio_put_file,                          &
             mio_get_file,                          &
             mio_store_instruction_set,             &
             mio_process_template_files,            &
             mio_delete_output_files,               &
             mio_write_model_input_files,           &
             mio_read_model_output_files,           &
             mio_get_message_string,                &
             mio_get_status,                        &
             mio_get_dimensions,                    &
             mio_finalise

contains


subroutine mio_get_message_string(ifail,mess_len,amessage_out)

      implicit none
      integer    :: mess_len
      integer, intent(out)       :: ifail
      integer :: i
      character(len=mess_len) :: amessage_out

      ifail=0
      amessage_out=amessage
      !do i=1,len(amessage)
      !    amessage_out(i:i) = amessage(i:i)
      !end do
      !write(*,*)amessage_out
      !write(*,*)trim(amessage)
      return

end subroutine mio_get_message_string


! jwhite - Nov 11 2015 - removed optional args from three subroutines - these can't be used for C calling




!subroutine mio_initialise(ifail,numin,numout,npar,nobs,precision,decpoint)
subroutine mio_initialise(ifail,numin,numout,npar,nobs)

! -- Subroutine MIO_INITIALISE initialises the model_input_output module.

      implicit none      
      
      integer, intent(out)                     :: ifail   ! indicates failure condition
      integer, intent(in)                      :: numin   ! number of model input files
      integer, intent(in)                      :: numout  ! number of model output files
      integer, intent(in)                      :: npar    ! number of parameters
      integer, intent(in)                      :: nobs    ! number of observations
      !character (len=*), intent(in), optional  :: precision !jwhite
      !character (len=*), intent(in), optional  :: decpoint !jwhite
      
      integer                             :: ierr
      character (len=10)                  :: atemp
      
      errsub='Error in subroutine MIO_INITIALISE:'
      ifail=0

      if((numin.le.0).or.(numout.le.0))then
        write(amessage,10) trim(errsub)
10      format(a,' variables NUMIN and NUMOUT must both be supplied as positive.')
        ifail=1
        return
      end if
      allocate(tempfile(numin),modinfile(numin),pardel(numin),stat=ierr)
      if(ierr.ne.0)then
        write(*,*) 'numin',numin
        write(amessage,20) trim(errsub)
        ifail=1
        return
      end if
      allocate(insfile(numout),modoutfile(numout),mrkdel(numout),stat=ierr)
      if(ierr.ne.0)then
        write(*,*) 'numout',numout
        write(amessage,20) trim(errsub)
20      format(a,' cannot allocate sufficient memory to store model interface filenames.')
        ifail=1
        return
      end if
      tempfile=' '        ! an array
      modinfile=' '       ! an array
      insfile=' '         ! an array
      modoutfile=' '      ! an array
      
      numinfile=numin
      numoutfile=numout
      ! jwhite - commented out
      !if(present(precision))then
      !  atemp=adjustl(precision)
      !  call mio_casetrans(atemp,'lo')
      !  if(atemp(1:6).eq.'double')then
      !    precis=1
      !  else 
      !    precis=0
      !  end if
      !else
      !  precis=0
      !end if
      !if(present(decpoint))then
      !  atemp=adjustl(decpoint)
      !  call mio_casetrans(atemp,'lo')
      !  if(atemp(1:7).eq.'nopoint')then
      !    nopnt=1
      !  else
      !    nopnt=0
      !  end if
      !else
      !  nopnt=0
      !end if
      precis=1 !jwhite
      nopnt=0 !jwhite
      mcall=0

! -- Various work arrays are allocated.

      if((npar.le.0).or.(nobs.le.0))then
        write(amessage,30) trim(errsub)
30      format(a,' variables NPAR and NOBS must both be supplied as positive.')
        ifail=1
        return
      end if

      allocate(obsn1(nobs),obsn2(nobs),iiobs(nobs),stat=ierr)
      if(ierr.ne.0) go to 9200
      allocate(nw(npar),pword(npar),stat=ierr)
      if(ierr.ne.0) go to 9200
      
      return

9200  write(amessage,9210) trim(errsub)
9210  format(a,' cannot allocate sufficient memory to store model interface work arrays.')
      ifail=1
      return
      
end subroutine mio_initialise



subroutine mio_put_file(ifail,itype,inum,filename)

! -- Subroutine MIO_PUT_FILE supplies a filename for storage in the model_input_output
!    module.

       implicit none
       
       integer, intent(out)            :: ifail      ! indicates error condition
       integer, intent(in)             :: itype      ! type of file
       integer, intent(in)             :: inum       ! file number
       character (len=*), intent(in)   :: filename   ! name of file
       !character*180  :: filename
       integer                         :: icount,i,j,k
       
       errsub = 'Error in subroutine MIO_PUT_FILE:'

       ifail=0
       if((itype.lt.1).or.(itype.gt.4))then
         write(amessage,10) trim(errsub)
10       format(a,' supplied second argument out of range.')
         ifail=1
         return
       end if
       if((itype.eq.1).or.(itype.eq.2))then
         if((inum.lt.1).or.(inum.gt.numinfile))go to 9000
         if(itype.eq.1)then
           tempfile(inum)=filename
         else
           modinfile(inum)=filename
         end if
       else if ((itype.eq.3).or.(itype.eq.4))then
         if((inum.lt.1).or.(inum.gt.numoutfile))go to 9000
         if(itype.eq.3)then
           insfile(inum)=filename
         else
           modoutfile(inum)=filename
           if(asize.ne.0)then
             icount=0
             do i=1,asize
               if(a(i).eq.achar(2))then
                 icount=icount+1
                 if(icount.eq.inum)then
                   k=0
                   do j=i+2,i+2+len(modoutfile(1))-1
                     k=k+1
                     if(k.gt.len(filename))exit
                     a(j)=filename(k:k)
                   end do
                   go to 8000
                 end if
               end if
             end do
           end if
         end if
       end if
8000   continue

       return

9000   write(amessage,9010) trim(errsub)
9010   format(a,' supplied third argument out of range.')
       ifail=1
       return

end subroutine mio_put_file



subroutine mio_get_file(ifail,itype,inum,filename)

! -- Subroutine MIO_GET_FILE retrieves template, model input, instruction and model output filenames.

       implicit none
       
       integer, intent(out)            :: ifail      ! indicates error condition
       integer, intent(in)             :: itype      ! type of file
       integer, intent(in)             :: inum       ! file number
       character (len=*), intent(out)   :: filename   ! name of file
       
       errsub = 'Error in subroutine MIO_GET_FILE:'

       ifail=0
       if((itype.lt.1).or.(itype.gt.4))then
         write(amessage,10) trim(errsub)
10       format(a,' supplied second argument out of range.')
         ifail=1
         return
       end if
       if((itype.eq.1).or.(itype.eq.2))then
         if((inum.lt.1).or.(inum.gt.numinfile))go to 9000
         if(itype.eq.1)then
           filename=tempfile(inum)
         else
           filename=modinfile(inum)
         end if
       else if ((itype.eq.3).or.(itype.eq.4))then
         if((inum.lt.1).or.(inum.gt.numoutfile))go to 9000
         if(itype.eq.3)then
           filename=insfile(inum)
         else
           filename=modoutfile(inum)
         end if
       end if

       return

9000   write(amessage,9010) trim(errsub)
9010   format(a,' supplied third argument out of range.')
       ifail=1
       return

end subroutine mio_get_file



subroutine mio_store_instruction_set(ifail)

! -- Subroutine MIO_STORE_INSTRUCTION_SET reads all instruction files, storing the
!    instructions contained therein for more efficient later access.

        implicit none
        
        integer, intent(out)          :: ifail
        
        integer                       :: i,nblbmx,iunit,ierr,nblc,j,ins,isum,numl

        ifail=0
        errsub='Error in user-supplied instruction set for reading model output files:'

        nblbmx=0
        asize=0
        numl=0
        ninstr=0
        do 200 i=1,numoutfile
          if(insfile(i).eq.' ')then
            write(amessage,10)
10          format('Error in subroutine MIO_STORE_INSTRUCTION_SET: names have not ',  &
            'been provided for all instruction files.')
            ifail=1
            return
          end if
          if(modoutfile(i).eq.' ')then
            write(amessage,11)
11          format('Error in subroutine MIO_STORE_INSTRUCTION_SET: names have not ', &
            'been provided for all model output files.')
            ifail=1
            return
          end if
          call mio_addquote(insfile(i),afile)
          iunit=mio_nextunit()
          open(unit=iunit,file=insfile(i),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,20) trim(afile)
20          format('Cannot open instruction file ',a,'.')
            ifail=1
            return
          end if
          read(iunit,'(A)',end=9400,err=9000) dline
          call mio_remchar(dline,achar(9))
          call mio_casetrans(dline,'lo')
          if((dline(1:3).ne.'pif').and.(dline(1:3).ne.'jif'))go to 9400
          mrkdel(i)=dline(5:5)
          if(mrkdel(i).eq.' ') go to 9400
50        read(iunit,'(a)',end=180,err=9000) dline
          call mio_remchar(dline,achar(9))
          if(index(dline,mrkdel(i)).eq.0) call mio_cmprss()
          nblc=len_trim(dline)
          if(nblc.eq.0) go to 50
          if(nblc.gt.nblbmx)nblbmx=nblc
          ninstr=ninstr+1
          do 60 j=1,nblc
            if(dline(j:j).ne.' ') then
              if((dline(j:j).eq.'L').or.(dline(j:j).eq.'l')) numl=numl+1
              go to 100
            end if
60        continue
100       asize=asize+nblc
          go to 50
180       close(unit=iunit)
200     continue
        nblbmx=nblbmx+1
        do 300 i=1,numoutfile
          asize=asize+2+len(modoutfile(i))
300     continue
        ninstr=ninstr+numoutfile
        
! -- Memory is allocated for storage of instructions.

    allocate(a(asize),ll(numl),lcins(ninstr),stat=ierr)
    if(ierr.ne.0)then
      write(amessage,310)
 310      format('Cannot allocate sufficient memory to store instruction set.')
          ifail=1
          return
        end if
        a=' ' ! a is an array

! -- The instruction set is now re-read and stored.

        ins=0
        isum=0
        do 400 i=1,numoutfile
          call mio_addquote(insfile(i),afile)
          iunit=mio_nextunit()
          open(unit=iunit,file=insfile(i),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,20) trim(afile)
            ifail=1
            return
          end if
          read(iunit,*,err=9000)
          ins=ins+1
          dline(1:1)=achar(2)
          dline(2:2)=' '
          dline(3:)=modoutfile(i)
          lcins(ins)=isum+1
          nblc=len(modoutfile(i))+2
          do 320 j=1,nblc
            a(j+isum)=dline(j:j)
320       continue
          isum=isum+nblc
350       read(iunit,322,end=181,err=9000) dline
322       format(a)
          call mio_remchar(dline,achar(9))
          if(index(dline,mrkdel(i)).eq.0) call mio_cmprss()
          nblc=len_trim(dline)
          if(nblc.eq.0) go to 350
          ins=ins+1
          lcins(ins)=isum+1
          do 370 j=1,nblc
            a(j+isum)=dline(j:j)
370       continue
          isum=isum+nblc
          go to 350
181       close(unit=iunit)
400     continue
        
        instruction_status=1
        return

9000    write(amessage,9010) trim(afile)
9010    format('Unable to read instruction file ',a,'.')
        ifail=1
        return
9400    write(amessage,9410) trim(afile)
9410    format('Header of "pif" or "jif" followed by space, followed by marker delimiter ', &
        'expected on first line of instruction file ',a,'.')
        ifail=1
        return
        
end subroutine mio_store_instruction_set


subroutine mio_process_template_files(ifail,npar,apar)

! -- Subroutine MIO_PROCESS_TEMPLATE_FILES does rudmentary checking of template files.
!    However its main role is to find the smallest character width to which each
!    parameter will be written.

        implicit none
        
        integer, intent(out)                :: ifail ! indicates error condition
        integer, intent(in)                 :: npar  ! number of parameters
        !character (len=*), dimension(npar)  :: apar  ! parameter names
        character (len=name_len)  :: apar(npar)
!! Note that parameter names must be supplied in lower case.
        
        integer                             :: ipar,i,ierr,iline,nblc,j2,j1,jfail,nnw,iunit
        integer tpar_len
        character (len=10)                  :: aline
        character (len=name_len)                  :: tpar
        type(c_ptr) cptr
        tpar_len = name_len
        errsub='Error in subroutine MIO_PROCESS_TEMPLATE_FILES:'

         call hash_add_table(apar, tpar_len, npar, cptr)
        ifail=0

        if(npar.le.0)then
          write(amessage,10) trim(errsub)
10        format(a,' value supplied for NPAR is zero or negative.')
          ifail=1
          call hash_free(cptr)
          return
        end if

        ipar=1
        do 400 i=1,npar
          nw(i)=1000
400     continue
        DO 500 I=1,numinfile
          if(tempfile(i).eq.' ')then
            write(amessage,401)
401         format('Error in subroutine MIO_PROCESS_TEMPLATE_FILES: names have not ',  &
            'been provided for all template files.')
            ifail=1
            call hash_free(cptr)
            return
          end if
          if(modinfile(i).eq.' ')then
            write(amessage,402)
402         format('Error in subroutine MIO_PROCESS_TEMPLATE_FILES: names have not ',  &
            'been provided for all model input files.')
            ifail=1
            call hash_free(cptr)
            return
          end if
          call mio_addquote(tempfile(i),afile)
          iunit=mio_nextunit()
          open(unit=iunit,file=tempfile(i),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,410) trim(afile)
410         format('Cannot open template file ',a,'.')
            ifail=1
            call hash_free(cptr)
            return
          end if
          read(iunit,'(a)',err=9000,end=9200) dline
          call mio_casetrans(dline(1:3),'lo')
          if((dline(1:3).ne.'ptf').and.   &
             (dline(1:3).ne.'jtf'))go to 9200
          pardel(i)=dline(5:5)
          if(pardel(i).eq.' ') go to 9200
          iline=1
520       iline=iline+1
          read(iunit,'(a)',err=9000,end=680) dline
          nblc=len_trim(dline)
          j2=0
550       if(j2.ge.nblc) go to 520
          j1=index(dline(j2+1:nblc),pardel(i))
          if(j1.eq.0) go to 520
          j1=j1+j2
          j2=index(dline(j1+1:nblc),pardel(i))
          if(j2.eq.0)then
            call mio_writint(aline,iline)
            write(amessage,555) trim(aline),trim(afile)
555         format('Unbalanced parameter delimiters at line ',a,' of template file ',a,'.')
            ifail=1
            call hash_free(cptr)
            return
          end if
          j2=j2+j1
          call mio_parnam(jfail,j1,j2,tpar)
          if(jfail.eq.1)then
            call mio_writint(aline,iline)
            write(amessage,556) trim(aline),trim(afile)
556         format('Parameter space less than three characters wide at line ',a,  &
            ' of file ',a,'.')
            ifail=1
            call hash_free(cptr)
            return
          else if (jfail.eq.2)then
            call mio_writint(aline,iline)
            write(amessage,557) trim(aline),trim(afile)
557         format('Blank parameter space at line ',a,' of file ',a,'.')
            ifail=1
            call hash_free(cptr)
            return
          end if 
          call hash_which1_lower(cptr, tpar, tpar_len, ipar, jfail)
!          call mio_which1(jfail,npar,ipar,apar,tpar)
          if(jfail.ne.0)then
            call mio_writint(aline,iline)
            write(amessage,558) trim(tpar),trim(aline),trim(afile)
558         format('Parameter "',a,'" cited on line ',a,' of template file ',a,   &
            ' has not been supplied with a value.')
            ifail=1
            call hash_free(cptr)
            return
          end if
          nnw=j2-j1+1
          if(nnw.lt.nw(ipar)) nw(ipar)=nnw
          go to 550
680       close(unit=iunit)
500     continue
        do 800 i=1,npar
          if(nw(i).eq.1000)then
            write(amessage,690) trim(apar(i))
690         format('Parameter "',a,'" is not cited on any template file.')
            ifail=1
            call hash_free(cptr)
            return
          end if
800     continue

        template_status=1
        call hash_free(cptr)
        return

9000    write(amessage,9010) trim(afile)
9010    format('Unable to read template file ',a,'.')
        ifail=1
        call hash_free(cptr)
        return
9200    write(amessage,9210) trim(afile)
9210    format('"ptf" or "jtf" header, followed by space, followed by parameter delimiter ',  &
        'expected on first line of template file ',a,'.')
        ifail=1
        call hash_free(cptr)
        return

end subroutine mio_process_template_files


subroutine mio_delete_output_files(ifail,asldir)

! -- Subroutine MIO_DELETE_OUTPUT_FILES deletes model output files.

        implicit none
        
        integer, intent(out)                :: ifail     ! indicates error condition
        character(*), intent(in), optional  :: asldir    ! slave directory

        logical                             :: lexist
        integer                             :: ierr,jerr,iunit,i,idir,n
        character*1                         :: aa
        character*200                       :: mofile

        idir=0
        if(present(asldir))then
          if(asldir.ne.' ')then
            n=len_trim(asldir)
            aa=asldir(n:n)
            if((aa.ne.achar(47)).and.(aa.ne.achar(92)))then
!            if((aa.ne.'\').and.(aa.ne.'/'))then
              write(amessage,5)
5             format('Error in subroutine MIO_DELETE_OUTPUT_FILES: final character ', &
              'of ASLDIR variable must be "\" (PC) or "/" (UNIX).')
              ifail=1
              return
            end if
            idir=1
          end if
        end if

        ifail=0
        do i=1,numoutfile
          if(idir.eq.0)then
            mofile=modoutfile(i)
          else
            mofile=trim(asldir)//trim(modoutfile(i))
          end if
          inquire(file=mofile,exist=lexist)
          if(lexist)then
            iunit=mio_nextunit()
            open(unit=iunit,file=mofile,status='old',iostat=ierr)
            if(ierr.eq.0)then
              close(unit=iunit,status='delete',iostat=jerr)
              if(jerr.ne.0)then
                call mio_addquote(mofile,afile)
                write(amessage,10) trim(afile)
10              format('Cannot delete model output file ',a,' prior to running model.')
                ifail=1
                return
              end if
            else
              call mio_addquote(mofile,afile)
              write(amessage,10) trim(afile)
              ifail=1
              return
            end if
          end if
        end do
        
        return
        
end subroutine mio_delete_output_files 

! Note that there seems to be a bug in lf95 in that if there is an error deleting
! the file there is a compiler error rather than ierr being given a nonzero value.
        



!subroutine mio_write_model_input_files(ifail,npar,apar,pval,asldir)
subroutine mio_write_model_input_files(ifail,npar,apar,pval)

! -- Subroutine MIO_WRITE_MODEL_INPUT_FILES writes model input files based on
!    current parameter values and a set of model input template files.

        implicit none
        
        integer, intent(out)                             :: ifail  ! error condition indicator
        integer, intent(in)                              :: npar   ! number of parameters
        !character (len=*), intent(in), dimension(npar)   :: apar   ! parameter names
        character*200  :: apar(npar)
        double precision, intent(inout), dimension(npar) :: pval   ! parameter values
        !character (len=*), intent(in), optional          :: asldir ! slave directory - jwhite
        character (len=1)                                :: asldir ! slave directory - jwhite
        

        integer                          :: ipar,ipp,jfail,ifile,iunit,iunit1,iline, &
                                            lc,j1,j2,j,ierr,idir,n
        integer tpar_len
        double precision                 :: tval
        character (len=1)                :: aa
        character (len=name_len)               :: tpar
        character (len=name_len)              :: mifile
        type(c_ptr) cptr
        tpar_len = name_len
        errsub='Error writing parameters to model input file(s):'
        ifail=0
        idir=0
        !jwhite
        asldir = ''
!        if(present(asldir))then
!          if(asldir.ne.' ')then
!            n=len_trim(asldir)
!            aa=asldir(n:n)
!            if((aa.ne.achar(47)).and.(aa.ne.achar(92)))then
!!            if((aa.ne.'\').and.(aa.ne.'/'))then
!              write(amessage,5)
!5             format('Error in subroutine MIO_WRITE_MODEL_INPUT_FILES: final character ', &
!              'of ASLDIR variable must be "\" (PC) or "/" (UNIX).')
!              ifail=1
!              return
!            end if
!            idir=1
!          end if
!        end if


        call hash_add_table(apar, tpar_len, npar, cptr)
! -- Each of the parameter words is filled.

        ipar=1
        do 100 ipp=1,npar 
          call mio_wrtsig(jfail,pval(ipp),pword(ipp),nw(ipp),precis,tval,nopnt)
          if(jfail.lt.0)then
            write(amessage,10) trim(apar(ipp))
10          format('Internal error condition has arisen while attempting to write ', &
            'current value of parameter "',a,'" to model input file.')
            go to 9900
          else if (jfail.eq.1)then
            write(amessage,11) trim(errsub),trim (apar(ipp))
11          format(a,' exponent of parameter "',a,'" is too large or too small for ', &
            'single precision protocol.')
            go to 9900
          else if (jfail.eq.2)then
            write(amessage,12) trim(errsub),trim (apar(ipp))
12          format(a,' exponent of parameter "',a,'" is too large or too small for ', &
            'double precision protocol.')
            go to 9900
          else if (jfail.eq.3)then
            write(amessage,13) trim(errsub),trim (apar(ipp))
13          format(a,' field width of parameter "',a,'" on at least one template file ', &
            'is too small to represent current parameter value. The number is too large ', &
            'to fit, or too small to be represented with any precision.')
            go to 9900
          end if
          pval(ipp)=tval
100     continue

! -- Next the substitutions in the template files are made.
        do 500 ifile=1,numinfile
          call mio_addquote(tempfile(ifile),afile)
          iunit=mio_nextunit()
          open(unit=iunit,file=tempfile(ifile),status='old',iostat=ierr)
          if(ierr.ne.0)then
            write(amessage,110) trim(errsub),trim(afile)
110         format(a,' cannot open template file ',a,'.')
            go to 9900
          end if
          iunit1=mio_nextunit()
          if(idir.eq.0)then
            mifile=modinfile(ifile)
          else
            mifile=trim(asldir)//trim(modinfile(ifile))
          end if
          open(unit=iunit1,file=mifile,iostat=ierr)
          if(ierr.ne.0)then
            call mio_addquote(mifile,afile)
            write(amessage,115) trim(errsub),trim(afile)
115         format(a,' cannot open model input file ',a,' to write updated parameter ',  &
            'values prior to running model.')
            go to 9900
          end if
          read(iunit,*,err=9000)
          iline=1
120       iline=iline+1
          read(iunit,22,end=400,err=9000) dline
22        format(a)
          lc=len_trim(dline)
          j2=0
150       if(j2.ge.lc) go to 300
          j1=index(dline(j2+1:lc),pardel(ifile))
          if(j1.eq.0) go to 300
          j1=j1+j2
          j2=index(dline(j1+1:lc),pardel(ifile))
          j2=j2+j1
          call mio_parnam(jfail,j1,j2,tpar)
          call hash_which1_lower(cptr, tpar, tpar_len, ipar, jfail)
!          call mio_which1(jfail,npar,ipar,apar,tpar)
!       The following works when space bigger than pword(:nblnk(pword))
!       dline(j1:j2)=pword(ipar)(:nblnk(pword(ipar)))
          do 160 j=j1,j2
            dline(j:j)=' '
160       continue
          j=len_trim(pword(ipar))
          dline(j2-j+1:j2)=pword(ipar)(1:j)
          go to 150
300       write(iunit1,22,err=320) trim(dline)
          go to 120
320       call mio_addquote(mifile,afile)
          write(amessage,321) trim(errsub),trim(afile)
321       format(a,' cannot write to model input file ',a,'.')
          go to 9900
400       close(unit=iunit)
          close(unit=iunit1,iostat=ierr)
          if(ierr.ne.0)then
            call mio_addquote(mifile,afile)
            write(amessage,321) trim(errsub),trim(afile)
            go to 9900
          end if
500     continue
        call hash_free(cptr)
        return

9000    write(amessage,9010) trim(afile)
9010    format('Unable to read template file ',a,'.')
        go to 9900
9900    ifail=1
        call hash_free(cptr)
        return        
        
end subroutine mio_write_model_input_files



!subroutine mio_read_model_output_files(ifail,nobs,aobs,obs,instruction,asldir) !jwhite
subroutine mio_read_model_output_files(ifail,nobs,aobs,obs,instruction)

! -- Subroutine MIO_READ_MODEL_OUTPUT_FILES reads model output files using an instruction
!    set.

!! Important note: if an error condition occurs and the INSTRUCTION variable is not blank,
!! then this variable contains the offending instruction. This should be reproduced with 
!! the error message.

        implicit none
        
        integer, intent(out)                            :: ifail   ! error indicator
        integer, intent(in)                             :: nobs    ! number of observations
        !character (len=*), intent(in), dimension(nobs)  :: aobs    ! observation names
        character*200 :: aobs(nobs)
        double precision, intent(out), dimension(nobs)  :: obs     ! observation values
        !character (len=*), intent(out)                  :: instruction ! instruction of error
        character (len=20000)  :: instruction
        !character (len=*), optional                     :: asldir   ! slave working directory
        character (len=1)                               :: asldir   ! slave working directory - jwhite
        
        
        integer             :: ifile,il,jobs,cil,iobs,begins,ins,nblb,i,n,    &
                               n1,n2,insnum,nblc,dumflg,marktyp,almark,iunit, &
                               nol,jfail,insfle,mrktyp,ierr,j2,j1,n3,num1,num2,j,idir, &
                               ilstart
        integer obsnam_len
        double precision    :: rtemp
        character (len=1)   :: aa,mkrdel
        character (len=name_len)  :: obsnam
        character (len=15)  :: fmt
        character (len=10)  :: anum
        character (len=200) :: flenme
        type(c_ptr) cptr
        instruction=' '
        obs=-1.1d270    ! an array
        obsnam_len = name_len
        errsub='Error reading model output file(s):'
        ifail=0
        idir=0
        asldir = '' !jwhite
        ! jwhite - commented out
!        if(present(asldir))then
!          if(asldir.ne.' ')then
!            n=len_trim(asldir)
!            aa=asldir(n:n)
!            if((aa.ne.achar(47)).and.(aa.ne.achar(92)))then
!!            if((aa.ne.'\').and.(aa.ne.'/'))then
!              write(amessage,5)
!5             format('Error in subroutine MIO_READ_MODEL_OUTPUT_FILES: final character ', &
!              'of ASLDIR variable must be "\" (PC) or "/" (UNIX).')
!              ifail=1
!              return
!            end if
!            idir=1
!          end if
!        end if

        call hash_add_table(aobs, obsnam_len, nobs, cptr)
        mcall=mcall+1
        ifile=0
        il=0
        jobs=0
        mkrdel=mrkdel(1)
        cil=0
        iobs=1
        begins=0

        ins=1
10      if(ins.lt.ninstr)then
          nblb=lcins(ins+1)-lcins(ins)
        else
          nblb=asize-lcins(ins)+1
        end if
        instruction=' '
        do 20 i=1,nblb
          instruction(i:i)=a(lcins(ins)+i-1)
20      continue
25      n2=0
        insnum=0

50      call mio_getint(jfail,instruction,n1,n2,nblb,mkrdel)
        if(jfail.ne.0)then
          write(amessage,49) trim(errsub)
49        format(a,' missing marker delimiter in user-supplied instruction.')
          go to 9995
        end if
51      if(n1.eq.0) go to 1000
        insnum=insnum+1
        if(insnum.eq.1)then
          if(instruction(n1:n1).ne.'&') then
            mrktyp=0
            almark=1
            begins=0
            ilstart=il
          else
            if(ins.eq.insfle+1) then
              write(amessage,52) trim(errsub)
52            format(a,' first instruction line in instruction file cannot start ', &
              'with continuation character.')
              go to 9995
            end if
            if(begins.eq.1)then
              ins=ins-1
              go to 10
            end if
          end if
        end if
        if(ichar(instruction(n1:n1)).eq.2)then
          if(ifile.ne.0) close(unit=iunit)
          do 60 i=n1+1,nblb
            if(instruction(i:i).ne.' ') go to 70
60        continue
70        flenme=instruction(i:nblb)
          if(idir.eq.1) flenme=trim(asldir)//trim(flenme)
          iunit=mio_nextunit()
          do i=1,4
            open(unit=iunit,file=flenme,status='old',iostat=ierr)
            if(ierr.eq.0) exit
            call mio_wait(100)
          end do
          call mio_addquote(flenme,afile)
          if(ierr.ne.0)then
            write(amessage,71) trim (errsub),trim(afile)
71          format(a,' cannot open model output file ',a,'.')
            instruction=' '
            go to 9995
          end if
          ifile=ifile+1
          cil=0
          mkrdel=mrkdel(ifile)
          insfle=ins
          go to 1000
        else if((instruction(n1:n1).eq.'l').or.(instruction(n1:n1).eq.'L'))then
          if(il.ne.ilstart)then
            write(amessage,72) trim(errsub)
72          format(a,' line advance item can only occur at the beginning of an instruction line.')
            go to 9995
          end if
          almark=0
          il=il+1
          if(mcall.eq.1)then
            if(n2.le.n1) go to 9050     ! put in pest
            write(fmt,150) n2-n1
150         format('(i',i4,')')
            read(instruction(n1+1:n2),fmt,err=9050) nol
            ll(il)=nol
          else
            nol=ll(il)
          end if
          if(nol.gt.1) then
            do 160 i=1,nol-1
              read(iunit,*,end=9100,err=9850)
              cil=cil+1
160         continue
          end if
          read(iunit,22,end=9100,err=9850) dline
22        format(a)
          if(index(dline,char(9)).ne.0) call mio_tabrep()
          cil=cil+1
          nblc=len_trim(dline)
          mrktyp=1
          j1=0
        else if(instruction(n1:n1).eq.mkrdel)then
          if(mrktyp.eq.0)then
200         read(iunit,22,end=9100,err=9850) dline
            if(index(dline,char(9)).ne.0) call mio_tabrep()
            cil=cil+1
            j1=index(dline,instruction(n1+1:n2-1))
            if(j1.eq.0) go to 200
            nblc=len_trim(dline)
            j1=j1+n2-n1-2
            mrktyp=1
          else
            if(j1.ge.nblc) then
              if(almark.eq.1) then
                begins=1
                go to 25
              end if
              go to 9200
            end if
            j2=index(dline(j1+1:nblc),instruction(n1+1:n2-1))
            if(j2.eq.0) then
              if(almark.eq.1) then
                begins=1
                go to 25
              end if
              go to 9200
            end if
            j1=j1+j2
            j1=j1+n2-n1-2
          end if
        else if(instruction(n1:n1).eq.'&')then
          if(insnum.ne.1) then
            write(amessage,201) trim(errsub)
201         format(a,' if present, continuation character must be first instruction on ', &
            'an instruction line.')
            go to 9995
          end if
        else if((instruction(n1:n1).EQ.'w').or.(instruction(n1:n1).eq.'W'))then
          almark=0
          if(j1.ge.nblc) go to 9400
          j2=index(dline(j1+1:nblc),' ')
          if(j2.eq.0) go to 9400
          j1=j1+j2
          do 210 i=j1,nblc
            if(dline(i:i).ne.' ') go to 220
210       continue
          i=nblc+1
220       j1=i-1
        else if((instruction(n1:n1).eq.'t').or.(instruction(n1:n1).eq.'T'))then
          almark=0
          if(n2.le.n1) go to 9000       ! put in PEST
          write(fmt,150) n2-n1
          read(instruction(n1+1:n2),fmt,err=9000) j2
          if(j2.lt.j1) then
            call mio_writint(anum,cil)
            write(amessage,221) trim(errsub),trim(anum),trim(afile)
221         format(a,' backwards move to tab position not allowed on line ',a,  &
            ' of model output file ',a,'.')
            go to 9995
          end if
          j1=j2
          if(j1.gt.nblc) then
            call mio_writint(anum,cil)
            write(amessage,222) trim(errsub),trim(anum),trim(afile)
222         format(a,' tab position beyond end of line at line ',a,' of ', &
            'model output file ',a,'.')
            go to 9995
          end if
        else if((instruction(n1:n1).eq.'[').or.(instruction(n1:n1).eq.'('))then
          almark=0
          aa=instruction(n1:n1)
          jobs=jobs+1
          if(mcall.eq.1)then
            if(aa.eq.'[')then
              n3=index(instruction(n1:n2),']')
            else
              n3=index(instruction(n1:n2),')')
            end if
            if(n3.eq.0)then
              call mio_writint(anum,cil)
              write(amessage,226) trim(errsub)
226           format(a,' missing "]" or ")" character in instruction.')
              go to 9995
            end if
            n3=n3+n1-1
            obsnam=instruction(n1+1:n3-1)
            call hash_which1_lower(cptr, obsnam, obsnam_len, iobs, jfail)
!            call mio_which1(jfail,nobs,iobs,aobs,obsnam)
            if(jfail.ne.0) go to 9700
            call mio_getnum(jfail,instruction,n3,n2,num1,num2,fmt)
            IF(jfail.ne.0) then
              write(amessage,223) trim(errsub)
223           format(a,' cannot interpret user-supplied instruction for reading model ', &
              'output file.')
              go to 9995
            end if
            obsn1(jobs)=num1
            obsn2(jobs)=num2
            iiobs(jobs)=iobs
          else
            num1=obsn1(jobs)
            num2=obsn2(jobs)
            iobs=iiobs(jobs)
          end if
          if(aa.eq.'(') then
            call mio_gettot(jfail,num1,num2,nblc)
            if(jfail.ne.0)then
              call mio_writint(anum,cil)
              write(amessage,224) trim (errsub),trim(aobs(iobs)),trim(anum),   &
              trim(afile)
224           format(a,' cannot find observation "',a,'" on line ',a,     &
              ' of model output file ',a,'.')
              go to 9995
            end if
          else
            if(num1.gt.nblc)then
              call mio_writint(anum,cil)
              write(amessage,224) trim(errsub),trim(aobs(iobs)),trim(anum),trim(afile)
              go to 9995
            end if
            if(num2.gt.nblc) num2=nblc
            if(dline(num1:num2).eq.' ')then
              call mio_writint(anum,cil)
              write(amessage,224) trim(errsub),trim(aobs(iobs)),trim(anum),trim(afile)
              go to 9995
            end if
          end if
          write(fmt,250) num2-num1+1
250       format('(f',i4,'.0)')
          if(obs(iobs).gt.-1.0d270) go to 9870
          read(dline(num1:num2),fmt,err=260) obs(iobs)
          j1=num2
          go to 50
260       continue
          call mio_writint(anum,cil)
          write(amessage,261) trim(errsub),trim(aobs(iobs)),trim(anum),trim(afile)
261       format(a,' cannot read observation "',a,'" from line ',a,     &
          ' of model output file ',a,'.')
          go to 9995
        else if(instruction(n1:n1).eq.'!') then
          almark=0
          call mio_casetrans(instruction(n1+1:n2-1),'lo')
          if((n2-n1.ne.4).or.(instruction(n1+1:n2-1).ne.'dum'))then
            jobs=jobs+1
            if(mcall.eq.1) then
              obsnam=instruction(n1+1:n2-1)
              call hash_which1_lower(cptr, obsnam, obsnam_len, iobs, jfail)
!              call mio_which1(jfail,nobs,iobs,aobs,obsnam)
              if(jfail.ne.0) go to 9700
              iiobs(jobs)=iobs
            else
              iobs=iiobs(jobs)
            end if
            dumflg=0
          else
            dumflg=1
          end if
          call mio_getnxt(jfail,j1,num1,num2,nblc)
          if(jfail.ne.0) then
            if(dumflg.eq.0) then
              call mio_writint(anum,cil)
              write(amessage,224) trim(errsub),trim(aobs(iobs)),trim(anum),trim(afile)
              go to 9995
            else
              call mio_writint(anum,cil)
              write(amessage,224) trim(errsub),'dum',trim(anum),trim(afile)
              go to 9995
            end if
          end if
          write(fmt,250) num2-num1+1
          read(dline(num1:num2),fmt,err=270) rtemp
          if(dumflg.eq.0)then
            if(obs(iobs).gt.-1.0d270) go to 9870
            obs(iobs)=rtemp
          end if
          j1=num2
          go to 50
270       call mio_getint(jfail,instruction,n1,n2,nblb,mkrdel)
          if(jfail.ne.0) then
            write(amessage,271) trim(errsub)
271         format(a,' missing marker delimiter in user-supplied instruction set.')
            go to 9995
          end if  
          if(n1.eq.0)then
            if(dumflg.eq.1) go to 9900
            go to 9800
          end if
          if(instruction(n1:n1).ne.mkrdel) then
            if(dumflg.eq.1) go to 9900
            go to 9800
          end if
          j2=index(dline(j1+1:nblc),instruction(n1+1:n2-1))
          if(j2.eq.0) then
            if(dumflg.eq.1) go to 9900
            go to 9800
          end if
          num2=j1+j2-1
          if(num2.lt.num1)then
            if(dumflg.eq.1) go to 9900
            go to 9800
          end if
          write(fmt,250) num2-num1+1
          if(dumflg.eq.1)then
            read(dline(num1:num2),fmt,err=9900) rtemp
          else
            if(obs(iobs).gt.-1.0d270) go to 9870
            read(dline(num1:num2),fmt,err=9800) obs(iobs)
          end if
          j1=num2
          go to 51
        else
          write(amessage,272) trim(errsub)
272       format(a,' cannot interpret user-supplied instruction for reading model ',  &
          'output file.')
          go to 9995
        end if
        go to 50
1000    ins=ins+1
        if(ins.le.ninstr) go to 10

        if(mcall.eq.1)then
          do 1100 i=1,nobs
          do 1050 j=1,jobs
          if(iiobs(j).eq.i) go to 1100
1050      continue
          write(amessage,1051) trim(errsub),trim(aobs(i))
1051      format(a,' observation "',a,'" not referenced in the user-supplied instruction set.')
          instruction=' '
          go to 9995
1100      continue
        end if

        close(unit=iunit)
        instruction=' '
        call hash_free(cptr)
        return

9000    write(amessage,9010) trim(errsub)
9010    format(a,' cannot read tab position from user-supplied instruction.')
        go to 9995
9050    write(amessage,9060) trim(errsub)
9060    format(a,' cannot read line advance item from user-supplied instruction.')
        go to 9995
9100    write(amessage,9110) trim(errsub),trim(afile)
9110    format(a,' unexpected end to model output file ',a,'.')
        go to 9995
9200    call mio_writint(anum,cil)
        write(amessage,9210) trim(errsub),trim(anum),trim(afile)
9210    format(a,' unable to find secondary marker on line ',a,   &
        ' of model output file ',a,'.')
        go to 9995
9400    call mio_writint(anum,cil)
        write(amessage,9410) trim(errsub),trim(anum),trim(afile)
9410    format(a,' unable to find requested whitespace, or whitespace ',  &
        'precedes end of line at line ',a,' of model output file ',a,'.')
        go to 9995
9700    write(amessage,9710) trim(errsub),trim(obsnam)
9710    format(a,' observation name "',a,'" from user-supplied instruction set ',&
        'is not cited in main program input file.')
        go to 9995
9800    call mio_writint(anum,cil)
        write(amessage,9810) trim(errsub),trim(aobs(iobs)),trim(anum),trim(afile)
9810    format(a,' cannot read observation "',a,'" from line ',a,   &
        ' of model output file ',a,'.')
!        instruction=' '
        go to 9995
9850    write(amessage,9860) trim(afile)
9860    format('Unable to read model output file ',a,'.')
        instruction=' '
        go to 9995
9870    write(amessage,9880) trim(errsub),trim(aobs(iobs))
9880    format(a,' observation "',a,'" already cited in instruction set.')
        go to 9995
9900    call mio_writint(anum,cil)
        write(amessage,9810) trim(errsub),'dum',trim(anum),trim(afile)
!        instruction=' '
        go to 9995

9995    ifail=1
        mcall=mcall-1
        if(mcall.lt.0) mcall=0
        close(unit=iunit)
        call hash_free(cptr)
        return

end subroutine mio_read_model_output_files



subroutine mio_finalise(ifail)

! -- Subroutine MIO_FINALISE de-allocates memory usage by the model_input_output
!    module.

      implicit none
      
      integer, intent(out)                :: ifail
      integer                             :: ierr

      ifail=0
      
      if(allocated(tempfile))deallocate(tempfile,stat=ierr)
      if(allocated(modinfile))deallocate(modinfile,stat=ierr)
      if(allocated(pardel))deallocate(pardel,stat=ierr)
      if(allocated(insfile))deallocate(insfile,stat=ierr)
      if(allocated(modoutfile))deallocate(modoutfile,stat=ierr)
      if(allocated(mrkdel))deallocate(mrkdel,stat=ierr)
      if(allocated(nw))deallocate(nw,stat=ierr)
      if(allocated(pword))deallocate(pword,stat=ierr)
      if(allocated(obsn1))deallocate(obsn1,stat=ierr)
      if(allocated(obsn2))deallocate(obsn2,stat=ierr)
      if(allocated(iiobs))deallocate(iiobs,stat=ierr)
      if(allocated(a))deallocate(a,stat=ierr)
      if(allocated(ll))deallocate(ll,stat=ierr)
      if(allocated(lcins))deallocate(lcins,stat=ierr)
      asize=0

      return      

end subroutine mio_finalise



subroutine mio_get_status(template_status_out,instruction_status_out)

      implicit none

      integer, intent(out)  :: template_status_out,instruction_status_out

      template_status_out=template_status
      instruction_status_out=instruction_status

      return

end subroutine mio_get_status


subroutine mio_get_dimensions(numinfile_out,numoutfile_out)

      implicit none

      integer, intent(out)  :: numinfile_out,numoutfile_out

      numinfile_out=numinfile
      numoutfile_out=numoutfile

      return

end subroutine mio_get_dimensions






subroutine mio_getnxt(ifail,j1,num1,num2,nblc)

! -- Subroutine MIO_GETNXT gets the next space-delimited word.

        implicit none        

        integer ifail
        integer j1,num1,num2,nblc,i

        ifail=0
        do 20 i=j1+1,nblc
        if(dline(i:i).ne.' ') go to 50
20      continue
        ifail=1
        return
50      num1=i
        i=index(dline(num1:nblc),' ')
        if(i.eq.0) then
          num2=nblc
        else
          num2=num1+i-2
        end if

        return

end subroutine mio_getnxt


subroutine mio_gettot(ifail,j1,j2,nblc)

! -- Subroutine MIO_GETTOT determines the exact position occupied by a number.

        implicit none
        integer ifail
        integer j1,j2,nblc,i
 
        ifail=0
        if(j1.gt.nblc)then
          ifail=1
          return
        end if
        if(j2.gt.nblc)j2=nblc
        if(dline(j2:j2).eq.' ') then
          do 10 i=j2,j1,-1
          if(dline(i:i).ne.' ')then
            j2=i
            go to 100
          end if
10        continue
          ifail=1
          return
        else
          if(j2.eq.nblc) go to 100
          do 20 i=j2,nblc
          if(dline(i:i).eq.' ') then
            j2=i-1
            go to 100
          end if
20        continue
          j2=nblc
        end if
100     if(j1.eq.1) go to 200
        do 120 i=j1,1,-1
        if(dline(i:i).eq.' ') then
          j1=i+1
          go to 200
        end if
120     continue
        j1=1
200     return

end subroutine mio_gettot


subroutine mio_getint(ifail,buf,n1,n2,nblb,mrkdel)

! -- Subroutine MIO_GETINT gets the next stored instruction for processing.

        integer n1,n2,nblb,i,ii
        integer ifail
        character mrkdel
        character*(*) buf

        ifail=0
        if(n2.ge.nblb) then
          n1=0
          return
        end if
        do 10 i=n2+1,nblb
        if((buf(i:i).ne.' ').and.(ichar(buf(i:i)).ne.9)) go to 50
10      continue
        n1=0
        return
50      n1=i
        if(buf(n1:n1).ne.mrkdel)then
          i=index(buf(n1:nblb),' ')
          ii=index(buf(n1:nblb),char(9))
          if((i.eq.0).and.(ii.eq.0))then
            i=0
          else if(i.eq.0)then
            i=ii
          else if(ii.eq.0) then
            i=i
          else
            i=min(i,ii)
          end if
          if(i.ne.0) then
            n2=n1+i-2
          else
            n2=nblb
          end if
        else
          if(n1.eq.nblb)then
            ifail=1
            return
          end if
          i=index(buf(n1+1:nblb),mrkdel)
          if(i.eq.0) then
            ifail=1
            return
          end if
          n2=n1+i
        end if

        return
 end subroutine mio_getint


subroutine mio_tabrep()

! -- Subroutine MIO_TABREP replaces a tab by blank space(s).

        integer llen,i,j,k,nblc

        llen=len(dline)
        do 10 i=llen,1,-1
        if(dline(i:i).ne.' ') go to 20
10      continue
        return
20      nblc=i

        i=0
30      i=i+1
        if(i.gt.nblc)return
        if(ichar(dline(i:i)).ne.9) go to 30
        j=((i-1)/8+1)*8-i
        if(j.eq.0) then
          dline(i:i)=' '
        else
          dline(i:i)=' '
          nblc=nblc+j
          if(nblc.gt.llen) nblc=llen
          do 50 k=nblc,((i-1)/8+1)*8,-1
          dline(k:k)=dline(k-j:k-j)
50        continue
          do 60 k=i+1,min(nblc,i+j)
          dline(k:k)=' '
60        continue
          i=i+j
        end if
        go to 30

end subroutine mio_tabrep



subroutine mio_parnam(ifail,j1,j2,tpar)

! -- Subroutine MIO_PARNAM extracts a parameter name from a string.

        implicit none
        integer, intent(out)          :: ifail   ! reports error condition
        integer, intent(in)           :: j1,j2   ! beginning and end of word
        character (len=name_len), intent(out):: tpar    ! the extracted parameter
        
        integer             :: i,j

        ifail=0
        tpar=' '
        if(j2-j1.le.1) then
          ifail=1
          return
        end if
        do 10 i=j1+1,j2-1
        if(dline(i:i).eq.' ') go to 10
        go to 30
10      continue
        ifail=2
        return
30      j=min(name_len,j2-i) !jwhite 29 July 2018
        tpar(1:j)=dline(i:i+j-1)
        return

end subroutine mio_parnam



integer function mio_nextunit()

! -- Function MIO_NEXTUNIT determines the lowest unit number available for
! -- opening.

        logical::lopen

        do mio_nextunit=10,100
          inquire(unit=mio_nextunit,opened=lopen)
          if(.not.lopen) return
        end do
        write(6,10)
10      format(' *** No more unit numbers to open files ***')
        stop

end function mio_nextunit


subroutine mio_cmprss()

! -- Subroutine MIO_CMPRSS compresses an instruction line by removing excess
! -- blank characters.

        implicit none
        
        integer nblc,j

        if(dline.eq.' ') return
10      nblc=len_trim(dline)
        j=index(dline(1:nblc),'  ')
        if(j.ne.0) then
          dline(j+1:)=adjustl(dline(j+1:))
          go to 10
        end if
        return

end subroutine mio_cmprss


subroutine mio_getnum(ifail,buf,n3,n2,num1,num2,fmt)

! -- Subroutine MIO_GETNUM retrieves character positions from fixed and
!    semi-fixed observation instructions.

        integer n3,num1,num2,i,n2
        integer ifail
        character*(*) buf
        character*(*) fmt

        ifail=0
        i=index(buf(n3+1:n2),':')
        if(i.eq.0) go to 100
        write(fmt,20) i-1
20      format('(i',i3,')')
        read(buf(n3+1:n3+i-1),fmt,err=100) num1
        n3=n3+i
        i=n2-n3
        if(i.lt.1) go to 100
        write(fmt,20) i
        read(buf(n3+1:n2),fmt,err=100) num2
        return
100     ifail=1
        return

end subroutine mio_getnum


subroutine mio_which1(ifail,npar,ipar,apar,tpar)

! -- Subroutine MIO_WHICH1 finds a string in an array of strings.

        implicit none

        integer, intent(out)                :: ifail    ! error indicator
        integer, intent(in)                 :: npar     ! number of parameters
        integer, intent(inout)              :: ipar     ! where to start the search
        character (len=*), intent(in), dimension(npar)  :: apar     ! parameter names
        character (len=name_len), intent(inout)    :: tpar     ! parameter name to look for

        integer                             :: i

        ifail=0
        if((ipar.lt.1).or.(ipar.gt.npar)) ipar=1
        call mio_casetrans(tpar,'lo')
        if(tpar.eq.apar(ipar)) return
        if(ipar.ne.npar)then
          do 20 i=ipar+1,npar
            if(tpar.eq.apar(i))then
              ipar=i
              return
            end if
20        continue
        end if
        if(ipar.ne.1)then
          do 40 i=ipar-1,1,-1
          if(tpar.eq.apar(i)) then
            ipar=i
            return
          end if
40        continue
        end if
        ifail=1
        return

end subroutine mio_which1


subroutine mio_casetrans(string,hi_or_lo)

! -- Subroutine MIO_CASETRANS converts a string to upper or lower case.

        implicit none

    character (len=*), intent(inout)        :: string
    character (len=*), intent(in)           :: hi_or_lo
    character                               :: alo, ahi
    integer                                 :: inc,i

    if(hi_or_lo.eq.'lo') then
      alo='A'; ahi='Z'; inc=iachar('a')-iachar('A')
    else if(hi_or_lo.eq.'hi') then
      alo='a'; ahi='z'; inc=iachar('A')-iachar('a')
    else
          write(6,*) ' *** Illegal call to subroutine CASETRANS ***'
          stop
    endif

    do i=1,len_trim(string)
      if((string(i:i).ge.alo).and.(string(i:i).le.ahi)) &
      string(i:i)=achar(iachar(string(i:i))+inc)
    end do

    return

end subroutine mio_casetrans


subroutine mio_wait(nsec)

! -- Subroutine MIO_WAIT hangs around for NSECS hundredths of a second.

        implicit none

        integer ddate(8),iticks,iticks1,nsec

        call date_and_time(values=ddate)
        iticks=ddate(5)*360000+ddate(6)*6000+ddate(7)*100+ddate(8)/10
10      call date_and_time(values=ddate)
        iticks1=ddate(5)*360000+ddate(6)*6000+ddate(7)*100+ddate(8)/10
        if(iticks1.lt.iticks) iticks1=iticks1+8640000
        if(iticks1.lt.iticks+nsec) go to 10

        return

  end subroutine mio_wait


  subroutine mio_writint(atemp,ival)

!       Subroutine MIO_WRITINT writes an integer to a character variable.

        integer ival
        character*6 afmt
        character*(*) atemp

        afmt='(i   )'
        write(afmt(3:5),'(i3)') len(atemp)
        write(atemp,afmt)ival
        atemp=adjustl(atemp)
        return

 end subroutine mio_writint


 subroutine mio_addquote(afile,aqfile)

! -- Subroutine MIO_ADDQUOTE adds quotes to a filename if it has a space in it.

        implicit none

        character (len=*), intent(in)   :: afile
        character (len=*), intent(out)  :: aqfile
        integer nbb

        if(index(trim(afile),' ').eq.0)then
          aqfile=afile
        else
          aqfile(1:1)='"'
          aqfile(2:)=trim(afile)
          nbb=len_trim(aqfile)+1
          aqfile(nbb:nbb)='"'
        end if

        return

end subroutine mio_addquote


subroutine mio_remchar(astring,ach)

       implicit none

       character*(*), intent(inout) :: astring
       character*(*), intent(in)    :: ach

       integer ll,ii,icount

       icount=0
       ll=len_trim(ach)

10     ii=index(astring,ach)
       if(ii.eq.0) then
         if(icount.eq.0)return
         go to 20
       end if
       icount=icount+1
       astring(ii:ii-1+ll)=' '
       go to 10

20     astring=adjustl(astring)
       return

end subroutine mio_remchar



subroutine mio_wrtsig(ifail,val,word,nw,precis,tval,nopnt)
! --
! -- Subroutine WRTSIG writes a number into a confined space with maximum
! -- precision.
! --

!       failure criteria:
!           ifail= 1 ...... number too large or small for single precision type
!           ifail= 2 ...... number too large or small for double precision type
!           ifail= 3 ...... field width too small to represent number
!           ifail=-1 ...... internal error type 1
!           ifail=-2 ...... internal error type 2
!           ifail=-3 ...... internal error type 3

        integer precis,lw,pos,inc,d,p,w,j,jj,k,jexp,n,jfail,nw,epos,pp,nopnt,kexp,iflag,lexp
        integer ifail,itemp
        double precision val,tval
        character*200 tword,ttword,fmt*14
        character*200 word
        character*200 wtemp

!       The following line overcomes what appears to be a bug in the LF90
!       compiler

#ifdef LAHEY
        if(abs(val).lt.1.0d-300) val=0.0d0
#endif
        !jwhite 7/12/101 - check for nan value
        if (val.ne.val) then
            ifail = -1
            return
        endif
        
        lexp=0
        iflag=0
        word=' '
        wtemp=' '
        pos=1
        if(val.lt.0.0d0)pos=0
#ifdef USE_D_FORMAT
        write(tword,'(1PD23.15D3)') val
#else
        write(tword,'(1PE23.15E3)') val
#endif
        call mio_casetrans(tword,'hi')
        read(tword(20:23),'(i4)') jexp
        epos=1
        if(jexp.lt.0)epos=0

        jfail=0
        ifail=0
        if(precis.eq.0)then
          lw=min(15,nw)
        else
          lw=min(23,nw)
          !lw = nw
        end if

        n=0
        if(nopnt.eq.1)n=n+1
        if(pos.eq.1)n=n+1
        if(precis.eq.0)then
          if(abs(jexp).gt.38)then
            ifail=1
            return
          end if
          if(pos.eq.1) then
            if(lw.ge.13) then
              write(word,'(1pe13.7)',err=80) val
              ! jwhite 27/11/2018 - hack for padding with leading zeros
              if ((lw.ne.13).and.(nw.ge.lw)) then !ge to account for leading single space
                wtemp = adjustl(word)
                do itemp=lw,nw
                  wtemp = '0' // wtemp
                end do
                word = trim(wtemp)
              end if
              go to 200
            end if
          else
            if(lw.ge.14)then
              write(word,'(1pe14.7)',err=80) val
              ! jwhite 27/11/2018 - hack for padding with leading zeros
              if (nw.gt.lw) then
                wtemp = adjustl(word)
                do itemp=lw,nw
                  wtemp = '0' // wtemp
                end do
                word = trim(wtemp)
              end if
              go to 200
            end if
          end if
          if(lw.ge.14-n) then
            lw=14-n
            go to 80
          end if
        else
          if(abs(jexp).gt.275)then
            ifail=2
            return
          end if
          if(pos.eq.1) then
            if(lw.ge.22) then
              
#ifdef USE_D_FORMAT
              write(word,'(1PD22.15D3)',err=80) val
#else
              write(word,'(1PE22.15E3)',err=80) val
#endif        
              ! jwhite 27/11/2018 - hack for padding with leading zeros
              if ((lw.ne.22).and.(nw.ge.lw)) then !ge to account for leading single space
                wtemp = adjustl(word)
                do itemp=lw,nw
                  wtemp = '0' // wtemp
                end do
                word = trim(wtemp)
              end if
              go to 200
            end if
          else
            if(lw.ge.23) then
#ifdef USE_D_FORMAT
              write(word,'(1PD23.15D3)',err=80) val
#else
              write(word,'(1PE23.15E3)',err=80) val
#endif
              !jwhite - leading zero hack
              if (nw.gt.lw) then
#ifdef USE_D_FORMAT
              write(word,'(1PD22.15D3)',err=80) abs(val)
#else
              write(word,'(1PE22.15E3)',err=80) abs(val)
#endif
                wtemp = adjustl(word)
                do itemp=lw,nw-1
                  wtemp = '0' // wtemp
                end do
                wtemp = '-' // wtemp
                word = trim(wtemp)
              end if
              go to 200
            end if
          end if
          if(lw.ge.23-n)then
            lw=23-n
            go to 80
          end if
        end if

        if(nopnt.eq.1)then
          if((jexp.eq.lw-2+pos).or.(jexp.eq.lw-3+pos))then
            write(fmt,15)lw+1
15          format('(f',i2,'.0)')
            write(word,fmt,err=19) val
            if(index(word,'*').ne.0) go to 19
            if(word(1:1).eq.' ') go to 19
            word(lw+1:lw+1)=' '
            go to 200
          end if
        end if
19      d=min(lw-2+pos,lw-jexp-3+pos)
20      if(d.lt.0) go to 80
        write(fmt,30) lw,d
30      format('(f',i2,'.',i2,')')
        write(word,fmt,err=80) val
        if(index(word,'*').ne.0) then
          d=d-1
          go to 20
        end if
        k=index(word,'.')
        if(k.eq.0)then
          ifail=-1
          return
        end if
        if((k.eq.1).or.((pos.eq.0).and.(k.eq.2)))then
          do 70 j=1,3
          if(k+j.gt.lw) go to 75
          if(word(k+j:k+j).ne.'0') go to 200
70        continue
          go to 80
75        ifail=3
          return
        end if
        go to 200

80      word=' '
        if(nopnt.eq.0)then
          d=lw-7
          if(pos.eq.1) d=d+1
          if(epos.eq.1) d=d+1
          if(abs(jexp).lt.100) d=d+1
          if(abs(jexp).lt.10) d=d+1
          if((jexp.ge.100).and.(jexp-(d-1).lt.100))then
            p=1+(jexp-99)
            d=d+1
            lexp=99
          else if((jexp.ge.10).and.(jexp-(d-1).lt.10))then
            p=1+(jexp-9)
            d=d+1
            lexp=9
          else if((jexp.eq.-10).or.(jexp.eq.-100)) then
            iflag=1
            d=d+1
          else
            p=1
          end if
          inc=0
85        if(d.le.0) go to 300
          if(iflag.eq.0)then
            write(fmt,100,err=300) p,d+7,d-1
          else
            write(fmt,100,err=300) 0,d+8,d
          end if
          write(tword,fmt) val
          call mio_casetrans(tword,'hi')
          if(iflag.eq.1) go to 87
          read(tword(d+4:d+7),'(i4)',err=500) kexp
          if(((kexp.eq.10).and.((jexp.eq.9).or.(lexp.eq.9))).or.     &
          ((kexp.eq.100).and.((jexp.eq.99).or.lexp.eq.99))) then
            if(inc.eq.0)then
              if(lexp.eq.0)then
                if(d-1.eq.0) then
                  d=d-1
                else
                  p=p+1
                end if
              else if(lexp.eq.9)then
                if(jexp-(d-2).lt.10) then
                  p=p+1
                else
                  d=d-1
                end if
              else if(lexp.eq.99)then
                if(jexp-(d-2).lt.100)then
                  p=p+1
                else
                  d=d-1
                end if
              end if
              inc=inc+1
              go to 85
            end if
          end if
#ifdef USE_D_FORMAT
87        j=index(tword,'D')
#else
87        j=index(tword,'E')
#endif
          go to 151
        end if
        inc=0
        p=lw-2
        pp=jexp-(p-1)
        if(pp.ge.10)then
          p=p-1
          if(pp.ge.100)p=p-1
        else if(pp.lt.0)then
          p=p-1
          if(pp.le.-10)then
            p=p-1
            if(pp.le.-100)p=p-1
          end if
        end if
        if(pos.eq.0)p=p-1
90      continue
        d=p-1
        w=d+8
        write(fmt,100) p,w,d
        if(d.lt.0)then
          if(jfail.eq.1) go to 300
          jfail=1
          p=p+1
          go to 90
        end if
#ifdef USE_D_FORMAT
100     format('(',I2,'pD',I2,'.',I2,'D3)')
#else
100     format('(',I2,'pE',I2,'.',I2,'E3)')
#endif
        write(tword,fmt) val
        call mio_casetrans(tword,'hi')
#ifdef USE_D_FORMAT
        j=index(tword,'D')
#else
        j=index(tword,'E')
#endif
        if(tword(j-1:j-1).ne.'.')then
          ifail=-1
          return
        end if
        n=1
        if(tword(j+1:j+1).eq.'-') n=n+1
        if(tword(j+2:j+2).ne.'0') then
          n=n+2
          go to 120
        end if
        if(tword(j+3:j+3).ne.'0') n=n+1
120     n=n+1
        if(j+n-2-pos.lt.lw)then
          if(inc.eq.-1) go to 150
          ttword=tword
          p=p+1
          inc=1
          go to 90
        else if(j+n-2-pos.eq.lw) then
          go to 150
        else
          if(inc.eq.1)then
            tword=ttword
            go to 150
          end if
          if(jfail.eq.1) go to 300
          p=p-1
          inc=-1
          go to 90
        end if

150     j=index(tword,'.')
151     if(pos.eq.0)then
          k=1
        else
         k=2
        end if
        word(1:j-k)=tword(k:j-1)
        jj=j
        j=j-k+1
        !if(precis.eq.0)then
        !  word(j:j)='E'
        !else
        !  word(j:j)='E'
        !end if
#ifdef USE_D_FORMAT
        word(j:j) = 'D'
#else
        word(j:j) = 'E'
#endif
       
        jj=jj+2
        if(nopnt.eq.0) jj=jj-1
        if(tword(jj:jj).eq.'-')then
          j=j+1
          word(j:j)='-'
        end if
        if(tword(jj+1:jj+1).ne.'0')then
          j=j+2
          word(j-1:j)=tword(jj+1:jj+2)
          go to 180
        end if
        if(tword(jj+2:jj+2).ne.'0')then
          j=j+1
          word(j:j)=tword(jj+2:jj+2)
        end if
180     j=j+1
        word(j:j)=tword(jj+3:jj+3)
        if(iflag.eq.1)then
          if(pos.eq.1)then
            jj=1
          else
            jj=2
          end if
          n=len_trim(word)
          do 190 j=jj,n-1
190       word(j:j)=word(j+1:j+1)
          word(n:n)=' '
        end if

200     if(len_trim(word).gt.nw)then
          ifail=-2
          return
        end if
        write(fmt,30) nw,0
        read(word,fmt,err=400) tval
        return
300     ifail=3
        return
400     ifail=-3
        return
500     ifail=-2
        return

end subroutine mio_wrtsig

end module model_input_output_interface
    
! define wrapper functions to call model_input_output_interface from C and C++
    subroutine mio_initialise_w(ifail,numin,numout,npar,nobs)
        use model_input_output_interface
        implicit none        
        integer, intent(out)                     :: ifail   ! indicates failure condition
        integer, intent(in)                      :: numin   ! number of model input files
        integer, intent(in)                      :: numout  ! number of model output files
        integer, intent(in)                      :: npar    ! number of parameters
        integer, intent(in)                      :: nobs    ! number of observations
        call mio_initialise(ifail,numin,numout,npar,nobs)
        return
    end subroutine mio_initialise_w
    
    subroutine mio_put_file_w(ifail,itype,inum,filename)
        use model_input_output_interface
        implicit none
        integer, intent(out)            :: ifail      ! indicates error condition
        integer, intent(in)             :: itype      ! type of file
        integer, intent(in)             :: inum       ! file number
        character (len=*), intent(in)   :: filename   ! name of file
        call mio_put_file(ifail,itype,inum,filename)
        return
    end subroutine mio_put_file_w
    
    subroutine mio_get_file_w(ifail,itype,inum,filename)
       use model_input_output_interface
       implicit none
       integer, intent(out)            :: ifail      ! indicates error condition
       integer, intent(in)             :: itype      ! type of file
       integer, intent(in)             :: inum       ! file number
       character (len=*), intent(out)   :: filename   ! name of file
       call mio_get_file(ifail,itype,inum,filename)
       return
    end subroutine mio_get_file_w
       
    subroutine mio_store_instruction_set_w(ifail)
        use model_input_output_interface
        implicit none
        integer, intent(out)          :: ifail
        call mio_store_instruction_set(ifail)
        return 
    end subroutine mio_store_instruction_set_w
        
    subroutine mio_process_template_files_w(ifail,npar,apar)
        use model_input_output_interface
        implicit none
        integer, intent(out)                :: ifail ! indicates error condition
        integer, intent(in)                 :: npar  ! number of parameters
        !character (len=*), dimension(npar)  :: apar  ! parameter names
        character*200  :: apar(npar)
        call mio_process_template_files(ifail,npar,apar)
        return
    end subroutine mio_process_template_files_w
        
    subroutine mio_delete_output_files_w(ifail,asldir)
        use model_input_output_interface
        implicit none
        integer, intent(out)                :: ifail     ! indicates error condition
        character(*), intent(in), optional  :: asldir    ! slave directory
        call mio_delete_output_files(ifail,asldir)
        return
    end subroutine mio_delete_output_files_w
        
    subroutine mio_write_model_input_files_w(ifail,npar,apar,pval)
        use model_input_output_interface
        implicit none
        integer, intent(out)                             :: ifail  ! error condition indicator
        integer, intent(in)                              :: npar   ! number of parameters
        !character (len=*), intent(in), dimension(npar)   :: apar   ! parameter names
        character*200  :: apar(npar)
        double precision, intent(inout), dimension(npar) :: pval   ! parameter values
        call mio_write_model_input_files(ifail,npar,apar,pval)
        return 
    end subroutine mio_write_model_input_files_w
    
    subroutine mio_read_model_output_files_w(ifail,nobs,aobs,obs)
        use model_input_output_interface
        implicit none
        integer, intent(out)                            :: ifail   ! error indicator
        integer, intent(in)                             :: nobs    ! number of observations
        !character (len=*), intent(in), dimension(nobs)  :: aobs    ! observation names
        character*200 :: aobs(nobs)
        double precision, intent(out), dimension(nobs)  :: obs     ! observation values
        !character (len=*), intent(out)                  :: instruction ! instruction of error
        character (len=20000)  :: instruction
        call mio_read_model_output_files(ifail,nobs,aobs,obs,instruction)
        return
    end subroutine mio_read_model_output_files_w
            
    subroutine mio_get_message_string_w(ifail,mess_len,amessage_out)
        use model_input_output_interface
        implicit none
        integer    :: mess_len
        integer, intent(out)       :: ifail
        integer :: i
        character(len=mess_len) :: amessage_out
        call mio_get_message_string(ifail,mess_len,amessage_out)
        return
    end subroutine mio_get_message_string_w
      
    subroutine mio_get_status_w(template_status_out,instruction_status_out)
        use model_input_output_interface
        implicit none
        integer, intent(out)  :: template_status_out,instruction_status_out
        call mio_get_status(template_status_out,instruction_status_out)
        return
    end subroutine mio_get_status_w
             
    subroutine mio_get_dimensions_w(numinfile_out,numoutfile_out)
        use model_input_output_interface
        implicit none
        integer, intent(out)  :: numinfile_out,numoutfile_out
        call mio_get_dimensions(numinfile_out,numoutfile_out)
        return
    end subroutine mio_get_dimensions_w
    
    subroutine mio_finalise_w(ifail)
        use model_input_output_interface
        implicit none
        integer, intent(out)                :: ifail
        call mio_finalise(ifail)
        return 
    end subroutine mio_finalise_w

! Notes:-

! -- As presently programmed, this module allows a parameter to be unrepresented in all model input
!      template files. However all observations must be cited in instruction files. This can be easily
!      altered of course.

