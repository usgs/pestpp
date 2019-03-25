subroutine getwrd(line,wrd,wrdnum)
!========================================================================================
!==== This subroutine returns the wrdnum-th trimmed character substring of a long    ====
!==== character string.                                                              ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------
  integer,intent(in)::wrdnum
  integer::i,j,k
  
  character(len=100),intent(in)::line
  character(len=100),intent(out)::wrd
!----------------------------------------------------------------------------------------

  wrd = ' '
  !
  i = 0
  j = 1
  k = 0
  !
  do
    !
    i = i + 1
    !
    if (line(i:i) /= ' ') then
      !
      if (j == wrdnum) then
        !
        do
          !
          k = k + 1
          !
          wrd(k:k) = line(i:i)
          !
          i = i + 1
          !
          if (line(i:i) == ' ') exit
          !
        end do
        !
      else
        !
        j = j + 1
        !
        do
          !
          i = i + 1
          !
          if (line(i:i) == ' ') exit
          !
        end do
        !        
      end if
      !
    end if
    !
    if (k > 0) exit
    !
  end do
  
end subroutine getwrd