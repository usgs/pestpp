subroutine instantrm(port)
!========================================================================================
!==== This subroutine instantiates the PEST++ run manager                            ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  external rmif_create_panther
  
  integer::rmif_create_panther,err
  
  character(len=20),intent(in)::port
  character(len=20)::storfile,rmi_info_file
!----------------------------------------------------------------------------------------

  err = 0
  !
  storfile      = 'tmp_run_data.bin    '
  rmi_info_file = 'run_manager_info.txt'
  !
  err = rmif_create_panther(storfile, 20,         &
                         port,     20,         &
                         rmi_info_file, 20,  1,&
                         1.00d+02, 1.00d+02, 1.00d+30)
  !
  if (err /= 0) then
    write(*,'(A,I0)')'PANTHER failed to instantiate --> err = ',err
    stop
  end if
  
end subroutine instantrm



subroutine initialrm(re)
!========================================================================================
!==== This subroutine initializes the PEST++ run manager                             ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  external rmif_initialize
  external rmif_reinitialize
  
  integer,intent(in)::re
  integer::rmif_initialize,rmif_reinitialize,err
!----------------------------------------------------------------------------------------

  err = 0
  !
  if (re == 0) then
    !
!   initialize run manager - allocate memory initialize parameter and observation names
    err = rmif_initialize(parnme,15,npar,obsnme,15,nobs)
    !
  else if (re == 1) then
    !
!   reinitialize run manager for a new set of runs
    err = rmif_reinitialize()   
    !
  end if
  !
  if (err /= 0) then
    write(*,'(A,I0)')'PANTHER failed to initialize --> err = ',err
    stop
  end if
  
end subroutine initialrm



subroutine modelrm(sr,irun,fail)
!========================================================================================
!==== This subroutine sends and receives runs from the PEST++ run manager            ====
!====    by Adam Siade                                                               ====
!========================================================================================
!========================================================================================

  use psodat
  
  implicit none
  
! specifications:
!----------------------------------------------------------------------------------------  
  external rmif_add_run
  external rmif_get_run
  
  integer,intent(in)::sr
  integer,intent(inout)::irun,fail
  integer::rmif_add_run,rmif_get_run,err
!----------------------------------------------------------------------------------------

  err  = 0
  fail = 0
  !
  if (sr == 1) then
    !
    err = rmif_add_run(parval,npar,irun)
    !
    modeval = modeval + 1
    !
    if (err /= 0) then
      write(*,'(A,I0)')'Model run was not added to the queue properly --> err = ',err
      stop
    end if
    !
  else if (sr == 0) then
    !
    err = rmif_get_run(irun,parval,npar,mobsval,nobs)
    !
    if (err /= 0) then
      fail = 1
    end if
    !
  end if
  
end subroutine modelrm
      

































