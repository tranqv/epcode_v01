!-----------------------------------------------------------------------
!
      program test_ep_gstate_qp_real16 
!
      use epcode_mod 
      implicit none 
!
!!    integer, parameter :: rk = 8 
      integer, parameter :: rk = 16 
!
!     LOCAL VARIABLE:
!
      real(rk), allocatable :: e(:), g(:,:)
      real(rk) :: tole 
      integer  :: nome, npar, ng 
      integer(8) :: neiv
      character :: prefix*200
      character :: fname*90, fexec*90
      real(8) :: dt  
      integer :: istat  
!
!
!     Default configuration file: 
!!    fname = 'epconf_default.txt'
      fname = 'ep_conf.txt'
!
      if ( command_argument_count() < 1 ) then 
         call get_command_argument( 0, fexec )
         write(*,'(A)') "USAGE: "//trim(fexec)//" epconf.txt"
      else 
         call get_command_argument( 1, fname )
      endif 
!
      write(*,'(/,A)') "Configuration file: "// trim(fname) 
!
      call ep_readconf ( & 
           fname, nome, npar, ng, g, e, neiv, tole, prefix, istat )
!
      if ( istat .ne. 0 ) stop 
!
      call evaltime ( 0, dt )
      call ep_gstate ( &
           nome, npar, e, ng, g, neiv, prefix, tole ) 
      call evaltime ( 1, dt )
!
      write(*,'(/,a,1pe14.7)') "CPU-time (sec):", dt 
!
      if ( allocated(e) ) deallocate ( e )
      if ( allocated(g) ) deallocate ( g ) 
!
      STOP 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      contains 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine evaltime ( istat, dt )
      implicit none 
!
!     ARGUMENTS: 
      integer :: istat 
      real(kind=8) :: dt 
!
!     DEPENDENCES: 
      intrinsic :: system_clock, dble
!
!     LOCAL VARIABLES: 
      integer(kind=8) :: clock_rate, clock_max, t  
!
!     EXECUTABLE STATEMENTS: 
!
      if ( istat .eq. 0 ) then 
         call system_clock ( t, clock_rate, clock_max )
         dt = dble(t) 
      else 
         call system_clock ( t, clock_rate, clock_max )
         dt =  ( dble(t) - dt )/ dble(clock_rate)
      endif 
      return 
      end subroutine evaltime
!***********************************************************************
      end program test_ep_gstate_qp_real16 
!***********************************************************************
