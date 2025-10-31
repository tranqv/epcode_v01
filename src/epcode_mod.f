!-----------------------------------------------------------------------
      module epcode_mod  
!
      implicit none
      private 
!-----------------------------------------------------------------------
      public :: ep_gstate
      interface ep_gstate
!-----
!     real(kind=8):
!
      subroutine ep_gstate_1a_dp ( 
     $           nome, npar, e, g, neiv, prefix, tole ) 
      implicit none 
      include 'epcode_inc_dp.h'
      integer, intent(in) :: nome, npar 
      real(wrp), intent(in) :: e(nome), g, tole
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_1a_dp 
!
      subroutine ep_gstate_1b_dp ( 
     $           nome, npar, e, g, neiv, prefix ) 
      implicit none 
      include 'epcode_inc_dp.h'
      integer, intent(in) :: nome, npar 
      real(wrp), intent(in) :: e(nome), g
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_1b_dp 
!
      subroutine ep_gstate_2a_dp ( 
     $           nome, npar, e, g, neiv, prefix, tole ) 
      implicit none 
      include 'epcode_inc_dp.h'
      integer, intent(in) :: nome, npar 
      real(wrp), intent(in) :: e(nome), g(nome,nome), tole
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_2a_dp 
!
      subroutine ep_gstate_2b_dp ( 
     $           nome, npar, e, g, neiv, prefix ) 
      implicit none 
      include 'epcode_inc_dp.h'
      integer, intent(in) :: nome, npar 
      real(wrp), intent(in) :: e(nome), g(nome,nome)
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_2b_dp 
!
      subroutine ep_gstate_3a_dp ( 
     $           nome, npar, e, ng, g, neiv, prefix, tole ) 
      implicit none 
      include 'epcode_inc_dp.h'
      integer, intent(in) :: nome, npar, ng 
      real(wrp), intent(in) :: e(nome), g(ng,ng), tole
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_3a_dp 
!
      subroutine ep_gstate_3b_dp ( 
     $           nome, npar, e, ng, g, neiv, prefix ) 
      implicit none 
      include 'epcode_inc_dp.h'
      integer, intent(in) :: nome, npar, ng
      real(wrp), intent(in) :: e(nome), g(ng,ng)
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_3b_dp 
!-----
!     real(kind=16):
!
      subroutine ep_gstate_1a_qp ( 
     $           nome, npar, e, g, neiv, prefix, tole ) 
      implicit none 
      include 'epcode_inc_qp.h'
      integer, intent(in) :: nome, npar 
      real(wrp), intent(in) :: e(nome), g, tole
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_1a_qp 
!
      subroutine ep_gstate_1b_qp ( 
     $           nome, npar, e, g, neiv, prefix ) 
      implicit none 
      include 'epcode_inc_qp.h'
      integer, intent(in) :: nome, npar 
      real(wrp), intent(in) :: e(nome), g
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_1b_qp 
!
      subroutine ep_gstate_2a_qp ( 
     $           nome, npar, e, g, neiv, prefix, tole ) 
      implicit none 
      include 'epcode_inc_qp.h'
      integer, intent(in) :: nome, npar 
      real(wrp), intent(in) :: e(nome), g(nome,nome), tole
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_2a_qp 
!
      subroutine ep_gstate_2b_qp ( 
     $           nome, npar, e, g, neiv, prefix ) 
      implicit none 
      include 'epcode_inc_qp.h'
      integer, intent(in) :: nome, npar 
      real(wrp), intent(in) :: e(nome), g(nome,nome)
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_2b_qp 
!
      subroutine ep_gstate_3a_qp ( 
     $           nome, npar, e, ng, g, neiv, prefix, tole ) 
      implicit none 
      include 'epcode_inc_qp.h'
      integer, intent(in) :: nome, npar, ng 
      real(wrp), intent(in) :: e(nome), g(ng,ng), tole
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_3a_qp 
!
      subroutine ep_gstate_3b_qp ( 
     $           nome, npar, e, ng, g, neiv, prefix ) 
      implicit none 
      include 'epcode_inc_qp.h'
      integer, intent(in) :: nome, npar, ng
      real(wrp), intent(in) :: e(nome), g(ng,ng)
      integer(wip), intent(in) :: neiv 
      character(len=*), intent(inout) :: prefix
      end subroutine ep_gstate_3b_qp 
!-----
      end interface
!-----------------------------------------------------------------------
      public :: ep_readconf
      interface ep_readconf
         module procedure ep_readconf_dp
         module procedure ep_readconf_qp
      end interface 
!-----------------------------------------------------------------------
      contains 
!-----------------------------------------------------------------------
      subroutine ep_readconf_dp ( 
     $           fname, nome, npar, ng, g, e, neiv, tole,  
     $           prefix, istat )
      implicit none
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      character(len=*), intent(in) :: fname
      integer, intent(out) :: nome, npar, ng 
      real(wrp), allocatable, intent(out) :: e(:), g(:,:)
      integer(wip), intent(out) :: neiv
      real(wrp), intent(out) :: tole 
      integer, intent(out) :: istat
      character(len=*), intent(out) :: prefix
!
!     LOCAL VARIABLES:
      integer,parameter :: icha_bs  = 32           ! White Space
      integer,parameter :: icha_tab = 9            ! Horizontal Tab
      integer,parameter :: icha_c1  = 33           ! '!' 
      integer,parameter :: icha_c2  = 35           ! '#'
      integer,parameter :: icha_eq  = iachar('=')  ! '=' 
      integer :: 
     $   i, j, k, l, io, n, ios, i0, j0 
      character :: 
     $   charg*550, chnam*50, chval*500, chn*3, chfmt*50, ch*1 
      real(wrp) :: tmp1 
!
!     DEPENDENCES: 
      intrinsic :: iachar, len_trim, trim, adjustl 
!
!     SHARED PARAMETERS:
      logical 
     $   flag_prep, 
     $   flag_wrtinp 
      common /epcode_shared_flags/ 
     $   flag_prep, 
     $   flag_wrtinp 
!
!     EXECUTABLE STATEMENTS: 
!
!     Read for STATES firstly:
      call ep_readconf_a_dp ( 
     $     fname, nome, ng, g, e,  
     $     istat )
!
!     Then, read for other infomation
      istat = 0 
      open(newunit=io, file=fname, status='old', iostat=istat)
      if ( istat .ne. 0 ) then
         write(*,'(t15,a)') 'ERROR: could not open file '//fname 
         return
      endif
      prefix = 'ep_dp'
      neiv = 1 
      tole = 0 
      n = 0
      do
         charg = '' 
         read(unit=io, fmt='(a)', end=101, err=102) charg
         n = n + 1 
         charg = adjustl(charg)
         l = len_trim(charg)
         if ( l .eq. 0 ) cycle
         if ( charg(1:l) .eq. '__EOF__' ) goto 101 
         j = iachar(charg(1:1))
         if (j .eq. icha_c1 .or. j .eq. icha_c2) cycle 
         i = 0  
         do k = 1, l 
            if ( iachar(charg(k:k)) .ne. icha_eq ) cycle 
            i = 1  
            exit 
         enddo
         if ( i .eq. 0 ) cycle
         chnam = charg(:k-1)
         chval = charg(k+1:l)
         select case ( uppercase(trim(adjustl(chnam))) )
!!       case('NOME','nome','Nome')
         case('NOME')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = '' 
               endif 
            enddo   
            read(chval, fmt=*, err=102) nome 
            if ( flag_wrtinp )
     $         write(*,10) n, 'NOME =', nome
            if ( allocated(e) ) deallocate(e)
            allocate ( e(nome) )
            do k = 1, nome 
               e(k) = k
            enddo 
!!          write(*,13) 'Initialize', 'E(:) =', e
!!       case('NPAR','npar','Npar')
         case('NPAR')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo 
            read(chval, fmt=*, err=102) npar 
            if ( flag_wrtinp )
     $         write(*,10) n, 'NPAR =', npar 
!!       case('G','g')
         case('G')
            chval = trim(adjustl(chval))
            if ( len_trim(chval) .ne. 0 )  then 
!!             g = constant 
               ng = 1 
               k = 1
               if ( flag_wrtinp )
     $            write(*,10) n, 'NG =', ng  
               if ( allocated(g) ) deallocate(g) 
               allocate( g(ng,ng) )          
               read(chval, fmt=*, err=102) g
               if ( flag_wrtinp )
     $            write(*,11) n, 'G =', g(1,1)
            else 
!!             g = matrix  
               ng = nome 
               if ( flag_wrtinp )
     $            write(*,10) n, 'NG =', ng  
               if ( allocated(g) ) deallocate(g) 
               allocate( g(ng,ng) )          
               do j = 1, ng  
               do k = j, ng   
                  do 
                     read(unit=io, fmt='(a)', end=101, err=102) 
     $                  chval 
                     n = n + 1  
                     if ( len_trim(chval) .ne. 0 ) then 
                        chval = adjustl(chval)
                        j0 = iachar(chval(1:1))
                        if ( j0 .eq. icha_c1 .or. 
     $                       j0 .eq. icha_c2 ) chval = ''  
                     endif 
                     if ( len_trim(chval) .ne. 0 ) exit  
                  enddo
                  chval = trim(adjustl(chval))
                  if ( flag_wrtinp )
     $               write(*,12) 
     $               n, 'j, k, G(j,k):', '' 
                  read(chval, fmt=*, err=102) i0, j0, tmp1 
                  g(i0,j0) = tmp1 
                  g(j0,i0) = tmp1 
                  if ( flag_wrtinp )
     $               write(*,*) i0, j0, tmp1  
               enddo
               enddo 
               if ( flag_wrtinp ) then 
                  chn = '' 
                  write(chn,fmt='(I3)') nome 
                  chfmt = '(I5, ":",2x,'//
     $               trim(adjustl(chn))//'(1X,ES9.2))'
                  write(*,12) n, 'G(:,:) =', '' 
                  do k = 1, ng 
                     write(*,fmt=chfmt) k, g(k,:)
                  enddo 
               endif 
            endif 
!!       case('PREFIX','prefix','Prefix')
         case('PREFIX')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1 
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo 
            prefix = trim(adjustl(chval))
            if ( flag_wrtinp ) 
     $         write(*,12) n, ' PREFIX =', trim(prefix)
!!       case('MAXNEIV','maxneiv','Maxneiv','NEIV','neiv','Neiv')
         case('MAXNEIV','NEIV')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if ( j .eq. icha_c1 .or. 
     $                 j .eq. icha_c2) chval = ''  
               endif 
            enddo 
            read(chval, fmt=*, err=102) neiv
            if ( flag_wrtinp )
     $         write(*,10) n, 'NEIV =', neiv
!!       case('E','e')
         case('E')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1 
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo
            if ( allocated(e) ) deallocate(e)
            allocate ( e(nome) )
            read(chval, fmt=*, iostat=k) e
            if ( k .eq. 0 ) then
               if ( flag_wrtinp )
     $            write(*,12) 
     $            n, 'Reading E(:)', '' 
            else  
               backspace(unit=io) 
               do j = 1, nome 
                  do 
                     read(unit=io, fmt='(a)', end=101, err=102) 
     $                  chval 
                     n = n + 1  
                     if ( len_trim(chval) .ne. 0 ) then 
                        chval = adjustl(chval)
                        j0 = iachar(chval(1:1))
                        if ( j0 .eq. icha_c1 .or. 
     $                       j0 .eq. icha_c2 ) chval = ''  
                     endif 
                     if ( len_trim(chval) .ne. 0 ) exit  
                  enddo
                  chval = trim(adjustl(chval))
                  if ( flag_wrtinp )
     $               write(*,12) 
     $               n, 'j, E(j):', '' 
                  read(chval, fmt=*, err=102) tmp1 
                  e(j) = tmp1 
                  if ( flag_wrtinp )
     $               write(*,*) j, tmp1  
               enddo 
            endif
            if ( flag_wrtinp )
     $         write(*,15) n, 'E(:) =', e
!!       case('TOLE','tole','Tole','TOL','tol','Tol')
         case('TOLE')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo
            read(chval, fmt=*, err=102) tole
            if ( flag_wrtinp )
     $         write(*,11) n, 'TOLE =', tole
!!       case('PREP','prep','Prep')
         case('PREP')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo
            read(chval, fmt=*, err=102) flag_prep
            if ( flag_wrtinp )
     $         write(*,16) n, 'flag_prep =', flag_prep
         case default 
            cycle 
         end select          
      enddo    
 101  continue
      if ( flag_wrtinp )
     $   write(*,14) n, 'Finish reading the configuration file.'
      close(io)
      return
 102  continue
      istat = 100  
      write(*,12) 
     $   n, 'ERROR as reading record.', 
     $   'The configuration process is terminated!'
      return  
 10   format( '>>> line ',i3,": ", 1x, a30, 1x, i5 ) 
 11   format( '>>> line ',i3,": ", 1x, a30, 1x, e14.7 ) 
 12   format( '>>> line ',i3,": ", 1x, a30, 1x, a ) 
 13   format(  14x,                1x, a30, /, 5(3x,f7.3) ) 
 14   format( '>>> line ',i3,": ", a, // ) 
 15   format( '>>> line ',i3,": ", 1x, a30, /, 5(3x,f7.3) ) 
 16   format( '>>> line ',i3,": ", 1x, a30, 1x, L ) 
 17   format( '>>> line ',i3,": ", 1x, a30, 1x, 2(i5,1x,i5) ) 
      end subroutine ep_readconf_dp
!-----------------------------------------------------------------------
!
!     Read the setting parameters defined by the colon ':'.
!     +  if ng=1, skip reading G and use g=constant defined 
!        out of this scope.
!
      subroutine ep_readconf_a_dp ( 
     $           fname, nome, ng, g, e,  
     $           istat )
      implicit none
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      character(len=*), intent(in) :: fname
      integer, intent(inout) :: ng 
      integer, intent(out) :: nome 
      real(wrp), allocatable, intent(out) :: e(:), g(:,:)
      integer, intent(out) :: istat
!
!     LOCAL VARIABLES:
      integer,parameter :: icha_bs  = 32              ! White Space
      character,parameter :: char_bs= achar(icha_bs)  ! White Space
      integer,parameter :: icha_tab = 9               ! Horizontal Tab
      integer,parameter :: icha_c1  = 33              ! '!' 
      integer,parameter :: icha_c2  = 35              ! '#'
      integer,parameter :: icha_eq  = iachar('=')     ! '=' 
      integer,parameter :: icha_colon = iachar(':')   ! ':' 
      integer,parameter :: icha_0   = iachar('0')     ! '0'
      integer,parameter :: icha_9   = iachar('9')     ! '9'
      integer :: 
     $   i, j, k, l, io, n, i0, i1, i2, j0, j1, j2, 
     $   k0, i3, i4 
      character :: 
     $   charg*550, chnam*50, chval*500, chn*3, chfmt*50, ch*1
      character(len=:), allocatable :: ch1, ch2
      integer :: nsta
      integer, allocatable :: ista(:) 
      real(wrp), allocatable :: esta(:), vsta(:,:)
      real(wrp) :: rtmp 
!
!     DEPENDENCES: 
      intrinsic :: iachar, len_trim, trim, adjustl, real  
!
!     SHARED PARAMETERS:
      logical 
     $   flag_prep, 
     $   flag_wrtinp 
      common /epcode_shared_flags/ 
     $   flag_prep, 
     $   flag_wrtinp 
!
!     EXECUTABLE STATEMENTS: 
!
      istat = 0 
      open(newunit=io, file=fname, status='old', iostat=istat)
      if ( istat .ne. 0 ) then
         write(*,'(t15,a)') 'ERROR: could not open file '//fname 
         return
      endif
      n = 0
      do
         charg = '' 
         read(unit=io, fmt='(a)', end=101, err=102) charg
         n = n + 1 
         charg = adjustl(charg)
         l = len_trim(charg)
         if ( l .eq. 0 ) cycle
         if ( charg(1:l) .eq. '__EOF__' ) goto 101 
         j = iachar(charg(1:1))
         if (j .eq. icha_c1 .or. j .eq. icha_c2) cycle 
         i = 0  
         do k = 1, l 
            if ( iachar(charg(k:k)) .ne. icha_colon ) cycle 
            i = 1  
            exit 
         enddo
         if ( i .eq. 0 ) cycle
         chnam = charg(:k-1)
         chval = trim(adjustl(charg(k+1:l)))
         select case ( uppercase(trim(adjustl(chnam))) )
!!       case('STATES','states','States')
         case('STATES')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo
            if ( flag_wrtinp )
     $         write(*,12) n, 'States:', trim(chval) 
            chval = trim(adjustl(chval)) // char_bs
            nsta = 0 
            do i = 1, len(chval) - 1 
               if ( iachar(chval(i:i)) .eq. icha_eq ) 
     $            nsta = nsta + 1 
            enddo   
            if ( nsta .eq. 0 ) cycle 
            if ( flag_wrtinp )
     $         write(*,10) n, 'N_STATES =', nsta
            if ( allocated(ista) ) deallocate(ista) 
            allocate ( ista(nsta) ) 
            i1 = 1
            i2 = 0 
            k = len(chval)
            i = i2 
            do while ( i .le. k )
               i = i + 1  
               if ( iachar(chval(i:i)) .ne. icha_bs .and. 
     $              iachar(chval(i:i)) .ne. icha_tab ) cycle
               i2 = i 
               ch1 = chval(i1:i2-1) 
               do while ( iachar(chval(i2:i2)) .eq. icha_bs .or.
     $                    iachar(chval(i2:i2)) .eq. icha_tab ) 
                  i2 = i2 + 1 
               enddo  
               ch2 = chval(i2:k)
               i1 = i2
               i = i2 
               do i3 = 1, len(ch1)
                  i4 = iachar(ch1(i3:i3))
                  if ( i4 .lt. icha_0 .or. i4 .gt. icha_9 ) 
     $               ch1(i3:i3) = char_bs
               enddo 
               read( ch1, * ) i0, j0, k0   
               ista(i0) = (j0 + 1)/2
            enddo 
            nome = 0 
            do i0 = 1, nsta
               nome = nome + ista(i0) 
            enddo
            if ( flag_wrtinp ) then 
               chn = '' 
               write(chn,fmt='(I3)') nsta 
               chfmt = '('//trim(adjustl(chn))//'(I5,1x))'
               write(charg,fmt=chfmt) ista 
               write(*,12) n, 'I_States =', trim(charg) 
               write(*,10) n, 'NOME =', nome
            endif 
!!       case('V0','v0')
         case('V0')
            if ( ng .eq. 1 ) cycle 
            ng = nome  
            if ( flag_wrtinp )
     $         write(*,10) n, 'NG =', ng  
            if ( allocated(vsta) ) deallocate(vsta)
            allocate ( vsta(nsta,nsta) )
            if ( allocated(g) ) deallocate(g) 
            allocate( g(ng,ng) )         
            do j = 1, nsta
            do k = j, nsta 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1   
               do while ( len_trim(chval) .eq. 0 ) 
                  read(unit=io, fmt='(a)', end=101, err=102) chval 
                  n = n + 1  
                  if ( len_trim(chval) .ne. 0 ) then 
                     chval = adjustl(chval)
                     j0 = iachar(chval(1:1))
                     if ( j0 .eq. icha_c1 .or. 
     $                    j0 .eq. icha_c2 ) chval = ''  
                  endif 
               enddo
               chval = trim(adjustl(chval))
               if ( flag_wrtinp )
     $            write(*,12) 
     $            n, 'j, k, V0_states(j,k):', '' 
               read(chval, fmt=*, err=102) i0, j0, rtmp 
               vsta(i0,j0) = rtmp 
               vsta(j0,i0) = rtmp 
               if ( flag_wrtinp )
     $            write(*,*) i0, j0, rtmp  
            enddo
            enddo 
            if ( flag_wrtinp ) then 
               write(*,12) n, 'V0_states(:,:) =', '' 
               chn = '' 
               write(chn,fmt='(I3)') nsta 
               chfmt = '(I5, ":",2x,'//trim(adjustl(chn))//'(1X,ES9.2))'
               do k = 1, nsta
                  write(*,fmt=chfmt) k, vsta(k,:)
               enddo
           endif 
            do j0 = 1, nsta
               do i0 = 1, nsta 
                  vsta(i0,j0) = 2 * vsta(i0,j0) / 
     $            sqrt( real(2*ista(i0)*2*ista(j0),kind=wrp) ) 
               enddo 
            enddo 
            i2 = 0 
            do i0 = 1, nsta
               i1 = ista(i0)
               j2 = 0 
               do j0 = 1, nsta 
                  j1 = ista(j0)
                  g( i2+1:i2+i1, j2+1:j2+j1 ) = vsta(i0,j0)
                  j2 = j2+j1 
               enddo 
               i2 = i2+i1 
            enddo
            if ( flag_wrtinp ) then 
               chn = '' 
               write(chn,fmt='(I3)') nome 
               chfmt = '(I5, ":",2x,'//
     $            trim(adjustl(chn))//'(1X,ES9.2))'
               write(*,12) n, 'G(:,:) =', '' 
               do k = 1, ng 
                  write(*,fmt=chfmt) k, g(k,:)
               enddo 
            endif 
!!       case('G','g')
         case('G')
            if ( ng .eq. 1 ) cycle 
            ng = nome 
            if ( flag_wrtinp )
     $         write(*,10) n, 'NG =', ng  
            if ( allocated(vsta) ) deallocate(vsta)
            allocate ( vsta(nsta,nsta) )
            if ( allocated(g) ) deallocate(g) 
            allocate( g(ng,ng) )         
            do j = 1, nsta
            do k = j, nsta 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1   
               do while ( len_trim(chval) .eq. 0 ) 
                  read(unit=io, fmt='(a)', end=101, err=102) chval 
                  n = n + 1  
                  if ( len_trim(chval) .ne. 0 ) then 
                     chval = adjustl(chval)
                     j0 = iachar(chval(1:1))
                     if ( j0 .eq. icha_c1 .or. 
     $                    j0 .eq. icha_c2 ) chval = ''  
                  endif 
               enddo
               chval = trim(adjustl(chval))
               if ( flag_wrtinp )
     $            write(*,12) 
     $            n, 'j, k, G_states(j,k):', '' 
               read(chval, fmt=*, err=102) i0, j0, rtmp 
               vsta(i0,j0) = rtmp 
               vsta(j0,i0) = rtmp 
               if ( flag_wrtinp )
     $            write(*,*) i0, j0, rtmp  
            enddo
            enddo 
            if ( flag_wrtinp ) then 
               write(*,12) n, 'G_states(:,:) =', ''
               chn = '' 
               write(chn,fmt='(I3)') nsta 
               chfmt = '(I5, ":",2x,'//
     $            trim(adjustl(chn))//'(1X,ES9.2))'
               do k = 1, nsta
                  write(*,fmt=chfmt) k, vsta(k,:)
               enddo    
            endif 
            i2 = 0 
            do i0 = 1, nsta
               i1 = ista(i0)
               j2 = 0 
               do j0 = 1, nsta 
                  j1 = ista(j0)
                  g( i2+1:i2+i1, j2+1:j2+j1 ) = vsta(i0,j0)
                  j2 = j2+j1 
               enddo 
               i2 = i2+i1 
            enddo 
            if ( flag_wrtinp ) then 
               chn = '' 
               write(chn,fmt='(I3)') nome 
               chfmt = '(I5, ":",2x,'//
     $            trim(adjustl(chn))//'(1X,ES9.2))'
               write(*,12) n, 'G(:,:) =', '' 
               do k = 1, ng 
                  write(*,fmt=chfmt) k, g(k,:)
               enddo 
            endif 
!!       case('E','e')
         case('E')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1 
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if ( j .eq. icha_c1 .or. 
     $                 j .eq. icha_c2) chval = ''  
               endif 
            enddo
            if ( nsta .eq. 0 ) cycle 
            if ( allocated(esta) ) deallocate(esta)
            allocate ( esta(nsta) )
            read(chval, fmt=*, iostat=i1) esta
            if ( i1 .ne. 0 ) then 
               backspace(unit=io) 
               do j = 1, nsta
                  read(unit=io, fmt=*, err=102) esta(j)
               enddo 
            endif
            if ( flag_wrtinp ) 
     $         write(*,15) n, 'E_states(:) =', esta
            if ( allocated(e) ) deallocate(e)
            allocate ( e(nome) )
            i2 = 0 
            do i0 = 1, nsta
               i1 = ista(i0)
               e(i2+1:i2+i1) = esta(i0)
               i2 = i2+i1 
            enddo 
            if ( flag_wrtinp )
     $         write(*,15) n, 'E(:) =', e
         case default 
            cycle 
         end select          
      enddo    
 101  continue
      if ( flag_wrtinp ) 
     $   write(*,14) 
     $   n, 'Finish reading the configuration file for STATES.'
      close(io)
      return
 102  continue
      close(io)
      istat = 100  
      write(*,12) 
     $   n, 'ERROR as reading record.', 
     $   'The configuration process is terminated!'
      return  
 10   format( '>>> line ',i3,": ", 1x, a30, 1x, i5 ) 
 11   format( '>>> line ',i3,": ", 1x, a30, 1x, e14.7 ) 
 12   format( '>>> line ',i3,": ", 1x, a30, 1x, a ) 
 13   format(  14x,                1x, a30, /, 5(3x,f7.3) ) 
 14   format( '>>> line ',i3,": ", a, // ) 
 15   format( '>>> line ',i3,": ", 1x, a30, /, 5(3x,f7.3) ) 
      end subroutine ep_readconf_a_dp
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine ep_readconf_qp ( 
     $           fname, nome, npar, ng, g, e, neiv, tole,  
     $           prefix, istat )
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      character(len=*), intent(in) :: fname
      integer, intent(out) :: nome, npar, ng 
      real(wrp), allocatable, intent(out) :: e(:), g(:,:)
      integer(wip), intent(out) :: neiv
      real(wrp), intent(out) :: tole 
      integer, intent(out) :: istat
      character(len=*), intent(out) :: prefix
!
!     LOCAL VARIABLES:
      integer,parameter :: icha_bs  = 32           ! White Space
      integer,parameter :: icha_tab = 9            ! Horizontal Tab
      integer,parameter :: icha_c1  = 33           ! '!' 
      integer,parameter :: icha_c2  = 35           ! '#'
      integer,parameter :: icha_eq  = iachar('=')  ! '=' 
      integer :: 
     $   i, j, k, l, io, n, ios, i0, j0 
      character :: 
     $   charg*550, chnam*50, chval*500, chn*3, chfmt*50, ch*1 
      real(wrp) :: tmp1 
!
!     DEPENDENCES: 
      intrinsic :: iachar, len_trim, trim, adjustl 
!
!     SHARED PARAMETERS:
      logical 
     $   flag_prep, 
     $   flag_wrtinp 
      common /epcode_shared_flags/ 
     $   flag_prep, 
     $   flag_wrtinp 
!
!     EXECUTABLE STATEMENTS: 
!
!     Read for STATES firstly:
      call ep_readconf_a_qp ( 
     $     fname, nome, ng, g, e,  
     $     istat )
!
!     Then, read for other infomation
      istat = 0 
      open(newunit=io, file=fname, status='old', iostat=istat)
      if ( istat .ne. 0 ) then
         write(*,'(t15,a)') 'ERROR: could not open file '//fname 
         return
      endif
      prefix = 'ep_dp'
      neiv = 1 
      tole = 0 
      n = 0
      do
         charg = '' 
         read(unit=io, fmt='(a)', end=101, err=102) charg
         n = n + 1 
         charg = adjustl(charg)
         l = len_trim(charg)
         if ( l .eq. 0 ) cycle
         if ( charg(1:l) .eq. '__EOF__' ) goto 101 
         j = iachar(charg(1:1))
         if (j .eq. icha_c1 .or. j .eq. icha_c2) cycle 
         i = 0  
         do k = 1, l 
            if ( iachar(charg(k:k)) .ne. icha_eq ) cycle 
            i = 1  
            exit 
         enddo
         if ( i .eq. 0 ) cycle
         chnam = charg(:k-1)
         chval = charg(k+1:l)
         select case ( uppercase(trim(adjustl(chnam))) )
!!       case('NOME','nome','Nome')
         case('NOME')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = '' 
               endif 
            enddo   
            read(chval, fmt=*, err=102) nome 
            if ( flag_wrtinp )
     $         write(*,10) n, 'NOME =', nome
            if ( allocated(e) ) deallocate(e)
            allocate ( e(nome) )
            do k = 1, nome 
               e(k) = k
            enddo 
!!          write(*,13) 'Initialize', 'E(:) =', e
!!       case('NPAR','npar','Npar')
         case('NPAR')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo 
            read(chval, fmt=*, err=102) npar 
            if ( flag_wrtinp )
     $         write(*,10) n, 'NPAR =', npar 
!!       case('G','g')
         case('G')
            chval = trim(adjustl(chval))
            if ( len_trim(chval) .ne. 0 )  then 
!!             g = constant 
               ng = 1 
               k = 1
               if ( flag_wrtinp )
     $            write(*,10) n, 'NG =', ng  
               if ( allocated(g) ) deallocate(g) 
               allocate( g(ng,ng) )          
               read(chval, fmt=*, err=102) g
               if ( flag_wrtinp )
     $            write(*,11) n, 'G =', g(1,1)
            else 
!!             g = matrix  
               ng = nome 
               if ( flag_wrtinp )
     $            write(*,10) n, 'NG =', ng  
               if ( allocated(g) ) deallocate(g) 
               allocate( g(ng,ng) )          
               do j = 1, ng  
               do k = j, ng   
                  do 
                     read(unit=io, fmt='(a)', end=101, err=102) 
     $                  chval 
                     n = n + 1  
                     if ( len_trim(chval) .ne. 0 ) then 
                        chval = adjustl(chval)
                        j0 = iachar(chval(1:1))
                        if ( j0 .eq. icha_c1 .or. 
     $                       j0 .eq. icha_c2 ) chval = ''  
                     endif 
                     if ( len_trim(chval) .ne. 0 ) exit  
                  enddo
                  chval = trim(adjustl(chval))
                  if ( flag_wrtinp )
     $               write(*,12) 
     $               n, 'j, k, G(j,k):', '' 
                  read(chval, fmt=*, err=102) i0, j0, tmp1 
                  g(i0,j0) = tmp1 
                  g(j0,i0) = tmp1 
                  if ( flag_wrtinp )
     $               write(*,*) i0, j0, tmp1  
               enddo
               enddo 
               if ( flag_wrtinp ) then 
                  chn = '' 
                  write(chn,fmt='(I3)') nome 
                  chfmt = '(I5, ":",2x,'//
     $               trim(adjustl(chn))//'(1X,ES9.2))'
                  write(*,12) n, 'G(:,:) =', '' 
                  do k = 1, ng 
                     write(*,fmt=chfmt) k, g(k,:)
                  enddo 
               endif 
            endif 
!!       case('PREFIX','prefix','Prefix')
         case('PREFIX')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1 
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo 
            prefix = trim(adjustl(chval))
            if ( flag_wrtinp ) 
     $         write(*,12) n, ' PREFIX =', trim(prefix)
!!       case('MAXNEIV','maxneiv','Maxneiv','NEIV','neiv','Neiv')
         case('MAXNEIV','NEIV')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if ( j .eq. icha_c1 .or. 
     $                 j .eq. icha_c2) chval = ''  
               endif 
            enddo 
            read(chval, fmt=*, err=102) neiv
            if ( flag_wrtinp )
     $         write(*,10) n, 'NEIV =', neiv
!!       case('E','e')
         case('E')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1 
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo
            if ( allocated(e) ) deallocate(e)
            allocate ( e(nome) )
            read(chval, fmt=*, iostat=k) e
            if ( k .eq. 0 ) then
               if ( flag_wrtinp )
     $            write(*,12) 
     $            n, 'Reading E(:)', '' 
            else  
               backspace(unit=io) 
               do j = 1, nome 
                  do 
                     read(unit=io, fmt='(a)', end=101, err=102) 
     $                  chval 
                     n = n + 1  
                     if ( len_trim(chval) .ne. 0 ) then 
                        chval = adjustl(chval)
                        j0 = iachar(chval(1:1))
                        if ( j0 .eq. icha_c1 .or. 
     $                       j0 .eq. icha_c2 ) chval = ''  
                     endif 
                     if ( len_trim(chval) .ne. 0 ) exit  
                  enddo
                  chval = trim(adjustl(chval))
                  if ( flag_wrtinp )
     $               write(*,12) 
     $               n, 'j, E(j):', '' 
                  read(chval, fmt=*, err=102) tmp1 
                  e(j) = tmp1 
                  if ( flag_wrtinp )
     $               write(*,*) j, tmp1  
               enddo 
            endif
            if ( flag_wrtinp )
     $         write(*,15) n, 'E(:) =', e
!!       case('TOLE','tole','Tole','TOL','tol','Tol')
         case('TOLE')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo
            read(chval, fmt=*, err=102) tole
            if ( flag_wrtinp )
     $         write(*,11) n, 'TOLE =', tole
!!       case('PREP','prep','Prep')
         case('PREP')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo
            read(chval, fmt=*, err=102) flag_prep
            if ( flag_wrtinp )
     $         write(*,16) n, 'flag_prep =', flag_prep
         case default 
            cycle 
         end select          
      enddo    
 101  continue
      if ( flag_wrtinp )
     $   write(*,14) n, 'Finish reading the configuration file.'
      close(io)
      return
 102  continue
      istat = 100  
      write(*,12) 
     $   n, 'ERROR as reading record.', 
     $   'The configuration process is terminated!'
      return  
 10   format( '>>> line ',i3,": ", 1x, a30, 1x, i5 ) 
 11   format( '>>> line ',i3,": ", 1x, a30, 1x, e14.7 ) 
 12   format( '>>> line ',i3,": ", 1x, a30, 1x, a ) 
 13   format(  14x,                1x, a30, /, 5(3x,f7.3) ) 
 14   format( '>>> line ',i3,": ", a, // ) 
 15   format( '>>> line ',i3,": ", 1x, a30, /, 5(3x,f7.3) ) 
 16   format( '>>> line ',i3,": ", 1x, a30, 1x, L ) 
 17   format( '>>> line ',i3,": ", 1x, a30, 1x, 2(i5,1x,i5) ) 
      end subroutine ep_readconf_qp
!-----------------------------------------------------------------------
!
!     Read the setting parameters defined by the colon ':'.
!
      subroutine ep_readconf_a_qp ( 
     $           fname, nome, ng, g, e,  
     $           istat )
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      character(len=*), intent(in) :: fname
      integer, intent(inout) :: ng 
      integer, intent(out) :: nome 
      real(wrp), allocatable, intent(out) :: e(:), g(:,:)
      integer, intent(out) :: istat
!
!     LOCAL VARIABLES:
      integer,parameter :: icha_bs  = 32              ! White Space
      character,parameter :: char_bs= achar(icha_bs)  ! White Space
      integer,parameter :: icha_tab = 9               ! Horizontal Tab
      integer,parameter :: icha_c1  = 33              ! '!' 
      integer,parameter :: icha_c2  = 35              ! '#'
      integer,parameter :: icha_eq  = iachar('=')     ! '=' 
      integer,parameter :: icha_colon = iachar(':')   ! ':' 
      integer,parameter :: icha_0   = iachar('0')     ! '0'
      integer,parameter :: icha_9   = iachar('9')     ! '9'
      integer :: 
     $   i, j, k, l, io, n, i0, i1, i2, j0, j1, j2, 
     $   k0, i3, i4 
      character :: 
     $   charg*550, chnam*50, chval*500, chn*3, chfmt*50, ch*1
      character(len=:), allocatable :: ch1, ch2
      integer :: nsta
      integer, allocatable :: ista(:) 
      real(wrp), allocatable :: esta(:), vsta(:,:)
      real(wrp) :: rtmp 
!
!     DEPENDENCES: 
      intrinsic :: iachar, len_trim, trim, adjustl, real  
!
!     SHARED PARAMETERS:
      logical 
     $   flag_prep, 
     $   flag_wrtinp 
      common /epcode_shared_flags/ 
     $   flag_prep, 
     $   flag_wrtinp 
!
!     EXECUTABLE STATEMENTS: 
!
      istat = 0 
      open(newunit=io, file=fname, status='old', iostat=istat)
      if ( istat .ne. 0 ) then
         write(*,'(t15,a)') 'ERROR: could not open file '//fname 
         return
      endif
      n = 0
      do
         charg = '' 
         read(unit=io, fmt='(a)', end=101, err=102) charg
         n = n + 1 
         charg = adjustl(charg)
         l = len_trim(charg)
         if ( l .eq. 0 ) cycle
         if ( charg(1:l) .eq. '__EOF__' ) goto 101 
         j = iachar(charg(1:1))
         if (j .eq. icha_c1 .or. j .eq. icha_c2) cycle 
         i = 0  
         do k = 1, l 
            if ( iachar(charg(k:k)) .ne. icha_colon ) cycle 
            i = 1  
            exit 
         enddo
         if ( i .eq. 0 ) cycle
         chnam = charg(:k-1)
         chval = trim(adjustl(charg(k+1:l)))
         select case ( uppercase(trim(adjustl(chnam))) )
!!       case('STATES','states','States')
         case('STATES')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1  
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if (j .eq. icha_c1 .or. j .eq. icha_c2) chval = ''  
               endif 
            enddo
            if ( flag_wrtinp )
     $         write(*,12) n, 'States:', trim(chval) 
            chval = trim(adjustl(chval)) // char_bs
            nsta = 0 
            do i = 1, len(chval) - 1 
               if ( iachar(chval(i:i)) .eq. icha_eq ) 
     $            nsta = nsta + 1 
            enddo   
            if ( nsta .eq. 0 ) cycle 
            if ( flag_wrtinp )
     $         write(*,10) n, 'N_STATES =', nsta
            if ( allocated(ista) ) deallocate(ista) 
            allocate ( ista(nsta) ) 
            i1 = 1
            i2 = 0 
            k = len(chval)
            i = i2 
            do while ( i .le. k )
               i = i + 1  
               if ( iachar(chval(i:i)) .ne. icha_bs .and. 
     $              iachar(chval(i:i)) .ne. icha_tab ) cycle
               i2 = i 
               ch1 = chval(i1:i2-1) 
               do while ( iachar(chval(i2:i2)) .eq. icha_bs .or.
     $                    iachar(chval(i2:i2)) .eq. icha_tab ) 
                  i2 = i2 + 1 
               enddo  
               ch2 = chval(i2:k)
               i1 = i2
               i = i2 
               do i3 = 1, len(ch1)
                  i4 = iachar(ch1(i3:i3))
                  if ( i4 .lt. icha_0 .or. i4 .gt. icha_9 ) 
     $               ch1(i3:i3) = char_bs
               enddo 
               read( ch1, * ) i0, j0, k0   
               ista(i0) = (j0 + 1)/2
            enddo 
            nome = 0 
            do i0 = 1, nsta
               nome = nome + ista(i0) 
            enddo
            if ( flag_wrtinp ) then 
               chn = '' 
               write(chn,fmt='(I3)') nsta 
               chfmt = '('//trim(adjustl(chn))//'(I5,1x))'
               write(charg,fmt=chfmt) ista 
               write(*,12) n, 'I_States =', trim(charg) 
               write(*,10) n, 'NOME =', nome
            endif 
!!       case('V0','v0')
         case('V0')
            if ( ng .eq. 1 ) cycle 
            ng = nome  
            if ( flag_wrtinp )
     $         write(*,10) n, 'NG =', ng  
            if ( allocated(vsta) ) deallocate(vsta)
            allocate ( vsta(nsta,nsta) )
            if ( allocated(g) ) deallocate(g) 
            allocate( g(ng,ng) )         
            do j = 1, nsta
            do k = j, nsta 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1   
               do while ( len_trim(chval) .eq. 0 ) 
                  read(unit=io, fmt='(a)', end=101, err=102) chval 
                  n = n + 1  
                  if ( len_trim(chval) .ne. 0 ) then 
                     chval = adjustl(chval)
                     j0 = iachar(chval(1:1))
                     if ( j0 .eq. icha_c1 .or. 
     $                    j0 .eq. icha_c2 ) chval = ''  
                  endif 
               enddo
               chval = trim(adjustl(chval))
               if ( flag_wrtinp )
     $            write(*,12) 
     $            n, 'j, k, V0_states(j,k):', '' 
               read(chval, fmt=*, err=102) i0, j0, rtmp 
               vsta(i0,j0) = rtmp 
               vsta(j0,i0) = rtmp 
               if ( flag_wrtinp )
     $            write(*,*) i0, j0, rtmp  
            enddo
            enddo 
            if ( flag_wrtinp ) then 
               write(*,12) n, 'V0_states(:,:) =', '' 
               chn = '' 
               write(chn,fmt='(I3)') nsta 
               chfmt = '(I5, ":",2x,'//trim(adjustl(chn))//'(1X,ES9.2))'
               do k = 1, nsta
                  write(*,fmt=chfmt) k, vsta(k,:)
               enddo
           endif 
            do j0 = 1, nsta
               do i0 = 1, nsta 
                  vsta(i0,j0) = 2 * vsta(i0,j0) / 
     $            sqrt( real(2*ista(i0)*2*ista(j0),kind=wrp) ) 
               enddo 
            enddo 
            i2 = 0 
            do i0 = 1, nsta
               i1 = ista(i0)
               j2 = 0 
               do j0 = 1, nsta 
                  j1 = ista(j0)
                  g( i2+1:i2+i1, j2+1:j2+j1 ) = vsta(i0,j0)
                  j2 = j2+j1 
               enddo 
               i2 = i2+i1 
            enddo
            if ( flag_wrtinp ) then 
               chn = '' 
               write(chn,fmt='(I3)') nome 
               chfmt = '(I5, ":",2x,'//
     $            trim(adjustl(chn))//'(1X,ES9.2))'
               write(*,12) n, 'G(:,:) =', '' 
               do k = 1, ng 
                  write(*,fmt=chfmt) k, g(k,:)
               enddo 
            endif 
!!       case('G','g')
         case('G')
            if ( ng .eq. 1 ) cycle 
            ng = nome 
            if ( flag_wrtinp )
     $         write(*,10) n, 'NG =', ng  
            if ( allocated(vsta) ) deallocate(vsta)
            allocate ( vsta(nsta,nsta) )
            if ( allocated(g) ) deallocate(g) 
            allocate( g(ng,ng) )         
            do j = 1, nsta
            do k = j, nsta 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1   
               do while ( len_trim(chval) .eq. 0 ) 
                  read(unit=io, fmt='(a)', end=101, err=102) chval 
                  n = n + 1  
                  if ( len_trim(chval) .ne. 0 ) then 
                     chval = adjustl(chval)
                     j0 = iachar(chval(1:1))
                     if ( j0 .eq. icha_c1 .or. 
     $                    j0 .eq. icha_c2 ) chval = ''  
                  endif 
               enddo
               chval = trim(adjustl(chval))
               if ( flag_wrtinp )
     $            write(*,12) 
     $            n, 'j, k, G_states(j,k):', '' 
               read(chval, fmt=*, err=102) i0, j0, rtmp 
               vsta(i0,j0) = rtmp 
               vsta(j0,i0) = rtmp 
               if ( flag_wrtinp )
     $            write(*,*) i0, j0, rtmp  
            enddo
            enddo 
            if ( flag_wrtinp ) then 
               write(*,12) n, 'G_states(:,:) =', ''
               chn = '' 
               write(chn,fmt='(I3)') nsta 
               chfmt = '(I5, ":",2x,'//
     $            trim(adjustl(chn))//'(1X,ES9.2))'
               do k = 1, nsta
                  write(*,fmt=chfmt) k, vsta(k,:)
               enddo    
            endif 
            i2 = 0 
            do i0 = 1, nsta
               i1 = ista(i0)
               j2 = 0 
               do j0 = 1, nsta 
                  j1 = ista(j0)
                  g( i2+1:i2+i1, j2+1:j2+j1 ) = vsta(i0,j0)
                  j2 = j2+j1 
               enddo 
               i2 = i2+i1 
            enddo 
            if ( flag_wrtinp ) then 
               chn = '' 
               write(chn,fmt='(I3)') nome 
               chfmt = '(I5, ":",2x,'//
     $            trim(adjustl(chn))//'(1X,ES9.2))'
               write(*,12) n, 'G(:,:) =', '' 
               do k = 1, ng 
                  write(*,fmt=chfmt) k, g(k,:)
               enddo 
            endif 
!!       case('E','e')
         case('E')
            do while ( len_trim(chval) .eq. 0 ) 
               read(unit=io, fmt='(a)', end=101, err=102) chval 
               n = n + 1 
               if ( len_trim(chval) .ne. 0 ) then 
                  chval = adjustl(chval)
                  j = iachar(chval(1:1))
                  if ( j .eq. icha_c1 .or. 
     $                 j .eq. icha_c2) chval = ''  
               endif 
            enddo
            if ( nsta .eq. 0 ) cycle 
            if ( allocated(esta) ) deallocate(esta)
            allocate ( esta(nsta) )
            read(chval, fmt=*, iostat=i1) esta
            if ( i1 .ne. 0 ) then 
               backspace(unit=io) 
               do j = 1, nsta
                  read(unit=io, fmt=*, err=102) esta(j)
               enddo 
            endif
            if ( flag_wrtinp ) 
     $         write(*,15) n, 'E_states(:) =', esta
            if ( allocated(e) ) deallocate(e)
            allocate ( e(nome) )
            i2 = 0 
            do i0 = 1, nsta
               i1 = ista(i0)
               e(i2+1:i2+i1) = esta(i0)
               i2 = i2+i1 
            enddo 
            if ( flag_wrtinp )
     $         write(*,15) n, 'E(:) =', e
         case default 
            cycle 
         end select          
      enddo    
 101  continue
      if ( flag_wrtinp ) 
     $   write(*,14) 
     $   n, 'Finish reading the configuration file for STATES.'
      close(io)
      return
 102  continue
      close(io)
      istat = 100  
      write(*,12) 
     $   n, 'ERROR as reading record.', 
     $   'The configuration process is terminated!'
      return  
 10   format( '>>> line ',i3,": ", 1x, a30, 1x, i5 ) 
 11   format( '>>> line ',i3,": ", 1x, a30, 1x, e14.7 ) 
 12   format( '>>> line ',i3,": ", 1x, a30, 1x, a ) 
 13   format(  14x,                1x, a30, /, 5(3x,f7.3) ) 
 14   format( '>>> line ',i3,": ", a, // ) 
 15   format( '>>> line ',i3,": ", 1x, a30, /, 5(3x,f7.3) ) 
      end subroutine ep_readconf_a_qp 
!-----------------------------------------------------------------------
      elemental function uppercase (str_i) result(str_o)
      implicit none
!
!     ARGUMENTS: 
      character(len=*), intent(in) :: str_i
      character(len=len(str_i)) :: str_o
!
!     LOCAL VARIABLES:
      integer :: i, j
!
!     DEPENDENCES: 
      intrinsic :: achar, iachar, len
!
!     EXECUTABLE STATEMENTS: 
!
      do i = 1, len(str_i)
         j = iachar(str_i(i:i))
         if (j .ge. iachar("a") .and. j .le. iachar("z") ) then
            str_o(i:i) = achar(iachar(str_i(i:i)) - 32 )
         else
            str_o(i:i) = str_i(i:i)
         endif
      enddo
      return 
      end function uppercase
!-----------------------------------------------------------------------
      end module epcode_mod 
!-----------------------------------------------------------------------
!
!
!
