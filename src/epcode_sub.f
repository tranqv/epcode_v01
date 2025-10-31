!-----------------------------------------------------------------------
!***********************************************************************
!     real(kind=8)
!-----
      subroutine ep_gstate_1a_dp ( 
     $           nome, npar, e, g, neiv, 
     $           prefix, tole ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar 
      real(wrp), intent(in) :: 
     $   e(nome), g, tole
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     LOCAL VARIABLE:
      real(wrp)
     $   g_sca(1,1)
!
!     DEPENDENCES:
      external 
     $   ep_gstate_main_dp  
!
!     EXECUTABLE STATEMENTS: 
!
      g_sca(1,1) = g 
      call ep_gstate_main_dp ( 
     $     nome, npar, e, 1, g_sca, neiv, 
     $     prefix, tole ) 
      return 
      end subroutine ep_gstate_1a_dp
!-----------------------------------------------------------------------
      subroutine ep_gstate_1b_dp ( 
     $           nome, npar, e, g, neiv, 
     $           prefix ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar 
      real(wrp), intent(in) :: 
     $   e(nome), g
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_1a_dp 
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_1a_dp ( 
     $     nome, npar, e, g, neiv, 
     $     prefix, 0.0_wrp ) 
      return 
      end subroutine ep_gstate_1b_dp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_2a_dp ( 
     $           nome, npar, e, g, neiv, 
     $           prefix, tole ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar 
      real(wrp), intent(in) :: 
     $   e(nome), g(nome,nome), tole
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_main_dp  
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_main_dp ( 
     $     nome, npar, e, nome, g, neiv, 
     $     prefix, tole ) 
      return 
      end subroutine ep_gstate_2a_dp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_2b_dp ( 
     $           nome, npar, e, g, neiv, 
     $           prefix ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar 
      real(wrp), intent(in) :: 
     $   e(nome), g(nome,nome)
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_2a_dp 
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_2a_dp ( 
     $     nome, npar, e, g, neiv, 
     $     prefix, 0.0_wrp ) 
      return 
      end subroutine ep_gstate_2b_dp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_3a_dp ( 
     $           nome, npar, e, ng, g, neiv, 
     $           prefix, tole ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar, ng 
      real(wrp), intent(in) :: 
     $   e(nome), g(ng,ng), tole
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_main_dp  
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_main_dp ( 
     $     nome, npar, e, ng, g, neiv, 
     $     prefix, tole ) 
      return 
      end subroutine ep_gstate_3a_dp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_3b_dp ( 
     $           nome, npar, e, ng, g, neiv, 
     $           prefix ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar, ng 
      real(wrp), intent(in) :: 
     $   e(nome), g(ng,ng)
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_3a_dp 
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_3a_dp ( 
     $     nome, npar, e, ng, g, neiv, 
     $     prefix, 0.0_wrp ) 
      return 
      end subroutine ep_gstate_3b_dp 
!-----------------------------------------------------------------------
!
!     The core of EP (real*8)
!
      subroutine ep_gstate_main_dp ( 
     $           nome, npar, e, ng, g, neiv, 
     $           prefix, tole ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar, ng 
      real(wrp), intent(in) :: 
     $   e(nome), g(ng,ng), tole
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     LOCAL VARIABLES AND CONSTANTS:
      real(wrp) 
     $   g_diag(63), g_matr(63,63)
      integer(wip)
     $   binom_tab(0:63,0:63)
      common /ep_gstate_main_shared/ binom_tab
      logical, parameter ::
     $   flag_screen = .true.,   
     $   flag_wrtres = .true.,  
!!   $   flag_wrtout = .true.,   
     $   flag_wrtout = .false.,  
!!   $   flag_wrtsiz = .true.,  
     $   flag_wrtsiz = .false., 
     $   flag_memest = .true., 
     $   flag_wrtnsj = .false.,  
     $   flag_wrtmat = .false., 
     $   flag_wrtevv = .false., 
     $   flag_wrtocc = .false.
      integer(wip)
     $   lsj, meiv 
      integer(i1b), allocatable :: 
     $   njjb(:), nsjb(:), nb1(:), nb2(:)
      integer(wip), allocatable :: 
     $   nsj(:), njj(:), mjp(:)
      integer(wip), allocatable :: 
     $   jah(:)
      real(wrp), allocatable :: 
     $   aah(:)
      real(wrp), allocatable :: 
     $   occ(:,:), ecc(:,:)
      real(wrp)
     $   tmp 
      integer
     $   mjp_len, ljp, nps, nmns, lb1, i, i1, i2, isum, 
     $   k, l, nt0, iox, ns 
      integer(wip)
     $   nrow, ncol, ir, ic, isj, ijj, ljj,
     $   nnz, iele, nele, jos, state_i, state_f,  
     $   njir, njic, ntmp, nde
      integer(wip)      
     $   ncv, lawork, lworkl, 
     $   io_resid, io_workd, io_eivec, io_eival, io_workl
      real(wrp), allocatable :: 
     $   awork(:)
      logical, allocatable :: 
     $   choose(:)
      character
     $   chn*5, chpack*30, chfmt*50, 
     $   fort00*90, fort01*90, fort02*90, fort03*90, 
     $   fort04*90, fort05*90, header*200
      integer ::
     $   iou00, iou01, iou02, iou03, iou04, iou05 
      integer(i1b)
     $   k1, k2  
      logical 
     $   g_is_matrix, use_arpack  
      integer(wip)      
     $   num_of_bytes_in_total, num_of_bytes_in_arrays,
     $   num_of_bytes_in_locvar
      integer(8) 
     $   iram 
!
!     DEPENDENCES:
      intrinsic :: 
     $   mod, min, max, sum, real, int, trim, adjustl, len_trim 
      integer(wip), external ::
     $   binom
      integer(wip), external ::
     $   npairing 
      integer(wip), external ::
     $   hashs
      external ::
     $   prt_nsj_njj, 
     $   deccon, 
     $   gens, genp, gind 
      character(len=80), external ::
     $   cparams0 
      external ::
     $   arpack_exec_dp, 
     $   lapack_exec_dp,  
     $   copy_s_to_d_symmsr_dp,
     $   higham_sum_dp 
      integer(wip), external ::
     $   wilaenv
      external 
     $   get_ram_in_bytes,
     $   ep_getmem_dp
!
!     SHARED PARAMETERS:
      logical 
     $   flag_prep, 
     $   flag_wrtinp 
      common /epcode_shared_flags/ 
     $   flag_prep, 
     $   flag_wrtinp 
!
!     DEFAULT VALUES:
      data                      
     $   flag_prep   /.false./,
!!   $   flag_wrtinp /.true./  
     $   flag_wrtinp /.false./ 
!
!     EXECUTABLE STATEMENTS: 
!
      if ( nome.lt.1 .or. nome.gt.63 ) 
     $   stop "ERROR: NOME out of [1,63]"
      if ( npar.lt.1 .or. npar.gt.nome*2 ) 
     $   stop "ERROR: NPAR out of [1,2*NOME]"
      if ( .not. (ng.eq.1 .or. ng.eq.nome) ) 
     $   stop "ERROR: NG /= 1 and NG /= NOME"
      g_is_matrix = ( ng .ne. 1 ) 
      if ( g_is_matrix ) then 
         g_diag = 0
         g_matr = 0 
         do k = 1, nome
            do l = 1, nome
               g_matr(l,k) = -g(l,k) 
            enddo
            g_diag(k) = -g(k,k)
         enddo
      else 
         g_diag(1) = -g(1,1)
      endif
      binom_tab = 0
      do k = 0, nome 
      do i = 0, k
         binom_tab ( k, i ) = binom ( k, i )
      enddo
      enddo 
      ns   = mod(npar,2)
      nmns = nome - ns 
      nps  = npar/2 - ns/2 
      mjp_len = nps * ( nmns - nps )
      lsj  = binom_tab ( nome, ns )
      if ( lsj .lt. 1 ) stop "ERROR: lsj < 1" 
      ljj  = binom_tab ( nmns, nps )
      if ( ljj .lt. 1 ) stop "ERROR: ljj < 1" 
      meiv = min( ljj, neiv )
      if ( meiv .lt. 1 ) stop "ERROR: no. of eigenvalues < 1" 
      allocate ( njj(ljj) ) 
      njj = 0 
      call gens ( nmns, nps, ljj, njj )
      ncol = ljj 
      nrow = ncol 
      allocate ( occ(nome,meiv), ecc(nome,meiv) )
      occ = 0.0_wrp 
      ecc = 0.0_wrp 
!!    if ( flag_wrtout ) 
!!   $   write(*,'(1x,a15,1x,i20)') 
!!   $   'Ncol = Nrow =', nrow 
      iele = 0 
      do ir = 1_wip, ljj
         iele = iele + npairing ( njj(ir), nmns )
      enddo 
      nele = iele + nrow
      nnz  = 2*( nele - nrow ) + nrow
!!    if ( flag_wrtout ) 
!!   $   write(*,
!!   $   '(/,a,/,1x,a15,1x,i20,/,3(1x,a15,1x,i20,1x,a,/))') 
!!   $   'MATRIX INFORMATION:',
!!   $   'Ncol = Nrow =', nrow, 
!!   $   'Nele =', nele, 
!!   $   '(No. Non-Zero elements of the half of H)', 
!!   $   'NNZ =',  nnz, 
!!   $   '(No. Non-Zero elements of H)',
!!   $   'LJJ**2 - NNZ =', ljj*int(ljj,kind=wip) - nnz, 
!!   $   '(No. Zeros of H)'
      nele = nele + 1 
      if ( flag_prep ) goto 201  
      allocate ( jah(nele) )
      if ( g_is_matrix ) then 
         allocate ( aah(nele) )
      else 
         allocate ( aah(nrow) )
      endif 
      aah = 0.0_wrp  
      jah = 0_wip  
      iele = ljj + 1  
      allocate ( mjp ( mjp_len ) )
      if ( g_is_matrix ) then 
         do ir = 1, ljj
            call genp ( nmns, nps, njj(ir), mjp_len, mjp, ljp ) 
            do k = 1, ljp
               ic = hashs ( nmns, nps, mjp(k) ) 
               iele = iele + 1 
               if ( k .eq. 1 ) jah(ir) = iele 
               jah(iele) = ic 
               call gind ( njj(ir), njj(ic), k1, k2 ) 
               aah(iele) = g_matr ( k1, k2 )
            enddo 
         enddo 
      else 
         do ir = 1, ljj
            call genp ( nmns, nps, njj(ir), mjp_len, mjp, ljp ) 
            do k = 1, ljp
               ic = hashs ( nmns, nps, mjp(k) ) 
               iele = iele + 1 
               if ( k .eq. 1 ) jah(ir) = iele 
               jah(iele) = ic 
            enddo 
         enddo 
      endif 
      deallocate ( mjp )
      ic = ljj 
      do while ( jah(ic) .eq. 0 ) 
         ic = ic - 1 
      enddo 
      if ( ic .eq. ljj ) then 
         jah(ic+1) = jah(ic)
      else 
         jah(ic+1) = jah(ic) + 1
         do ir = ic+2, ljj+1
            jah(ir) = jah(ir-1) 
         enddo 
      endif 
!
!               PREFIX
!               |
!     header = "EP10_n10p11g0600s01" for NOME=10, NPAR=11, G=0.6, NS=1
!                    |  |  |    |        
!              NOME=10  |  G=0.6|    
!                 NPAR=11       |
!                               NS=1
!
      header = ''
      if ( len_trim(adjustl(prefix)) .gt. 0 )
     $   header = trim(adjustl(header)) // trim(adjustl(prefix))
      header = trim(adjustl(header)) //
     $         trim(cparams0 ( nome, npar, real(g(1,1)), ns) )
!
      if ( flag_wrtres ) then 
         fort00 = trim(adjustl(header))//'_result.txt' 
         open( newunit=iou00, file=trim(adjustl(fort00)), 
     $         status='replace' )
         rewind(iou00)
      endif 
      if ( flag_wrtmat ) then 
         fort01 = trim(adjustl(header))//'_iajaaa.bin'
         open( newunit=iou01, file=trim(adjustl(fort01)), 
     $   status='replace', form='unformatted' )
         rewind(iou01)
      endif 
      if ( flag_wrtnsj ) then 
         fort02 = trim(adjustl(header))//'_nsjnjj.txt'
         open( newunit=iou02, file=trim(adjustl(fort02)), 
     $         status='replace' )
         rewind(iou02)
      endif 
      if ( flag_wrtsiz ) then 
         fort03 = trim(adjustl(header))//'_dtsize.txt' 
         open( newunit=iou03, file=trim(adjustl(fort03)), 
     $         status='replace' )
         rewind(iou03)
      endif 
      if ( flag_wrtevv ) then  
         fort04 = trim(adjustl(header))//'_eivave.bin'
         open( newunit=iou04, file=trim(adjustl(fort04)), 
     $         status='replace', form='unformatted' )
         rewind(iou04)
      endif 
      if ( flag_wrtocc ) then 
         fort05 = trim(adjustl(header))//'_eivocc.bin' 
         open( newunit=iou05, file=trim(adjustl(fort05)), 
     $      status='replace', form='unformatted' )
         rewind(iou05)
         prefix = trim(adjustl(fort05))
      endif 
!
      allocate ( 
     $   njjb(nome), nsjb(nome), nb1(nome), nb2(nome) )
      njjb = 0
      nsjb = 0
      nb1  = 0
      nb2  = 0 
      if ( flag_wrtocc ) then 
         write(iou05) ns, nome, meiv
      endif 
 201  continue 
      use_arpack = ( ljj .gt. 100 ) 
      if ( use_arpack ) then 
         chpack = 'ARPACK'
         ncv    = max( meiv*2, meiv + 10 )    
         lworkl = ncv * ( ncv + 8 )
         lawork = ncv*2 + nrow*ncv + nrow + 3*nrow + lworkl
         io_eival = 0 
         io_eivec = io_eival + ncv*2
         io_resid = io_eivec + ncv*nrow 
         io_workd = io_resid + nrow 
         io_workl = io_workd + nrow*3
      else
         chpack = 'LAPACK'
         ncv = 0 
         lworkl = wilaenv ( 
     $      1, 'WRSYTRD', 'U', nrow, -1_wip, -1_wip, -1_wip )
         lworkl = max( 1_wip, ( lworkl + 2 )*nrow )
         io_eivec = 0
         io_eival = io_eivec + nrow*nrow 
         io_workl = io_eival + nrow 
         lawork = nrow*nrow + nrow + lworkl  
      endif 
      if ( flag_prep ) goto 202 
      if ( use_arpack ) then 
         allocate ( awork ( lawork ) ) 
         allocate ( choose ( ncv ) )
         awork = 0.0_wrp  
         choose = .false.
      else
         allocate ( awork ( lawork ) ) 
         awork = 0.0_wrp 
         call copy_s_to_d_symmsr_dp  ( 
     $        nrow, aah, jah, g_diag(1), g_is_matrix, 
     $        awork(io_eivec+1) )
      endif 
 202  continue 
      if ( flag_wrtsiz ) then 
         iox = iou03
      else
         iox = 123  
         goto 506 
      endif 
 505  continue 
      write(iox,'(72("-"))') 
!!    write(iox,'( 6(/,1x,a15,1x,i20) )') 
!!   $   'NOME =', nome, 
!!   $   'NPAR =', npar,
!!   $   'NS   =', ns,
!!   $   'LSJ  =', lsj, 
!!   $   'LJJ  =', ljj, 
!!   $   'NEIV =', meiv 
      write(iox, '( 3(/,1x,a15,1x,i20))') 
     $   'NOME =', nome, 
     $   'NPAR =', npar,
     $   'NEIV =', meiv
      write(iox, '( 1x,a15,1x,e20.13 )') 
     $   'TOLE =', tole  
      write(iox, '( 1x,a15,1x, a )') 
     $   'PREFIX =', trim(prefix)
      if ( ng .eq. 1 ) then 
         write(iox, '( 1x,a15,1x,e20.13 )') 
     $   'G =', g(1,1)  
      else
         write(iox,'( 1x, a15 )') 'G ='
         chn = '' 
         write(chn,fmt='(I3)') nome 
         chfmt = '(I5, ":",2x,'//
     $      trim(adjustl(chn))//'(1X,ES9.2))'
         do k = 1, ng 
            write(iox,fmt=chfmt) k, g(k,:)
         enddo 
      endif
      write(iox,'( 1x, a15 )') 'SPE ='
      chfmt = '(I5, ":", 3x, ES9.2 )'
      do k = 1, nome 
         write(iox,fmt=chfmt) k, e(k)
      enddo 
      write(iox,
     $   '(/,a,/,1x,a15,1x,i20,/,3(1x,a15,1x,i20,1x,a,/))') 
     $   'MATRIX INFORMATION:',
     $   'Ncol = Nrow =', nrow, 
     $   'Nele =', nele-1, 
     $   '(No. Non-Zero elements of the half of H)', 
     $   'NNZ =',  nnz, 
     $   '(No. Non-Zero elements of H)',
     $   'LJJ**2 - NNZ =', ljj*int(ljj,kind=wip) - nnz, 
     $   '(No. Zeros of H)'
  506 continue 
      if (iox .ne. 6 .and. flag_screen ) then 
         iox = 6 
         goto 505 
      endif
      call ep_getmem_dp (
     $     nome, npar, ng, neiv, prefix, 
     $     ljj, nele, mjp_len, lsj, lawork, ncv, 
     $     num_of_bytes_in_arrays, 
     $     num_of_bytes_in_locvar, 
     $     num_of_bytes_in_total ) 
      if ( flag_wrtsiz ) then 
         iox = iou03
      else
         iox = 123  
         goto 406 
      endif 
  405 continue 
      write(iox,'(/,A)') 'COMPUTER MEMORY REQUIREMENT:' 
      write(iox,'( 3x, I20, A10, 2(E17.10, A) )') 
     $    num_of_bytes_in_total, ' (bytes) =', 
     $    num_of_bytes_in_total * dble(2)**(-20), ' (MB) =',
     $    num_of_bytes_in_total * dble(2)**(-30), ' (GB)'
  406 continue 
      if (iox .ne. 6 .and. flag_screen ) then 
         iox = 6 
         goto 405 
      endif
  400 continue
      if ( flag_prep ) goto 800  
      isj = 1 
      nsjb = 0 
      if (mod(npar,2) .ne. 0) nsjb ( (npar-1)/2 + 1 ) = 1 
      k = 0 
      do i = 1, nome
         if ( nsjb(i) .ne. 0 ) cycle 
         k = k + 1
         nb1(k) = i
      enddo
      lb1 = k 
      if ( g_is_matrix ) then 
         do iele = 1_wip, ljj
            njjb = nsjb 
            call deccon ( njj(iele), nome, nb2 )
            do i = 1, lb1  
               njjb( nb1(i) ) = nb2(i)*2 
            enddo
            aah(iele) = sum ( e*njjb + 0.25_wrp*g_diag(1:nome)*
     $                        (njjb - nsjb)*(4 - njjb - nsjb) ) 
         enddo
      else 
         do iele = 1_wip, ljj
            njjb = nsjb 
            call deccon ( njj(iele), nome, nb2 )
            do i = 1, lb1  
               njjb( nb1(i) ) = nb2(i)*2 
            enddo
            aah(iele) = sum ( e*njjb + 0.25_wrp*g_diag(1)*
     $                        (njjb - nsjb)*(4 - njjb - nsjb) ) 
         enddo
      endif 
      if ( flag_wrtmat ) then  
         write(iou01) nrow, nele 
         write(iou01) jah
         write(iou01) aah
      endif 
      if ( use_arpack ) then 
         awork = 0.0_wrp  
         choose = .false.
         call arpack_exec_dp  ( 
     $        nrow, meiv, ncv, lworkl,
     $        awork(io_eival+1), awork(io_eivec+1), 
     $        awork(io_resid+1), awork(io_workd+1), 
     $        awork(io_workl+1), choose, 
     $        nele, aah, jah, g_diag(1), g_is_matrix, 
     $        flag_wrtout, tole )
      else
         jos = io_eivec - nrow 
         do ir = 1_wip, nrow 
            awork ( jos + ir*(1+nrow) ) = aah (ir) 
         enddo 
         call lapack_exec_dp  ( 
     $        nrow, 
     $        awork(io_eivec+1), awork(io_eival+1),  
     $        lworkl, awork(io_workl+1),  
     $        nele, aah, jah, flag_wrtout )
      endif 
      if ( flag_wrtevv ) then  
         write(iou04) meiv, ljj
         jos = io_eivec 
         do ic = 1_wip, meiv
            write(iou04) ic, awork(io_eival+ic)
            write(iou04) awork(jos+1:jos+ljj)
            jos = jos + ljj
         enddo  
      endif 
      nde = 2**ns 
      occ = 0.0_wrp 
      ecc = 0.0_wrp 
      do ir = 1_wip, ljj 
         njjb = nsjb 
         call deccon ( njj(ir), nome, nb2 )
         do i = 1, lb1  
            njjb( nb1(i) ) = nb2(i)*2 
         enddo
         jos = io_eivec + ir
         do ic = 1_wip, meiv
            tmp = awork ( jos )**2 * 0.5_wrp 
            jos = jos + nrow 
            do k = 1, nome
               call higham_sum_dp ( 
     $              occ(k,ic), tmp*njjb(k), ecc(k,ic) )  
            enddo 
         enddo 
      enddo
      chn = '' 
      write(chn,fmt='(I5)') nome 
      if ( flag_wrtres ) then
         select case ( wrp )
         case(kind(1.0e0))
            write(iou00,
     $      '("#",A9,1X,A22,1X,A22)') 
     $      'i', 'Eigenvalue(i)',  
     $      '( Occ(j,i), j=1,'//trim(adjustl(chn))//' )'
            chfmt = '(I10,1X,ES22.8,'//
     $            trim(adjustl(chn))//'(1X,ES22.8))'
         case(kind(1.0d0))
            write(iou00,
     $      '("#",A9,1X,A22,1X,A22)') 
     $      'i', 'Eigenvalue(i)',  
     $      '( Occ(j,i), j=1,'//trim(adjustl(chn))//' )'
            chfmt = '(I10,1X,ES22.15,'//
     $            trim(adjustl(chn))//'(1X,ES22.15))'
         case(kind(1.0d0)*2)
            write(iou00,
     $      '("#",A9,1X,A40,1X,A40)') 
     $      'i', 'Eigenvalue(i)',  
     $      '( Occ(j,i), j=1,'//trim(adjustl(chn))//' )'
            chfmt = '(I10,1X,ES40.32,'//
     $            trim(adjustl(chn))//'(1X,ES40.32))'
         end select 
         k = len_trim(chfmt)
         do ic = 1_wip, meiv 
            write(iou00,fmt=chfmt(1:k))
     $         ic, awork(io_eival+ic), occ(1:nome,ic)
         enddo
      endif 
      if ( flag_wrtocc ) then 
         write(iou05) awork(io_eival+1:io_eival+meiv), occ
      endif 
      if ( flag_wrtnsj ) then  
         allocate ( nsj(lsj) )  
         nsj = 0
         call gens ( nome, ns, lsj, nsj )
         call prt_nsj_njj (iou02, nome, lsj, nsj, ljj, njj)
         deallocate ( nsj )
      endif 
      if ( flag_screen ) then  
         write(*,'(/,a)') 'CALCULATIONS COMPLETED' 
         if (  flag_wrtres .or. 
     $         flag_wrtmat .or.
     $         flag_wrtnsj .or. 
     $         flag_wrtsiz .or. 
     $         flag_wrtevv .or. 
     $         flag_wrtocc )   
     $      write(*,'(a)') 'Check:' 
      endif 
      if ( flag_wrtres .and. flag_screen )   
     $   write(*,11) '', adjustl(fort00), '(ascii)'
      if ( flag_wrtmat .and. flag_screen )   
     $   write(*,11) '', adjustl(fort01), '(binary)'
      if ( flag_wrtnsj .and. flag_screen )   
     $   write(*,11) '', adjustl(fort02), '(ascii)'
      if ( flag_wrtsiz .and. flag_screen )   
     $   write(*,11) '', adjustl(fort03), '(ascii)'
      if ( flag_wrtevv .and. flag_screen )   
     $   write(*,11) '', adjustl(fort04), '(binary)' 
      if ( flag_wrtocc .and. flag_screen )   
     $   write(*,11) '', adjustl(fort05), '(binary)'    
      if ( flag_wrtres ) close (iou00)
      if ( flag_wrtmat ) close (iou01)
      if ( flag_wrtnsj ) close (iou02)
      if ( flag_wrtsiz ) close (iou03)
      if ( flag_wrtevv ) close (iou04)
      if ( flag_wrtocc ) close (iou05)
      if ( flag_wrtout .and. flag_screen ) then 
         write(*,'(/,A)') 
     $   'Eigenvalue and Occupation numbers:' 
         write(chn,'(I5)') 1 + nome
         select case ( wrp )
         case(kind(1.0e0))
            chfmt = '('//trim(adjustl(chn))//'ES16.8)'
         case(kind(1.0d0))
            chfmt = '('//trim(adjustl(chn))//'ES23.15)'
         case(kind(1.0d0)*2)
            chfmt = '('//trim(adjustl(chn))//'ES40.32)'
         end select 
         k = len_trim ( chfmt )
         do ic = 1_wip, meiv 
            write(*,fmt=chfmt(1:k)) 
     $      awork(io_eival+ic), occ(1:nome,ic)
         enddo
      endif 
  800 continue  
      if (allocated(nsj    )) deallocate (nsj    )
      if (allocated(njj    )) deallocate (njj    )
      if (allocated(njjb   )) deallocate (njjb   )
      if (allocated(nsjb   )) deallocate (nsjb   )
      if (allocated(nb1    )) deallocate (nb1    )
      if (allocated(nb2    )) deallocate (nb2    )
      if (allocated(jah    )) deallocate (jah    )
      if (allocated(aah    )) deallocate (aah    )
      if (allocated(occ    )) deallocate (occ    )
      if (allocated(ecc    )) deallocate (ecc    )
      if (allocated(mjp    )) deallocate (mjp    )
      if (allocated(awork  )) deallocate (awork  ) 
      if (allocated(choose )) deallocate (choose )
      return 
 11   format ( A10, 1x, A50, 1x, A10 )
      end subroutine ep_gstate_main_dp 
!-----------------------------------------------------------------------
!
!     This is a copied version of the program published in ARPACK: 
!        arpack96/EXAMPLES/SYM/dsdrv1.f
!     Simple program to illustrate the idea of reverse communication
!     in regular mode for a standard symmetric eigenvalue problem.
!     Information of the orginal version of this program:
!     +  Author:
!           Richard Lehoucq
!           Danny Sorensen
!           Chao Yang
!           Dept. of Computational &
!           Applied Mathematics
!           Rice University
!           Houston, Texas
!     +  FILE: sdrv1.F   SID: 2.4   DATE OF SID: 4/22/96   RELEASE: 2
!
!-----
      subroutine arpack_exec_dp  ( 
     $           nrow, neiv, ncv, lworkl,
     $           d, v, resid, workd, workl, choose,
     $           nm, am, jm, g_cons, g_is_matrix, 
     $           flag_wrtout, tole )
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   nrow, neiv, ncv, lworkl
      real(wrp)
     $   d(ncv,2), v(nrow,ncv), resid(nrow), 
     $   workd(3*nrow), workl(lworkl), 
     $   tole
      logical
     $   choose(ncv)
      integer(wip)
     $   nm
      real(wrp)
     $   am(*)
      integer(wip)
     $   jm(*)
      real(wrp)
     $   g_cons  
      logical 
     $   g_is_matrix, flag_wrtout 
!
!     LOCAL VARIABLES: 
      integer(wip)      
     $   iparam(11), ipntr(11)
      integer      
     $   ido, icount 
      integer(wip) 
     $   ldv, info, j, ierr, nconv, maxitr, ishfts, mode
      real(wrp)
     $   tol, sigma
      logical
     $   first, rvec
      character         
     $   bmat*1, which*2
!
!     CONSTANTS: 
      real(wrp), parameter :: 
     &   zero = 0.0_wrp 
!
!     DEPENDENCES: 
!
!     +  BLAS & LAPACK routines:
      real(wrp) 
     &   dpnrm2
      external          
     $   dpnrm2, 
     $   dpaxpy
!
!     +  ARPACK: 
      external 
     $   dpsaupd,
     $   dpseupd
!
      external 
     $   matvec_1_dp, 
     $   matvec_2_dp  
!
      intrinsic         
     $   abs
!
!     EXECUTABLE STATEMENTS: 
!
      ldv   = nrow 
      bmat  = 'I'
      which = 'SA'    
      tol   = tole
      ido   = 0
      info  = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 500
      mode   = 1
!
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
      icount = 0 
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call dpsaupd ( 
     $        ido, bmat, nrow, which, neiv, tol, resid, 
     $        ncv, v, ldv, iparam, ipntr, workd, workl,
     $        lworkl, info )

!
      if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
!
         if ( g_is_matrix ) then  
            call matvec_2_dp (
     $           nrow, workd(ipntr(1)), am, jm, 
     $           workd(ipntr(2)) )
         else 
            call matvec_1_dp (
     $           nrow, workd(ipntr(1)), am, jm, g_cons, 
     $           workd(ipntr(2)) )
         endif 
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
         icount = icount + 1 
         goto 10
!
      endif 
! 
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | ERROR message, check the |
!        | documentation in DNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check the documentation of _saupd'
         print *, ' icount =', icount 
         print *, ' '
!
      else 
!
!        %-------------------------------------------%
!        | No fatal ERRORs occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call dpseupd ( 
     $        rvec, 'All', choose, d, v, ldv, sigma, bmat,
     $        nrow, which, neiv, tol, resid, ncv, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0 ) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DSEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd. '
            print *, ' icount =', icount 
            print *, ' '
!
         else 
!
            first  = .true.
            nconv  = iparam(5)
!             
            if ( flag_wrtout )
     $         write(*,*) 'nconv  = iparam(5) =', nconv 
!             
            do j = 1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A@x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
!              Use WORKD instead of AX.
!
!              workd := A@x   
               if ( g_is_matrix ) then  
                  call matvec_2_dp ( 
     $                 nrow, v(1,j), am, jm, workd(1) )
               else 
                  call matvec_1_dp ( 
     $                 nrow, v(1,j), am, jm, g_cons, workd(1) )
               endif 
!
!              workd := workd - lambda*x = A@x - lambda*x  
               call dpaxpy ( 
     $              nrow, -d(j,1), v(1,j), 1_wip, workd(1), 1_wip )
!
!              d(j,2) := ||workd|| = ||A@x - lambda*x|| 
               d(j,2) = dpnrm2 ( 
     $                  nrow, workd(1), 1_wip )
!
!              d(j,2) := ||A@x - lambda*x|| / |lambda|
               d(j,2) = d(j,2) / abs( d(j,1) )
!               
            enddo 
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
            if ( flag_wrtout ) then 
               write(*,'(/)')
               write(*,100) 
     $         ' ', 'Ritz values', 'Rel. Residuals'
               write(*,100) 
     $         'k', 'lam(k)', '||H@v(k) - lam(k)*v(k)|| / |lam(k)|'
               do j = 1, nconv
                  write(*,200) j, d(j,1), d(j,2) 
               enddo
            endif 
!
         endif
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit' //
     $               ' Arnoldi update, try increasing NCV.'
            print *, ' '
         endif      
!
         if ( flag_wrtout ) then 
            print *, ' '
            print *, ' _SDRV1 '
            print *, ' ====== '
            print *, ' '
            print *, ' The order of the matrix is ', nrow
            print *, ' The number of Ritz values requested is ',neiv
            print *, ' The number of Arnoldi vectors generated',
     $               ' (NCV) is ', ncv
            print *, ' What portion of the spectrum: ', which
            print *, ' The number of converged Ritz values is ', 
     $                 nconv 
            print *, ' The number of Implicit Arnoldi update',
     $               ' iterations taken is ', iparam(3)
            print *, ' The number of OP*x is ', iparam(9)
            print *, ' The convergence criterion is ', tol
            print *, ' The number of calls y=A@x: icount =', icount 
            print *, ' '
         endif 
!
      endif
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue
!
      if ( flag_wrtout )
     $   write(*,'("Done. Exit ARPACK.")')
      return 
 100  format( A9, 1x, a15, 2x, a36 )
 200  format( i9, ':', E15.5, 2x, E36.5 )
      end subroutine arpack_exec_dp 
!-----------------------------------------------------------------------
      subroutine matvec_1_dp ( 
     $           n, x, am, jm, g, y ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      integer(wip), intent(in) ::
     $   n
      real(wrp), intent(in) ::
     $   x(*)
      real(wrp), intent(in) ::
     $   am(*)
      integer(wip), intent(in) ::
     $   jm(*)
      real(wrp), intent(in) :: 
     $   g
      real(wrp), intent(out) ::
     $   y(*) 
!
!     LOCAL VARIABLES:
      integer(wip)
     $   i, k, j 
!
!     EXECUTABLE STATEMENTS: 
!
      do i = 1, n
         y(i) = am(i)*x(i)
      enddo
      do i = 1, n
         do k = jm(i), jm(i+1)-1
            j = jm(k) 
            y(i) = y(i) + g*x(j)
            y(j) = y(j) + g*x(i)
         enddo 
      enddo
      return 
      end subroutine matvec_1_dp 
!-----------------------------------------------------------------------
      subroutine matvec_2_dp ( 
     $           n, x, am, jm, y ) 
!
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      integer(wip), intent(in) ::
     $   n
      real(wrp), intent(in) ::
     $   x(*)
      real(wrp), intent(in) ::
     $   am(*)
      integer(wip), intent(in) ::
     $   jm(*)
      real(wrp), intent(out) ::
     $   y(*) 
!
!     LOCAL VARIABLES:
      real(wrp) 
     $   h_ij  
      integer(wip)
     $   i, k, j 
!
!     EXECUTABLE STATEMENTS: 
!
      do i = 1, n
         y(i) = am(i)*x(i)
      enddo
      do i = 1, n
         do k = jm(i), jm(i+1)-1
            j = jm(k) 
            h_ij = am(k)
            y(i) = y(i) + h_ij*x(j)
            y(j) = y(j) + h_ij*x(i)
         enddo 
      enddo
      return 
      end subroutine matvec_2_dp 
!-----------------------------------------------------------------------
      subroutine copy_s_to_d_symmsr_dp ( 
     $           n, am, jm, g, g_is_mat, bm ) 
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   n
      real(wrp)
     $   am(*) 
      integer(wip)
     $   jm(*)
      real(wrp)  
     $   g
      logical
     $   g_is_mat 
      real(wrp)
     $   bm(n,n) 
!
!     DEPENDENCES: 
      external
     $   copy_s_to_d_symmsr_1_dp, 
     $   copy_s_to_d_symmsr_2_dp
!
!     EXECUTABLE STATEMENTS: 
!
      if ( g_is_mat ) then 
         call copy_s_to_d_symmsr_2_dp ( 
     $        n, am, jm, bm ) 
      else 
         call copy_s_to_d_symmsr_1_dp ( 
     $        n, am, jm, g, bm ) 
      endif 
      return
      end subroutine copy_s_to_d_symmsr_dp
!-----------------------------------------------------------------------
      subroutine copy_s_to_d_symmsr_1_dp ( 
     $           n, am, jm, g, bm ) 
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   n
      real(wrp)
     $   am(*) 
      integer(wip)
     $   jm(*)
      real(wrp)  
     $   g
      real(wrp)
     $   bm(n,n) 
!
!     LOCAL VARIABLES:
      integer(wip)
     $   ir, k, jc
!
!     EXECUTABLE STATEMENTS: 
!
      bm = 0.0_wrp 
      do ir = 1, n
         bm(ir,ir) = am(ir)
      enddo 
      do ir = 1, n
         do k  = jm(ir), jm(ir+1)-1 
            jc = jm(k)
            bm(ir,jc) = g
            bm(jc,ir) = g
         enddo 
      enddo 
      return
      end subroutine copy_s_to_d_symmsr_1_dp
!-----------------------------------------------------------------------
      subroutine copy_s_to_d_symmsr_2_dp ( 
     $           n, am, jm, bm ) 
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   n
      real(wrp)
     $   am(*) 
      integer(wip)
     $   jm(*)
      real(wrp)  
     $   g
      real(wrp)
     $   bm(n,n) 
!
!     LOCAL VARIABLES:
      integer(wip)
     $   ir, k, jc
      real(wrp)
     $   h_ij 
!
!     EXECUTABLE STATEMENTS: 
!
      bm = 0.0_wrp 
      do ir = 1, n
         bm(ir,ir) = am(ir)
      enddo 
      do ir = 1, n
         do k  = jm(ir), jm(ir+1)-1 
            jc = jm(k)
            h_ij = am(k)
            bm(ir,jc) = h_ij 
            bm(jc,ir) = h_ij 
         enddo 
      enddo
      return
      end subroutine copy_s_to_d_symmsr_2_dp
!-----------------------------------------------------------------------
      subroutine lapack_exec_dp ( 
     $           nrow, bmat, eival, lwork, work, 
     $           nm, am, jm, flag_wrtout )
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   nrow
      real(wrp)
     $   bmat ( nrow, nrow )
      real(wrp)
     $   eival ( nrow )
      integer(wip)
     $   lwork
      real(wrp)
     $   work ( lwork )
      integer(wip)
     $   nm
      real(wrp)
     $   am ( nm )
      integer(wip)
     $   jm ( nm )
      logical
     $   flag_wrtout
!
!     DEPENDENCES: 
      real(wrp) 
     &   dpnrm2
      external          
     $   dpsyev,
     $   dpnrm2, 
     $   dpaxpy,
     $   dpsymv 
!
!     LOCAL VARIABLES: 
      integer
     $   info
      real(wrp), allocatable :: 
     $   amat(:,:), workd(:), d(:)
      integer(wip) 
     $   j
!
!     EXECUTABLE STATEMENTS: 
!
      if ( flag_wrtout )
     $   write(*,'("Inside LAPACK, LWORK =", i12)') lwork 
      allocate ( amat(nrow,nrow), workd(nrow), d(nrow) ) 
      amat = bmat 
      call dpsyev ( 
     $     'V','U', nrow, bmat(1,1), nrow, eival(1), 
     $     work(1), lwork, info )
!
!     Compute the residual norm ||A@x - lambda*x||/|lambda|: 
      do j = 1, nrow 
!
!        workd := A@x
         call dpsymv ( 
     $        'U', nrow, 1.0_wrp, amat, nrow, bmat(1,j), 1_wip,
     $        0.0_wrp, workd(1), 1_wip )
!
!        workd := workd - lambda*x = A@x - lambda*x  
         call dpaxpy ( 
     $        nrow, -eival(j), bmat(1,j), 1_wip, workd(1), 1_wip )
!
!        d(j) := ||workd|| = ||A@x - lambda*x|| 
         d(j) = dpnrm2 ( 
     $          nrow, workd(1), 1_wip )
!
!        d(j) := ||A@x - lambda*x|| / |lambda|
         d(j) = d(j) / abs( eival(j) )
!               
      enddo 
      if ( flag_wrtout ) then 
         write(*,'(/)')
         write(*,10) 
     $      ' ', 'Eigenvalues', 'Rel. Residuals'
         write(*,10) 
     $      'k', 'lam(k)', '||H@v(k) - lam(k)*v(k)|| / |lam(k)|'
         do j = 1, nrow
            write(*,20) j, eival(j), d(j) 
         enddo
      endif 
      deallocate ( amat, workd, d ) 
      if ( flag_wrtout )     
     $   write(*,'(/,"Done. Exit LAPACK.")')
      return 
 10   format( A9, 1x, a15, 2x, a36 )
 20   format( i9, ':', E15.5, 2x, E36.5 )
      end subroutine lapack_exec_dp 
!-----------------------------------------------------------------------
!  Ref:
!     Higham, Nicholas J. "The accuracy of floating point summation." 
!     SIAM Journal on Scientific Computing 14, no. 4 (1993): 783-799.
!
      subroutine higham_sum_dp ( s, x, e )
      implicit none 
      include 'epcode_inc_dp.h'     
!     
!     ARGUMENTS:  
      real(wrp)      
     $   s, x, e 
!
!     LOCAL VARIABLES: 
      real(wrp)      
     $   t, y, r 
!
!     EXECUTABLE STATEMENTS:  
!
      t = s            ! t:   _____t1____.____t2___
      y = x + e        ! y:              .____y1___ ____y2___
      s = t + y        ! s:   _____t1____.__t2+y1__
      r = t - s        ! r:              .___-y1___
      e = r + y        ! e:              .          ____y2___
      return 
      end subroutine higham_sum_dp 
!-----------------------------------------------------------------------
!***********************************************************************
!     real(kind=16)
!-----
      subroutine ep_gstate_1a_qp ( 
     $           nome, npar, e, g, neiv, 
     $           prefix, tole ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar 
      real(wrp), intent(in) :: 
     $   e(nome), g, tole
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     LOCAL VARIABLE:
      real(wrp)
     $   g_sca(1,1)
!
!     DEPENDENCES:
      external 
     $   ep_gstate_main_qp  
!
!     EXECUTABLE STATEMENTS: 
!
      g_sca(1,1) = g 
      call ep_gstate_main_qp ( 
     $     nome, npar, e, 1, g_sca, neiv, 
     $     prefix, tole ) 
      return 
      end subroutine ep_gstate_1a_qp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_1b_qp ( 
     $           nome, npar, e, g, neiv, 
     $           prefix ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar 
      real(wrp), intent(in) :: 
     $   e(nome), g
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_1a_qp 
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_1a_qp ( 
     $     nome, npar, e, g, neiv, 
     $     prefix, 0.0_wrp ) 
      return 
      end subroutine ep_gstate_1b_qp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_2a_qp ( 
     $           nome, npar, e, g, neiv, 
     $           prefix, tole ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar 
      real(wrp), intent(in) :: 
     $   e(nome), g(nome,nome), tole
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_main_qp  
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_main_qp ( 
     $     nome, npar, e, nome, g, neiv, 
     $     prefix, tole ) 
      return 
      end subroutine ep_gstate_2a_qp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_2b_qp ( 
     $           nome, npar, e, g, neiv, 
     $           prefix ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar 
      real(wrp), intent(in) :: 
     $   e(nome), g(nome,nome)
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_2a_qp 
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_2a_qp ( 
     $     nome, npar, e, g, neiv, 
     $     prefix, 0.0_wrp ) 
      return 
      end subroutine ep_gstate_2b_qp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_3a_qp ( 
     $           nome, npar, e, ng, g, neiv, 
     $           prefix, tole ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar, ng 
      real(wrp), intent(in) :: 
     $   e(nome), g(ng,ng), tole
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_main_qp  
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_main_qp ( 
     $     nome, npar, e, ng, g, neiv, 
     $     prefix, tole ) 
      return 
      end subroutine ep_gstate_3a_qp 
!-----------------------------------------------------------------------
      subroutine ep_gstate_3b_qp ( 
     $           nome, npar, e, ng, g, neiv, 
     $           prefix ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar, ng 
      real(wrp), intent(in) :: 
     $   e(nome), g(ng,ng)
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     DEPENDENCES:
      external 
     $   ep_gstate_3a_qp 
!
!     EXECUTABLE STATEMENTS: 
!
      call ep_gstate_3a_qp ( 
     $     nome, npar, e, ng, g, neiv, 
     $     prefix, 0.0_wrp ) 
      return 
      end subroutine ep_gstate_3b_qp 
!-----------------------------------------------------------------------
!
!     The core of EP (real*16) 
!
      subroutine ep_gstate_main_qp ( 
     $           nome, npar, e, ng, g, neiv, 
     $           prefix, tole ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer, intent(in) :: 
     $   nome, npar, ng 
      real(wrp), intent(in) :: 
     $   e(nome), g(ng,ng), tole
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(inout) :: 
     $   prefix
!
!     LOCAL VARIABLES AND CONSTANTS:
      real(wrp) 
     $   g_diag(63), g_matr(63,63)
      integer(wip)
     $   binom_tab(0:63,0:63)
      common /ep_gstate_main_shared/ binom_tab
      logical, parameter ::
     $   flag_screen = .true.,   
     $   flag_wrtres = .true.,  
!!   $   flag_wrtout = .true.,   
     $   flag_wrtout = .false.,  
!!   $   flag_wrtsiz = .true.,  
     $   flag_wrtsiz = .false., 
     $   flag_memest = .true., 
     $   flag_wrtnsj = .false.,  
     $   flag_wrtmat = .false., 
     $   flag_wrtevv = .false., 
     $   flag_wrtocc = .false.
      integer(wip)
     $   lsj, meiv 
      integer(i1b), allocatable :: 
     $   njjb(:), nsjb(:), nb1(:), nb2(:)
      integer(wip), allocatable :: 
     $   nsj(:), njj(:), mjp(:)
      integer(wip), allocatable :: 
     $   jah(:)
      real(wrp), allocatable :: 
     $   aah(:)
      real(wrp), allocatable :: 
     $   occ(:,:), ecc(:,:)
      real(wrp)
     $   tmp 
      integer
     $   mjp_len, ljp, nps, nmns, lb1, i, i1, i2, isum, 
     $   k, l, nt0, iox, ns 
      integer(wip)
     $   nrow, ncol, ir, ic, isj, ijj, ljj,
     $   nnz, iele, nele, jos, state_i, state_f,  
     $   njir, njic, ntmp, nde
      integer(wip)      
     $   ncv, lawork, lworkl, 
     $   io_resid, io_workd, io_eivec, io_eival, io_workl
      real(wrp), allocatable :: 
     $   awork(:)
      logical, allocatable :: 
     $   choose(:)
      character
     $   chn*5, chpack*30, chfmt*50, 
     $   fort00*90, fort01*90, fort02*90, fort03*90, 
     $   fort04*90, fort05*90, header*200
      integer ::
     $   iou00, iou01, iou02, iou03, iou04, iou05 
      integer(i1b)
     $   k1, k2  
      logical 
     $   g_is_matrix, use_arpack  
      integer(wip)      
     $   num_of_bytes_in_total, num_of_bytes_in_arrays,
     $   num_of_bytes_in_locvar
      integer(8) 
     $   iram 
!
!     DEPENDENCES:
      intrinsic :: 
     $   mod, min, max, sum, real, int, trim, adjustl, len_trim 
      integer(wip), external ::
     $   binom
      integer(wip), external ::
     $   npairing 
      integer(wip), external ::
     $   hashs
      external ::
     $   prt_nsj_njj, 
     $   deccon, 
     $   gens, genp, gind 
      character(len=80), external ::
     $   cparams0 
      external ::
     $   arpack_exec_qp, 
     $   lapack_exec_qp,  
     $   copy_s_to_d_symmsr_qp,
     $   higham_sum_qp 
      integer(wip), external ::
     $   wilaenv
      external 
     $   get_ram_in_bytes,
     $   ep_getmem_qp
!
!     SHARED PARAMETERS:
      logical 
     $   flag_prep, 
     $   flag_wrtinp 
      common /epcode_shared_flags/ 
     $   flag_prep, 
     $   flag_wrtinp 
!
!     DEFAULT VALUES:
      data                      
     $   flag_prep   /.false./,
!!   $   flag_wrtinp /.true./  
     $   flag_wrtinp /.false./ 
!
!     EXECUTABLE STATEMENTS: 
!
      if ( nome.lt.1 .or. nome.gt.63 ) 
     $   stop "ERROR: NOME out of [1,63]"
      if ( npar.lt.1 .or. npar.gt.nome*2 ) 
     $   stop "ERROR: NPAR out of [1,2*NOME]"
      if ( .not. (ng.eq.1 .or. ng.eq.nome) ) 
     $   stop "ERROR: NG /= 1 and NG /= NOME"
      g_is_matrix = ( ng .ne. 1 ) 
      if ( g_is_matrix ) then 
         g_diag = 0
         g_matr = 0 
         do k = 1, nome
            do l = 1, nome
               g_matr(l,k) = -g(l,k) 
            enddo
            g_diag(k) = -g(k,k)
         enddo
      else 
         g_diag(1) = -g(1,1)
      endif
      binom_tab = 0
      do k = 0, nome 
      do i = 0, k
         binom_tab ( k, i ) = binom ( k, i )
      enddo
      enddo 
      ns   = mod(npar,2)
      nmns = nome - ns 
      nps  = npar/2 - ns/2 
      mjp_len = nps * ( nmns - nps )
      lsj  = binom_tab ( nome, ns )
      if ( lsj .lt. 1 ) stop "ERROR: lsj < 1" 
      ljj  = binom_tab ( nmns, nps )
      if ( ljj .lt. 1 ) stop "ERROR: ljj < 1" 
      meiv = min( ljj, neiv )
      if ( meiv .lt. 1 ) stop "ERROR: no. of eigenvalues < 1" 
      allocate ( njj(ljj) ) 
      njj = 0 
      call gens ( nmns, nps, ljj, njj )
      ncol = ljj 
      nrow = ncol 
      allocate ( occ(nome,meiv), ecc(nome,meiv) )
      occ = 0.0_wrp 
      ecc = 0.0_wrp 
!!    if ( flag_wrtout ) 
!!   $   write(*,'(1x,a15,1x,i20)') 
!!   $   'Ncol = Nrow =', nrow 
      iele = 0 
      do ir = 1_wip, ljj
         iele = iele + npairing ( njj(ir), nmns )
      enddo 
      nele = iele + nrow
      nnz  = 2*( nele - nrow ) + nrow
!!    if ( flag_wrtout ) 
!!   $   write(*,
!!   $   '(/,a,/,1x,a15,1x,i20,/,3(1x,a15,1x,i20,1x,a,/))') 
!!   $   'MATRIX INFORMATION:',
!!   $   'Ncol = Nrow =', nrow, 
!!   $   'Nele =', nele, 
!!   $   '(No. Non-Zero elements of the half of H)', 
!!   $   'NNZ =',  nnz, 
!!   $   '(No. Non-Zero elements of H)',
!!   $   'LJJ**2 - NNZ =', ljj*int(ljj,kind=wip) - nnz, 
!!   $   '(No. Zeros of H)'
      nele = nele + 1 
      if ( flag_prep ) goto 201  
      allocate ( jah(nele) )
      if ( g_is_matrix ) then 
         allocate ( aah(nele) )
      else 
         allocate ( aah(nrow) )
      endif 
      aah = 0.0_wrp  
      jah = 0_wip  
      iele = ljj + 1  
      allocate ( mjp ( mjp_len ) )
      if ( g_is_matrix ) then 
         do ir = 1, ljj
            call genp ( nmns, nps, njj(ir), mjp_len, mjp, ljp ) 
            do k = 1, ljp
               ic = hashs ( nmns, nps, mjp(k) ) 
               iele = iele + 1 
               if ( k .eq. 1 ) jah(ir) = iele 
               jah(iele) = ic 
               call gind ( njj(ir), njj(ic), k1, k2 ) 
               aah(iele) = g_matr ( k1, k2 )
            enddo 
         enddo 
      else 
         do ir = 1, ljj
            call genp ( nmns, nps, njj(ir), mjp_len, mjp, ljp ) 
            do k = 1, ljp
               ic = hashs ( nmns, nps, mjp(k) ) 
               iele = iele + 1 
               if ( k .eq. 1 ) jah(ir) = iele 
               jah(iele) = ic 
            enddo 
         enddo 
      endif 
      deallocate ( mjp )
      ic = ljj 
      do while ( jah(ic) .eq. 0 ) 
         ic = ic - 1 
      enddo 
      if ( ic .eq. ljj ) then 
         jah(ic+1) = jah(ic)
      else 
         jah(ic+1) = jah(ic) + 1
         do ir = ic+2, ljj+1
            jah(ir) = jah(ir-1) 
         enddo 
      endif 
!
!               PREFIX
!               |
!     header = "EP10_n10p11g0600s01" for NOME=10, NPAR=11, G=0.6, NS=1
!                    |  |  |    |        
!              NOME=10  |  G=0.6|    
!                 NPAR=11       |
!                               NS=1
!
      header = ''
      if ( len_trim(adjustl(prefix)) .gt. 0 )
     $   header = trim(adjustl(header)) // trim(adjustl(prefix))
      header = trim(adjustl(header)) //
     $         trim(cparams0 ( nome, npar, real(g(1,1)), ns) )
!
      if ( flag_wrtres ) then 
         fort00 = trim(adjustl(header))//'_result.txt' 
         open( newunit=iou00, file=trim(adjustl(fort00)), 
     $         status='replace' )
         rewind(iou00)
      endif 
      if ( flag_wrtmat ) then 
         fort01 = trim(adjustl(header))//'_iajaaa.bin'
         open( newunit=iou01, file=trim(adjustl(fort01)), 
     $   status='replace', form='unformatted' )
         rewind(iou01)
      endif 
      if ( flag_wrtnsj ) then 
         fort02 = trim(adjustl(header))//'_nsjnjj.txt'
         open( newunit=iou02, file=trim(adjustl(fort02)), 
     $         status='replace' )
         rewind(iou02)
      endif 
      if ( flag_wrtsiz ) then 
         fort03 = trim(adjustl(header))//'_dtsize.txt' 
         open( newunit=iou03, file=trim(adjustl(fort03)), 
     $         status='replace' )
         rewind(iou03)
      endif 
      if ( flag_wrtevv ) then  
         fort04 = trim(adjustl(header))//'_eivave.bin'
         open( newunit=iou04, file=trim(adjustl(fort04)), 
     $         status='replace', form='unformatted' )
         rewind(iou04)
      endif 
      if ( flag_wrtocc ) then 
         fort05 = trim(adjustl(header))//'_eivocc.bin' 
         open( newunit=iou05, file=trim(adjustl(fort05)), 
     $      status='replace', form='unformatted' )
         rewind(iou05)
         prefix = trim(adjustl(fort05))
      endif 
!
      allocate ( 
     $   njjb(nome), nsjb(nome), nb1(nome), nb2(nome) )
      njjb = 0
      nsjb = 0
      nb1  = 0
      nb2  = 0 
      if ( flag_wrtocc ) then 
         write(iou05) ns, nome, meiv
      endif 
 201  continue 
      use_arpack = ( ljj .gt. 100 ) 
      if ( use_arpack ) then 
         chpack = 'ARPACK'
         ncv    = max( meiv*2, meiv + 10 )    
         lworkl = ncv * ( ncv + 8 )
         lawork = ncv*2 + nrow*ncv + nrow + 3*nrow + lworkl
         io_eival = 0 
         io_eivec = io_eival + ncv*2
         io_resid = io_eivec + ncv*nrow 
         io_workd = io_resid + nrow 
         io_workl = io_workd + nrow*3
      else
         chpack = 'LAPACK'
         ncv = 0 
         lworkl = wilaenv ( 
     $      1, 'WRSYTRD', 'U', nrow, -1_wip, -1_wip, -1_wip )
         lworkl = max( 1_wip, ( lworkl + 2 )*nrow )
         io_eivec = 0
         io_eival = io_eivec + nrow*nrow 
         io_workl = io_eival + nrow 
         lawork = nrow*nrow + nrow + lworkl  
      endif 
      if ( flag_prep ) goto 202 
      if ( use_arpack ) then 
         allocate ( awork ( lawork ) ) 
         allocate ( choose ( ncv ) )
         awork = 0.0_wrp  
         choose = .false.
      else
         allocate ( awork ( lawork ) ) 
         awork = 0.0_wrp 
         call copy_s_to_d_symmsr_qp  ( 
     $        nrow, aah, jah, g_diag(1), g_is_matrix, 
     $        awork(io_eivec+1) )
      endif 
 202  continue 
      if ( flag_wrtsiz ) then 
         iox = iou03
      else
         iox = 123  
         goto 506 
      endif 
 505  continue 
      write(iox,'(72("-"))') 
!!    write(iox,'( 6(/,1x,a15,1x,i20) )') 
!!   $   'NOME =', nome, 
!!   $   'NPAR =', npar,
!!   $   'NS   =', ns,
!!   $   'LSJ  =', lsj, 
!!   $   'LJJ  =', ljj, 
!!   $   'NEIV =', meiv 
      write(iox, '( 3(/,1x,a15,1x,i20))') 
     $   'NOME =', nome, 
     $   'NPAR =', npar,
     $   'NEIV =', meiv
      write(iox, '( 1x,a15,1x,e20.13 )') 
     $   'TOLE =', tole  
      write(iox, '( 1x,a15,1x, a )') 
     $   'PREFIX =', trim(prefix)
      if ( ng .eq. 1 ) then 
         write(iox, '( 1x,a15,1x,e20.13 )') 
     $   'G =', g(1,1)  
      else
         write(iox,'( 1x, a15 )') 'G ='
         chn = '' 
         write(chn,fmt='(I3)') nome 
         chfmt = '(I5, ":",2x,'//
     $      trim(adjustl(chn))//'(1X,ES9.2))'
         do k = 1, ng 
            write(iox,fmt=chfmt) k, g(k,:)
         enddo 
      endif
      write(iox,'( 1x, a15 )') 'SPE ='
      chfmt = '(I5, ":", 3x, ES9.2 )'
      do k = 1, nome 
         write(iox,fmt=chfmt) k, e(k)
      enddo 
      write(iox,
     $   '(/,a,/,1x,a15,1x,i20,/,3(1x,a15,1x,i20,1x,a,/))') 
     $   'MATRIX INFORMATION:',
     $   'Ncol = Nrow =', nrow, 
     $   'Nele =', nele-1, 
     $   '(No. Non-Zero elements of the half of H)', 
     $   'NNZ =',  nnz, 
     $   '(No. Non-Zero elements of H)',
     $   'LJJ**2 - NNZ =', ljj*int(ljj,kind=wip) - nnz, 
     $   '(No. Zeros of H)'
  506 continue 
      if (iox .ne. 6 .and. flag_screen ) then 
         iox = 6 
         goto 505 
      endif
      call ep_getmem_qp (
     $     nome, npar, ng, neiv, prefix, 
     $     ljj, nele, mjp_len, lsj, lawork, ncv, 
     $     num_of_bytes_in_arrays, 
     $     num_of_bytes_in_locvar, 
     $     num_of_bytes_in_total ) 
      if ( flag_wrtsiz ) then 
         iox = iou03
      else
         iox = 123  
         goto 406 
      endif 
  405 continue 
      write(iox,'(/,A)') 'COMPUTER MEMORY REQUIREMENT:' 
      write(iox,'( 3x, I20, A10, 2(E17.10, A) )') 
     $    num_of_bytes_in_total, ' (bytes) =', 
     $    num_of_bytes_in_total * dble(2)**(-20), ' (MB) =',
     $    num_of_bytes_in_total * dble(2)**(-30), ' (GB)'
  406 continue 
      if (iox .ne. 6 .and. flag_screen ) then 
         iox = 6 
         goto 405 
      endif
  400 continue
      if ( flag_prep ) goto 800  
      isj = 1 
      nsjb = 0 
      if (mod(npar,2) .ne. 0) nsjb ( (npar-1)/2 + 1 ) = 1 
      k = 0 
      do i = 1, nome
         if ( nsjb(i) .ne. 0 ) cycle 
         k = k + 1
         nb1(k) = i
      enddo
      lb1 = k 
      if ( g_is_matrix ) then 
         do iele = 1_wip, ljj
            njjb = nsjb 
            call deccon ( njj(iele), nome, nb2 )
            do i = 1, lb1  
               njjb( nb1(i) ) = nb2(i)*2 
            enddo
            aah(iele) = sum ( e*njjb + 0.25_wrp*g_diag(1:nome)*
     $                        (njjb - nsjb)*(4 - njjb - nsjb) ) 
         enddo
      else 
         do iele = 1_wip, ljj
            njjb = nsjb 
            call deccon ( njj(iele), nome, nb2 )
            do i = 1, lb1  
               njjb( nb1(i) ) = nb2(i)*2 
            enddo
            aah(iele) = sum ( e*njjb + 0.25_wrp*g_diag(1)*
     $                        (njjb - nsjb)*(4 - njjb - nsjb) ) 
         enddo
      endif 
      if ( flag_wrtmat ) then  
         write(iou01) nrow, nele 
         write(iou01) jah
         write(iou01) aah
      endif 
      if ( use_arpack ) then 
         awork = 0.0_wrp  
         choose = .false.
         call arpack_exec_qp  ( 
     $        nrow, meiv, ncv, lworkl,
     $        awork(io_eival+1), awork(io_eivec+1), 
     $        awork(io_resid+1), awork(io_workd+1), 
     $        awork(io_workl+1), choose, 
     $        nele, aah, jah, g_diag(1), g_is_matrix, 
     $        flag_wrtout, tole )
      else
         jos = io_eivec - nrow 
         do ir = 1_wip, nrow 
            awork ( jos + ir*(1+nrow) ) = aah (ir) 
         enddo 
         call lapack_exec_qp  ( 
     $        nrow, 
     $        awork(io_eivec+1), awork(io_eival+1),  
     $        lworkl, awork(io_workl+1),  
     $        nele, aah, jah, flag_wrtout )
      endif 
      if ( flag_wrtevv ) then  
         write(iou04) meiv, ljj
         jos = io_eivec 
         do ic = 1_wip, meiv
            write(iou04) ic, awork(io_eival+ic)
            write(iou04) awork(jos+1:jos+ljj)
            jos = jos + ljj
         enddo  
      endif 
      nde = 2**ns 
      occ = 0.0_wrp 
      ecc = 0.0_wrp 
      do ir = 1_wip, ljj 
         njjb = nsjb 
         call deccon ( njj(ir), nome, nb2 )
         do i = 1, lb1  
            njjb( nb1(i) ) = nb2(i)*2 
         enddo
         jos = io_eivec + ir
         do ic = 1_wip, meiv
            tmp = awork ( jos )**2 * 0.5_wrp 
            jos = jos + nrow 
            do k = 1, nome
               call higham_sum_qp ( 
     $              occ(k,ic), tmp*njjb(k), ecc(k,ic) )  
            enddo 
         enddo 
      enddo
      chn = '' 
      write(chn,fmt='(I5)') nome 
      if ( flag_wrtres ) then
         select case ( wrp )
         case(kind(1.0e0))
            write(iou00,
     $      '("#",A9,1X,A22,1X,A22)') 
     $      'i', 'Eigenvalue(i)',  
     $      '( Occ(j,i), j=1,'//trim(adjustl(chn))//' )'
            chfmt = '(I10,1X,ES22.8,'//
     $            trim(adjustl(chn))//'(1X,ES22.8))'
         case(kind(1.0d0))
            write(iou00,
     $      '("#",A9,1X,A22,1X,A22)') 
     $      'i', 'Eigenvalue(i)',  
     $      '( Occ(j,i), j=1,'//trim(adjustl(chn))//' )'
            chfmt = '(I10,1X,ES22.15,'//
     $            trim(adjustl(chn))//'(1X,ES22.15))'
         case(kind(1.0d0)*2)
            write(iou00,
     $      '("#",A9,1X,A40,1X,A40)') 
     $      'i', 'Eigenvalue(i)',  
     $      '( Occ(j,i), j=1,'//trim(adjustl(chn))//' )'
            chfmt = '(I10,1X,ES40.32,'//
     $            trim(adjustl(chn))//'(1X,ES40.32))'
         end select 
         k = len_trim(chfmt)
         do ic = 1_wip, meiv 
            write(iou00,fmt=chfmt(1:k))
     $         ic, awork(io_eival+ic), occ(1:nome,ic)
         enddo
      endif 
      if ( flag_wrtocc ) then 
         write(iou05) awork(io_eival+1:io_eival+meiv), occ
      endif 
      if ( flag_wrtnsj ) then  
         allocate ( nsj(lsj) )  
         nsj = 0
         call gens ( nome, ns, lsj, nsj )
         call prt_nsj_njj (iou02, nome, lsj, nsj, ljj, njj)
         deallocate ( nsj )
      endif 
      if ( flag_screen ) then  
         write(*,'(/,a)') 'CALCULATIONS COMPLETED' 
         if (  flag_wrtres .or. 
     $         flag_wrtmat .or.
     $         flag_wrtnsj .or. 
     $         flag_wrtsiz .or. 
     $         flag_wrtevv .or. 
     $         flag_wrtocc )   
     $      write(*,'(a)') 'Check:' 
      endif 
      if ( flag_wrtres .and. flag_screen )   
     $   write(*,11) '', adjustl(fort00), '(ascii)'
      if ( flag_wrtmat .and. flag_screen )   
     $   write(*,11) '', adjustl(fort01), '(binary)'
      if ( flag_wrtnsj .and. flag_screen )   
     $   write(*,11) '', adjustl(fort02), '(ascii)'
      if ( flag_wrtsiz .and. flag_screen )   
     $   write(*,11) '', adjustl(fort03), '(ascii)'
      if ( flag_wrtevv .and. flag_screen )   
     $   write(*,11) '', adjustl(fort04), '(binary)' 
      if ( flag_wrtocc .and. flag_screen )   
     $   write(*,11) '', adjustl(fort05), '(binary)'    
      if ( flag_wrtres ) close (iou00)
      if ( flag_wrtmat ) close (iou01)
      if ( flag_wrtnsj ) close (iou02)
      if ( flag_wrtsiz ) close (iou03)
      if ( flag_wrtevv ) close (iou04)
      if ( flag_wrtocc ) close (iou05)
      if ( flag_wrtout .and. flag_screen ) then 
         write(*,'(/,A)') 
     $   'Eigenvalue and Occupation numbers:' 
         write(chn,'(I5)') 1 + nome
         select case ( wrp )
         case(kind(1.0e0))
            chfmt = '('//trim(adjustl(chn))//'ES16.8)'
         case(kind(1.0d0))
            chfmt = '('//trim(adjustl(chn))//'ES23.15)'
         case(kind(1.0d0)*2)
            chfmt = '('//trim(adjustl(chn))//'ES40.32)'
         end select 
         k = len_trim ( chfmt )
         do ic = 1_wip, meiv 
            write(*,fmt=chfmt(1:k)) 
     $      awork(io_eival+ic), occ(1:nome,ic)
         enddo
      endif 
  800 continue  
      if (allocated(nsj    )) deallocate (nsj    )
      if (allocated(njj    )) deallocate (njj    )
      if (allocated(njjb   )) deallocate (njjb   )
      if (allocated(nsjb   )) deallocate (nsjb   )
      if (allocated(nb1    )) deallocate (nb1    )
      if (allocated(nb2    )) deallocate (nb2    )
      if (allocated(jah    )) deallocate (jah    )
      if (allocated(aah    )) deallocate (aah    )
      if (allocated(occ    )) deallocate (occ    )
      if (allocated(ecc    )) deallocate (ecc    )
      if (allocated(mjp    )) deallocate (mjp    )
      if (allocated(awork  )) deallocate (awork  ) 
      if (allocated(choose )) deallocate (choose )
      return 
 11   format ( A10, 1x, A50, 1x, A10 )
      end subroutine ep_gstate_main_qp 
!-----------------------------------------------------------------------
!
!     This is a copied version of the program published in ARPACK: 
!        arpack96/EXAMPLES/SYM/dsdrv1.f
!     Simple program to illustrate the idea of reverse communication
!     in regular mode for a standard symmetric eigenvalue problem.
!     Information of the orginal version of this program:
!     +  Author:
!           Richard Lehoucq
!           Danny Sorensen
!           Chao Yang
!           Dept. of Computational &
!           Applied Mathematics
!           Rice University
!           Houston, Texas
!     +  FILE: sdrv1.F   SID: 2.4   DATE OF SID: 4/22/96   RELEASE: 2
!
!-----
      subroutine arpack_exec_qp  ( 
     $           nrow, neiv, ncv, lworkl,
     $           d, v, resid, workd, workl, choose,
     $           nm, am, jm, g_cons, g_is_matrix, 
     $           flag_wrtout, tole )
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   nrow, neiv, ncv, lworkl
      real(wrp)
     $   d(ncv,2), v(nrow,ncv), resid(nrow), 
     $   workd(3*nrow), workl(lworkl), 
     $   tole
      logical
     $   choose(ncv)
      integer(wip)
     $   nm
      real(wrp)
     $   am(*)
      integer(wip)
     $   jm(*)
      real(wrp)
     $   g_cons  
      logical 
     $   g_is_matrix, flag_wrtout 
!
!     LOCAL VARIABLES: 
      integer(wip)      
     $   iparam(11), ipntr(11)
      integer      
     $   ido, icount 
      integer(wip) 
     $   ldv, info, j, ierr, nconv, maxitr, ishfts, mode
      real(wrp)
     $   tol, sigma
      logical
     $   first, rvec
      character         
     $   bmat*1, which*2
!
!     CONSTANTS: 
      real(wrp), parameter :: 
     &   zero = 0.0_wrp 
!
!     DEPENDENCES: 
!
!     +  BLAS & LAPACK routines:
      real(wrp) 
     &   qpnrm2
      external          
     $   qpnrm2, 
     $   qpaxpy
!
!     +  ARPACK: 
      external 
     $   qpsaupd,
     $   qpseupd
!
      external 
     $   matvec_1_qp, 
     $   matvec_2_qp  
!
      intrinsic         
     $   abs
!
!     EXECUTABLE STATEMENTS: 
!
      ldv   = nrow 
      bmat  = 'I'
      which = 'SA'    
      tol   = tole
      ido   = 0
      info  = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DNAUPD is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | DNAUPD.                                           |
!     %---------------------------------------------------%
!
      ishfts = 1
      maxitr = 500
      mode   = 1
!
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
!
      icount = 0 
!
 10   continue
!
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take | 
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
         call qpsaupd ( 
     $        ido, bmat, nrow, which, neiv, tol, resid, 
     $        ncv, v, ldv, iparam, ipntr, workd, workl,
     $        lworkl, info )

!
      if (ido .eq. -1 .or. ido .eq. 1) then
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | The user should supply his/her own        |
!           | matrix vector multiplication routine here |
!           | that takes workd(ipntr(1)) as the input   |
!           | vector, and return the matrix vector      |
!           | product to workd(ipntr(2)).               | 
!           %-------------------------------------------%
!
         if ( g_is_matrix ) then  
            call matvec_2_qp (
     $           nrow, workd(ipntr(1)), am, jm, 
     $           workd(ipntr(2)) )
         else 
            call matvec_1_qp (
     $           nrow, workd(ipntr(1)), am, jm, g_cons, 
     $           workd(ipntr(2)) )
         endif 
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call DSAUPD again. |
!           %-----------------------------------------%
!
         icount = icount + 1 
         goto 10
!
      endif 
! 
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | ERROR message, check the |
!        | documentation in DNAUPD. |
!        %--------------------------%
!
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check the documentation of _saupd'
         print *, ' icount =', icount 
         print *, ' '
!
      else 
!
!        %-------------------------------------------%
!        | No fatal ERRORs occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
         rvec = .true.
!
         call qpseupd ( 
     $        rvec, 'All', choose, d, v, ldv, sigma, bmat,
     $        nrow, which, neiv, tol, resid, ncv, v, ldv, 
     &        iparam, ipntr, workd, workl, lworkl, ierr )
!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
         if ( ierr .ne. 0 ) then
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of DSEUPD. |
!           %------------------------------------%
!
            print *, ' '
            print *, ' Error with _seupd, info = ', ierr
            print *, ' Check the documentation of _seupd. '
            print *, ' icount =', icount 
            print *, ' '
!
         else 
!
            first  = .true.
            nconv  = iparam(5)
!             
            if ( flag_wrtout )
     $         write(*,*) 'nconv  = iparam(5) =', nconv 
!             
            do j = 1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A@x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
!              Use WORKD instead of AX.
!
!              workd := A@x   
               if ( g_is_matrix ) then  
                  call matvec_2_qp ( 
     $                 nrow, v(1,j), am, jm, workd(1) )
               else 
                  call matvec_1_qp ( 
     $                 nrow, v(1,j), am, jm, g_cons, workd(1) )
               endif 
!
!              workd := workd - lambda*x = A@x - lambda*x  
               call qpaxpy ( 
     $              nrow, -d(j,1), v(1,j), 1_wip, workd(1), 1_wip )
!
!              d(j,2) := ||workd|| = ||A@x - lambda*x|| 
               d(j,2) = qpnrm2 ( 
     $                  nrow, workd(1), 1_wip )
!
!              d(j,2) := ||A@x - lambda*x|| / |lambda|
               d(j,2) = d(j,2) / abs( d(j,1) )
!               
            enddo 
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
            if ( flag_wrtout ) then 
               write(*,'(/)')
               write(*,100) 
     $         ' ', 'Ritz values', 'Rel. Residuals'
               write(*,100) 
     $         'k', 'lam(k)', '||H@v(k) - lam(k)*v(k)|| / |lam(k)|'
               do j = 1, nconv
                  write(*,200) j, d(j,1), d(j,2) 
               enddo
            endif
!
         endif
!
!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%
!
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit' //
     $               ' Arnoldi update, try increasing NCV.'
            print *, ' '
         endif      
!
         if ( flag_wrtout ) then 
            print *, ' '
            print *, ' _SDRV1 '
            print *, ' ====== '
            print *, ' '
            print *, ' Size of the matrix is ', nrow
            print *, ' The number of Ritz values requested is ',neiv
            print *, ' The number of Arnoldi vectors generated',
     $               ' (NCV) is ', ncv
            print *, ' What portion of the spectrum: ', which
            print *, ' The number of converged Ritz values is ', 
     $                 nconv 
            print *, ' The number of Implicit Arnoldi update',
     $               ' iterations taken is ', iparam(3)
            print *, ' The number of OP*x is ', iparam(9)
            print *, ' The convergence criterion is ', tol
            print *, ' The number of calls y=A@x: icount =', icount 
            print *, ' '
         endif 
!
      endif
!
!     %---------------------------%
!     | Done with program dsdrv1. |
!     %---------------------------%
!
 9000 continue
!
      if ( flag_wrtout )
     $   write(*,'("Done. Exit ARPACK.")')
      return 
 100  format( A9, 1x, a15, 2x, a36 )
 200  format( i9, ':', E15.5, 2x, E36.5 )
      end subroutine arpack_exec_qp 
!-----------------------------------------------------------------------
      subroutine matvec_1_qp ( n, x, am, jm, g, y ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip), intent(in) ::
     $   n
      real(wrp), intent(in) ::
     $   x(*)
      real(wrp), intent(in) ::
     $   am(*)
      integer(wip), intent(in) ::
     $   jm(*)
      real(wrp), intent(in) :: 
     $   g
      real(wrp), intent(out) ::
     $   y(*) 
!
!     LOCAL VARIABLES:
      integer(wip)
     $   i, k, j 
!
!     EXECUTABLE STATEMENTS: 
!
      do i = 1, n
         y(i) = am(i)*x(i)
      enddo
      do i = 1, n
         do k = jm(i), jm(i+1)-1
            j = jm(k) 
            y(i) = y(i) + g*x(j)
            y(j) = y(j) + g*x(i)
         enddo 
      enddo
      return 
      end subroutine matvec_1_qp
!-----------------------------------------------------------------------
      subroutine matvec_2_qp ( n, x, am, jm, y ) 
!
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip), intent(in) ::
     $   n
      real(wrp), intent(in) ::
     $   x(*)
      real(wrp), intent(in) ::
     $   am(*)
      integer(wip), intent(in) ::
     $   jm(*)
      real(wrp), intent(out) ::
     $   y(*) 
!
!     LOCAL VARIABLES:
      real(wrp) 
     $   h_ij  
      integer(wip)
     $   i, k, j 
!
!     EXECUTABLE STATEMENTS: 
!
      do i = 1, n
         y(i) = am(i)*x(i)
      enddo
      do i = 1, n
         do k = jm(i), jm(i+1)-1
            j = jm(k) 
            h_ij = am(k)
            y(i) = y(i) + h_ij*x(j)
            y(j) = y(j) + h_ij*x(i)
         enddo 
      enddo
      return 
      end subroutine matvec_2_qp
!-----------------------------------------------------------------------
      subroutine copy_s_to_d_symmsr_qp ( 
     $           n, am, jm, g, g_is_mat, bm ) 
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   n
      real(wrp)
     $   am(*) 
      integer(wip)
     $   jm(*)
      real(wrp)  
     $   g
      logical
     $   g_is_mat 
      real(wrp)
     $   bm(n,n) 
!
!     DEPENDENCES: 
      external
     $   copy_s_to_d_symmsr_1_qp, 
     $   copy_s_to_d_symmsr_2_qp
!
!     EXECUTABLE STATEMENTS: 
!
      if ( g_is_mat ) then 
         call copy_s_to_d_symmsr_2_qp ( 
     $        n, am, jm, bm ) 
      else 
         call copy_s_to_d_symmsr_1_qp ( 
     $        n, am, jm, g, bm ) 
      endif 
      return
      end subroutine copy_s_to_d_symmsr_qp
!-----------------------------------------------------------------------
      subroutine copy_s_to_d_symmsr_1_qp ( 
     $           n, am, jm, g, bm ) 
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   n
      real(wrp)
     $   am(*) 
      integer(wip)
     $   jm(*)
      real(wrp)  
     $   g
      real(wrp)
     $   bm(n,n) 
!
!     LOCAL VARIABLES:
      integer(wip)
     $   ir, k, jc
!
!     EXECUTABLE STATEMENTS: 
!
      bm = 0.0_wrp 
      do ir = 1, n
         bm(ir,ir) = am(ir)
      enddo 
      do ir = 1, n
         do k  = jm(ir), jm(ir+1)-1 
            jc = jm(k)
            bm(ir,jc) = g
            bm(jc,ir) = g
         enddo 
      enddo 
      return
      end subroutine copy_s_to_d_symmsr_1_qp
!-----------------------------------------------------------------------
      subroutine copy_s_to_d_symmsr_2_qp ( 
     $           n, am, jm, bm ) 
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   n
      real(wrp)
     $   am(*) 
      integer(wip)
     $   jm(*)
      real(wrp)  
     $   g
      real(wrp)
     $   bm(n,n) 
!
!     LOCAL VARIABLES:
      integer(wip)
     $   ir, k, jc
      real(wrp)
     $   h_ij 
!
!     EXECUTABLE STATEMENTS: 
!
      bm = 0.0_wrp 
      do ir = 1, n
         bm(ir,ir) = am(ir)
      enddo 
      do ir = 1, n
         do k  = jm(ir), jm(ir+1)-1 
            jc = jm(k)
            h_ij = am(k)
            bm(ir,jc) = h_ij 
            bm(jc,ir) = h_ij 
         enddo 
      enddo
      return
      end subroutine copy_s_to_d_symmsr_2_qp
!-----------------------------------------------------------------------
      subroutine lapack_exec_qp ( 
     $           nrow, bmat, eival, lwork, work, 
     $           nm, am, jm, flag_wrtout )
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   nrow
      real(wrp)
     $   bmat ( nrow, nrow )
      real(wrp)
     $   eival ( nrow )
      integer(wip)
     $   lwork
      real(wrp)
     $   work ( lwork )
      integer(wip)
     $   nm
      real(wrp)
     $   am ( nm )
      integer(wip)
     $   jm ( nm )
      logical
     $   flag_wrtout
!
!     DEPENDENCES: 
      real(wrp) 
     &   qpnrm2
      external          
     $   qpsyev,
     $   qpnrm2, 
     $   qpaxpy,
     $   qpsymv 
!
!     LOCAL VARIABLES: 
      integer
     $   info
      real(wrp), allocatable :: 
     $   amat(:,:), workd(:), d(:)
      integer(wip) 
     $   j
!
!     EXECUTABLE STATEMENTS: 
!
      if ( flag_wrtout )
     $   write(*,'("Inside LAPACK, LWORK =", i12)') lwork 
      allocate ( amat(nrow,nrow), workd(nrow), d(nrow) ) 
      amat = bmat 
      call qpsyev ( 
     $     'V','U', nrow, bmat(1,1), nrow, eival(1), 
     $     work(1), lwork, info )
!
!     Compute the residual norm ||A@x - lambda*x||/|lambda|: 
      do j = 1, nrow 
!
!        workd := A@x
         call qpsymv ( 
     $        'U', nrow, 1.0_wrp, amat, nrow, bmat(1,j), 1_wip,
     $        0.0_wrp, workd(1), 1_wip )
!
!        workd := workd - lambda*x = A@x - lambda*x  
         call qpaxpy ( 
     $        nrow, -eival(j), bmat(1,j), 1_wip, workd(1), 1_wip )
!
!        d(j) := ||workd|| = ||A@x - lambda*x|| 
         d(j) = qpnrm2 ( 
     $          nrow, workd(1), 1_wip )
!
!        d(j) := ||A@x - lambda*x|| / |lambda|
         d(j) = d(j) / abs( eival(j) )
!               
      enddo 
      if ( flag_wrtout ) then 
         write(*,'(/)')
         write(*,10) ' ', 'Eigenvalues', 'Rel. Residuals'
         write(*,10) 'k', 'lam(k)',      '||H@v - lam(k)*v|| / |lam(k)|'
         do j = 1, nrow
            write(*,20) j, eival(j), d(j) 
         enddo
      endif 
      deallocate ( amat, workd, d ) 
      if ( flag_wrtout )     
     $   write(*,'(/,"Done. Exit LAPACK.")')
      return 
 10   format( A9, 1x, a15, 2x, a30 )
 20   format( i9, ':', E15.5, 2x, E30.5 )
      end subroutine lapack_exec_qp 
!-----------------------------------------------------------------------
!  Ref:
!     Higham, Nicholas J. "The accuracy of floating point summation." 
!     SIAM Journal on Scientific Computing 14, no. 4 (1993): 783-799.
!
      subroutine higham_sum_qp ( s, x, e )
      implicit none 
      include 'epcode_inc_qp.h'     
!     
!     ARGUMENTS:  
      real(wrp)      
     $   s, x, e 
!
!     LOCAL VARIABLES: 
      real(wrp)      
     $   t, y, r 
!
!     EXECUTABLE STATEMENTS:  
!
      t = s            ! t:   _____t1____.____t2___
      y = x + e        ! y:              .____y1___ ____y2___
      s = t + y        ! s:   _____t1____.__t2+y1__
      r = t - s        ! r:              .___-y1___
      e = r + y        ! e:              .          ____y2___
      return 
      end subroutine higham_sum_qp 
!-----------------------------------------------------------------------
!***********************************************************************
      subroutine prt_nsj_njj ( 
     $           iou, nome, lsj, nsj, ljj, njj )
      implicit none 
      include 'epcode_inc_dp.h' 
!
!     ARGUMENTS:
      integer
     $   iou
      integer
     $   nome
      integer(wip)
     $   lsj, ljj 
      integer(wip)
     $   nsj(lsj), njj(ljj)
!
!     DEPENDENCES:
      external 
     $   deccon
!
!     LOCAL VARIABLES:
      integer(wip)
     $   i, j
      integer
     $   k, lb1 
       integer(i1b), allocatable ::
     $   njjb(:), nsjb(:), nb1(:), nb2(:)    
      character
     $   chnm1*2, chnp1*2 
      character(len=:), allocatable ::
     $   chafmt, chbfmt 
!
!     EXECUTABLE STATEMENTS: 
!
      chnm1 = '' 
      write(chnm1,fmt='(I2)') nome - 1
      chnp1 = '' 
      write(chnp1,fmt='(I2)') nome + 1
      chbfmt = 'B'//trim(adjustl(chnp1))//'.'//trim(adjustl(chnp1)) 
      chafmt = '(a5, 2(I11,1x), '//chbfmt//',1x,"(",I1,' //
     $          trim(adjustl(chnm1)) // 'I2,")")'
      allocate ( 
     $   njjb(nome), nsjb(nome), nb1(nome), nb2(nome) ) 
      njjb = 0 
      nsjb = 0
      nb1  = 0
      nb2  = 0 
      do i = 1_wip, lsj 
         call deccon ( nsj(i), nome, nsjb )
         lb1 = 0 
         do k = 1, nome
            if ( nsjb(k) .ne. 0 ) cycle 
            lb1 = lb1 + 1
            nb1(lb1) = k
         enddo
         do j = 1_wip, ljj
            njjb = nsjb 
            call deccon ( njj(j), nome, nb2 )
            do k = 1, lb1  
               njjb( nb1(k) ) = nb2(k)*2 
            enddo
            write ( iou, fmt=chafmt )
     $         'nsj:', i, nsj(i), nsj(i), (nsjb(k), k=1,nome)
            write ( iou, fmt=chafmt )
     $         'njj:', j, njj(j), njj(j), (njjb(k), k=1,nome)
         enddo 
      enddo 
      deallocate ( njjb, nsjb, nb1, nb2 )
      return 
      end subroutine prt_nsj_njj 
!----------------------------------------------------------------------
      function binom ( n, k ) 
      implicit none 
      include 'epcode_inc_dp.h' 
!
      integer(wip)
     $   binom 
!
!     ARGUMENTS:
      integer
     $   n, k
!
!     LOCAL VARIABLES:
      integer    
     $   i 
      real(kind=dpr)   
     $   c
!
!     EXECUTABLE STATEMENTS:
!
      binom = 0 
      if ( n.lt.0 .or. k.lt.0 .or. n.lt.k ) return 
      c = 1
      do i = 1, k
         c = c * (n-k+i) / i
      enddo
      binom = int(c,kind=wip)
      return 
      end function binom
!-----------------------------------------------------------------------
      subroutine deccon ( j, n, nk )
      implicit none 
      include 'epcode_inc_dp.h' 
!       
!     ARGUMENTS: 
      integer(wip)
     $   j
      integer
     $   n
      integer(i1b)
     $   nk(n)
!
!     DEPENDENCES: 
      intrinsic ::
     $   int, iand, ishft
!
!     LOCAL VARIABLES:
      integer
     $   i, l 
      integer(wip)
     $   k
!
!     EXECUTABLE STATEMENTS: 
!
      k = j 
      do i = 1, n 
         nk(i) = int(iand(k,1_wip),kind=i1b) 
         k = ishft(k,-1_wip)                
         if ( k .eq. 0_wip ) exit  
      enddo 
      do l = i+1, n 
         nk(l) = 0
      enddo 
      return 
      end subroutine deccon
!-----------------------------------------------------------------------
      subroutine gens ( n, p, l, s )   
      implicit none 
      include 'epcode_inc_dp.h' 
!
!     ARGUMENTS: 
      integer,intent(in) :: n, p
      integer(wip),intent(in) :: l 
      integer(wip),dimension(l),intent(out) :: s
!
!     DEPENDENCES: 
      intrinsic :: ior, iand, popcnt  
!
!     LOCAL VARIABLES:
      integer(kind=wip) :: k, j0, j1, j2
      integer(kind=i1b) :: k0 
!
!     EXECUTABLE STATEMENTS: 
!
      s(1) = 2**p - 1 
      do k = 2, l
         j0 = s(k-1) 
         j1 = ior( j0, j0-1 ) 
         j2 = j0 - iand( j1+1, j0 )
         k0 =-1 
         do while ( j2 .ne. 0 )
            k0 = k0 + 1 
            j2 = iand( j2, j2-1 )
         enddo
         s(k) = j1 + 2**(k0) 
      enddo 
      return 
      end subroutine gens
!-----------------------------------------------------------------------
      function hashs ( n, p, s )  result ( js )
      implicit none
      include 'epcode_inc_dp.h' 
!
      integer(wip) :: js 
!
!     ARGUMENTS: 
      integer,intent(in) :: n, p 
      integer(wip),intent(in) :: s
!
!     DEPENDENCES: 
      integer(wip),external :: binom 
      intrinsic :: iand, ishft 
      integer(wip)
     $   binom_tab(0:63,0:63)
      common /ep_gstate_main_shared/ binom_tab
!
!     LOCAL VARIABLES:
      integer(wip) :: t
      integer :: il, ih  
!
!     EXECUTABLE STATEMENTS: 
!
      t = s  
      js = 1
      ih = 0
      il = 0 
      do while ( t .ne. 0 )
         if ( iand(t,1_wip) .ne. 0 ) then 
            ih = ih + 1
            js = js + binom_tab ( il, ih )
            if ( ih .eq. p ) exit 
         endif
         t = ishft ( t, -1 )
         il = il + 1 
      enddo 
      return 
      end function hashs
!-----------------------------------------------------------------------
      subroutine genp ( n, p, s, lr, r, k ) 
      implicit none 
      include 'epcode_inc_dp.h' 
!
!     ARGUMENTS: 
      integer,intent(in) :: n, p  
      integer(wip),intent(in) :: s 
      integer,intent(in) :: lr 
      integer(wip),intent(inout) :: r(lr)
      integer,intent(out) :: k 
!
!     DEPENDENCES: 
      intrinsic :: iand, ishft, ibset, ibclr
!
!     LOCAL VARIABLES:
      integer(kind=i1b) :: j, i, l, a(n)
      integer(wip) :: t 
!
!     EXECUTABLE STATEMENTS: 
!
      t = s 
      k = 0
      j = 0
      a = 0
      i = 0 
      do while ( i < n ) 
         if ( iand(t,1_wip) .ne. 0 ) then 
            j = j + 1
            a(j) = i
         else 
            do l = j, 1, -1
               k = k + 1 
               r(k) = s
               r(k) = ibset( r(k), i )
               r(k) = ibclr( r(k), a(l) )
            enddo 
         endif
         t = ishft( t, -1 )
         i = i + 1  
      enddo
      return  
      end subroutine genp
!-----------------------------------------------------------------------
      function npairing ( s, n ) result ( k ) 
      implicit none 
      include 'epcode_inc_dp.h' 
!
      integer(wip) :: k  
!
!     ARGUMENTS: 
      integer(wip),intent(in) :: s 
      integer,intent(in) :: n  
!
!     DEPENDENCES: 
      intrinsic :: ishft, iand 
!
!     LOCAL VARIABLES:
      integer(kind=i1b) :: j, i 
      integer(wip) :: t 
!
!     EXECUTABLE STATEMENTS: 
!
      t = s 
      j = 0 
      k = 0
      i = 0 
      do while ( t .ne. 0 ) 
         if ( iand(t,1_wip) .ne. 0 ) then  
            j = j + 1 
         else 
            k = k + j 
         endif 
         t = ishft( t, -1 )
         i = i + 1 
      enddo
      k = k + j*(n-i)
      return 
      end function npairing
!-----------------------------------------------------------------------
      subroutine gind ( s1, s2, k1, k2 )
      implicit none
      include 'epcode_inc_dp.h' 
      integer(wip), intent(in) :: s1, s2 
      integer(i1b), intent(out) :: k1, k2 
!
!     LOCAL VARIABLES:
      integer(wip) :: x, x1, x2 
!
!     DEPENDENCES: 
      intrinsic :: ieor, iand, ior, popcnt
!
!     EXECUTABLE STATEMENTS: 
!
      x  = ieor ( s1, s2 )
      x2 = iand ( x, x-1 )
      x1 = x - x2
      k2 = popcnt ( ior( x2, x2-1 ) )  
      k1 = popcnt ( ior( x1, x1-1 ) )
      return 
      end subroutine gind
!-----------------------------------------------------------------------
      function cparams0 ( n, npar, g, ns ) result ( st ) 
      implicit none
! 
      character(len=80)
     $   st 
!
!     ARGUMENTS: 
      integer, intent(in) :: 
     $   n, npar
      real, intent(in) ::
     $   g 
      integer, intent(in) :: 
     $   ns
!
!     LOCAL VARIABLES: 
      character(len=5)
     $   ch 
!
!     EXECUTABLE STATEMENTS: 
!
      st = '_'
      write(ch,fmt='(I0.2)') N 
      st = trim(adjustl(st))//'n'//trim(adjustl(ch))
      write(ch,fmt='(I0.2)') npar 
      st = trim(adjustl(st))//'p'//trim(adjustl(ch))
      write(ch,fmt='(f0.3)') g 
      if      ( ch(1:1) .eq. '.' )  then 
         ch(1:1) = '0'
      else if ( ch(1:1) .eq. '1' )  then 
         ch(2:2) = '1'
         ch(1:1) = ' '
      endif 
      st = trim(adjustl(st))//'g'//trim(adjustl(ch))
      write(ch,fmt='(I0.2)') ns
      st = trim(adjustl(st))//'s'//trim(adjustl(ch))
      return 
      end function cparams0
!-----------------------------------------------------------------------
      subroutine get_ram_in_bytes ( iram )
      implicit none 
!
!     ARGUMENTS: 
      integer(8), intent(out) :: iram 
!
!     LOCAL VARIABLES: 
      character(len=:), allocatable :: fname
      character(len=500) :: cline, ctmp 
      integer :: istat, io  
!
!     EXECUTABLE STATEMENTS: 
!
      fname = '/proc/meminfo'
      istat = 0 
      open(newunit=io, file=fname, status='old', iostat=istat)
      iram = 0 
      if ( istat .ne. 0 ) then 
         write(*,*) 'WARNING: could not open the file '// fname 
         return 
      endif 
      read(unit=io, fmt='(a)') cline 
      close(io)
      read(cline,*) ctmp, iram 
      iram = iram * int(2,kind=8)**10 
      return 
      end subroutine 
!-----------------------------------------------------------------------
      subroutine ep_getmem_dp (
     $           nome, npar, ng, neiv, prefix, 
     $           ljj, nele, mjp_len, lsj, lawork, ncv, 
     $           num_of_bytes_in_arrays, 
     $           num_of_bytes_in_locvar, 
     $           num_of_bytes_in_total ) 
      implicit none 
      include 'epcode_inc_dp.h'
!
!     ARGUMENTS: 
      integer, intent(in) :: 
     $   nome, npar, ng 
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(in) :: 
     $   prefix
      integer(wip), intent(out) :: 
     $   num_of_bytes_in_total, num_of_bytes_in_arrays,
     $   num_of_bytes_in_locvar
!
!     LOCAL VARIABLES:
      integer(wip)
     $   binom_tab(0:63,0:63)
      common /ep_gstate_main_shared/ binom_tab
      logical, parameter ::
     $   flag_screen = .true.,   
     $   flag_wrtres = .true.,  
     $   flag_wrtout = .false.,  
     $   flag_wrtsiz = .false., 
     $   flag_memest = .true., 
     $   flag_wrtnsj = .false.,  
     $   flag_wrtmat = .false., 
     $   flag_wrtevv = .false., 
     $   flag_wrtocc = .false.
      integer(wip)
     $   lsj, meiv 
      real(wrp) :: 
     $   tole, tmp 
      integer
     $   mjp_len, ljp, nps, nmns, lb1, i, i1, i2, isum, 
     $   k, l, nt0, iox, ns 
      integer(wip)
     $   nrow, ncol, ir, ic, isj, ijj, ljj,
     $   nnz, iele, nele, jos, state_i, state_f,  
     $   njir, njic, ntmp, nde
      integer(wip)      
     $   ncv, lawork, lworkl, 
     $   io_resid, io_workd, io_eivec, io_eival, io_workl
      character
     $   chn*5, chpack*30, chfmt*50, 
     $   fort00*90, fort01*90, fort02*90, fort03*90, 
     $   fort04*90, fort05*90, header*200
      integer ::
     $   iou00, iou01, iou02, iou03, iou04, iou05 
      integer(i1b)
     $   k1, k2  
      logical 
     $   g_is_matrix, use_arpack  
      integer(8) 
     $   iram 
      integer(wip)      
     $   sizeof_e, sizeof_g, sizeof_njjb, sizeof_nsjb,
     $   sizeof_nb1, sizeof_nb2, sizeof_njj,
     $   sizeof_jah, sizeof_mjp, sizeof_aah, 
     $   sizeof_occ, sizeof_ecc, sizeof_nsj, 
     $   sizeof_awork, sizeof_choose, sizeof_g_diag,
     $   sizeof_g_matr, sizeof_binom_tab
      logical 
     $   flag_prep, 
     $   flag_wrtinp 
!
!     EXECUTABLE STATEMENTS: 
!
      sizeof_e = nome * wrp
      sizeof_g = ng*ng * wrp
      sizeof_g_diag = 63 * wrp
      sizeof_g_matr = 63 * 63 * wrp
      sizeof_binom_tab = 64 * 64 * wip
      sizeof_njjb = nome * i1b
      sizeof_nsjb = nome * i1b
      sizeof_nb1 = nome * i1b
      sizeof_nb2 = nome * i1b
      sizeof_njj = ljj * wip  
      sizeof_nsj = lsj * wip  
      sizeof_mjp = mjp_len * wip  
      sizeof_jah = nele * wip  
      if ( g_is_matrix ) then 
         sizeof_aah = nele * wrp 
      else 
         sizeof_aah = ljj * wrp 
      endif 
      sizeof_occ = nome * meiv * wrp  
      sizeof_ecc = nome * meiv * wrp  
      sizeof_awork = lawork * wrp  
      sizeof_choose = ncv * sizeof(flag_screen)   
      num_of_bytes_in_arrays = 
     $   sizeof_e + sizeof_g + 
     $   sizeof_njjb + sizeof_nsjb + 
     $   sizeof_nb1  + sizeof_nb2  +
     $   sizeof_njj  + sizeof_jah  + 
     $   sizeof_mjp  + sizeof_aah  + 
     $   sizeof_occ  + sizeof_ecc  + 
     $   sizeof_nsj  + sizeof_mjp  + 
     $   sizeof_awork + sizeof_g_diag +
     $   sizeof_g_matr + sizeof_binom_tab +  
     $   sizeof_choose
      num_of_bytes_in_locvar =
     $   sizeof(nome) + sizeof(npar) + 
     $   sizeof(ng) + sizeof(tole) + 
     $   sizeof(neiv) + sizeof(prefix) + 
     $   sizeof(flag_screen) + 
     $   sizeof(flag_wrtres) +
     $   sizeof(flag_wrtout) +
     $   sizeof(flag_wrtsiz) +
     $   sizeof(flag_memest) +
     $   sizeof(flag_wrtnsj) +
     $   sizeof(flag_wrtmat) +
     $   sizeof(flag_wrtevv) +
     $   sizeof(flag_wrtocc) +
     $   sizeof(flag_prep) + 
     $   sizeof(flag_wrtinp) +
     $   sizeof(g_is_matrix) + 
     $   sizeof(use_arpack) + 
     $   sizeof(lsj) + sizeof(ns) + 
     $   sizeof(neiv) + sizeof(tmp) + 
     $   sizeof(mjp_len) + sizeof(ljp) + 
     $   sizeof(nps) + sizeof(nmns) + 
     $   sizeof(lb1) + sizeof(i) + 
     $   sizeof(i1) + sizeof(i2) + 
     $   sizeof(isum) + sizeof(k) + 
     $   sizeof(l) + sizeof(nt0) + 
     $   sizeof(iox) + sizeof(nrow) + 
     $   sizeof(ncol) + sizeof(ir) + 
     $   sizeof(isj) + sizeof(ijj) + 
     $   sizeof(ljj) + sizeof(nnz) + 
     $   sizeof(iele) + sizeof(nele) + 
     $   sizeof(jos) + sizeof(state_i) + 
     $   sizeof(state_f) + sizeof(njir) + 
     $   sizeof(njic) + sizeof(ntmp) + 
     $   sizeof(nde) + sizeof(ncv) + 
     $   sizeof(lawork) + sizeof(lworkl) + 
     $   sizeof(io_resid) + sizeof(io_workd) + 
     $   sizeof(io_eivec) + sizeof(io_eival) + 
     $   sizeof(io_workl) + sizeof(ic) +  
     $   sizeof(chn) + sizeof(chpack) + 
     $   sizeof(fort00) + sizeof(fort01) + 
     $   sizeof(fort02) + sizeof(fort03) + 
     $   sizeof(fort04) + sizeof(fort05) + 
     $   sizeof(header) + sizeof(iou00) + 
     $   sizeof(iou01) + sizeof(iou02) + 
     $   sizeof(iou03) + sizeof(iou04) + 
     $   sizeof(iou05)  
      num_of_bytes_in_total = 
     $   num_of_bytes_in_arrays + num_of_bytes_in_locvar
!!    write(*,'(/,A)') 'COMPUTER MEMORY REQUIREMENT:' 
!!    write(*,'( 3x, I20, A10, 2(E17.10, A) )') 
!!   $    num_of_bytes_in_total, ' (bytes) =', 
!!   $    num_of_bytes_in_total * dble(2)**(-20), ' (MB) =',
!!   $    num_of_bytes_in_total * dble(2)**(-30), ' (GB)'
      return 
      end subroutine 
!-----------------------------------------------------------------------
      subroutine ep_getmem_qp (
     $           nome, npar, ng, neiv, prefix, 
     $           ljj, nele, mjp_len, lsj, lawork, ncv, 
     $           num_of_bytes_in_arrays, 
     $           num_of_bytes_in_locvar, 
     $           num_of_bytes_in_total ) 
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer, intent(in) :: 
     $   nome, npar, ng 
      integer(wip), intent(in) :: 
     $   neiv  
      character(len=*), intent(in) :: 
     $   prefix
      integer(wip), intent(out) :: 
     $   num_of_bytes_in_total, num_of_bytes_in_arrays,
     $   num_of_bytes_in_locvar
!
!     LOCAL VARIABLES:
      integer(wip)
     $   binom_tab(0:63,0:63)
      common /ep_gstate_main_shared/ binom_tab
      logical, parameter ::
     $   flag_screen = .true.,   
     $   flag_wrtres = .true.,  
     $   flag_wrtout = .false.,  
     $   flag_wrtsiz = .false., 
     $   flag_memest = .true., 
     $   flag_wrtnsj = .false.,  
     $   flag_wrtmat = .false., 
     $   flag_wrtevv = .false., 
     $   flag_wrtocc = .false.
      integer(wip)
     $   lsj, meiv 
      real(wrp) :: 
     $   tole, tmp 
      integer
     $   mjp_len, ljp, nps, nmns, lb1, i, i1, i2, isum, 
     $   k, l, nt0, iox, ns 
      integer(wip)
     $   nrow, ncol, ir, ic, isj, ijj, ljj,
     $   nnz, iele, nele, jos, state_i, state_f,  
     $   njir, njic, ntmp, nde
      integer(wip)      
     $   ncv, lawork, lworkl, 
     $   io_resid, io_workd, io_eivec, io_eival, io_workl
      character
     $   chn*5, chpack*30, chfmt*50, 
     $   fort00*90, fort01*90, fort02*90, fort03*90, 
     $   fort04*90, fort05*90, header*200
      integer ::
     $   iou00, iou01, iou02, iou03, iou04, iou05 
      integer(i1b)
     $   k1, k2  
      logical 
     $   g_is_matrix, use_arpack  
      integer(8) 
     $   iram 
      integer(wip)      
     $   sizeof_e, sizeof_g, sizeof_njjb, sizeof_nsjb,
     $   sizeof_nb1, sizeof_nb2, sizeof_njj,
     $   sizeof_jah, sizeof_mjp, sizeof_aah, 
     $   sizeof_occ, sizeof_ecc, sizeof_nsj, 
     $   sizeof_awork, sizeof_choose, sizeof_g_diag,
     $   sizeof_g_matr, sizeof_binom_tab
      logical 
     $   flag_prep, 
     $   flag_wrtinp 
!
!     EXECUTABLE STATEMENTS: 
!
      sizeof_e = nome * wrp
      sizeof_g = ng*ng * wrp
      sizeof_g_diag = 63 * wrp
      sizeof_g_matr = 63 * 63 * wrp
      sizeof_binom_tab = 64 * 64 * wip
      sizeof_njjb = nome * i1b
      sizeof_nsjb = nome * i1b
      sizeof_nb1 = nome * i1b
      sizeof_nb2 = nome * i1b
      sizeof_njj = ljj * wip  
      sizeof_nsj = lsj * wip  
      sizeof_mjp = mjp_len * wip  
      sizeof_jah = nele * wip  
      if ( g_is_matrix ) then 
         sizeof_aah = nele * wrp 
      else 
         sizeof_aah = ljj * wrp 
      endif 
      sizeof_occ = nome * meiv * wrp  
      sizeof_ecc = nome * meiv * wrp  
      sizeof_awork = lawork * wrp  
      sizeof_choose = ncv * sizeof(flag_screen)   
      num_of_bytes_in_arrays = 
     $   sizeof_e + sizeof_g + 
     $   sizeof_njjb + sizeof_nsjb + 
     $   sizeof_nb1  + sizeof_nb2  +
     $   sizeof_njj  + sizeof_jah  + 
     $   sizeof_mjp  + sizeof_aah  + 
     $   sizeof_occ  + sizeof_ecc  + 
     $   sizeof_nsj  + sizeof_mjp  + 
     $   sizeof_awork + sizeof_g_diag +
     $   sizeof_g_matr + sizeof_binom_tab +  
     $   sizeof_choose
      num_of_bytes_in_locvar =
     $   sizeof(nome) + sizeof(npar) + 
     $   sizeof(ng) + sizeof(tole) + 
     $   sizeof(neiv) + sizeof(prefix) + 
     $   sizeof(flag_screen) + 
     $   sizeof(flag_wrtres) +
     $   sizeof(flag_wrtout) +
     $   sizeof(flag_wrtsiz) +
     $   sizeof(flag_memest) +
     $   sizeof(flag_wrtnsj) +
     $   sizeof(flag_wrtmat) +
     $   sizeof(flag_wrtevv) +
     $   sizeof(flag_wrtocc) +
     $   sizeof(flag_prep) + 
     $   sizeof(flag_wrtinp) +
     $   sizeof(g_is_matrix) + 
     $   sizeof(use_arpack) + 
     $   sizeof(lsj) + sizeof(ns) + 
     $   sizeof(neiv) + sizeof(tmp) + 
     $   sizeof(mjp_len) + sizeof(ljp) + 
     $   sizeof(nps) + sizeof(nmns) + 
     $   sizeof(lb1) + sizeof(i) + 
     $   sizeof(i1) + sizeof(i2) + 
     $   sizeof(isum) + sizeof(k) + 
     $   sizeof(l) + sizeof(nt0) + 
     $   sizeof(iox) + sizeof(nrow) + 
     $   sizeof(ncol) + sizeof(ir) + 
     $   sizeof(isj) + sizeof(ijj) + 
     $   sizeof(ljj) + sizeof(nnz) + 
     $   sizeof(iele) + sizeof(nele) + 
     $   sizeof(jos) + sizeof(state_i) + 
     $   sizeof(state_f) + sizeof(njir) + 
     $   sizeof(njic) + sizeof(ntmp) + 
     $   sizeof(nde) + sizeof(ncv) + 
     $   sizeof(lawork) + sizeof(lworkl) + 
     $   sizeof(io_resid) + sizeof(io_workd) + 
     $   sizeof(io_eivec) + sizeof(io_eival) + 
     $   sizeof(io_workl) + sizeof(ic) +  
     $   sizeof(chn) + sizeof(chpack) + 
     $   sizeof(fort00) + sizeof(fort01) + 
     $   sizeof(fort02) + sizeof(fort03) + 
     $   sizeof(fort04) + sizeof(fort05) + 
     $   sizeof(header) + sizeof(iou00) + 
     $   sizeof(iou01) + sizeof(iou02) + 
     $   sizeof(iou03) + sizeof(iou04) + 
     $   sizeof(iou05)  
      num_of_bytes_in_total = 
     $   num_of_bytes_in_arrays + num_of_bytes_in_locvar
!!    write(*,'(/,A)') 'COMPUTER MEMORY REQUIREMENT:' 
!!    write(*,'( 3x, I20, A10, 2(E17.10, A) )') 
!!   $    num_of_bytes_in_total, ' (bytes) =', 
!!   $    num_of_bytes_in_total * dble(2)**(-20), ' (MB) =',
!!   $    num_of_bytes_in_total * dble(2)**(-30), ' (GB)'
      return 
      end subroutine 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***********************************************************************
