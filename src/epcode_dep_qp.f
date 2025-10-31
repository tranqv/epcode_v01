!***********************************************************************
!-----------------------------------------------------------------------
!
!@@@  ARPACK::
!
!-----------------------------------------------------------------------
!
!\BeginDoc
!
!\Name: qpsaupd
!
!\Description: 
!
!  Reverse communication interface for the Implicitly Restarted Arnoldi 
!  Iteration. For symmetric problems this reduces to a variant of the Lanczos 
!  method. This method has been designed to compute approximations to a 
!  few eigenpairs of a linear operator OP that is real and symmetric 
!  with respect to a real positive semi-definite symmetric matrix B, 
!  i.e.
!                   
!       B*OP = (OP')*B.  
!
!  Another way to express this condition is 
!
!       < x,OPy > = < OPx,y >  where < z,w > = z'Bw 
!
!  where A' denotes transpose of A.
!  
!  In the standard eigenproblem B is the identity matrix.  
!
!  The computed approximate eigenvalues are called Ritz values and
!  the corresponding approximate eigenvectors are called Ritz vectors.
!
!  DSAUPD is usually called iteratively to solve one of the 
!  following problems:
!
!  Mode 1:  A*x = lambda*x, A symmetric 
!           ===> OP = A  and  B = I.
!
!  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
!           ===> OP = inv[M]*A  and  B = M.
!           ===> (If M can be factored see remark 3 below)
!
!  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
!           ===> OP = (inv[K - sigma*M])*M  and  B = M. 
!           ===> Shift-and-Invert mode
!
!  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite, 
!           KG symmetric indefinite
!           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
!           ===> Buckling mode
!
!  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
!           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
!           ===> Cayley transformed mode
!
!  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
!        should be accomplished either by a direct method
!        using a sparse matrix factorization and solving
!
!           [A - sigma*M]*w = v  or M*w = v,
!
!        or through an iterative method for solving these
!        systems.  If an iterative method is used, the
!        convergence test must be more stringent than
!        the accuracy requirements for the eigenvalue
!        approximations.
!
!\Usage:
!  call qpsaupd 
!     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first 
!          call to qpsaupd.  IDO will be set internally to
!          indicate the type of operation to be performed.  Control is
!          then given back to the calling routine which has the
!          responsibility to carry out the requested operation and call
!          qpsaupd with the result.  The operand is given in
!          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
!          (If Mode = 2 see remark 5 below)
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    In mode 3,4 and 5, the vector B * X is already
!                    available in WORKD(ipntr(3)).  It does not
!                    need to be recomputed in forming OP * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO =  3: compute the IPARAM(8) shifts where
!                    IPNTR(11) is the pointer into WORKL for
!                    placing the shifts. See remark 6 below.
!          IDO = 99: done
!          -------------------------------------------------------------
!             
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  WHICH   Character*2.  (INPUT)
!          Specify which of the Ritz values of OP to compute.
!
!          'LA' - compute the NEV largest (algebraic) eigenvalues.
!          'SA' - compute the NEV smallest (algebraic) eigenvalues.
!          'LM' - compute the NEV largest (in magnitude) eigenvalues.
!          'SM' - compute the NEV smallest (in magnitude) eigenvalues. 
!          'BE' - compute NEV eigenvalues, half from each end of the
!                 spectrum.  When NEV is odd, compute one more from the
!                 high end than from the low end.
!           (see remark 1 below)
!
!  NEV     Integer.  (INPUT)
!          Number of eigenvalues of OP to be computed. 0 < NEV < N.
!
!  TOL     real(wrp) scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value 
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
!          If TOL .LE. 0. is passed a default is set:
!          DEFAULT = DLAMCH('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine DLAMCH).
!
!  RESID   real(wrp) array of length N.  (INPUT/OUTPUT)
!          On INPUT: 
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector. 
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V (less than or equal to N).
!          This will indicate how many Lanczos vectors are generated 
!          at each iteration.  After the startup phase in which NEV 
!          Lanczos vectors are generated, the algorithm generates 
!          NCV-NEV Lanczos vectors at each subsequent update iteration.
!          Most of the cost in generating each Lanczos vector is in the 
!          matrix-vector product OP*x. (See remark 4 below).
!
!  V       real(wrp) N by NCV array.  (OUTPUT)
!          The NCV columns of V contain the Lanczos basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          -------------------------------------------------------------
!          ISHIFT = 0: the shifts are provided by the user via
!                      reverse communication.  The NCV eigenvalues of
!                      the current tridiagonal matrix T are returned in
!                      the part of WORKL array corresponding to RITZ.
!                      See remark 6 below.
!          ISHIFT = 1: exact shifts with respect to the reduced 
!                      tridiagonal matrix T.  This is equivalent to 
!                      restarting the iteration with a starting vector 
!                      that is a linear combination of Ritz vectors 
!                      associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = LEVEC
!          No longer referenced. See remark 2 below.
!
!          IPARAM(3) = MXITER
!          On INPUT:  maximum number of Arnoldi update iterations allowed. 
!          On OUTPUT: actual number of Arnoldi update iterations taken. 
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" Ritz values.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used. 
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4,5; See under \Description of qpsaupd for the 
!          five modes available.
!
!          IPARAM(8) = NP
!          When ido = 3 and the user provides shifts through reverse
!          communication (IPARAM(1)=0), qpsaupd returns NP, the number
!          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
!          6 below.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*x operations,
!                  NUMOPB = total number of B*x operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.        
!
!  IPNTR   Integer array of length 11.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD and WORKL
!          arrays for matrices/vectors used by the Lanczos iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X in WORKD.
!          IPNTR(2): pointer to the current result vector Y in WORKD.
!          IPNTR(3): pointer to the vector B * X in WORKD when used in 
!                    the shift-and-invert mode.
!          IPNTR(4): pointer to the next available location in WORKL
!                    that is untouched by the program.
!          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
!          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
!          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
!                    with the Ritz values located in RITZ in WORKL.
!          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
!
!          Note: IPNTR(8:10) is only referenced by qpseupd. See Remark 2.
!          IPNTR(8): pointer to the NCV RITZ values of the original system.
!          IPNTR(9): pointer to the NCV corresponding error bounds.
!          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
!                     of the tridiagonal matrix T. Only referenced by
!                     qpseupd if RVEC = .TRUE. See Remarks.
!          -------------------------------------------------------------
!          
!  WORKD   real(wrp) work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD 
!          as temporary workspace during the iteration. Upon termination
!          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
!          subroutine qpseupd uses this output.
!          See Data Distribution Note below.  
!
!  WORKL   real(wrp) work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  See Data Distribution Note below.
!
!  LWORKL  Integer.  (INPUT)
!          LWORKL must be at least NCV**2 + 8*NCV .
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)  
!                returns the number of wanted converged Ritz values.
!          =  2: No longer an informational error. Deprecated starting
!                with release 2 of ARPACK.
!          =  3: No shifts could be applied during a cycle of the 
!                Implicitly restarted Arnoldi iteration. One possibility 
!                is to increase the size of NCV relative to NEV. 
!                See remark 4 below.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV must be greater than NEV and less than or equal to N.
!          = -4: The maximum number of Arnoldi update iterations allowed
!                must be greater than zero.
!          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work array WORKL is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Informatinal error from LAPACK routine qpsteqr.
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -13: NEV and WHICH = 'BE' are incompatable.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization. The user is advised to check that
!                   enough workspace and array storage has been allocated.
!
!
!\Remarks
!  1. The converged Ritz values are always returned in ascending 
!     algebraic order.  The computed Ritz values are approximate
!     eigenvalues of OP.  The selection of WHICH should be made
!     with this in mind when Mode = 3,4,5.  After convergence, 
!     approximate eigenvalues of the original problem may be obtained 
!     with the ARPACK subroutine qpseupd. 
!
!  2. If the Ritz vectors corresponding to the converged Ritz values
!     are needed, the user must call qpseupd immediately following completion
!     of qpsaupd. This is new starting with version 2.1 of ARPACK.
!
!  3. If M can be factored into a Cholesky factorization M = LL'
!     then Mode = 2 should not be selected.  Instead one should use
!     Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular 
!     linear systems should be solved with L and L' rather
!     than computing inverses.  After convergence, an approximate
!     eigenvector z of the original problem is recovered by solving
!     L'z = x  where x is a Ritz vector of OP.
!
!  4. At present there is no a-priori analysis to guide the selection
!     of NCV relative to NEV.  The only formal requrement is that NCV > NEV.
!     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
!     the same type are to be solved, one should experiment with increasing
!     NCV while keeping NEV fixed for a given test problem.  This will 
!     usually decrease the required number of OP*x operations but it
!     also increases the work and storage required to maintain the orthogonal
!     basis vectors.   The optimal "cross-over" with respect to CPU time
!     is problem dependent and must be determined empirically.
!
!  5. If IPARAM(7) = 2 then in the Reverse commuication interface the user
!     must do the following. When IDO = 1, Y = OP * X is to be computed.
!     When IPARAM(7) = 2 OP = inv(B)*A. After computing A*X the user
!     must overwrite X with A*X. Y is then the solution to the linear set
!     of equations B*Y = A*X.
!
!  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
!     NP = IPARAM(8) shifts in locations: 
!     1   WORKL(IPNTR(11))           
!     2   WORKL(IPNTR(11)+1)         
!                        .           
!                        .           
!                        .      
!     NP  WORKL(IPNTR(11)+NP-1). 
!
!     The eigenvalues of the current tridiagonal matrix are located in 
!     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
!     order defined by WHICH. The associated Ritz estimates are located in
!     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
!
!-----------------------------------------------------------------------
!
!\Data Distribution Note:
!
!  Fortran-D syntax:
!  ================
!  REAL       RESID(N), V(LDV,NCV), WORKD(3*N), WORKL(LWORKL)
!  DECOMPOSE  D1(N), D2(N,NCV)
!  ALIGN      RESID(I) with D1(I)
!  ALIGN      V(I,J)   with D2(I,J)
!  ALIGN      WORKD(I) with D1(I)     range (1:N)
!  ALIGN      WORKD(I) with D1(I-N)   range (N+1:2*N)
!  ALIGN      WORKD(I) with D1(I-2*N) range (2*N+1:3*N)
!  DISTRIBUTE D1(BLOCK), D2(BLOCK,:)
!  REPLICATED WORKL(LWORKL)
!
!  Cray MPP syntax:
!  ===============
!  REAL       RESID(N), V(LDV,NCV), WORKD(N,3), WORKL(LWORKL)
!  SHARED     RESID(BLOCK), V(BLOCK,:), WORKD(BLOCK,:)
!  REPLICATED WORKL(LWORKL)
!  
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!  8. R.B. Lehoucq, D.C. Sorensen, "Implementation of Some Spectral
!     Transformations in a k-Step Arnoldi Method". In Preparation.
!
!\Routines called:
!     qpsaup2  ARPACK routine that implements the Implicitly Restarted
!             Arnoldi Iteration.
!     qpstats  ARPACK routine that initialize timing and other statistics
!             variables.
!     wivout   ARPACK utility routine that prints integers.
!     qpsecond  ARPACK utility routine for timing.
!     qpvout   ARPACK utility routine that prints vectors.
!     qplamch  LAPACK routine that determines machine constants.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
! 
!\Revision history:
!     12/15/93: Version ' 2.4'
!
!\SCCS Information: @(#) 
! FILE: saupd.F   SID: 2.7   DATE OF SID: 8/27/96   RELEASE: 2 
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpsaupd ( 
     $           ido, bmat, n, which, nev, tol, resid, 
     $           ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, 
     $           info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer    
     $   ido
      character(len=1) 
     $   bmat
      character(len=2)
     $   which
      integer(wip) 
     $   n, nev, ncv, ldv, lworkl,
     $   iparam(11), ipntr(11),
     $   info
      real(wrp)
     $   tol, resid(n), v(ldv,ncv),
     $   workd(3*n), workl(lworkl)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real(wrp)
     $   one, zero
      parameter (
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer(wip)
     $   np, nev0, j, ldh, ldq, ritz, bounds, iq, iw, next, ih,
     $   nb, mxiter, ishift, mode, iupd 
      integer    
     $   ierr, msglvl
!
      save       
     $   bounds, ierr, ih, iq, ishift, iupd, iw,
     $   ldh, ldq, msglvl, mxiter, mode, nb,
     $   nev0, next, np, ritz
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   
     $   qpsaup2, 
     $   qpvout, 
     $   wivout, 
     $   qpsecond, 
     $   qpstats
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
!
!     lapack function. use FUNCTION WFLAMCH ( CMACH ) instead
!
      real(wrp)
     $   qplamch
      external   
     $   qplamch
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido .eq. 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call qpstats
         call qpsecond (t0)
         msglvl = msaupd
!
         ierr   = 0
         ishift = iparam(1)
         mxiter = iparam(3)
         nb     = iparam(4)
!
!        %--------------------------------------------%
!        | Revision 2 performs only implicit restart. |
!        %--------------------------------------------%
!
         iupd   = 1
         mode   = iparam(7)
!
!        %----------------%
!        | Error checking |
!        %----------------%
!
         if (n .le. 0) then
            ierr = -1
         else if (nev .le. 0) then
            ierr = -2
         else if (ncv .le. nev .or.  ncv .gt. n) then
            ierr = -3
         endif
!
!!       %----------------------------------------------%
!!       | NP is the number of additional steps to      |
!!       | extend the length NEV Lanczos factorization. |
!!       %----------------------------------------------%
!!       np     = ncv - nev
! 
         if (mxiter .le. 0)                     
     $      ierr = -4
         if (which .ne. 'LM' .and.
     $       which .ne. 'SM' .and.
     $       which .ne. 'LA' .and.
     $       which .ne. 'SA' .and.
     $       which .ne. 'BE' )                  
     $      ierr = -5
         if (bmat .ne. 'I' .and. bmat .ne. 'G') 
     $      ierr = -6
!
         if (lworkl .lt. ncv**2 + 8*ncv)        
     $      ierr = -7
         if (mode .lt. 1 .or. mode .gt. 5) then
            ierr = -10
         else if (mode .eq. 1 .and. bmat .eq. 'G') then
            ierr = -11
         else if (ishift .lt. 0 .or. ishift .gt. 1) then
            ierr = -12
         else if (nev .eq. 1 .and. which .eq. 'BE') then
            ierr = -13
         endif
! 
!        %------------%
!        | Error Exit |
!        %------------%
!
         if (ierr .ne. 0) then
            info = ierr
            ido  = 99
            goto 9000
         endif
! 
!        %------------------------%
!        | Set default parameters |
!        %------------------------%
!
         if (nb .le. 0) nb = 1
         if (tol .le. zero)  
     $      tol = qplamch (
     $            'EpsMach')
!
!        %----------------------------------------------%
!        | NP is the number of additional steps to      |
!        | extend the length NEV Lanczos factorization. |
!        | NEV0 is the local variable designating the   |
!        | size of the invariant subspace desired.      |
!        %----------------------------------------------%
!
         np     = ncv - nev
         nev0   = nev 
! 
!        %-----------------------------%
!        | Zero out internal workspace |
!        %-----------------------------%
!
         do j = 1, ncv**2 + 8*ncv
            workl(j) = zero
         enddo 
! 
!        %-------------------------------------------------------%
!        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
!        | etc... and the remaining workspace.                   |
!        | Also update pointer to be used on output.             |
!        | Memory is laid out as follows:                        |
!        | workl(1:2*ncv) := generated tridiagonal matrix        |
!        | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
!        | workl(3*ncv+1:3*ncv+ncv) := computed error bounds     |
!        | workl(4*ncv+1:4*ncv+ncv*ncv) := rotation matrix Q     |
!        | workl(4*ncv+ncv*ncv+1:7*ncv+ncv*ncv) := workspace     |
!        %-------------------------------------------------------%
!
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritz   = ih     + 2*ldh
         bounds = ritz   + ncv
         iq     = bounds + ncv
         iw     = iq     + ncv**2
         next   = iw     + 3*ncv
!
         ipntr( 4) = next
         ipntr( 5) = ih
         ipntr( 6) = ritz
         ipntr( 7) = bounds
         ipntr(11) = iw
      endif
!
!     %-------------------------------------------------------%
!     | Carry out the Implicitly restarted Lanczos Iteration. |
!     %-------------------------------------------------------%
!
!!    write(*,'( A30, 3x, "ido =", i7)') "(before qpsaup2)",  ido 
!
      call qpsaup2 ( 
     $     ido, bmat, n, which, nev0, np, tol, resid, mode, iupd,
     $     ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritz),
     $     workl(bounds), workl(iq), ldq, workl(iw), ipntr, workd,
     $     info )
!
!!    write(*,'( A30, 3x, "ido =", i7)') "(after  qpsaup2)",  ido 
!
!     %--------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication |
!     | to compute operations involving OP or shifts.    |
!     %--------------------------------------------------%
!
      if (ido .eq. 3 ) iparam(8) = np
!
!     Converged. Exit!
      if (ido .ne. 99) goto 9000
! 
      iparam( 3) = mxiter
      iparam( 5) = np
      iparam( 9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
!
!     %------------------------------------%
!     | Exit if there was an informational |
!     | error within qpsaup2.               |
!     %------------------------------------%
!
      if (info .lt. 0) goto 9000
      if (info .eq. 2) info = 3
!
      if (msglvl .gt. 0) then
         call wivout (
     $        logfil, 1_wip, (/mxiter/), ndigit,
     $        '_saupd: number of update iterations taken' )
         call wivout (
     $        logfil, 1_wip, (/np/), ndigit,
     $        '_saupd: number of "converged" Ritz values' )
         call qpvout (
     $        logfil, np, workl(ritz), ndigit, 
     $        '_saupd: final Ritz values' )
         call qpvout (
     $        logfil, np, workl(bounds), ndigit, 
     $        '_saupd: corresponding error bounds' )
      endif 
!
      call qpsecond (t1)
      tsaupd = t1 - t0
! 
      if (msglvl .gt. 0) then
!
!        %--------------------------------------------------------%
!        | Version Number & Version Date are defined in version.h |
!        %--------------------------------------------------------%
!
         write (6,1000)
         write (6,1100) 
     $      mxiter, nopx, nbx, nrorth, nitref, nrstrt,
     $      tmvopx, tmvbx, tsaupd, tsaup2, tsaitr, titref,
     $      tgetv0, tseigt, tsgets, tsapps, tsconv
 1000    format (//,
     $      5x, '==========================================',/
     $      5x, '= Symmetric implicit Arnoldi update code =',/
     $      5x, '= Version Number:', ' 2.4', 19x, ' =',/
     $      5x, '= Version Date:  ', ' 07/31/96', 14x, ' =',/
     $      5x, '==========================================',/
     $      5x, '= Summary of timing statistics           =',/
     $      5x, '==========================================',//)
 1100    format (
     $      5x, 'Total number update iterations             = ', i5,/
     $      5x, 'Total number of OP*x operations            = ', i5,/
     $      5x, 'Total number of B*x operations             = ', i5,/
     $      5x, 'Total number of reorthogonalization steps  = ', i5,/
     $      5x, 'Total number of iterative refinement steps = ', i5,/
     $      5x, 'Total number of restart steps              = ', i5,/
     $      5x, 'Total time in user OP*x operation          = ', f12.6,/
     $      5x, 'Total time in user B*x operation           = ', f12.6,/
     $      5x, 'Total time in Arnoldi update routine       = ', f12.6,/
     $      5x, 'Total time in saup2 routine                = ', f12.6,/
     $      5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/
     $      5x, 'Total time in reorthogonalization phase    = ', f12.6,/
     $      5x, 'Total time in (re)start vector generation  = ', f12.6,/
     $      5x, 'Total time in trid eigenvalue subproblem   = ', f12.6,/
     $      5x, 'Total time in getting the shifts           = ', f12.6,/
     $      5x, 'Total time in applying the shifts          = ', f12.6,/
     $      5x, 'Total time in convergence testing          = ', f12.6)
      endif
! 
 9000 continue
      return
      end subroutine qpsaupd
!
!-----------------------------------------------------------------------
!
!\BeginDoc
!
!\Name: qpseupd
!
!\Description: 
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) the corresponding approximate eigenvectors,
!
!      (2) an orthonormal (Lanczos) basis for the associated approximate
!          invariant subspace,
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  (Lanczos) basis is always computed.  There is an additional storage cost 
!  of n*nev if both are requested (in this case a separate array Z must be 
!  supplied).
!
!  These quantities are obtained from the Lanczos factorization computed
!  by DSAUPD for the linear operator OP prescribed by the MODE selection
!  (see IPARAM(7) in DSAUPD documentation.)  DSAUPD must be called before
!  this routine is called. These approximate eigenvalues and vectors are 
!  commonly called Ritz values and Ritz vectors respectively.  They are 
!  referred to as such in the comments that follow.   The computed orthonormal 
!  basis for the invariant subspace corresponding to these Ritz values is 
!  referred to as a Lanczos basis.
!
!  See documentation in the header of the subroutine DSAUPD for a definition 
!  of OP as well as other terms and the relation of computed Ritz values 
!  and vectors of OP with respect to the given problem  A*z = lambda*B*z.  
!
!  The approximate eigenvalues of the original problem are returned in
!  ascending algebraic order.  The user may elect to call this routine
!  once for each desired Ritz vector and store it peripherally if desired.
!  There is also the option of computing a selected set of these vectors
!  with a single call.
!
!\Usage:
!  call qpseupd 
!     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
!       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!  RVEC    LOGICAL  (INPUT) 
!          Specifies whether Ritz vectors corresponding to the Ritz value 
!          approximations to the eigenproblem A*z = lambda*B*z are computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute Ritz vectors.
!
!  HOWMNY  Character*1  (INPUT) 
!          Specifies how many Ritz vectors are wanted and the form of Z
!          the matrix of Ritz vectors. See remark 1 below.
!          = 'A': compute NEV Ritz vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NEV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
!          If HOWMNY = 'A' , SELECT is not referenced.
!
!  D       real(wrp) array of dimension NEV.  (OUTPUT)
!          On exit, D contains the Ritz value approximations to the
!          eigenvalues of A*z = lambda*B*z. The values are returned
!          in ascending order. If IPARAM(7) = 3,4,5 then D represents
!          the Ritz values of OP computed by qpsaupd transformed to
!          those of the original eigensystem A*z = lambda*B*z. If 
!          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
!          as the those of A*z = lambda*B*z.
!
!  Z       real(wrp) N by NEV array if HOWMNY = 'A'.  (OUTPUT)
!          On exit, Z contains the B-orthonormal Ritz vectors of the
!          eigensystem A*z = lambda*B*z corresponding to the Ritz
!          value approximations.
!          If  RVEC = .FALSE. then Z is not referenced.
!          NOTE: The array Z may be set equal to first NEV columns of the 
!          Arnoldi/Lanczos basis array V computed by DSAUPD.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
!
!  SIGMA   real(wrp)  (INPUT)
!          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
!          IPARAM(7) = 1 or 2.
!
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to DNAUPD that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, INFO
!
!         must be passed directly to DSEUPD following the last call
!         to DSAUPD.  These arguments MUST NOT BE MODIFIED between
!         the the last call to DSAUPD and the call to DSEUPD.
!
!  Two of these parameters (WORKL, INFO) are also output parameters:
!
!  WORKL   real(wrp) work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:4*ncv) contains information obtained in
!          qpsaupd.  They are not changed by qpseupd.
!          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
!          untransformed Ritz values, the computed error estimates,
!          and the associated eigenvector matrix of H.
!
!          Note: IPNTR(8:10) contains the pointer into WORKL for addresses
!          of the above information computed by qpseupd.
!          -------------------------------------------------------------
!          IPNTR(8): pointer to the NCV RITZ values of the original system.
!          IPNTR(9): pointer to the NCV corresponding error bounds.
!          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
!                     of the tridiagonal matrix T. Only referenced by
!                     qpseupd if RVEC = .TRUE. See Remarks.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV must be greater than NEV and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Information error from LAPACK routine qpsteqr.
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: NEV and WHICH = 'BE' are incompatible.
!          = -14: DSAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
!          = -16: HOWMNY = 'S' not yet implemented
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!
!\Remarks
!  1. The converged Ritz values are always returned in increasing 
!     (algebraic) order.
!
!  2. Currently only HOWMNY = 'A' is implemented. It is included at this
!     stage for the user who wants to incorporate it. 
!
!\Routines called:
!     qpsesrt  ARPACK routine that sorts an array X, and applies the
!             corresponding permutation to a matrix A.
!     qpsortr  qpsortr  ARPACK sorting routine.
!     wivout   ARPACK utility routine that prints integers.
!     qpvout   ARPACK utility routine that prints vectors.
!     qpgeqr2  LAPACK routine that computes the QR factorization of
!             a matrix.
!     qplacpy  LAPACK matrix copy routine.
!     qplamch  LAPACK routine that determines machine constants.
!     qporm2r  LAPACK routine that applies an orthogonal matrix in
!             factored form.
!     qpsteqr  LAPACK routine that computes eigenvalues and eigenvectors
!             of a tridiagonal matrix.
!     qpger    Level 2 BLAS rank one update to a matrix.
!     qpcopy   Level 1 BLAS that copies one vector to another .
!     qpnrm2   Level 1 BLAS that computes the norm of a vector.
!     qpscal   Level 1 BLAS that scales a vector.
!     qpswap   Level 1 BLAS that swaps the contents of two vectors.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Chao Yang                    Houston, Texas
!     Dept. of Computational & 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
! 
!\Revision history:
!     12/15/93: Version ' 2.1'
!
!\SCCS Information: @(#) 
! FILE: seupd.F   SID: 2.7   DATE OF SID: 8/27/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpseupd (
     $           rvec, howmny, lselect, d, z, ldz, sigma, bmat,
     $           n, which, nev, tol, resid, ncv, v, ldv, iparam,
     $           ipntr, workd, workl, lworkl, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      logical    
     $   rvec
      character  
     $   howmny, bmat
      character(len=2)
     $   which
      integer(wip) 
     $   ncv, nev, n, ldz, ldv, lworkl, iparam(7), ipntr(11),
     $   info
      logical    
     $   lselect(ncv)
      real(wrp)
     $   d(nev), z(ldz,nev), sigma, tol, resid(n), v(ldv,ncv),
     $   workd(2*n), workl(lworkl)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real(wrp)
     $   one, zero
      parameter (
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character(len=6) 
     $   ctyp
      integer(wip)
     $   j, k, ih, ritz, bounds, ldh, ldq, ihd, ihb, 
     $   iq, iw, next, irz, ibd, ism, ilg, ktrord,
     $   leftptr, rghtptr, nconv, mode, ierr 
      integer    
     $   msglvl
      real(wrp)
     $   bnorm2, rnorm, temp, thres1, thres2, tempbnd, eps23
      logical    
     $   reord
!
!     %--------------%
!     | Local Arrays |
!     %--------------%
!
      real(wrp) 
     $   kv(2)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   
     $   qpcopy, 
     $   qpger, 
     $   qpgeqr2, 
     $   qplacpy, 
     $   qporm2r, 
     $   qpscal, 
     $   qpsesrt, 
     $   qpsteqr, 
     $   qpswap, 
     $   qpvout, 
     $   wivout, 
     $   qpsortr
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      real(wrp)
     $   qpnrm2, 
     $   qplamch
      external   
     $   qpnrm2, 
     $   qplamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    
     $   min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
! 
!     %------------------------%
!     | Set default parameters |
!     %------------------------%
!
      msglvl = mseupd
      mode   = iparam(7)
      nconv  = iparam(5)
      info = 0
!
!     %--------------%
!     | Quick return |
!     %--------------%
!
      if (nconv .eq. 0) goto 9000
      ierr = 0
!
      if (nconv .le. 0)                        
     $   ierr = -14 
      if (n .le. 0)                            
     $   ierr = -1
      if (nev .le. 0)                          
     $   ierr = -2
      if (ncv .le. nev .or.  ncv .gt. n)       
     $   ierr = -3
      if (which .ne. 'LM' .and.
     $    which .ne. 'SM' .and.
     $    which .ne. 'LA' .and.
     $    which .ne. 'SA' .and.
     $    which .ne. 'BE')                     
     $   ierr = -5
      if (bmat .ne. 'I' .and. bmat .ne. 'G')   
     $   ierr = -6
      if ( (howmny .ne. 'A' .and.
     $      howmny .ne. 'P' .and.
     $      howmny .ne. 'S') .and. rvec ) 
     $   ierr = -15
      if (rvec .and. howmny .eq. 'S')  
     $   ierr = -16
!
      if (rvec .and. lworkl .lt. ncv**2+8*ncv) 
     $   ierr = -7
!     
      if (mode .eq. 1 .or. mode .eq. 2) then
         ctyp = 'REGULR'
      else if (mode .eq. 3 ) then
         ctyp = 'SHIFTI'
      else if (mode .eq. 4 ) then
         ctyp = 'BUCKLE'
      else if (mode .eq. 5 ) then
         ctyp = 'CAYLEY'
      else 
         ierr = -10
      endif
      if (mode .eq. 1 .and. bmat .eq. 'G')     
     $   ierr = -11
      if (nev .eq. 1 .and. which .eq. 'BE')    
     $   ierr = -12
!
!     %------------%
!     | Error Exit |
!     %------------%
!
      if (ierr .ne. 0) then
         info = ierr
         goto 9000
      endif
!     
!     %-------------------------------------------------------%
!     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
!     | etc... and the remaining workspace.                   |
!     | Also update pointer to be used on output.             |
!     | Memory is laid out as follows:                        |
!     | workl(1:2*ncv) := generated tridiagonal matrix H      |
!     |       The subdiagonal is stored in workl(2:ncv).      |
!     |       The dead spot is workl(1) but upon exiting      |
!     |       qpsaupd stores the B-norm of the last residual   |
!     |       vector in workl(1). We use this !!!             |
!     | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
!     |       The wanted values are in the first NCONV spots. |
!     | workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates   |
!     |       The wanted values are in the first NCONV spots. |
!     | NOTE: workl(1:4*ncv) is set by qpsaupd and is not      |
!     |       modified by qpseupd.                             |
!     %-------------------------------------------------------%
!
!     %-------------------------------------------------------%
!     | The following is used and set by qpseupd.              |
!     | workl(4*ncv+1:4*ncv+ncv) := used as workspace during  |
!     |       computation of the eigenvectors of H. Stores    |
!     |       the diagonal of H. Upon EXIT contains the NCV   |
!     |       Ritz values of the original system. The first   |
!     |       NCONV spots have the wanted values. If MODE =   |
!     |       1 or 2 then will equal workl(2*ncv+1:3*ncv).    |
!     | workl(5*ncv+1:5*ncv+ncv) := used as workspace during  |
!     |       computation of the eigenvectors of H. Stores    |
!     |       the subdiagonal of H. Upon EXIT contains the    |
!     |       NCV corresponding Ritz estimates of the         |
!     |       original system. The first NCONV spots have the |
!     |       wanted values. If MODE = 1,2 then will equal    |
!     |       workl(3*ncv+1:4*ncv).                           |
!     | workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is  |
!     |       the eigenvector matrix for H as returned by     |
!     |       qpsteqr. Not referenced if RVEC = .False.        |
!     |       Ordering follows that of workl(4*ncv+1:5*ncv)   |
!     | workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=         |
!     |       Workspace. Needed by qpsteqr and by qpseupd.      |
!     | GRAND total of NCV*(NCV+8) locations.                 |
!     %-------------------------------------------------------%
!
!
      ih     = ipntr(5)
      ritz   = ipntr(6)
      bounds = ipntr(7)
      ldh    = ncv
      ldq    = ncv
      ihd    = bounds + ldh
      ihb    = ihd    + ldh
      iq     = ihb    + ldh
      iw     = iq     + ldh*ncv
      next   = iw     + 2*ncv
      ipntr( 4) = next
      ipntr( 8) = ihd
      ipntr( 9) = ihb
      ipntr(10) = iq
!
!     %----------------------------------------%
!     | irz points to the Ritz values computed |
!     |     by _seigt before exiting _saup2.   |
!     | ibd points to the Ritz estimates       |
!     |     computed by _seigt before exiting  |
!     |     _saup2.                            |
!     %----------------------------------------%
!
      irz = ipntr(11) + ncv
      ibd = irz + ncv
!
!
!     %---------------------------------%
!     | Set machine dependent constant. |
!     %---------------------------------%
!
      eps23 = qplamch (
     $        'Epsilon-Machine' ) 
      eps23 = eps23**( 2.0_wrp / 3.0_wrp )
!
!     %---------------------------------------%
!     | RNORM is B-norm of the RESID(1:N).    |
!     | BNORM2 is the 2 norm of B*RESID(1:N). |
!     | Upon exit of qpsaupd WORKD(1:N) has    |
!     | B*RESID(1:N).                         |
!     %---------------------------------------%
!
      rnorm = workl(ih)
      if      (bmat .eq. 'I') then
         bnorm2 = rnorm
      else if (bmat .eq. 'G') then
         bnorm2 = qpnrm2 ( 
     $            n, workd, 1_wip )
      endif
!
      if (rvec) then
!
!        %------------------------------------------------%
!        | Get the converged Ritz value on the boundary.  |
!        | This value will be used to dermine whether we  |
!        | need to reorder the eigenvalues and            |
!        | eigenvectors comupted by _steqr, and is        |
!        | referred to as the "threshold" value.          |
!        |                                                |
!        | A Ritz value gamma is said to be a wanted      |
!        | one, if                                        |
!        | abs(gamma) .ge. threshold, when WHICH = 'LM';  |
!        | abs(gamma) .le. threshold, when WHICH = 'SM';  |
!        | gamma      .ge. threshold, when WHICH = 'LA';  |
!        | gamma      .le. threshold, when WHICH = 'SA';  |
!        | gamma .le. thres1 .or. gamma .ge. thres2       |
!        |                            when WHICH = 'BE';  |
!        |                                                |
!        | Note: converged Ritz values and associated     |
!        | Ritz estimates have been placed in the first   |
!        | NCONV locations in workl(ritz) and             |
!        | workl(bounds) respectively. They have been     |
!        | sorted (in _saup2) according to the WHICH      |
!        | lselection criterion. (Except in the case       |
!        | WHICH = 'BE', they are sorted in an increasing |
!        | order.)                                        |
!        %------------------------------------------------%
!
         if ( which .eq. 'LM' .or. 
     $        which .eq. 'SM' .or. 
     $        which .eq. 'LA' .or. 
     $        which .eq. 'SA' ) then
!
             thres1 = workl(ritz)
!
             if (msglvl .gt. 2) then
                call qpvout ( 
     $               logfil, 1_wip, (/thres1/), ndigit,
     $         '_seupd: Threshold eigenvalue used for re-ordering' )
             endif
!
         else if (which .eq. 'BE') then
!
!            %------------------------------------------------%
!            | Ritz values returned from _saup2 have been     |
!            | sorted in increasing order.  Thus two          |
!            | "threshold" values (one for the small end, one |
!            | for the large end) are in the middle.          |
!            %------------------------------------------------%
!
             ism = max(nev,nconv) / 2
             ilg = ism + 1
             thres1 = workl(ism)
             thres2 = workl(ilg) 
!
             if (msglvl .gt. 2) then
                kv(1) = thres1
                kv(2) = thres2
                call qpvout ( 
     $               logfil, 2_wip, kv, ndigit,
     $         '_seupd: Threshold eigenvalues used for re-ordering' )
             endif
!
         endif
!
!        %----------------------------------------------------------%
!        | Check to see if all converged Ritz values appear within  |
!        | the first NCONV diagonal elements returned from _seigt.  |
!        | This is done in the following way:                       |
!        |                                                          |
!        | 1) For each Ritz value obtained from _seigt, compare it  |
!        |    with the threshold Ritz value computed above to       |
!        |    determine whether it is a wanted one.                 |
!        |                                                          |
!        | 2) If it is wanted, then check the corresponding Ritz    |
!        |    estimate to see if it has converged.  If it has, set  |
!        |    correponding entry in the logical array SELECT to     |
!        |    .TRUE..                                               |
!        |                                                          |
!        | If SELECT(j) = .TRUE. and j > NCONV, then there is a     |
!        | converged Ritz value that does not appear at the top of  |
!        | the diagonal matrix computed by _seigt in _saup2.        |
!        | Reordering is needed.                                    |
!        %----------------------------------------------------------%
!
         reord = .false.
         ktrord = 0
         do j = 0, ncv-1
            lselect(j+1) = .false.
            if (which .eq. 'LM') then
               if (abs(workl(irz+j)) .ge. abs(thres1)) then
                   tempbnd = max( eps23, abs(workl(irz+j)) )
                   if (workl(ibd+j) .le. tol*tempbnd) then
                      lselect(j+1) = .true.
                   endif
               endif
            else if (which .eq. 'SM') then
               if (abs(workl(irz+j)) .le. abs(thres1)) then
                   tempbnd = max( eps23, abs(workl(irz+j)) )
                   if (workl(ibd+j) .le. tol*tempbnd) then
                      lselect(j+1) = .true.
                   endif
               endif
            else if (which .eq. 'LA') then
               if (workl(irz+j) .ge. thres1) then
                  tempbnd = max( eps23, abs(workl(irz+j)) )
                  if (workl(ibd+j) .le. tol*tempbnd) then
                     lselect(j+1) = .true.
                  endif
               endif
            else if (which .eq. 'SA') then
               if (workl(irz+j) .le. thres1) then
                  tempbnd = max( eps23, abs(workl(irz+j)) )
                  if (workl(ibd+j) .le. tol*tempbnd) then
                     lselect(j+1) = .true.
                  endif
               endif
            else if (which .eq. 'BE') then
               if ( workl(irz+j) .le. thres1 .or.
     $              workl(irz+j) .ge. thres2 ) then
                  tempbnd = max( eps23, abs(workl(irz+j)) )
                  if ( workl(ibd+j) .le. tol*tempbnd ) then
                     lselect(j+1) = .true.
                  endif
               endif
            endif
            if (j+1 .gt. nconv ) reord = lselect(j+1) .or. reord
            if (lselect(j+1)) ktrord = ktrord + 1
         enddo 

!        %-------------------------------------------%
!        | If KTRORD .ne. NCONV, something is qpong. |
!        %-------------------------------------------%
!
         if (msglvl .gt. 2) then
             call wivout ( 
     $            logfil, 1_wip, (/ktrord/), ndigit,
     $            '_seupd: Number of specified eigenvalues' )
             call wivout ( 
     $            logfil, 1_wip, (/nconv/), ndigit,
     $            '_seupd: Number of "converged" eigenvalues' )
         endif
!
!        %-----------------------------------------------------------%
!        | Call LAPACK routine _steqr to compute the eigenvalues and |
!        | eigenvectors of the final symmetric tridiagonal matrix H. |
!        | Initialize the eigenvector matrix Q to the identity.      |
!        %-----------------------------------------------------------%
!
         call qpcopy (
     $        ncv-1, workl(ih+1), 1_wip, workl(ihb), 1_wip )
         call qpcopy (
     $        ncv, workl(ih+ldh), 1_wip, workl(ihd), 1_wip )
!
         call qpsteqr (
     $        'Identity', ncv, workl(ihd), workl(ihb),
     $        workl(iq), ldq, workl(iw), ierr )
!
         if (ierr .ne. 0) then
            info = -8
            goto 9000
         endif
!
         if (msglvl .gt. 1) then
            call qpcopy (
     $           ncv, workl(iq+ncv-1), ldq, workl(iw), 1_wip )
            call qpvout (
     $           logfil, ncv, workl(ihd), ndigit,
     $          '_seupd: NCV Ritz values of the final H matrix' )
            call qpvout (
     $           logfil, ncv, workl(iw), ndigit,
     $           '_seupd: last row of the eigenvector matrix for H' )
         endif
!
         if (reord) then
!
!           %---------------------------------------------%
!           | Reordered the eigenvalues and eigenvectors  |
!           | computed by _steqr so that the "converged"  |
!           | eigenvalues appear in the first NCONV       |
!           | positions of workl(ihd), and the associated |
!           | eigenvectors appear in the first NCONV      |
!           | columns.                                    |
!           %---------------------------------------------%
!
            leftptr = 1
            rghtptr = ncv
!
            if (ncv .eq. 1) goto 30
!
 20         if (lselect(leftptr)) then
!
!              %-------------------------------------------%
!              | Search, from the left, for the first Ritz |
!              | value that has not converged.             |
!              %-------------------------------------------%
!
               leftptr = leftptr + 1
!
            else if ( .not. lselect(rghtptr)) then
!
!              %----------------------------------------------%
!              | Search, from the right, the first Ritz value |
!              | that has converged.                          |
!              %----------------------------------------------%
!
               rghtptr = rghtptr - 1
!
            else
!
!              %----------------------------------------------%
!              | Swap the Ritz value on the left that has not |
!              | converged with the Ritz value on the right   |
!              | that has converged.  Swap the associated     |
!              | eigenvector of the tridiagonal matrix H as   |
!              | well.                                        |
!              %----------------------------------------------%
!
               temp = workl(ihd+leftptr-1)
               workl(ihd+leftptr-1) = workl(ihd+rghtptr-1)
               workl(ihd+rghtptr-1) = temp
               call qpcopy ( 
     $              ncv, workl(iq+ncv*(leftptr-1)), 1_wip,
     $              workl(iw), 1_wip )
               call qpcopy (
     $              ncv, workl(iq+ncv*(rghtptr-1)), 1_wip,
     $              workl(iq+ncv*(leftptr-1)), 1_wip )
               call qpcopy ( 
     $              ncv, workl(iw), 1_wip,
     $              workl(iq+ncv*(rghtptr-1)), 1_wip )
               leftptr = leftptr + 1
               rghtptr = rghtptr - 1
!
            endif
!
            if (leftptr .lt. rghtptr) goto 20
!
 30      endif
!
         if (msglvl .gt. 2) then
             call qpvout (
     $            logfil, ncv, workl(ihd), ndigit,
     $            '_seupd: The eigenvalues of H--reordered' )
         endif
!
!        %----------------------------------------%
!        | Load the converged Ritz values into D. |
!        %----------------------------------------%
!
         call qpcopy ( 
     $        nconv, workl(ihd), 1_wip, d, 1_wip )
!
      else
!
!        %-----------------------------------------------------%
!        | Ritz vectors not required. Load Ritz values into D. |
!        %-----------------------------------------------------%
!
         call qpcopy (
     $        nconv, workl(ritz), 1_wip, d, 1_wip )
         call qpcopy (
     $        ncv, workl(ritz), 1_wip, workl(ihd), 1_wip )
!
      endif
!
!     %------------------------------------------------------------------%
!     | Transform the Ritz values and possibly vectors and corresponding |
!     | Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values |
!     | (and corresponding data) are returned in ascending order.        |
!     %------------------------------------------------------------------%
!
      if (ctyp .eq. 'REGULR') then
!
!        %---------------------------------------------------------%
!        | Ascending sort of wanted Ritz values, vectors and error |
!        | bounds. Not necessary if only Ritz values are desired.  |
!        %---------------------------------------------------------%
!
         if (rvec) then
            call qpsesrt (
     $           'LA', rvec , nconv, d, ncv, workl(iq), ldq )
         else
            call qpcopy (
     $           ncv, workl(bounds), 1_wip, workl(ihb), 1_wip )
         endif
!
      else 
! 
!        %-------------------------------------------------------------%
!        | *  Make a copy of all the Ritz values.                      |
!        | *  Transform the Ritz values back to the original system.   |
!        |    For TYPE = 'SHIFTI' the transformation is                |
!        |             lambda = 1/theta + sigma                        |
!        |    For TYPE = 'BUCKLE' the transformation is                |
!        |             lambda = sigma * theta / ( theta - 1 )          |
!        |    For TYPE = 'CAYLEY' the transformation is                |
!        |             lambda = sigma * (theta + 1) / (theta - 1 )     |
!        |    where the theta are the Ritz values returned by qpsaupd.  |
!        | NOTES:                                                      |
!        | *The Ritz vectors are not affected by the transformation.   |
!        |  They are only reordered.                                   |
!        %-------------------------------------------------------------%
!
         call qpcopy (
     $        ncv, workl(ihd), 1_wip, workl(iw), 1_wip )
         if (ctyp .eq. 'SHIFTI') then 
            do k = 1,ncv
               workl(ihd+k-1) = one / workl(ihd+k-1) + sigma
            enddo 
         else if (ctyp .eq. 'BUCKLE') then
            do k = 1,ncv
               workl(ihd+k-1) = sigma * 
     $                          workl(ihd+k-1) / 
     $                          ( workl(ihd+k-1) - one )
            enddo 
         else if (ctyp .eq. 'CAYLEY') then
            do k = 1,ncv
               workl(ihd+k-1) = sigma * 
     $                          ( workl(ihd+k-1) + one ) /
     $                          ( workl(ihd+k-1) - one )
            enddo 
         endif
! 
!        %-------------------------------------------------------------%
!        | *  Store the wanted NCONV lambda values into D.             |
!        | *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)   |
!        |    into ascending order and apply sort to the NCONV theta   |
!        |    values in the transformed system. We'll need this to     |
!        |    compute Ritz estimates in the original system.           |
!        | *  Finally sort the lambda's into ascending order and apply |
!        |    to Ritz vectors if wanted. Else just sort lambda's into  |
!        |    ascending order.                                         |
!        | NOTES:                                                      |
!        | *workl(iw:iw+ncv-1) contain the theta ordered so that they  |
!        |  match the ordering of the lambda. We'll use them again for |
!        |  Ritz vector purification.                                  |
!        %-------------------------------------------------------------%
!
         call qpcopy (
     $        nconv, workl(ihd), 1_wip, d, 1_wip )
         call qpsortr (
     $        'LA', .true., nconv, workl(ihd), workl(iw) )
         if (rvec) then
            call qpsesrt ( 
     $           'LA', rvec , nconv, d, ncv, workl(iq), ldq )
         else
            call qpcopy (
     $           ncv, workl(bounds), 1_wip, workl(ihb), 1_wip )
            call qpscal (
     $           ncv, bnorm2/rnorm, workl(ihb), 1_wip )
            call qpsortr (
     $           'LA', .true., nconv, d, workl(ihb) )
         endif
!
      endif 
! 
!     %------------------------------------------------%
!     | Compute the Ritz vectors. Transform the wanted |
!     | eigenvectors of the symmetric tridiagonal H by |
!     | the Lanczos basis matrix V.                    |
!     %------------------------------------------------%
!
      if (rvec .and. howmny .eq. 'A') then
!    
!        %----------------------------------------------------------%
!        | Compute the QR factorization of the matrix representing  |
!        | the wanted invariant subspace located in the first NCONV |
!        | columns of workl(iq,ldq).                                |
!        %----------------------------------------------------------%
!     
         call qpgeqr2 (
     $        ncv, nconv, workl(iq), ldq, workl(iw+ncv), 
     $        workl(ihb), ierr )
!
!     
!        %--------------------------------------------------------%
!        | * Postmultiply V by Q.                                 |   
!        | * Copy the first NCONV columns of VQ into Z.           |
!        | The N by NCONV matrix Z is now a matrix representation |
!        | of the approximate invariant subspace associated with  |
!        | the Ritz values in workl(ihd).                         |
!        %--------------------------------------------------------%
!     
         call qporm2r (
     $        'Right', 'Notranspose', n, ncv, nconv, workl(iq),
     $        ldq, workl(iw+ncv), v, ldv, workd(n+1), ierr )
         call qplacpy (
     $        'All', n, nconv, v, ldv, z, ldz )
!
!        %-----------------------------------------------------%
!        | In order to compute the Ritz estimates for the Ritz |
!        | values in both systems, need the last row of the    |
!        | eigenvector matrix. Remember, it's in factored form |
!        %-----------------------------------------------------%
!
         do j = 1, ncv-1
            workl(ihb+j-1) = zero 
         enddo 
         workl(ihb+ncv-1) = one
         call qporm2r (
     $        'Left', 'Transpose', ncv, 1_wip, nconv, workl(iq),
     $        ldq, workl(iw+ncv), workl(ihb), ncv, (/temp/), ierr )
!
      else if (rvec .and. howmny .eq. 'S') then
!
!     Not yet implemented. See remark 2 above.
!
      endif
!
      if (ctyp .eq. 'REGULR' .and. rvec) then
!
         do j = 1, ncv
            workl(ihb+j-1) = rnorm * abs( workl(ihb+j-1) )
         enddo 
!
      else if (ctyp .ne. 'REGULR' .and. rvec) then
!
!        %-------------------------------------------------%
!        | *  Determine Ritz estimates of the theta.       |
!        |    If RVEC = .true. then compute Ritz estimates |
!        |               of the theta.                     |
!        |    If RVEC = .false. then copy Ritz estimates   |
!        |              as computed by qpsaupd.             |
!        | *  Determine Ritz estimates of the lambda.      |
!        %-------------------------------------------------%
!
         call qpscal (
     $        ncv, bnorm2, workl(ihb), 1_wip )
         if (ctyp .eq. 'SHIFTI') then 
!
            do k = 1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) ) / 
     $                          workl(iw+k-1)**2
            enddo 
!
         else if (ctyp .eq. 'BUCKLE') then
!
            do k = 1, ncv
               workl(ihb+k-1) = sigma * 
     $                          abs(workl(ihb+k-1)) / 
     $                          ( workl(iw+k-1) - one )**2
            enddo 
!
         else if (ctyp .eq. 'CAYLEY') then
!
            do k = 1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) / 
     $                               workl(iw+k-1) * 
     $                               ( workl(iw+k-1) - one) )
            enddo 
!
         endif
!
      endif
!
      if ( ctyp .ne. 'REGULR' .and. 
     $     msglvl .gt. 1 ) then
         call qpvout (
     $        logfil, nconv, d, ndigit,
     $        '_seupd: Untransformed converged Ritz values')
         call qpvout (
     $        logfil, nconv, workl(ihb), ndigit, 
     $        '_seupd: Ritz estimates of the untransformed Ritz values')
      else if (msglvl .gt. 1) then
         call qpvout (
     $        logfil, nconv, d, ndigit,
     $        '_seupd: Converged Ritz values' )
         call qpvout (
     $        logfil, nconv, workl(ihb), ndigit, 
     $        '_seupd: Associated Ritz estimates')
      endif
! 
!     %-------------------------------------------------%
!     | Ritz vector purification step. Formally perform |
!     | one of inverse subspace iteration. Only used    |
!     | for MODE = 3,4,5. See reference 7               |
!     %-------------------------------------------------%
!
      if ( rvec .and. 
     $     ( ctyp .eq. 'SHIFTI' .or. 
     $       ctyp .eq. 'CAYLEY' )  ) then
!
         do k = 0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1) / workl(iw+k)
         enddo 
!
      else if ( rvec .and. 
     $          ctyp .eq. 'BUCKLE' ) then
!
         do k = 0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1) / (workl(iw+k)-one)
         enddo 
!
      endif 
!
      if (ctyp .ne. 'REGULR')
     $   call qpger (
     $        n, nconv, one, resid, 1_wip, workl(iw), 1_wip, z, ldz )
!
 9000 continue
      return
      end subroutine qpseupd
!
!***********************************************************************
!
!     --> qpsaupd:
!
!     +  qpsaup2  ARPACK routine that implements the Implicitly Restarted
!                Arnoldi Iteration.
!     +  qpstats  ARPACK routine that initialize timing and other statistics
!                variables.
!     +  wivout   ARPACK utility routine that prints integers.
!     +  qpsecond  ARPACK utility routine for timing.
!     +  qpvout   ARPACK utility routine that prints vectors.
!*    +  qplamch  LAPACK routine that determines machine constants.
!
!     --> qpseupd:
!
!     +  qpsesrt  ARPACK routine that sorts an array X, and applies the
!                corresponding permutation to a matrix A.
!     +  qpsortr  ARPACK sorting routine.
!!    +  wivout   ARPACK utility routine that prints integers.
!!    +  qpvout   ARPACK utility routine that prints vectors.
!?    +  qpgeqr2  LAPACK routine that computes the QR factorization of
!                a matrix.
!?    +  qplacpy  LAPACK matrix copy routine.
!*    +  qplamch  LAPACK routine that determines machine constants.
!?    +  qporm2r  LAPACK routine that applies an orthogonal matrix in
!                factored form.
!*    +  qpsteqr  LAPACK routine that computes eigenvalues and eigenvectors
!                of a tridiagonal matrix.
!*    +  qpger    Level 2 BLAS rank one update to a matrix.
!*    +  qpcopy   Level 1 BLAS that copies one vector to another .
!*    +  qpnrm2   Level 1 BLAS that computes the norm of a vector.
!*    +  qpscal   Level 1 BLAS that scales a vector.
!*    +  qpswap   Level 1 BLAS that swaps the contents of two vectors.
!
!-----------------------------------------------------------------------
!
!\BeginDoc
!
!\Name: qpsaup2
!
!\Description: 
!  Intermediate level interface called by qpsaupd.
!
!\Usage:
!  call qpsaup2 
!     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
!       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, 
!       IPNTR, WORKD, INFO )
!
!\Arguments
!
!  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in qpsaupd.
!  MODE, ISHIFT, MXITER: see the definition of IPARAM in qpsaupd.
!  
!  NP      Integer.  (INPUT/OUTPUT)
!          Contains the number of implicit shifts to apply during 
!          each Arnoldi/Lanczos iteration.  
!          If ISHIFT=1, NP is adjusted dynamically at each iteration 
!          to accelerate convergence and prevent stagnation.
!          This is also roughly equal to the number of matrix-vector 
!          products (involving the operator OP) per Arnoldi iteration.
!          The logic for adjusting is contained within the current
!          subroutine.
!          If ISHIFT=0, NP is the number of shifts the user needs
!          to provide via reverse comunication. 0 < NP < NCV-NEV.
!          NP may be less than NCV-NEV since a leading block of the current
!          upper Tridiagonal matrix has split off and contains "unwanted"
!          Ritz values.
!          Upon termination of the IRA iteration, NP contains the number 
!          of "converged" wanted Ritz values.
!
!  IUPD    Integer.  (INPUT)
!          IUPD .EQ. 0: use explicit restart instead implicit update.
!          IUPD .NE. 0: use implicit update.
!
!  V       real(wrp) N by (NEV+NP) array.  (INPUT/OUTPUT)
!          The Lanczos basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling 
!          program.
!
!  H       real(wrp) (NEV+NP) by 2 array.  (OUTPUT)
!          H is used to store the generated symmetric tridiagonal matrix
!          The subdiagonal is stored in the first column of H starting 
!          at H(2,1).  The main diagonal is stored in the qpsecond column
!          of H starting at H(1,2). If qpsaup2 converges store the 
!          B-norm of the final residual vector in H(1,1).
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling 
!          program.
!
!  RITZ    real(wrp) array of length NEV+NP.  (OUTPUT)
!          RITZ(1:NEV) contains the computed Ritz values of OP.
!
!  BOUNDS  real(wrp) array of length NEV+NP.  (OUTPUT)
!          BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.
!
!  Q       real(wrp) (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
!          Private (replicated) work array used to accumulate the 
!          rotation in the shift application step.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!          
!  WORKL   real(wrp) array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  It is used in the computation of the 
!          tridiagonal eigenvalue problem, the calculation and
!          application of the shifts and convergence checking.
!          If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations
!          of WORKL are used in reverse communication to hold the user 
!          supplied shifts.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD for 
!          vectors used by the Lanczos iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in one of  
!                    the spectral transformation modes.  X is the current
!                    operand.
!          -------------------------------------------------------------
!          
!  WORKD   real(wrp) work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Lanczos iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration !!!!!!!!!!
!          See Data Distribution Note in qpsaupd.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =     0: Normal return.
!          =     1: All possible eigenvalues of OP has been found.  
!                   NP returns the size of the invariant subspace
!                   spanning the operator OP. 
!          =     2: No shifts could be applied.
!          =    -8: Error return from trid. eigenvalue calculation;
!                   This should never happen.
!          =    -9: Starting vector is zero.
!          = -9999: Could not build an Lanczos factorization.
!                   Size that was built in returned in NP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!
!\Routines called:
!     qpgetv0  ARPACK initial vector generation routine. 
!     qpsaitr  ARPACK Lanczos factorization routine.
!     qpsapps  ARPACK application of implicit shifts routine.
!     qpsconv  ARPACK convergence of Ritz values routine.
!     qpseigt  ARPACK compute Ritz values and error bounds routine.
!     qpsgets  ARPACK reorder Ritz values and error bounds routine.
!     qpsortr  ARPACK sorting routine.
!     wivout   ARPACK utility routine that prints integers.
!     qpsecond  ARPACK utility routine for timing.
!     qpvout   ARPACK utility routine that prints vectors.
!     qplamch  LAPACK routine that determines machine constants.
!     qpcopy   Level 1 BLAS that copies one vector to another.
!     qpdot    Level 1 BLAS that computes the scalar product of two vectors. 
!     qpnrm2   Level 1 BLAS that computes the norm of a vector.
!     qpscal   Level 1 BLAS that scales a vector.
!     qpswap   Level 1 BLAS that swaps two vectors.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
! 
!\Revision history:
!     12/15/93: Version ' 2.4'
!     xx/xx/95: Version ' 2.4'.  (R.B. Lehoucq)
!
!\SCCS Information: @(#) 
! FILE: saup2.F   SID: 2.6   DATE OF SID: 8/16/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpsaup2 (
     $           ido, bmat, n, which, nev, np, tol, resid, mode, 
     $           iupd, ishift, mxiter, v, ldv, h, ldh, ritz, bounds,
     $           q, ldq, workl, ipntr, workd, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer    
     $   ido
      character(len=1) 
     $   bmat
      character(len=2)
     $   which
      integer(wip) 
     $   n, nev, np, mode, iupd, ishift, 
     $   ldh, ldq, ldv, mxiter, ipntr(3)
      real(wrp)
     $   tol, resid(n), bounds(nev+np), h(ldh,2), 
     $   q(ldq,nev+np), ritz(nev+np), v(ldv,nev+np), 
     $   workd(3*n), workl(3*(nev+np))
      integer(wip) 
     $   info
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real(wrp)
     $   one, zero
      parameter (
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character(len=2) 
     $   wprime
      logical    
     $   cnorm, getv0, initv, update, ushift
      integer(wip) 
     $   iter, j, kplusp, nconv, nevbef, nev0,
     $   np0, nptemp, nevd2, nevm2, kp(3), ierr 
      integer    
     $   msglvl
      real(wrp)
     $   rnorm, temp, eps23
      save       
     $   cnorm, getv0, initv, update, ushift,
     $   iter, kplusp, msglvl, nconv, nev0, np0,
     $   rnorm, eps23
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   
     $   qpcopy, 
     $   qpgetv0, 
     $   qpsaitr, 
     $   qpscal, 
     $   qpsconv, 
     $   qpseigt, 
     $   qpsgets, 
     $   qpsapps, 
     $   qpsortr, 
     $   qpvout, 
     $   wivout, 
     $   qpsecond, 
     $   qpswap
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      real(wrp)
     $   qpdot, 
     $   qpnrm2, 
     $   qplamch
      external   
     $   qpdot, 
     $   qpnrm2, 
     $   qplamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    
     $   min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido .eq. 0) then
! 
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call qpsecond (t0)
         msglvl = msaup2
!
!        %---------------------------------%
!        | Set machine dependent constant. |
!        %---------------------------------%
!
         eps23 = qplamch ( 
     $           'Epsilon-Machine' )
         eps23 = eps23**( 2.0_wrp / 3.0_wrp )
!
!        %-------------------------------------%
!        | nev0 and np0 are integer variables  |
!        | hold the initial values of NEV & NP |
!        %-------------------------------------%
!
         nev0   = nev
         np0    = np
!
!        %-------------------------------------%
!        | kplusp is the bound on the largest  |
!        |        Lanczos factorization built. |
!        | nconv is the current number of      |
!        |        "converged" eigenvlues.      |
!        | iter is the counter on the current  |
!        |      iteration step.                |
!        %-------------------------------------%
!
         kplusp = nev0 + np0
         nconv  = 0
         iter   = 0
! 
!        %--------------------------------------------%
!        | Set flags for computing the first NEV steps |
!        | of the Lanczos factorization.              |
!        %--------------------------------------------%
!
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
!
         if (info .ne. 0) then
!
!        %--------------------------------------------%
!        | User provides the initial residual vector. |
!        %--------------------------------------------%
!
            initv = .true.
            info  = 0
         else
            initv = .false.
         endif
      endif
! 
!     %---------------------------------------------%
!     | Get a possibly random starting vector and   |
!     | force it into the range of the operator OP. |
!     %---------------------------------------------%
!
   10 continue
!
      if (getv0) then
         call qpgetv0 (
     $        ido, bmat, 1_wip, initv, n, 1_wip, 
     $        v, ldv, resid, rnorm,
     $        ipntr, workd, info )
!
         if (ido .ne. 99) goto 9000
!
         if (rnorm .eq. zero) then
!
!           %-----------------------------------------%
!           | The initial vector is zero. Error exit. | 
!           %-----------------------------------------%
!
            info = -9
            goto 1200
         endif
         getv0 = .false.
         ido  = 0
      endif
! 
!     %------------------------------------------------------------%
!     | Back from reverse communication: continue with update step |
!     %------------------------------------------------------------%
!
      if ( update ) goto 20
!
!     %-------------------------------------------%
!     | Back from computing user specified shifts |
!     %-------------------------------------------%
!
      if ( ushift ) goto 50
!
!     %-------------------------------------%
!     | Back from computing residual norm   |
!     | at the end of the current iteration |
!     %-------------------------------------%
!
      if ( cnorm )  goto 100
! 
!     %----------------------------------------------------------%
!     | Compute the first NEV steps of the Lanczos factorization |
!     %----------------------------------------------------------%
!
      call qpsaitr (
     $     ido, bmat, n, 0_wip, nev0, mode, resid, rnorm, 
     $     v, ldv, h, ldh, ipntr, workd, info )
! 
!     %---------------------------------------------------%
!     | ido .ne. 99 implies use of reverse communication  |
!     | to compute operations involving OP and possibly B |
!     %---------------------------------------------------%
!
      if (ido .ne. 99) goto 9000
!
      if (info .gt. 0) then
!
!        %-----------------------------------------------------%
!        | qpsaitr was unable to build an Lanczos factorization |
!        | of length NEV0. INFO is returned with the size of   |
!        | the factorization built. Exit main loop.            |
!        %-----------------------------------------------------%
!
         np   = info
         mxiter = iter
         info = -9999
         goto 1200
      endif
! 
!     %--------------------------------------------------------------%
!     |                                                              |
!     |           M A I N  LANCZOS  I T E R A T I O N  L O O P       |
!     |           Each iteration implicitly restarts the Lanczos     |
!     |           factorization in place.                            |
!     |                                                              |
!     %--------------------------------------------------------------%
! 
 1000 continue
!
         iter = iter + 1
!
         if (msglvl .gt. 0) then
            call wivout (
     $           logfil, 1_wip, (/iter/), ndigit, 
     $      '_saup2: **** Start of major iteration number ****' )
         endif
         if (msglvl .gt. 1) then
            call wivout (
     $           logfil, 1_wip, (/nev/), ndigit, 
     $     '_saup2: The length of the current Lanczos factorization' )
            call wivout (
     $           logfil, 1_wip, (/np/), ndigit, 
     $      '_saup2: Extend the Lanczos factorization by' )
         endif
! 
!        %------------------------------------------------------------%
!        | Compute NP additional steps of the Lanczos factorization. |
!        %------------------------------------------------------------%
!
         ido = 0
   20    continue
         update = .true.
!
         call qpsaitr (
     $        ido, bmat, n, nev, np, mode, resid, rnorm, v, 
     $        ldv, h, ldh, ipntr, workd, info )
! 
!        %---------------------------------------------------%
!        | ido .ne. 99 implies use of reverse communication  |
!        | to compute operations involving OP and possibly B |
!        %---------------------------------------------------%
!
         if (ido .ne. 99) goto 9000
!
         if (info .gt. 0) then
!
!           %-----------------------------------------------------%
!           | qpsaitr was unable to build an Lanczos factorization |
!           | of length NEV0+NP0. INFO is returned with the size  |  
!           | of the factorization built. Exit main loop.         |
!           %-----------------------------------------------------%
!
            np = info
            mxiter = iter
            info = -9999
            goto 1200
         endif
         update = .false.
!
         if (msglvl .gt. 1) then
            call qpvout (
     $           logfil, 1_wip, (/rnorm/), ndigit, 
     $     '_saup2: Current B-norm of residual for factorization' )
         endif
! 
!        %--------------------------------------------------------%
!        | Compute the eigenvalues and corresponding error bounds |
!        | of the current symmetric tridiagonal matrix.           |
!        %--------------------------------------------------------%
!
         call qpseigt (
     $        rnorm, kplusp, h, ldh, ritz, bounds, workl, ierr )
!
         if (ierr .ne. 0) then
            info = -8
            goto 1200
         endif
!
!        %----------------------------------------------------%
!        | Make a copy of eigenvalues and corresponding error |
!        | bounds obtained from _seigt.                       |
!        %----------------------------------------------------%
!
         call qpcopy ( 
     $        kplusp, ritz, 1_wip, workl(kplusp+1), 1_wip )
         call qpcopy ( 
     $        kplusp, bounds, 1_wip, workl(2*kplusp+1), 1_wip )
!
!        %---------------------------------------------------%
!        | Select the wanted Ritz values and their bounds    |
!        | to be used in the convergence test.               |
!        | The selection is based on the requested number of |
!        | eigenvalues instead of the current NEV and NP to  |
!        | prevent possible misconvergence.                  |
!        | * Wanted Ritz values := RITZ(NP+1:NEV+NP)         |
!        | * Shifts := RITZ(1:NP) := WORKL(1:NP)             |
!        %---------------------------------------------------%
!
         nev = nev0
         np = np0
         call qpsgets (
     $        ishift, which, nev, np, ritz, bounds, workl )
! 
!        %-------------------%
!        | Convergence test. |
!        %-------------------%
!
         call qpcopy (
     $        nev, bounds(np+1), 1_wip, workl(np+1), 1_wip )
         call qpsconv (
     $        nev, ritz(np+1), workl(np+1), tol, nconv )
!
         if (msglvl .gt. 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = nconv
            call wivout (
     $           logfil, 3_wip, kp, ndigit,
     $      '_saup2: NEV, NP, NCONV are' )
            call qpvout (
     $           logfil, kplusp, ritz, ndigit,
     $      '_saup2: The eigenvalues of H' )
            call qpvout (
     $           logfil, kplusp, bounds, ndigit,
     $      '_saup2: Ritz estimates of the current NCV Ritz values' )
         endif
!
!        %---------------------------------------------------------%
!        | Count the number of unwanted Ritz values that have zero |
!        | Ritz estimates. If any Ritz estimates are equal to zero |
!        | then a leading block of H of order equal to at least    |
!        | the number of Ritz values with zero Ritz estimates has  |
!        | split off. None of these Ritz values may be removed by  |
!        | shifting. Decrease NP the number of shifts to apply. If |
!        | no shifts may be applied, then prepare to exit          |
!        %---------------------------------------------------------%
!
         nptemp = np
         do j = 1, nptemp
            if (bounds(j) .eq. zero) then
               np = np - 1
               nev = nev + 1
            endif
         enddo 
! 
         if ( (nconv .ge. nev0) .or. 
     $        (iter .gt. mxiter) .or.
     $        (np .eq. 0) ) then
!     
!           %------------------------------------------------%
!           | Prepare to exit. Put the converged Ritz values |
!           | and corresponding bounds in RITZ(1:NCONV) and  |
!           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
!           | careful when NCONV > NP since we don't want to |
!           | swap overlapping locations.                    |
!           %------------------------------------------------%
!
            if (which .eq. 'BE') then
!
!              %-----------------------------------------------------%
!              | Both ends of the spectrum are requested.            |
!              | Sort the eigenvalues into algebraically decreasing  |
!              | order first then swap low end of the spectrum next  |
!              | to high end in appropriate locations.               |
!              | NOTE: when np < floor(nev/2) be careful not to swap |
!              | overlapping locations.                              |
!              %-----------------------------------------------------%
!
               wprime = 'SA'
               call qpsortr (
     $              wprime, .true., kplusp, ritz, bounds)
               nevd2 = nev / 2
               nevm2 = nev - nevd2 
               if ( nev .gt. 1 ) then
                  call qpswap ( 
     $                 min(nevd2,np), ritz(nevm2+1), 1_wip,
     $                 ritz( max(kplusp-nevd2+1,kplusp-np+1) ), 1_wip)
                  call qpswap ( 
     $                 min(nevd2,np), bounds(nevm2+1), 1_wip,
     $                 bounds(max(kplusp-nevd2+1,kplusp-np)+1), 1_wip)
               endif
!
            else
!
!              %--------------------------------------------------%
!              | LM, SM, LA, SA case.                             |
!              | Sort the eigenvalues of H into the an order that |
!              | is opposite to WHICH, and apply the resulting    |
!              | order to BOUNDS.  The eigenvalues are sorted so  |
!              | that the wanted part are always within the first |
!              | NEV locations.                                   |
!              %--------------------------------------------------%
!
               if (which .eq. 'LM') wprime = 'SM'
               if (which .eq. 'SM') wprime = 'LM'
               if (which .eq. 'LA') wprime = 'SA'
               if (which .eq. 'SA') wprime = 'LA'
!
               call qpsortr (
     $              wprime, .true., kplusp, ritz, bounds )
!
            endif
!
!           %--------------------------------------------------%
!           | Scale the Ritz estimate of each Ritz value       |
!           | by 1 / max(eps23,magnitude of the Ritz value).   |
!           %--------------------------------------------------%
!
            do j = 1, nev0
               temp = max( eps23, abs(ritz(j)) )
               bounds(j) = bounds(j) / temp
            enddo 
!
!           %----------------------------------------------------%
!           | Sort the Ritz values according to the scaled Ritz  |
!           | esitmates.  This will push all the converged ones  |
!           | towards the front of ritzr, ritzi, bounds          |
!           | (in the case when NCONV < NEV.)                    |
!           %----------------------------------------------------%
!
            wprime = 'LA'
            call qpsortr ( 
     $           wprime, .true., nev0, bounds, ritz )
!
!           %----------------------------------------------%
!           | Scale the Ritz estimate back to its original |
!           | value.                                       |
!           %----------------------------------------------%
!
            do j = 1, nev0
                temp = max( eps23, abs(ritz(j)) )
                bounds(j) = bounds(j)*temp
            enddo 
!
!           %--------------------------------------------------%
!           | Sort the "converged" Ritz values again so that   |
!           | the "threshold" values and their associated Ritz |
!           | estimates appear at the appropriate position in  |
!           | ritz and bound.                                  |
!           %--------------------------------------------------%
!
            if (which .eq. 'BE') then
!
!              %------------------------------------------------%
!              | Sort the "converged" Ritz values in increasing |
!              | order.  The "threshold" values are in the      |
!              | middle.                                        |
!              %------------------------------------------------%
!
               wprime = 'LA'
               call qpsortr ( 
     $              wprime, .true., nconv, ritz, bounds )
!
            else
!
!              %----------------------------------------------%
!              | In LM, SM, LA, SA case, sort the "converged" |
!              | Ritz values according to WHICH so that the   |
!              | "threshold" value appears at the front of    |
!              | ritz.                                        |
!              %----------------------------------------------%

               call qpsortr ( 
     $              which, .true., nconv, ritz, bounds )
!
            endif
!
!           %------------------------------------------%
!           |  Use h( 1,1 ) as storage to communicate  |
!           |  rnorm to _seupd if needed               |
!           %------------------------------------------%
!
            h(1,1) = rnorm
!
            if (msglvl .gt. 1) then
               call qpvout (
     $              logfil, kplusp, ritz, ndigit,
     $              '_saup2: Sorted Ritz values.' )
               call qpvout (
     $              logfil, kplusp, bounds, ndigit,
     $              '_saup2: Sorted ritz estimates.' )
            endif
!
!           %------------------------------------%
!           | Max iterations have been exceeded. | 
!           %------------------------------------%
!
            if (iter .gt. mxiter .and. nconv .lt. nev) info = 1
!
!           %---------------------%
!           | No shifts to apply. | 
!           %---------------------%
!
            if (np .eq. 0 .and. nconv .lt. nev0) info = 2
!
            np = nconv
            goto 1100
!
         else if (nconv .lt. nev .and. ishift .eq. 1) then
!
!           %---------------------------------------------------%
!           | Do not have all the requested eigenvalues yet.    |
!           | To prevent possible stagnation, adjust the number |
!           | of Ritz values and the shifts.                    |
!           %---------------------------------------------------%
!
            nevbef = nev
            nev = nev + min (nconv, np/2)
            if (nev .eq. 1 .and. kplusp .ge. 6) then
               nev = kplusp / 2
            else if (nev .eq. 1 .and. kplusp .gt. 2) then
               nev = 2
            endif
            np  = kplusp - nev
!     
!           %---------------------------------------%
!           | If the size of NEV was just increased |
!           | resort the eigenvalues.               |
!           %---------------------------------------%
!     
            if (nevbef .lt. nev) 
     $         call qpsgets (
     $              ishift, which, nev, np, ritz, bounds, workl )
!
         endif
!
         if (msglvl .gt. 0) then
            call wivout (
     $           logfil, 1_wip, (/nconv/), ndigit,
     $      '_saup2: no. of "converged" Ritz values at this iter.' )
            if (msglvl .gt. 1) then
               kp(1) = nev
               kp(2) = np
               call wivout (
     $              logfil, 2_wip, kp, ndigit,
     $         '_saup2: NEV and NP are' )
               call qpvout (
     $              logfil, nev, ritz(np+1), ndigit,
     $         '_saup2: "wanted" Ritz values.' )
               call qpvout (
     $              logfil, nev, bounds(np+1), ndigit,
     $         '_saup2: Ritz estimates of the "wanted" values ' )
            endif
         endif

! 
         if (ishift .eq. 0) then
!
!           %-----------------------------------------------------%
!           | User specified shifts: reverse communication to     |
!           | compute the shifts. They are returned in the first  |
!           | NP locations of WORKL.                              |
!           %-----------------------------------------------------%
!
            ushift = .true.
            ido = 3
            goto 9000
         endif
!
   50    continue
!
!        %------------------------------------%
!        | Back from reverse communication;   |
!        | User specified shifts are returned |
!        | in WORKL(1:*NP)                   |
!        %------------------------------------%
!
         ushift = .false.
! 
! 
!        %---------------------------------------------------------%
!        | Move the NP shifts to the first NP locations of RITZ to |
!        | free up WORKL.  This is for the non-exact shift case;   |
!        | in the exact shift case, qpsgets already handles this.   |
!        %---------------------------------------------------------%
!
         if (ishift .eq. 0) 
     $      call qpcopy (
     $           np, workl, 1_wip, ritz, 1_wip )
!
         if (msglvl .gt. 2) then
            call wivout (
     $           logfil, 1_wip, (/np/), ndigit,
     $           '_saup2: The number of shifts to apply ' )
            call qpvout (
     $           logfil, np, workl, ndigit,
     $           '_saup2: shifts selected' )
            if (ishift .eq. 1) then
               call qpvout (
     $              logfil, np, bounds, ndigit,
     $              '_saup2: corresponding Ritz estimates' )
             endif
         endif
! 
!        %---------------------------------------------------------%
!        | Apply the NP0 implicit shifts by QR bulge chasing.      |
!        | Each shift is applied to the entire tridiagonal matrix. |
!        | The first 2*N locations of WORKD are used as workspace. |
!        | After qpsapps is done, we have a Lanczos                 |
!        | factorization of length NEV.                            |
!        %---------------------------------------------------------%
!
         call qpsapps (
     $        n, nev, np, ritz, v, ldv, h, ldh, 
     $        resid, q, ldq, workd )
!
!        %---------------------------------------------%
!        | Compute the B-norm of the updated residual. |
!        | Keep B*RESID in WORKD(1:N) to be used in    |
!        | the first step of the next call to qpsaitr.  |
!        %---------------------------------------------%
!
         cnorm = .true.
         call qpsecond (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call qpcopy (
     $           n, resid, 1_wip, workd(n+1), 1_wip )
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
! 
!           %----------------------------------%
!           | Exit in order to compute B*RESID |
!           %----------------------------------%
! 
            goto 9000
         else if (bmat .eq. 'I') then
            call qpcopy (
     $           n, resid, 1_wip, workd, 1_wip )
         endif
! 
  100    continue
! 
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(1:N) := B*RESID            |
!        %----------------------------------%
!
         if (bmat .eq. 'G') then
            call qpsecond (t3)
            tmvbx = tmvbx + (t3 - t2)
         endif
! 
         if (bmat .eq. 'G') then         
            rnorm = qpdot (
     $              n, resid, 1_wip, workd, 1_wip )
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = qpnrm2 (
     $              n, resid, 1_wip )
         endif
         cnorm = .false.
  130    continue
!
         if (msglvl .gt. 2) then
            call qpvout (
     $           logfil, 1_wip, (/rnorm/), ndigit, 
     $      '_saup2: B-norm of residual for NEV factorization' )
            call qpvout (
     $           logfil, nev, h(1,2), ndigit,
     $           '_saup2: main diagonal of compressed H matrix' )
            call qpvout (
     $           logfil, nev-1, h(2,1), ndigit,
     $           '_saup2: subdiagonal of compressed H matrix' )
         endif
! 
      goto 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
! 
 1100 continue
!
      mxiter = iter
      nev = nconv
! 
 1200 continue
      ido = 99
!
!     %------------%
!     | Error exit |
!     %------------%
!
      call qpsecond (t1)
      tsaup2 = t1 - t0
! 
 9000 continue
      return
      end subroutine qpsaup2
!
!-----------------------------------------------------------------------
!
!\SCCS Information: @(#) 
! FILE: stats.F   SID: 2.1   DATE OF SID: 4/19/96   RELEASE: 2
!     %---------------------------------------------%
!     | Initialize statistic and timing information |
!     | for symmetric Arnoldi code.                 |
!     %---------------------------------------------%
 
      subroutine qpstats

!     %--------------------------------%
!     | See stat.doc for documentation |
!     %--------------------------------%
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
 
      tsaupd = 0.0D+0
      tsaup2 = 0.0D+0
      tsaitr = 0.0D+0
      tseigt = 0.0D+0
      tsgets = 0.0D+0
      tsapps = 0.0D+0
      tsconv = 0.0D+0
      titref = 0.0D+0
      tgetv0 = 0.0D+0
      trvec  = 0.0D+0
 
!     %----------------------------------------------------%
!     | User time including reverse communication overhead |
!     %----------------------------------------------------%
      tmvopx = 0.0D+0
      tmvbx  = 0.0D+0
      return
      end subroutine qpstats
!
!-----------------------------------------------------------------------
!
      subroutine qpsecond ( t )
      implicit none 
!
      real       
     $   t
!
!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     July 26, 1991
!
!     Edit:
!     +  2023-03-25_10-17-53
!        tqviet edit this subroutine using the cpu_time intrinsic function
!
!  Purpose
!  =======
!
!  SECOND returns the user time for a process in qpseconds.
!  This version gets the time from the system function ETIME.
!
!     .. Local Scalars ..
!     REAL               T1
!     ..
!     .. Local Arrays ..
!     REAL               TARRAY( 2 )
!     ..
!     .. External Functions ..
!     REAL               ETIME
!     EXTERNAL           ETIME
!     ..
!     .. Executable Statements ..
!
!
!     T1 = ETIME( TARRAY )
!     T  = TARRAY( 1 )
!
      call cpu_time ( t ) 
!
      return
      end subroutine qpsecond 
!
!-----------------------------------------------------------------------
!
!  Routine:    DVOUT
!
!  Purpose:    Real vector output routine.
!
!  Usage:      CALL DVOUT (LOUT, N, SX, IDIGIT, IFMT)
!
!  Arguments
!     N      - Length of array SX.  (Input)
!     SX     - Real array to be printed.  (Input)
!     IFMT   - Format to be used in printing array SX.  (Input)
!     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
!              If IDIGIT .LT. 0, printing is done with 72 columns.
!              If IDIGIT .GT. 0, printing is done with 132 columns.
!
!-----------------------------------------------------------------------
!
      subroutine qpvout ( 
     $           lout, n, sx, idigit, ifmt )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer
     $   lout
      integer(wip)    
     $   n
      integer
     $   idigit
      real(wrp)
     $   sx(*) 
      character(len=*)
     $   ifmt
!
!     LOCAL VARIABLES:
      character(len=80)
     $   line
      integer            
     $   i, lll, ndigit
      integer(wip)
     $   k1, k2
!
!     DEPENDENCES: 
      intrinsic          
     $   len, min
!
!     EXECUTABLE STATEMENTS: 
!
      lll = min( len( ifmt ), 80 )
      do i = 1, lll
         line( i: i ) = '-'
      enddo 
!
      do i = lll + 1, 80
         line( i: i ) = ' '
      enddo 
!
      write( lout, fmt = 9999 ) ifmt, line( 1: lll )
 9999 format( / 1x, a, / 1x, a )
!
      if ( n.le.0 ) return
      ndigit = idigit
      if( idigit .eq. 0 ) ndigit = 4
!
!=======================================================================
!             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
!=======================================================================
!
      if( idigit .lt. 0 ) then
         ndigit = -idigit
         if( ndigit.le.4 ) then
            do k1 = 1, n, 5
               k2 = min( n, k1+4 )
               write( lout, fmt = 9998 ) k1, k2, ( sx( i ), i = k1, k2 )
            enddo 
         else if( ndigit.le.6 ) then
            do k1 = 1, n, 4
               k2 = min( n, k1+3 )
               write( lout, fmt = 9997 ) k1, k2, ( sx( i ), i = k1, k2 )
            enddo 
         else if( ndigit.le.10 ) then
            do k1 = 1, n, 3
               k2 = min( n, k1+2 )
               write( lout, fmt = 9996 ) k1, k2, ( sx( i ), i = k1, k2 )
            enddo 
         else
            do k1 = 1, n, 2
               k2 = min( n, k1+1 )
               write( lout, fmt = 9995 ) k1, k2, ( sx( i ), i = k1, k2 )
            enddo 
         endif
!
!=======================================================================
!             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
!=======================================================================
!
      else
         if( ndigit.le.4 ) then
            do k1 = 1, n, 10
               k2 = min( n, k1+9 )
               write( lout, fmt = 9998 ) k1, k2, ( sx( i ), i = k1, k2 )
            enddo 
         else if( ndigit.le.6 ) then
            do k1 = 1, n, 8
               k2 = min( n, k1+7 )
               write( lout, fmt = 9997 ) k1, k2, ( sx( i ), i = k1, k2 )
            enddo 
         else if( ndigit.le.10 ) then
            do k1 = 1, n, 6
               k2 = min( n, k1+5 )
               write( lout, fmt = 9996 ) k1, k2, ( sx( i ), i = k1, k2 )
            enddo 
         else
            do k1 = 1, n, 5
               k2 = min( n, k1+4 )
               write( lout, fmt = 9995 ) k1, k2, ( sx( i ), i = k1, k2 )
            enddo 
         endif
      endif
      write( lout, fmt = 9994 )
      return
 9998 format( 1x, i4, ' - ', i4, ':',     1p, 10e12.3 )
 9997 format( 1x, i4, ' - ', i4, ':', 1x, 1p,  8e14.5 )
 9996 format( 1x, i4, ' - ', i4, ':', 1x, 1p,  6e18.9 )
 9995 format( 1x, i4, ' - ', i4, ':', 1x, 1p,  5e24.13 )
 9994 format( 1x, ' ' )
      end subroutine qpvout  
!
!-----------------------------------------------------------------------
!======================================================================
!
!     Duplicated by WFLAMCH of inc_lapack.f 
!
!     DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!======================================================================
!-----------------------------------------------------------------------
!
!\BeginDoc
!
!\Name: qpsesrt
!
!\Description:
!  Sort the array X in the order specified by WHICH and optionally 
!  apply the permutation to the columns of the matrix A.
!
!\Usage:
!  call qpsesrt
!     ( WHICH, APPLY, N, X, NA, A, LDA)
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> X is sorted into increasing order of magnitude.
!          'SM' -> X is sorted into decreasing order of magnitude.
!          'LA' -> X is sorted into increasing order of algebraic.
!          'SA' -> X is sorted into decreasing order of algebraic.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to A.
!          APPLY = .FALSE. -> do not apply the sorted order to A.
!
!  N       Integer.  (INPUT)
!          Dimension of the array X.
!
!  X      real(wrp) array of length N.  (INPUT/OUTPUT)
!          The array to be sorted.
!
!  NA      Integer.  (INPUT)
!          Number of rows of the matrix A.
!
!  A      real(wrp) array of length NA by N.  (INPUT/OUTPUT)
!         
!  LDA     Integer.  (INPUT)
!          Leading dimension of A.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines
!     qpswap  Level 1 BLAS that swaps the contents of two vectors.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Dept. of Computational &     Houston, Texas 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     12/15/93: Version ' 2.1'.
!               Adapted from the sort routine in LANSO and 
!               the ARPACK code qpsortr
!
!\SCCS Information: @(#) 
! FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpsesrt (
     $           which, apply, n, x, na, a, lda )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character(len=2) 
     $   which
      logical    
     $   apply
      integer(wip)
     $   lda, n, na
      real(wrp)
     $   x(0:n-1), a(lda,0:n-1)
!
!     LOCAL VARIABLES:
      integer(wip) 
     $   i, igap, j
      real(wrp)
     $   temp
!
!     DEPENDENCES:
      external   
     $   qpswap
!
!     EXECUTABLE STATEMENTS: 
!
      igap = n / 2
! 
      if (which .eq. 'SA') then
!
!        X is sorted into decreasing order of algebraic.
!
   10    continue
         if (igap .eq. 0) goto 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
!
            if ( j .lt. 0 ) goto 30
!
            if ( x(j) .lt. x(j+igap) ) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) 
     $            call qpswap ( 
     $                 na, a(1,j), 1_wip, a(1,j+igap), 1_wip )
            else
               goto 30
            endif
            j = j-igap
            goto 20
   30    continue
         igap = igap / 2
         goto 10
!
      else if (which .eq. 'SM') then
!
!        X is sorted into decreasing order of magnitude.
!
   40    continue
         if (igap .eq. 0) goto 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
!
            if (j.lt.0) goto 60
!
            if ( abs(x(j)) .lt. abs(x(j+igap)) ) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) 
     $            call qpswap ( 
     $                 na, a(1,j), 1_wip, a(1,j+igap), 1_wip )
            else
               goto 60
            endif
            j = j-igap
            goto 50
   60    continue
         igap = igap / 2
         goto 40
!
      else if (which .eq. 'LA') then
!
!        X is sorted into increasing order of algebraic.
!
   70    continue
         if (igap .eq. 0) goto 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
!
            if (j.lt.0) goto 90
!           
            if (x(j).gt.x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) 
     $            call qpswap ( 
     $                 na, a(1,j), 1_wip, a(1,j+igap), 1_wip )
            else
               goto 90
            endif
            j = j-igap
            goto 80
   90    continue
         igap = igap / 2
         goto 70
! 
      else if (which .eq. 'LM') then
!
!        X is sorted into increasing order of magnitude.
!
  100    continue
         if (igap .eq. 0) goto 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
!
            if (j.lt.0) goto 120
!
            if (abs(x(j)).gt.abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) 
     $            call qpswap ( 
     $                 na, a(1,j), 1_wip, a(1,j+igap), 1_wip )
            else
               goto 120
            endif
            j = j-igap
            goto 110
  120    continue
         igap = igap / 2
         goto 100
      endif
!
 9000 continue
      return
      end subroutine qpsesrt
!
!-----------------------------------------------------------------------
!
!\BeginDoc
!
!\Name: qpsortr
!
!\Description:
!  Sort the array X1 in the order specified by WHICH and optionally 
!  applies the permutation to the array X2.
!
!\Usage:
!  call qpsortr
!     ( WHICH, APPLY, N, X1, X2 )
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> X1 is sorted into increasing order of magnitude.
!          'SM' -> X1 is sorted into decreasing order of magnitude.
!          'LA' -> X1 is sorted into increasing order of algebraic.
!          'SA' -> X1 is sorted into decreasing order of algebraic.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to X2.
!          APPLY = .FALSE. -> do not apply the sorted order to X2.
!
!  N       Integer.  (INPUT)
!          Size of the arrays.
!
!  X1      real(wrp) array of length N.  (INPUT/OUTPUT)
!          The array to be sorted.
!
!  X2      real(wrp) array of length N.  (INPUT/OUTPUT)
!          Only referenced if APPLY = .TRUE.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Dept. of Computational &     Houston, Texas 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     12/16/93: Version ' 2.1'.
!               Adapted from the sort routine in LANSO.
!
!\SCCS Information: @(#) 
! FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpsortr (
     $           which, apply, n, x1, x2 )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character(len=2) 
     $   which
      logical    
     $   apply
      integer(wip)
     $   n
      real(wrp)
     $   x1(0:n-1), x2(0:n-1)
!
!     LOCAL VARIABLES:
      integer(wip) 
     $   i, igap, j
      real(wrp)
     $   temp
!
!
      igap = n / 2
! 
      if (which .eq. 'SA') then
!
!        X1 is sorted into decreasing order of algebraic.
!
   10    continue
         if (igap .eq. 0) goto 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
!
            if (j.lt.0) goto 30
!
            if ( x1(j) .lt. x1(j+igap) ) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               endif
            else
               goto 30
            endif
            j = j-igap
            goto 20
   30    continue
         igap = igap / 2
         goto 10
!
      else if (which .eq. 'SM') then
!
!        X1 is sorted into decreasing order of magnitude.
!
   40    continue
         if (igap .eq. 0) goto 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
!
            if (j.lt.0) goto 60
!
            if ( abs(x1(j)) .lt. abs(x1(j+igap)) ) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               endif
            else
               goto 60
            endif
            j = j-igap
            goto 50
   60    continue
         igap = igap / 2
         goto 40
!
      else if (which .eq. 'LA') then
!
!        X1 is sorted into increasing order of algebraic.
!
   70    continue
         if (igap .eq. 0) goto 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
!
            if (j.lt.0) goto 90
!           
            if ( x1(j) .gt. x1(j+igap) ) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               endif
            else
               goto 90
            endif
            j = j-igap
            goto 80
   90    continue
         igap = igap / 2
         goto 70
! 
      else if (which .eq. 'LM') then
!
!        X1 is sorted into increasing order of magnitude.
!
  100    continue
         if (igap .eq. 0) goto 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
!
            if (j.lt.0) goto 120
!
            if ( abs(x1(j)) .gt. abs(x1(j+igap)) ) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               endif
            else
               goto 120
            endif
            j = j-igap
            goto 110
  120    continue
         igap = igap / 2
         goto 100
      endif
!
 9000 continue
      return
      end subroutine qpsortr
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: qpgetv0
!
!\Description: 
!  Generate a random initial residual vector for the Arnoldi process.
!  Force the residual vector to be in the range of the operator OP.  
!
!\Usage:
!  call qpgetv0
!     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM, 
!       IPNTR, WORKD, IERR )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to qpgetv0.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B in the (generalized)
!          eigenvalue problem A*x = lambda*B*x.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  ITRY    Integer.  (INPUT)
!          ITRY counts the number of times that qpgetv0 is called.  
!          It should be set to 1 on the initial call to qpgetv0.
!
!  INITV   Logical variable.  (INPUT)
!          .TRUE.  => the initial residual vector is given in RESID.
!          .FALSE. => generate a random initial residual vector.
!
!  N       Integer.  (INPUT)
!          Dimension of the problem.
!
!  J       Integer.  (INPUT)
!          Index of the residual vector to be generated, with respect to
!          the Arnoldi process.  J > 1 in case of a "restart".
!
!  V       real(wrp) N by J array.  (INPUT)
!          The first J-1 columns of V contain the current Arnoldi basis
!          if this is a "restart".
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling 
!          program.
!
!  RESID   real(wrp) array of length N.  (INPUT/OUTPUT)
!          Initial residual vector to be generated.  If RESID is 
!          provided, force RESID into the range of the operator OP.
!
!  RNORM   real(wrp) scalar.  (OUTPUT)
!          B-norm of the generated residual.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!
!  WORKD   real(wrp) work array of length 2*N.  (REVERSE COMMUNICATION).
!          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
!
!  IERR    Integer.  (OUTPUT)
!          =  0: Normal exit.
!          = -1: Cannot generate a nontrivial restarted residual vector
!                in the range of the operator OP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     qpsecond  ARPACK utility routine for timing.
!     qpvout   ARPACK utility routine for vector output.
!     qplarnv  LAPACK routine for generating a random vector.
!     qpgemv   Level 2 BLAS routine for matrix vector multiplication.
!     qpcopy   Level 1 BLAS that copies one vector to another.
!     qpdot    Level 1 BLAS that computes the scalar product of two vectors. 
!     qpnrm2   Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\SCCS Information: @(#) 
! FILE: getv0.F   SID: 2.6   DATE OF SID: 8/27/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpgetv0 ( 
     $           ido, bmat, itry, initv, n, j, v, ldv, 
     $           resid, rnorm, ipntr, workd, ierr )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer    
     $   ido
      character(len=1)  
     $   bmat
      logical    
     $   initv
      integer(wip) 
     $   itry, j, ldv, n, ipntr(3)
      real(wrp)
     $   rnorm, resid(n), v(ldv,j), workd(2*n)
      integer(wip) 
     $   ierr
! 
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real(wrp)
     $   one, zero
      parameter (
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      logical    
     $   first, inits, orth
      integer    
     $   idist, iseed(4), msglvl, iter 
      integer(wip) 
     $   jj
      real(wrp)
     $   rnorm0
      save 
     $   first, iseed, inits, iter, msglvl, orth, rnorm0
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external
     $   qplarnv, 
     $   qpvout, 
     $   qpcopy, 
     $   qpgemv, 
     $   qpsecond
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      real(wrp)
     $   qpdot, 
     $   qpnrm2
      external 
     $   qpdot, 
     $   qpnrm2
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic 
     $   abs, sqrt
!
!     %-----------------%
!     | Data Statements |
!     %-----------------%
!
      data
     $   inits /.true./
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!
!     %-----------------------------------%
!     | Initialize the seed of the LAPACK |
!     | random number generator           |
!     %-----------------------------------%
!
      if (inits) then
          iseed(1) = 1
          iseed(2) = 3
          iseed(3) = 5
          iseed(4) = 7
          inits = .false.
      endif
!
      if (ido .eq.  0) then
! 
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call qpsecond (t0)
         msglvl = mgetv0
! 
         ierr   = 0
         iter   = 0
         first  = .false.
         orth   = .false.
!
!        %-----------------------------------------------------%
!        | Possibly generate a random starting vector in RESID |
!        | Use a LAPACK random number generator used by the    |
!        | matrix generation routines.                         |
!        |    idist = 1: uniform (0,1)  distribution;          |
!        |    idist = 2: uniform (-1,1) distribution;          |
!        |    idist = 3: normal  (0,1)  distribution;          |
!        %-----------------------------------------------------%
!
         if (.not.initv) then
            idist = 2
            call qplarnv (
     $           idist, iseed, n, resid )
         endif
! 
!        %----------------------------------------------------------%
!        | Force the starting vector into the range of OP to handle |
!        | the generalized problem when B is possibly (singular).   |
!        %----------------------------------------------------------%
!
         call qpsecond (t2)
         if (bmat .eq. 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call qpcopy (
     $           n, resid, 1_wip, workd, 1_wip )
            ido = -1
            goto 9000
         endif
      endif
! 
!     %-----------------------------------------%
!     | Back from computing OP*(initial-vector) |
!     %-----------------------------------------%
!
      if (first) goto 20
!
!     %-----------------------------------------------%
!     | Back from computing B*(orthogonalized-vector) |
!     %-----------------------------------------------%
!
      if (orth)  goto 40
! 
      if (bmat .eq. 'G') then
         call qpsecond (t3)
         tmvopx = tmvopx + (t3 - t2)
      endif
! 
!     %------------------------------------------------------%
!     | Starting vector is now in the range of OP; r = OP*r; |
!     | Compute B-norm of starting vector.                   |
!     %------------------------------------------------------%
!
      call qpsecond (t2)
      first = .TRUE.
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call qpcopy (
     $        n, workd(n+1), 1_wip, resid, 1_wip )
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         goto 9000
      else if (bmat .eq. 'I') then
         call qpcopy (
     $        n, resid, 1_wip, workd, 1_wip )
      endif
! 
   20 continue
!
      if (bmat .eq. 'G') then
         call qpsecond (t3)
         tmvbx = tmvbx + (t3 - t2)
      endif
! 
      first = .false.
      if (bmat .eq. 'G') then
          rnorm0 = qpdot (
     $             n, resid, 1_wip, workd, 1_wip )
          rnorm0 = sqrt(abs(rnorm0))
      else if (bmat .eq. 'I') then
           rnorm0 = qpnrm2 (
     $              n, resid, 1_wip )
      endif
      rnorm  = rnorm0
!
!     %---------------------------------------------%
!     | Exit if this is the very first Arnoldi step |
!     %---------------------------------------------%
!
      if (j .eq. 1) goto 50
! 
!     %----------------------------------------------------------------
!     | Otherwise need to B-orthogonalize the starting vector against |
!     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
!     | This is the case where an invariant subspace is encountered   |
!     | in the middle of the Arnoldi factorization.                   |
!     |                                                               |
!     |       s = V^{T}*B*r;   r = r - V*s;                           |
!     |                                                               |
!     | Stopping criteria used for iter. ref. is discussed in         |
!     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
!     %---------------------------------------------------------------%
!
      orth = .TRUE.
   30 continue
!
      call qpgemv (
     $     'T', n, j-1, one, v, ldv, workd, 1_wip, 
     $     zero, workd(n+1), 1_wip )
      call qpgemv (
     $     'N', n, j-1, -one, v, ldv, workd(n+1), 1_wip, 
     $     one, resid, 1_wip )
! 
!     %----------------------------------------------------------%
!     | Compute the B-norm of the orthogonalized starting vector |
!     %----------------------------------------------------------%
!
      call qpsecond (t2)
      if (bmat .eq. 'G') then
         nbx = nbx + 1
         call qpcopy (
     $        n, resid, 1_wip, workd(n+1), 1_wip )
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         goto 9000
      else if (bmat .eq. 'I') then
         call qpcopy (
     $        n, resid, 1_wip, workd, 1_wip )
      endif
! 
   40 continue
!
      if (bmat .eq. 'G') then
         call qpsecond (t3)
         tmvbx = tmvbx + (t3 - t2)
      endif
! 
      if (bmat .eq. 'G') then
         rnorm = qpdot (
     $           n, resid, 1_wip, workd, 1_wip)
         rnorm = sqrt(abs(rnorm))
      else if (bmat .eq. 'I') then
         rnorm = qpnrm2 ( 
     $           n, resid, 1_wip )
      endif
!
!     %--------------------------------------%
!     | Check for further orthogonalization. |
!     %--------------------------------------%
!
      if (msglvl .gt. 2) then
          call qpvout (
     $         logfil, 1_wip, (/rnorm0/), ndigit, 
     $         '_getv0: re-orthonalization ; rnorm0 is' )
          call qpvout (
     $         logfil, 1_wip, (/rnorm/), ndigit, 
     $         '_getv0: re-orthonalization ; rnorm is' )
      endif
!
      if (rnorm .gt. 0.717_wrp*rnorm0) goto 50
! 
      iter = iter + 1
      if (iter .le. 1) then
!
!        %-----------------------------------%
!        | Perform iterative refinement step |
!        %-----------------------------------%
!
         rnorm0 = rnorm
         goto 30
      else
!
!        %------------------------------------%
!        | Iterative refinement step "failed" |
!        %------------------------------------%
!
         do jj = 1, n
            resid(jj) = zero
         enddo 
         rnorm = zero
         ierr = -1
      endif
! 
   50 continue
!
      if (msglvl .gt. 0) then
         call qpvout (
     $        logfil, 1_wip, (/rnorm/), ndigit,
     $   '_getv0: B-norm of initial / restarted starting vector')
      endif
      if (msglvl .gt. 2) then
         call qpvout (
     $        logfil, n, resid, ndigit,
     $   '_getv0: initial / restarted starting vector')
      endif
      ido = 99
! 
      call qpsecond (t1)
      tgetv0 = tgetv0 + (t1 - t0)
! 
 9000 continue
      return
      end subroutine qpgetv0 
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: qpsaitr
!
!\Description: 
!  Reverse communication interface for applying NP additional steps to 
!  a K step symmetric Arnoldi factorization.
!
!  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
!
!          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
!
!  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
!
!          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
!
!  where OP and B are as in qpsaupd.  The B-norm of r_{k+p} is also
!  computed and returned.
!
!\Usage:
!  call qpsaitr
!     ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH, 
!       IPNTR, WORKD, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!                    This is for the restart phase to force the new
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y,
!                    IPNTR(3) is the pointer into WORK for B * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!          When the routine is used in the "shift-and-invert" mode, the
!          vector B * Q is already available and does not need to be
!          recomputed in forming OP * Q.
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of matrix B that defines the
!          semi-inner product for the operator OP.  See qpsaupd.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  K       Integer.  (INPUT)
!          Current order of H and the number of columns of V.
!
!  NP      Integer.  (INPUT)
!          Number of additional Arnoldi steps to take.
!
!  MODE    Integer.  (INPUT)
!          Signifies which form for "OP". If MODE=2 then
!          a reduction in the number of B matrix vector multiplies
!          is possible since the B-norm of OP*x is equivalent to
!          the inv(B)-norm of A*x.
!
!  RESID   real(wrp) array of length N.  (INPUT/OUTPUT)
!          On INPUT:  RESID contains the residual vector r_{k}.
!          On OUTPUT: RESID contains the residual vector r_{k+p}.
!
!  RNORM   real(wrp) scalar.  (INPUT/OUTPUT)
!          On INPUT the B-norm of r_{k}.
!          On OUTPUT the B-norm of the updated residual r_{k+p}.
!
!  V       real(wrp) N by K+NP array.  (INPUT/OUTPUT)
!          On INPUT:  V contains the Arnoldi vectors in the first K 
!          columns.
!          On OUTPUT: V contains the new NP Arnoldi vectors in the next
!          NP columns.  The first K columns are unchanged.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling 
!          program.
!
!  H       real(wrp) (K+NP) by 2 array.  (INPUT/OUTPUT)
!          H is used to store the generated symmetric tridiagonal matrix
!          with the subdiagonal in the first column starting at H(2,1)
!          and the main diagonal in the qpsecond column.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling 
!          program.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORK for 
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the 
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!          
!  WORKD   real(wrp) work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The calling program should not 
!          use WORKD as temporary workspace during the iteration !!!!!!
!          On INPUT, WORKD(1:N) = B*RESID where RESID is associated
!          with the K step Arnoldi factorization. Used to save some 
!          computation at the first step. 
!          On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
!          with the K+NP step Arnoldi factorization.
!
!  INFO    Integer.  (OUTPUT)
!          = 0: Normal exit.
!          > 0: Size of an invariant subspace of OP is found that is
!               less than K + NP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     qpgetv0  ARPACK routine to generate the initial vector.
!     wivout   ARPACK utility routine that prints integers.
!     qpvout   ARPACK utility routine that prints vectors.
!     qplamch  LAPACK routine that determines machine constants.
!     qplascl  LAPACK routine for careful scaling of a matrix.
!     qpgemv   Level 2 BLAS routine for matrix vector multiplication.
!     qpaxpy   Level 1 BLAS that computes a vector triad.
!     qpscal   Level 1 BLAS that scales a vector.
!     qpcopy   Level 1 BLAS that copies one vector to another .
!     qpdot    Level 1 BLAS that computes the scalar product of two vectors. 
!     qpnrm2   Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
! 
!\Revision history:
!     xx/xx/93: Version ' 2.4'
!
!\SCCS Information: @(#) 
! FILE: saitr.F   SID: 2.6   DATE OF SID: 8/28/96   RELEASE: 2
!
!\Remarks
!  The algorithm implemented is:
!  
!  restart = .false.
!  Given V_{k} = [v_{1}, ..., v_{k}], r_{k}; 
!  r_{k} contains the initial residual vector even for k = 0;
!  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already 
!  computed by the calling program.
!
!  betaj = rnorm ; p_{k+1} = B*r_{k} ;
!  For  j = k+1, ..., k+np  Do
!     1) if ( betaj < tol ) stop or restart depending on j.
!        if ( restart ) generate a new starting vector.
!     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];  
!        p_{j} = p_{j}/betaj
!     3) r_{j} = OP*v_{j} where OP is defined as in qpsaupd
!        For shift-invert mode p_{j} = B*v_{j} is already available.
!        wnorm = || OP*v_{j} ||
!     4) Compute the j-th step residual vector.
!        w_{j} =  V_{j}^T * B * OP * v_{j}
!        r_{j} =  OP*v_{j} - V_{j} * w_{j}
!        alphaj <- j-th component of w_{j}
!        rnorm = || r_{j} ||
!        betaj+1 = rnorm
!        If (rnorm > 0.717D+0*wnorm) accept step and go back to 1)
!     5) Re-orthogonalization step:
!        s = V_{j}'*B*r_{j}
!        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
!        alphaj = alphaj + s_{j};   
!     6) Iterative refinement step:
!        If (rnorm1 > 0.717D+0*rnorm) then
!           rnorm = rnorm1
!           accept step and go back to 1)
!        Else
!           rnorm = rnorm1
!           If this is the first time in step 6), goto 5)
!           Else r_{j} lies in the span of V_{j} numerically.
!              Set r_{j} = 0 and rnorm = 0; goto 1)
!        EndIf 
!  End Do
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpsaitr (
     $           ido, bmat, n, k, np, mode, resid, rnorm, 
     $           v, ldv, h, ldh, ipntr, workd, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer
     $   ido
      character(len=1)  
     $   bmat
      integer(wip)
     $   k, ldh, ldv, n, mode, np, ipntr(3),
     $   info
      real(wrp)
     $   rnorm, h(ldh,2), resid(n), v(ldv,k+np), workd(3*n)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real(wrp)
     $   one, zero
      parameter (
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      logical    
     $   first, orth1, orth2, rstart, step3, step4
      integer(wip)  
     $   ierr, ipj, irj, ivj, i, iter, itry, j, infol, jj
      integer    
     $   msglvl 
      real(wrp)
     $   rnorm1, wnorm, safmin, temp1
      save 
     $   orth1, orth2, rstart, step3, step4,
     $   ierr, ipj, irj, ivj, iter, itry, j, msglvl,
     $   rnorm1, safmin, wnorm
!
!     %-----------------------%
!     | Local Array Arguments | 
!     %-----------------------%
!
      real(wrp)
     $   xtemp(2)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external 
     $   qpaxpy, 
     $   qpcopy, 
     $   qpscal, 
     $   qpgemv, 
     $   qpgetv0, 
     $   qpvout, 
     $   qplascl, 
     $   wivout, 
     $   qpsecond
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      real(wrp)
     $   qpdot, 
     $   qpnrm2, 
     $   qplamch
      external 
     $   qpdot, 
     $   qpnrm2, 
     $   qplamch
!
!     %-----------------%
!     | Data statements |
!     %-----------------%
!
      data      
     $   first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (first) then
         first = .false.
!
!        %--------------------------------%
!        | safmin = safe minimum is such  |
!        | that 1/sfmin does not overflow |
!        %--------------------------------%
!
         safmin = qplamch ( 
     $            'safmin' )
      endif
!
      if (ido .eq. 0) then
! 
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call qpsecond (t0)
         msglvl = msaitr
! 
!        %------------------------------%
!        | Initial call to this routine |
!        %------------------------------%
!
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
! 
!        %--------------------------------%
!        | Pointer to the current step of |
!        | the factorization to build     |
!        %--------------------------------%
!
         j      = k + 1
! 
!        %------------------------------------------%
!        | Pointers used for reverse communication  |
!        | when using WORKD.                        |
!        %------------------------------------------%
!
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      endif
! 
!     %-------------------------------------------------%
!     | When in reverse communication mode one of:      |
!     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
!     | will be .true.                                  |
!     | STEP3: return from computing OP*v_{j}.          |
!     | STEP4: return from computing B-norm of OP*v_{j} |
!     | ORTH1: return from computing B-norm of r_{j+1}  |
!     | ORTH2: return from computing B-norm of          |
!     |        correction to the residual vector.       |
!     | RSTART: return from OP computations needed by   |
!     |         qpgetv0.                                 |
!     %-------------------------------------------------%
!
      if (step3)  goto 50
      if (step4)  goto 60
      if (orth1)  goto 70
      if (orth2)  goto 90
      if (rstart) goto 30
!
!     %------------------------------%
!     | Else this is the first step. |
!     %------------------------------%
! 
!     %--------------------------------------------------------------%
!     |                                                              |
!     |        A R N O L D I     I T E R A T I O N     L O O P       |
!     |                                                              |
!     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
!     %--------------------------------------------------------------%
!
 1000 continue
!
         if (msglvl .gt. 2) then
            call wivout (
     $           logfil, 1_wip, (/j/), ndigit, 
     $      '_saitr: generating Arnoldi vector no.' )
            call qpvout (
     $           logfil, 1_wip, (/rnorm/), ndigit, 
     $      '_saitr: B-norm of the current residual =' )
         endif
! 
!        %---------------------------------------------------------%
!        | Check for exact zero. Equivalent to determing whether a |
!        | j-step Arnoldi factorization is present.                |
!        %---------------------------------------------------------%
!
         if (rnorm .gt. zero) goto 40
!
!           %---------------------------------------------------%
!           | Invariant subspace found, generate a new starting |
!           | vector which is orthogonal to the current Arnoldi |
!           | basis and continue the iteration.                 |
!           %---------------------------------------------------%
!
            if (msglvl .gt. 0) then
               call wivout (
     $              logfil, 1_wip, (/j/), ndigit,
     $              '_saitr: ****** restart at step ******' )
            endif
! 
!           %---------------------------------------------%
!           | ITRY is the loop variable that controls the |
!           | maximum amount of times that a restart is   |
!           | attempted. NRSTRT is used by stat.h         |
!           %---------------------------------------------%
!
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
!
!           %--------------------------------------%
!           | If in reverse communication mode and |
!           | RSTART = .true. flow returns here.   |
!           %--------------------------------------%
!
            call qpgetv0 (
     $           ido, bmat, itry, .false., n, j, v, ldv, 
     $           resid, rnorm, ipntr, workd, ierr )
            if (ido .ne. 99) goto 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) goto 20
!
!              %------------------------------------------------%
!              | Give up after several restart attempts.        |
!              | Set INFO to the size of the invariant subspace |
!              | which spans OP and exit.                       |
!              %------------------------------------------------%
!
               info = j - 1
               call qpsecond (t1)
               tsaitr = tsaitr + (t1 - t0)
               ido = 99
               goto 9000
            endif
! 
   40    continue
!
!        %---------------------------------------------------------%
!        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
!        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
!        | when reciprocating a small RNORM, test against lower    |
!        | machine bound.                                          |
!        %---------------------------------------------------------%
!
         call qpcopy (
     $        n, resid, 1_wip, v(1,j), 1_wip )
         if (rnorm .ge. safmin) then
             temp1 = one / rnorm
             call qpscal (
     $            n, temp1, v(1,j), 1_wip )
             call qpscal (
     $            n, temp1, workd(ipj), 1_wip )
         else
!
!            %-----------------------------------------%
!            | To scale both v_{j} and p_{j} carefully |
!            | use LAPACK routine SLASCL               |
!            %-----------------------------------------%
!
             call qplascl (
     $            'General', i, i, rnorm, one, n, 1_wip, 
     $            v(1,j), n, infol )
             call qplascl (
     $            'General', i, i, rnorm, one, n, 1_wip, 
     $            workd(ipj), n, infol )
         endif
! 
!        %------------------------------------------------------%
!        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
!        | Note that this is not quite yet r_{j}. See STEP 4    |
!        %------------------------------------------------------%
!
         step3 = .true.
         nopx  = nopx + 1
         call qpsecond (t2)
         call qpcopy (
     $        n, v(1,j), 1_wip, workd(ivj), 1_wip )
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
! 
!        %-----------------------------------%
!        | Exit in order to compute OP*v_{j} |
!        %-----------------------------------%
! 
         goto 9000
   50    continue
! 
!        %-----------------------------------%
!        | Back from reverse communication;  |
!        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
!        %-----------------------------------%
!
         call qpsecond (t3)
         tmvopx = tmvopx + (t3 - t2)
! 
         step3 = .false.
!
!        %------------------------------------------%
!        | Put another copy of OP*v_{j} into RESID. |
!        %------------------------------------------%
!
         call qpcopy (
     $        n, workd(irj), 1_wip, resid, 1_wip )
! 
!        %-------------------------------------------%
!        | STEP 4:  Finish extending the symmetric   |
!        |          Arnoldi to length j. If MODE = 2 |
!        |          then B*OP = B*inv(B)*A = A and   |
!        |          we don't need to compute B*OP.   |
!        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
!        | assumed to have A*v_{j}.                  |
!        %-------------------------------------------%
!
         if (mode .eq. 2) goto 65
         call qpsecond (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
! 
!           %-------------------------------------%
!           | Exit in order to compute B*OP*v_{j} |
!           %-------------------------------------%
! 
            goto 9000
         else if (bmat .eq. 'I') then
              call qpcopy ( 
     $             n, resid, 1_wip, workd(ipj), 1_wip )
         endif
   60    continue
! 
!        %-----------------------------------%
!        | Back from reverse communication;  |
!        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
!        %-----------------------------------%
!
         if (bmat .eq. 'G') then
            call qpsecond (t3)
            tmvbx = tmvbx + (t3 - t2)
         endif 
!
         step4 = .false.
!
!        %-------------------------------------%
!        | The following is needed for STEP 5. |
!        | Compute the B-norm of OP*v_{j}.     |
!        %-------------------------------------%
!
   65    continue
         if (mode .eq. 2) then
!
!           %----------------------------------%
!           | Note that the B-norm of OP*v_{j} |
!           | is the inv(B)-norm of A*v_{j}.   |
!           %----------------------------------%
!
            wnorm = qpdot (
     $              n, resid, 1_wip, workd(ivj), 1_wip )
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'G') then         
            wnorm = qpdot (
     $              n, resid, 1_wip, workd(ipj), 1_wip )
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = qpnrm2 (
     $              n, resid, 1_wip )
         endif
!
!        %-----------------------------------------%
!        | Compute the j-th residual corresponding |
!        | to the j step factorization.            |
!        | Use Classical Gram Schmidt and compute: |
!        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
!        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
!        %-----------------------------------------%
!
!
!        %------------------------------------------%
!        | Compute the j Fourier coefficients w_{j} |
!        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
!        %------------------------------------------%
!
         if (mode .ne. 2 ) then
            call qpgemv (
     $           'T', n, j, one, v, ldv, workd(ipj), 1_wip, zero, 
     $            workd(irj), 1_wip )
         else if (mode .eq. 2) then
            call qpgemv (
     $           'T', n, j, one, v, ldv, workd(ivj), 1_wip, zero, 
     $           workd(irj), 1_wip )
         endif
!
!        %--------------------------------------%
!        | Orthgonalize r_{j} against V_{j}.    |
!        | RESID contains OP*v_{j}. See STEP 3. | 
!        %--------------------------------------%
!
         call qpgemv ( 
     $        'N', n, j, -one, v, ldv, workd(irj), 1_wip, one, 
     $        resid, 1_wip )
!
!        %--------------------------------------%
!        | Extend H to have j rows and columns. |
!        %--------------------------------------%
!
         h(j,2) = workd(irj + j - 1)
         if (j .eq. 1 .or.  rstart) then
            h(j,1) = zero
         else
            h(j,1) = rnorm
         endif
         call qpsecond (t4)
! 
         orth1 = .true.
         iter  = 0
! 
         call qpsecond (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call qpcopy (
     $           n, resid, 1_wip, workd(irj), 1_wip )
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
! 
!           %----------------------------------%
!           | Exit in order to compute B*r_{j} |
!           %----------------------------------%
! 
            goto 9000
         else if (bmat .eq. 'I') then
            call qpcopy (
     $           n, resid, 1_wip, workd(ipj), 1_wip )
         endif
   70    continue
! 
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH1 = .true. |
!        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
!        %---------------------------------------------------%
!
         if (bmat .eq. 'G') then
            call qpsecond (t3)
            tmvbx = tmvbx + (t3 - t2)
         endif
! 
         orth1 = .false.
!
!        %------------------------------%
!        | Compute the B-norm of r_{j}. |
!        %------------------------------%
!
         if (bmat .eq. 'G') then         
            rnorm = qpdot (
     $              n, resid, 1_wip, workd(ipj), 1_wip )
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = qpnrm2 (
     $              n, resid, 1_wip )
         endif
!
!        %-----------------------------------------------------------%
!        | STEP 5: Re-orthogonalization / Iterative refinement phase |
!        | Maximum NITER_ITREF tries.                                |
!        |                                                           |
!        |          s      = V_{j}^T * B * r_{j}                     |
!        |          r_{j}  = r_{j} - V_{j}*s                         |
!        |          alphaj = alphaj + s_{j}                          |
!        |                                                           |
!        | The stopping criteria used for iterative refinement is    |
!        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
!        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
!        | Determine if we need to correct the residual. The goal is |
!        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
!        %-----------------------------------------------------------%
!
         if (rnorm .gt. 0.717_wrp*wnorm) goto 100
         nrorth = nrorth + 1
! 
!        %---------------------------------------------------%
!        | Enter the Iterative refinement phase. If further  |
!        | refinement is necessary, loop back here. The loop |
!        | variable is ITER. Perform a step of Classical     |
!        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
!        %---------------------------------------------------%
!
   80    continue
!
         if (msglvl .gt. 2) then
            xtemp(1) = wnorm
            xtemp(2) = rnorm
            call qpvout (
     $           logfil, 2_wip, xtemp, ndigit, 
     $      '_saitr: re-orthonalization ; wnorm and rnorm are' )
         endif
!
!        %----------------------------------------------------%
!        | Compute V_{j}^T * B * r_{j}.                       |
!        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
!        %----------------------------------------------------%
!
         call qpgemv (
     $        'T', n, j, one, v, ldv, workd(ipj), 1_wip, 
     $        zero, workd(irj), 1_wip )
!
!        %----------------------------------------------%
!        | Compute the correction to the residual:      |
!        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
!        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
!        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
!        | H(j,j) is updated.                           |
!        %----------------------------------------------%
!
         call qpgemv (
     $        'N', n, j, -one, v, ldv, workd(irj), 1_wip, 
     $        one, resid, 1_wip )
!
         if (j .eq. 1 .or. rstart) h(j,1) = zero
         h(j,2) = h(j,2) + workd(irj + j - 1)
! 
         orth2 = .true.
         call qpsecond (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call qpcopy (
     $           n, resid, 1_wip, workd(irj), 1_wip )
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
! 
!           %-----------------------------------%
!           | Exit in order to compute B*r_{j}. |
!           | r_{j} is the corrected residual.  |
!           %-----------------------------------%
! 
            goto 9000
         else if (bmat .eq. 'I') then
            call qpcopy (
     $           n, resid, 1_wip, workd(ipj), 1_wip )
         endif
   90    continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH2 = .true. |
!        %---------------------------------------------------%
!
         if (bmat .eq. 'G') then
            call qpsecond (t3)
            tmvbx = tmvbx + (t3 - t2)
         endif
!
!        %-----------------------------------------------------%
!        | Compute the B-norm of the corrected residual r_{j}. |
!        %-----------------------------------------------------%
! 
         if (bmat .eq. 'G') then         
             rnorm1 = qpdot (
     $                n, resid, 1_wip, workd(ipj), 1_wip )
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = qpnrm2 ( 
     $                n, resid, 1_wip )
         endif
!
         if (msglvl .gt. 0 .and. iter .gt. 0) then
            call wivout (
     $           logfil, 1_wip, (/j/), ndigit,
     $           '_saitr: Iterative refinement for Arnoldi residual' )
            if (msglvl .gt. 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call qpvout (
     $               logfil, 2_wip, xtemp, ndigit,
     $          '_saitr: iterative refinement ; rnorm and rnorm1 are' )
            endif
         endif
! 
!        %-----------------------------------------%
!        | Determine if we need to perform another |
!        | step of re-orthogonalization.           |
!        %-----------------------------------------%
!
         if (rnorm1 .gt. 0.717_wrp*rnorm) then
!
!           %--------------------------------%
!           | No need for further refinement |
!           %--------------------------------%
!
            rnorm = rnorm1
! 
         else
!
!           %-------------------------------------------%
!           | Another step of iterative refinement step |
!           | is required. NITREF is used by stat.h     |
!           %-------------------------------------------%
!
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) goto 80
!
!           %-------------------------------------------------%
!           | Otherwise RESID is numerically in the span of V |
!           %-------------------------------------------------%
!
            do jj = 1, n
               resid(jj) = zero
            enddo
            rnorm = zero
         endif
! 
!        %----------------------------------------------%
!        | Branch here directly if iterative refinement |
!        | wasn't necessary or after at most NITER_REF  |
!        | steps of iterative refinement.               |
!        %----------------------------------------------%
!
  100    continue
! 
         rstart = .false.
         orth2  = .false.
! 
         call qpsecond (t5)
         titref = titref + (t5 - t4)
! 
!        %----------------------------------------------------------%
!        | Make sure the last off-diagonal element is non negative  |
!        | If not perform a similarity transformation on H(1:j,1:j) |
!        | and scale v(:,j) by -1.                                  |
!        %----------------------------------------------------------%
!
         if (h(j,1) .lt. zero) then
            h(j,1) = -h(j,1)
            if ( j .lt. k+np) then 
               call qpscal ( 
     $              n, -one, v(1,j+1), 1_wip )
            else
               call qpscal (
     $              n, -one, resid, 1_wip )
            endif
         endif
! 
!        %------------------------------------%
!        | STEP 6: Update  j = j+1;  Continue |
!        %------------------------------------%
!
         j = j + 1
         if (j .gt. k+np) then
            call qpsecond (t1)
            tsaitr = tsaitr + (t1 - t0)
            ido = 99
!
            if (msglvl .gt. 1) then
               call qpvout (
     $              logfil, k+np, h(1,2), ndigit, 
     $         '_saitr: main diagonal of matrix H of step K+NP.' )
               if (k+np .gt. 1) then
               call qpvout (
     $              logfil, k+np-1, h(2,1), ndigit, 
     $         '_saitr: sub diagonal of matrix H of step K+NP.')
               endif
            endif
!
            goto 9000
         endif
!
!        %--------------------------------------------------------%
!        | Loop back to extend the factorization by another step. |
!        %--------------------------------------------------------%
!
      goto 1000
! 
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
 9000 continue
      return
      end subroutine qpsaitr
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: qpsapps
!
!\Description:
!  Given the Arnoldi factorization
!
!     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!
!  apply NP shifts implicitly resulting in
!
!     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!
!  where Q is an orthogonal matrix of order KEV+NP. Q is the product of 
!  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi 
!  factorization becomes:
!
!     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
!
!\Usage:
!  call qpsapps
!     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Problem size, i.e. dimension of matrix A.
!
!  KEV     Integer.  (INPUT)
!          INPUT: KEV+NP is the size of the input matrix H.
!          OUTPUT: KEV is the size of the updated matrix HNEW.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be applied.
!
!  SHIFT   real(wrp) array of length NP.  (INPUT)
!          The shifts to be applied.
!
!  V       real(wrp) N by (KEV+NP) array.  (INPUT/OUTPUT)
!          INPUT: V contains the current KEV+NP Arnoldi vectors.
!          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
!          are in the first KEV columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       real(wrp) (KEV+NP) by 2 array.  (INPUT/OUTPUT)
!          INPUT: H contains the symmetric tridiagonal matrix of the
!          Arnoldi factorization with the subdiagonal in the 1st column
!          starting at H(2,1) and the main diagonal in the 2nd column.
!          OUTPUT: H contains the updated tridiagonal matrix in the 
!          KEV leading submatrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RESID   real(wrp) array of length (N).  (INPUT/OUTPUT)
!          INPUT: RESID contains the the residual vector r_{k+p}.
!          OUTPUT: RESID is the updated residual vector rnew_{k}.
!
!  Q       real(wrp) KEV+NP by KEV+NP work array.  (WORKSPACE)
!          Work array used to accumulate the rotations during the bulge
!          chase sweep.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKD   real(wrp) work array of length 2*N.  (WORKSPACE)
!          Distributed array used in the application of the accumulated
!          orthogonal matrix Q.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     wivout   ARPACK utility routine that prints integers. 
!     qpsecond  ARPACK utility routine for timing.
!     qpvout   ARPACK utility routine that prints vectors.
!     qplamch  LAPACK routine that determines machine constants.
!     qplartg  LAPACK Givens rotation construction routine.
!     qplacpy  LAPACK matrix copy routine.
!     qplaset  LAPACK matrix initialization routine.
!     qpgemv   Level 2 BLAS routine for matrix vector multiplication.
!     qpaxpy   Level 1 BLAS that computes a vector triad.
!     qpcopy   Level 1 BLAS that copies one vector to another.
!     qpscal   Level 1 BLAS that scales a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     12/16/93: Version ' 2.1'
!
!\SCCS Information: @(#) 
! FILE: sapps.F   SID: 2.5   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!  1. In this version, each shift is applied to all the subblocks of
!     the tridiagonal matrix H and not just to the submatrix that it 
!     comes from. This routine assumes that the subdiagonal elements 
!     of H that are stored in h(1:kev+np,1) are nonegative upon input
!     and enforce this condition upon output. This version incorporates
!     deflation. See code for documentation.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpsapps (
     $           n, kev, np, shift, v, ldv, h, ldh, resid, 
     $           q, ldq, workd )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer(wip) 
     $   kev, ldh, ldq, ldv, n, np
      real(wrp)
     $   h(ldh,2), q(ldq,kev+np), resid(n), shift(np), 
     $   v(ldv,kev+np), workd(2*n)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real(wrp)
     $   one, zero
      parameter (
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer(wip) 
     $   i, iend, istart, itop, j, jj, kplusp
      integer    
     $   msglvl

      logical    
     $   first
      real(wrp)
     $   a1, a2, a3, a4, big, c, epsmch, f, g, r, s
      save 
     $   epsmch, first
!
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external
     $   qpaxpy, 
     $   qpcopy, 
     $   qpscal, 
     $   qplacpy, 
     $   qplartg, 
     $   qplaset, 
     $   qpvout, 
     $   wivout, 
     $   qpsecond, 
     $   qpgemv
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      real(wrp)
     $   qplamch
      external   
     $   qplamch
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  
     $   abs
!
!     %----------------%
!     | Data statments |
!     %----------------%
!
      data
     $   first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (first) then
         epsmch = qplamch ( 
     $            'Epsilon-Machine' )
         first = .false.
      endif
      itop = 1
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call qpsecond (t0)
      msglvl = msapps
! 
      kplusp = kev + np 
! 
!     %----------------------------------------------%
!     | Initialize Q to the identity matrix of order |
!     | kplusp used to accumulate the rotations.     |
!     %----------------------------------------------%
!
      call qplaset (
     $     'All', kplusp, kplusp, zero, one, q, ldq )
!
!     %----------------------------------------------%
!     | Quick return if there are no shifts to apply |
!     %----------------------------------------------%
!
      if (np .eq. 0) goto 9000
! 
!     %----------------------------------------------------------%
!     | Apply the np shifts implicitly. Apply each shift to the  |
!     | whole matrix and not just to the submatrix from which it |
!     | comes.                                                   |
!     %----------------------------------------------------------%
!
      do 90 jj = 1, np
! 
         istart = itop
!
!        %----------------------------------------------------------%
!        | Check for splitting and deflation. Currently we consider |
!        | an off-diagonal element h(i+1,1) negligible if           |
!        |         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )   |
!        | for i=1:KEV+NP-1.                                        |
!        | If above condition tests true then we set h(i+1,1) = 0.  |
!        | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
!        %----------------------------------------------------------%
!
   20    continue
!
!        %------------------------------------------------%
!        | The following loop exits early if we encounter |
!        | a negligible off diagonal element.             |
!        %------------------------------------------------%
!
         do i = istart, kplusp-1
            big   = abs(h(i,2)) + abs(h(i+1,2))
            if (h(i+1,1) .le. epsmch*big) then
               if (msglvl .gt. 0) then
                  call wivout (
     $                 logfil, 1_wip, (/i/), ndigit, 
     $            '_sapps: deflation at row/column no.' )
                  call wivout (
     $                 logfil, 1_wip, (/jj/), ndigit, 
     $            '_sapps: occured before shift number.')
                  call qpvout (
     $                 logfil, 1_wip, h(i+1,1), ndigit, 
     $            '_sapps: the corresponding off diagonal element')
               endif
               h(i+1,1) = zero
               iend = i
               goto 40
            endif
         enddo 
         iend = kplusp
   40    continue
!
         if (istart .lt. iend) then
! 
!           %--------------------------------------------------------%
!           | Construct the plane rotation G'(istart,istart+1,theta) |
!           | that attempts to drive h(istart+1,1) to zero.          |
!           %--------------------------------------------------------%
!
            f = h(istart,2) - shift(jj)
            g = h(istart+1,1)
            call qplartg (
     $           f, g, c, s, r )
! 
!           %-------------------------------------------------------%
!           | Apply rotation to the left and right of H;            |
!           | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
!           | This will create a "bulge".                           |
!           %-------------------------------------------------------%
!
            a1 = c*h(istart,2)   + s*h(istart+1,1)
            a2 = c*h(istart+1,1) + s*h(istart+1,2)
            a4 = c*h(istart+1,2) - s*h(istart+1,1)
            a3 = c*h(istart+1,1) - s*h(istart,2) 
            h(istart,2)   = c*a1 + s*a2
            h(istart+1,2) = c*a4 - s*a3
            h(istart+1,1) = c*a3 + s*a4
! 
!           %----------------------------------------------------%
!           | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!           %----------------------------------------------------%
!
            do j = 1, min(istart+jj,kplusp)
               a1            =   c*q(j,istart) + s*q(j,istart+1)
               q(j,istart+1) = - s*q(j,istart) + c*q(j,istart+1)
               q(j,istart)   = a1
            enddo 
               
!
!
!           %----------------------------------------------%
!           | The following loop chases the bulge created. |
!           | Note that the previous rotation may also be  |
!           | done within the following loop. But it is    |
!           | kept separate to make the distinction among  |
!           | the bulge chasing sweeps and the first plane |
!           | rotation designed to drive h(istart+1,1) to  |
!           | zero.                                        |
!           %----------------------------------------------%
!
            do i = istart+1, iend-1
! 
!              %----------------------------------------------%
!              | Construct the plane rotation G'(i,i+1,theta) |
!              | that zeros the i-th bulge that was created   |
!              | by G(i-1,i,theta). g represents the bulge.   |
!              %----------------------------------------------%
!
               f = h(i,1)
               g = s*h(i+1,1)
!
!              %----------------------------------%
!              | Final update with G(i-1,i,theta) |
!              %----------------------------------%
!
               h(i+1,1) = c*h(i+1,1)
               call qplartg (
     $              f, g, c, s, r )
!
!              %-------------------------------------------%
!              | The following ensures that h(1:iend-1,1), |
!              | the first iend-2 off diagonal of elements |
!              | H, remain non negative.                   |
!              %-------------------------------------------%
!
               if (r .lt. zero) then
                  r = -r
                  c = -c
                  s = -s
               endif
! 
!              %--------------------------------------------%
!              | Apply rotation to the left and right of H; |
!              | H <- G * H * G',  where G = G(i,i+1,theta) |
!              %--------------------------------------------%
!
               h(i,1) = r
! 
               a1 = c*h(i,2)   + s*h(i+1,1)
               a2 = c*h(i+1,1) + s*h(i+1,2)
               a3 = c*h(i+1,1) - s*h(i,2)
               a4 = c*h(i+1,2) - s*h(i+1,1)
! 
               h(i,2)   = c*a1 + s*a2
               h(i+1,2) = c*a4 - s*a3
               h(i+1,1) = c*a3 + s*a4
! 
!              %----------------------------------------------------%
!              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!              %----------------------------------------------------%
!
               do j = 1, min( j+jj, kplusp )
                  a1       =   c*q(j,i) + s*q(j,i+1)
                  q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                  q(j,i)   = a1
               enddo 
!
            enddo
!
         endif
!
!        %--------------------------%
!        | Update the block pointer |
!        %--------------------------%
!
         istart = iend + 1
!
!        %------------------------------------------%
!        | Make sure that h(iend,1) is non-negative |
!        | If not then set h(iend,1) <-- -h(iend,1) |
!        | and negate the last column of Q.         |
!        | We have effectively carried out a        |
!        | similarity on transformation H           |
!        %------------------------------------------%
!
         if (h(iend,1) .lt. zero) then
             h(iend,1) = -h(iend,1)
             call qpscal ( 
     $            kplusp, -one, q(1,iend), 1_wip )
         endif
!
!        %--------------------------------------------------------%
!        | Apply the same shift to the next block if there is any |
!        %--------------------------------------------------------%
!
         if (iend .lt. kplusp) goto 20
!
!        %-----------------------------------------------------%
!        | Check if we can increase the the start of the block |
!        %-----------------------------------------------------%
!
         do i = itop, kplusp-1
            if (h(i+1,1) .gt. zero) goto 90
            itop  = itop + 1
         enddo 
!
!        %-----------------------------------%
!        | Finished applying the jj-th shift |
!        %-----------------------------------%
!
   90 continue
!
!     %------------------------------------------%
!     | All shifts have been applied. Check for  |
!     | more possible deflation that might occur |
!     | after the last shift is applied.         |                               
!     %------------------------------------------%
!
      do i = itop, kplusp-1
         big   = abs(h(i,2)) + abs(h(i+1,2))
         if (h(i+1,1) .le. epsmch*big) then
            if (msglvl .gt. 0) then
               call wivout (
     $              logfil, 1_wip, (/i/), ndigit, 
     $         '_sapps: deflation at row/column no.' )
               call qpvout (
     $              logfil, 1_wip, h(i+1,1), ndigit, 
     $         '_sapps: the corresponding off diagonal element')
            endif
            h(i+1,1) = zero
         endif
      enddo 
!
!     %-------------------------------------------------%
!     | Compute the (kev+1)-st column of (V*Q) and      |
!     | temporarily store the result in WORKD(N+1:2*N). |
!     | This is not necessary if h(kev+1,1) = 0.         |
!     %-------------------------------------------------%
!
      if ( h(kev+1,1) .gt. zero ) 
     $   call qpgemv (
     $        'N', n, kplusp, one, v, ldv,
     $        q(1,kev+1), 1_wip, zero, workd(n+1), 1_wip )
! 
!     %-------------------------------------------------------%
!     | Compute column 1 to kev of (V*Q) in backward order    |
!     | taking advantage that Q is an upper triangular matrix |    
!     | with lower bandwidth np.                              |
!     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
!     %-------------------------------------------------------%
!
      do i = 1, kev
         call qpgemv (
     $        'N', n, kplusp-i+1, one, v, ldv,
     $        q(1,kev-i+1), 1_wip, zero, workd, 1_wip )
         call qpcopy (
     $        n, workd, 1_wip, v(1,kplusp-i+1), 1_wip )
      enddo 
!
!     %-------------------------------------------------%
!     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
!     %-------------------------------------------------%
!
      call qplacpy (
     $     'All', n, kev, v(1,np+1), ldv, v, ldv )
! 
!     %--------------------------------------------%
!     | Copy the (kev+1)-st column of (V*Q) in the |
!     | appropriate place if h(kev+1,1) .ne. zero. |
!     %--------------------------------------------%
!
      if ( h(kev+1,1) .gt. zero ) 
     $     call qpcopy (
     $          n, workd(n+1), 1_wip, v(1,kev+1), 1_wip )
! 
!     %-------------------------------------%
!     | Update the residual vector:         |
!     |    r <- sigmak*r + betak*v(:,kev+1) |
!     | where                               |
!     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
!     |    betak = e_{kev+1}'*H*e_{kev}     |
!     %-------------------------------------%
!
      call qpscal (
     $     n, q(kplusp,kev), resid, 1_wip )
      if (h(kev+1,1) .gt. zero) 
     $   call qpaxpy (
     $        n, h(kev+1,1), v(1,kev+1), 1_wip, resid, 1_wip )
!
      if (msglvl .gt. 1) then
         call qpvout (
     $        logfil, 1_wip, q(kplusp,kev), ndigit, 
     $   '_sapps: sigmak of the updated residual vector' )
         call qpvout (
     $        logfil, 1_wip, h(kev+1,1), ndigit, 
     $   '_sapps: betak of the updated residual vector' )
         call qpvout (
     $        logfil, kev, h(1,2), ndigit, 
     $   '_sapps: updated main diagonal of H for next iteration' )
         if (kev .gt. 1) then
            call qpvout (
     $           logfil, kev-1, h(2,1), ndigit, 
     $      '_sapps: updated sub diagonal of H for next iteration' )
         endif
      endif
!
      call qpsecond (t1)
      tsapps = tsapps + (t1 - t0)
! 
 9000 continue 
      return
      end subroutine qpsapps
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: qpsconv
!
!\Description: 
!  Convergence testing for the symmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call qpsconv
!     ( N, RITZ, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZ    real(wrp) array of length N.  (INPUT)
!          The Ritz values to be checked for convergence.
!
!  BOUNDS  real(wrp) array of length N.  (INPUT)
!          Ritz estimates associated with the Ritz values in RITZ.
!
!  TOL     real(wrp) scalar.  (INPUT)
!          Desired relative accuracy for a Ritz value to be considered
!          "converged".
!
!  NCONV   Integer scalar.  (OUTPUT)
!          Number of "converged" Ritz values.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines called:
!     qpsecond  ARPACK utility routine for timing.
!     qplamch  LAPACK routine that determines machine constants. 
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Dept. of Computational &     Houston, Texas 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\SCCS Information: @(#) 
! FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.4, this routine no longer uses the
!        Parlett strategy using the gap conditions. 
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpsconv (
     $           n, ritz, bounds, tol, nconv )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer(wip) 
     $   n, nconv
      real(wrp)
     $   tol, ritz(n), bounds(n)

!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer(wip)  
     $   i
      real(wrp)
     $   temp, eps23
!
!     %-------------------%
!     | External routines |
!     %-------------------%
!
      real(wrp)
     $   qplamch
      external   
     $   qplamch

!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    
     $   abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      call qpsecond (t0)
!
      eps23 = qplamch ( 
     $        'Epsilon-Machine' ) 
      eps23 = eps23**( 2.0_wrp / 3.0_wrp )
!
      nconv  = 0
      do i = 1, n
!
!        %-----------------------------------------------------%
!        | The i-th Ritz value is considered "converged"       |
!        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
!        %-----------------------------------------------------%
!
         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) .le. tol*temp ) nconv = nconv + 1
      enddo 
! 
      call qpsecond (t1)
      tsconv = tsconv + (t1 - t0)
! 
      return
      end subroutine qpsconv
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: qpseigt
!
!\Description: 
!  Compute the eigenvalues of the current symmetric tridiagonal matrix
!  and the corresponding error bounds given the current residual norm.
!
!\Usage:
!  call qpseigt
!     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
!
!\Arguments
!  RNORM   real(wrp) scalar.  (INPUT)
!          RNORM contains the residual norm corresponding to the current
!          symmetric tridiagonal matrix H.
!
!  N       Integer.  (INPUT)
!          Size of the symmetric tridiagonal matrix H.
!
!  H       real(wrp) N by 2 array.  (INPUT)
!          H contains the symmetric tridiagonal matrix with the 
!          subdiagonal in the first column starting at H(2,1) and the 
!          main diagonal in qpsecond column.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling 
!          program.
!
!  EIG     real(wrp) array of length N.  (OUTPUT)
!          On output, EIG contains the N eigenvalues of H possibly 
!          unsorted.  The BOUNDS arrays are returned in the
!          same sorted order as EIG.
!
!  BOUNDS  real(wrp) array of length N.  (OUTPUT)
!          On output, BOUNDS contains the error estimates corresponding
!          to the eigenvalues EIG.  This is equal to RNORM times the
!          last components of the eigenvectors corresponding to the
!          eigenvalues in EIG.
!
!  WORKL   real(wrp) work array of length 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from qpstqrb.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     qpstqrb  ARPACK routine that computes the eigenvalues and the
!             last components of the eigenvectors of a symmetric
!             and tridiagonal matrix.
!     qpsecond  ARPACK utility routine for timing.
!     qpvout   ARPACK utility routine that prints vectors.
!     qpcopy   Level 1 BLAS that copies one vector to another.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University 
!     Dept. of Computational &     Houston, Texas 
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!
!\SCCS Information: @(#) 
! FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpseigt ( 
     $           rnorm, n, h, ldh, eig, bounds, workl, ierr )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer(wip)  
     $   ldh, n
      real(wrp)
     $   rnorm, eig(n), bounds(n), h(ldh,2), workl(3*n)
      integer(wip)  
     $   ierr
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real(wrp)
     $   zero
      parameter (
     $   zero = 0.0_wrp )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer(wip)  
     $   k
      integer    
     $   msglvl


!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external 
     $   qpcopy, 
     $   qpstqrb, 
     $   qpvout, 
     $   qpsecond
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------% 
!
      call qpsecond (t0)
      msglvl = mseigt
!
      if (msglvl .gt. 0) then
         call qpvout (
     $        logfil, n, h(1,2), ndigit,
     $        '_seigt: main diagonal of matrix H' )
         if (n .gt. 1) then
         call qpvout (
     $        logfil, n-1, h(2,1), ndigit,
     $        '_seigt: sub diagonal of matrix H' )
         endif
      endif
!
      call qpcopy (
     $     n, h(1,2), 1_wip, eig, 1_wip )
      call qpcopy (
     $     n-1, h(2,1), 1_wip, workl, 1_wip )
      call qpstqrb (
     $     n, eig, workl, bounds, workl(n+1), ierr )
      if (ierr .ne. 0) goto 9000
      if (msglvl .gt. 1) then
         call qpvout (
     $        logfil, n, bounds, ndigit,
     $   '_seigt: last row of the eigenvector matrix for H' )
      endif
!
!     %-----------------------------------------------%
!     | Finally determine the error bounds associated |
!     | with the n Ritz values of H.                  |
!     %-----------------------------------------------%
!
      do k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
      enddo 
! 
      call qpsecond (t1)
      tseigt = tseigt + (t1 - t0)
!
 9000 continue
      return
      end subroutine qpseigt
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: qpsgets
!
!\Description: 
!  Given the eigenvalues of the symmetric tridiagonal matrix H,
!  computes the NP shifts AMU that are zeros of the polynomial of 
!  degree NP which filters out components of the unwanted eigenvectors 
!  corresponding to the AMU's based on some given criteria.
!
!  NOTE: This is called even in the case of user specified shifts in 
!  order to sort the eigenvalues, and error bounds of H for later use.
!
!\Usage:
!  call qpsgets
!     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
!
!\Arguments
!  ISHIFT  Integer.  (INPUT)
!          Method for selecting the implicit shifts at each iteration.
!          ISHIFT = 0: user specified shifts
!          ISHIFT = 1: exact shift with respect to the matrix H.
!
!  WHICH   Character*2.  (INPUT)
!          Shift selection criteria.
!          'LM' -> KEV eigenvalues of largest magnitude are retained.
!          'SM' -> KEV eigenvalues of smallest magnitude are retained.
!          'LA' -> KEV eigenvalues of largest value are retained.
!          'SA' -> KEV eigenvalues of smallest value are retained.
!          'BE' -> KEV eigenvalues, half from each end of the spectrum.
!                  If KEV is odd, compute one more from the high end.
!
!  KEV      Integer.  (INPUT)
!          KEV+NP is the size of the matrix H.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be computed.
!
!  RITZ    real(wrp) array of length KEV+NP.  (INPUT/OUTPUT)
!          On INPUT, RITZ contains the eigenvalues of H.
!          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues 
!          are in the first NP locations and the wanted part is in 
!          the last KEV locations.  When exact shifts are selected, the
!          unwanted part corresponds to the shifts to be applied.
!
!  BOUNDS  real(wrp) array of length KEV+NP.  (INPUT/OUTPUT)
!          Error bounds corresponding to the ordering in RITZ.
!
!  SHIFTS  real(wrp) array of length NP.  (INPUT/OUTPUT)
!          On INPUT:  contains the user specified shifts if ISHIFT = 0.
!          On OUTPUT: contains the shifts sorted into decreasing order 
!          of magnitude with respect to the Ritz estimates contained in
!          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     qpsortr  ARPACK utility sorting routine.
!     wivout   ARPACK utility routine that prints integers.
!     qpsecond  ARPACK utility routine for timing.
!     qpvout   ARPACK utility routine that prints vectors.
!     qpcopy   Level 1 BLAS that copies one vector to another.
!     qpswap   Level 1 BLAS that swaps the contents of two vectors.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\Revision history:
!     xx/xx/93: Version ' 2.1'
!
!\SCCS Information: @(#) 
! FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpsgets ( 
     $           ishift, which, kev, np, ritz, bounds, shifts )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character(len=2) 
     $   which
      integer(wip)  
     $   ishift, kev, np
      real(wrp)
     $   bounds(kev+np), ritz(kev+np), shifts(np)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
!     include   'debug.h'
      integer  logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
      common /debug/ 
     $         logfil, ndigit, mgetv0,
     $         msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd,
     $         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd,
     $         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd
!
!     include   'stat.h'
      real       t0, t1, t2, t3, t4, t5
      save       t0, t1, t2, t3, t4, t5
      integer    nopx, nbx, nrorth, nitref, nrstrt
      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     $           nopx, nbx, nrorth, nitref, nrstrt,
     $           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
     $           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
     $           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
     $           tmvopx, tmvbx, tgetv0, titref, trvec
!
!     %------------%
!     | Parameters |
!     %------------%
!
      real(wrp)
     $   one, zero
      parameter (
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer(wip) 
     $   kevd2
      integer
     $   msglvl

!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external 
     $   qpswap, 
     $   qpcopy, 
     $   qpsortr, 
     $   qpsecond
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic 
     $   max, min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
! 
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call qpsecond (t0)
      msglvl = msgets
! 
      if (which .eq. 'BE') then
!
!        %-----------------------------------------------------%
!        | Both ends of the spectrum are requested.            |
!        | Sort the eigenvalues into algebraically increasing  |
!        | order first then swap high end of the spectrum next |
!        | to low end in appropriate locations.                |
!        | NOTE: when np < floor(kev/2) be careful not to swap |
!        | overlapping locations.                              |
!        %-----------------------------------------------------%
!
         call qpsortr (
     $        'LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2 
         if ( kev .gt. 1 ) then
            call qpswap ( 
     $           min(kevd2,np), ritz, 1_wip, 
     $           ritz(max(kevd2,np) + 1), 1_wip )
            call qpswap ( 
     $           min(kevd2,np), bounds, 1_wip, 
     $           bounds(max(kevd2,np) + 1), 1_wip )
         endif
!
      else
!
!        %----------------------------------------------------%
!        | LM, SM, LA, SA case.                               |
!        | Sort the eigenvalues of H into the desired order   |
!        | and apply the resulting order to BOUNDS.           |
!        | The eigenvalues are sorted so that the wanted part |
!        | are always in the last KEV locations.               |
!        %----------------------------------------------------%
!
         call qpsortr (
     $        which, .true., kev+np, ritz, bounds )
      endif
!
      if (ishift .eq. 1 .and. np .gt. 0) then
!     
!        %-------------------------------------------------------%
!        | Sort the unwanted Ritz values used as shifts so that  |
!        | the ones with largest Ritz estimates are first.       |
!        | This will tend to minimize the effects of the         |
!        | forward instability of the iteration when the shifts  |
!        | are applied in subroutine qpsapps.                     |
!        %-------------------------------------------------------%
!     
         call qpsortr (
     $        'SM', .true., np, bounds, ritz )
         call qpcopy (
     $        np, ritz, 1_wip, shifts, 1_wip )
      endif
! 
      call qpsecond (t1)
      tsgets = tsgets + (t1 - t0)
!
      if (msglvl .gt. 0) then
         call wivout (
     $        logfil, 1_wip, (/kev/), ndigit, 
     $        '_sgets: KEV is' )
         call wivout (
     $        logfil, 1_wip, (/np/), ndigit, 
     $        '_sgets: NP is' )
         call qpvout (
     $        logfil, kev+np, ritz, ndigit,
     $        '_sgets: Eigenvalues of current H matrix' )
         call qpvout (
     $        logfil, kev+np, bounds, ndigit, 
     $        '_sgets: Associated Ritz estimates' )
      endif
! 
      return
      end subroutine qpsgets
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: qpstqrb
!
!\Description:
!  Computes all eigenvalues and the last component of the eigenvectors
!  of a symmetric tridiagonal matrix using the implicit QL or QR method.
!
!  This is mostly a modification of the LAPACK routine qpsteqr.
!  See Remarks.
!
!\Usage:
!  call qpstqrb
!     ( N, D, E, Z, WORK, INFO )
!
!\Arguments
!  N       Integer.  (INPUT)
!          The number of rows and columns in the matrix.  N >= 0.
!
!  D       real(wrp) array, dimension (N).  (INPUT/OUTPUT)
!          On entry, D contains the diagonal elements of the
!          tridiagonal matrix.
!          On exit, D contains the eigenvalues, in ascending order.
!          If an error exit is made, the eigenvalues are correct
!          for indices 1,2,...,INFO-1, but they are unordered and
!          may not be the smallest eigenvalues of the matrix.
!
!  E       real(wrp) array, dimension (N-1).  (INPUT/OUTPUT)
!          On entry, E contains the subdiagonal elements of the
!          tridiagonal matrix in positions 1 through N-1.
!          On exit, E has been destroyed.
!
!  Z       real(wrp) array, dimension (N).  (OUTPUT)
!          On exit, Z contains the last row of the orthonormal 
!          eigenvector matrix of the symmetric tridiagonal matrix.  
!          If an error exit is made, Z contains the last row of the
!          eigenvector matrix associated with the stored eigenvalues.
!
!  WORK    real(wrp) array, dimension (max(1,2*N-2)).  (WORKSPACE)
!          Workspace used in accumulating the transformation for 
!          computing the last components of the eigenvectors.
!
!  INFO    Integer.  (OUTPUT)
!          = 0:  normal return.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = +i, the i-th eigenvalue has not converged
!                              after a total of  30*N  iterations.
!
!\Remarks
!  1. None.
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     qpaxpy   Level 1 BLAS that computes a vector triad.
!     qpcopy   Level 1 BLAS that copies one vector to another.
!     qpswap   Level 1 BLAS that swaps the contents of two vectors.
!     wlsame   LAPACK character comparison routine.
!     qplae2   LAPACK routine that computes the eigenvalues of a 2-by-2 
!             symmetric matrix.
!     qplaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric 
!             matrix.
!     qplamch  LAPACK routine that determines machine constants.
!     qplanst  LAPACK routine that computes the norm of a matrix.
!     qplapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     qplartg  LAPACK Givens rotation construction routine.
!     qplascl  LAPACK routine for careful scaling of a matrix.
!     qplaset  LAPACK matrix initialization routine.
!     qplasr   LAPACK routine that applies an orthogonal transformation to 
!             a matrix.
!     qplasrt  LAPACK sorting routine.
!     qpsteqr  LAPACK routine that computes eigenvalues and eigenvectors
!             of a symmetric tridiagonal matrix.
!     wxerbla  LAPACK error handler routine.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University           
!     Houston, Texas            
!
!\SCCS Information: @(#) 
! FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.5, this routine is a modified version
!        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
!        only commeted out and new lines inserted.
!        All lines commented out have "c$$$" at the beginning.
!        Note that the LAPACK version 1.0 subroutine SSTEQR contained
!        bugs. 
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine qpstqrb ( 
     $           n, d, e, z, work, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer(wip) 
     $   n
      real(wrp)
     $   d(n), e(n-1), z(n), work(2*n-2)
      integer(wip) 
     $   info
!
!     LOCAL VARIABLES:
      real(wrp)               
     $   zero, one, two, three
      parameter ( 
     $   zero  = 0.0_wrp, 
     $   one   = 1.0_wrp, 
     $   two   = 2.0_wrp, 
     $   three = 3.0_wrp )
      integer 
     $   maxit
      parameter ( 
     $   maxit = 30 )
      integer(wip) 
     $   i, ii, j, jtot, k, l, l1, lend,
     $   lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1,
     $   nm1, nmaxit
      integer
     $   iscale, icompz
      real(wrp)               
     $   anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2,
     $   s, safmax, safmin, ssfmax, ssfmin, tst
!     ..
!     .. external functions ..
      logical 
     $   wlsame
      real(wrp)
     $   qplamch, 
     $   qplanst, 
     $   qplapy2
      external 
     $   wlsame, 
     $   qplamch, 
     $   qplanst, 
     $   qplapy2
!     ..
!     .. external subroutines ..
      external
     $   qplae2, 
     $   qplaev2, 
     $   qplartg, 
     $   qplascl, 
     $   qplaset, 
     $   qplasr,
     $   qplasrt, 
     $   qpswap, 
     $   xerbla
!     ..
!     .. intrinsic functions ..
      intrinsic
     $   abs, max, sign, sqrt
!     ..
!     .. executable statements ..
!
!     test the input parameters.
!
      info = 0
!
!$$$      IF( LSAME( COMPZ, 'N' ) ) THEN
!$$$         ICOMPZ = 0
!$$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
!$$$         ICOMPZ = 1
!$$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
!$$$         ICOMPZ = 2
!$$$      ELSE
!$$$         ICOMPZ = -1
!$$$      END IF
!$$$      IF( ICOMPZ.LT.0 ) THEN
!$$$         INFO = -1
!$$$      ELSE IF( N.LT.0 ) THEN
!$$$         INFO = -2
!$$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
!$$$     $         N ) ) ) THEN
!$$$         INFO = -6
!$$$      END IF
!$$$      IF( INFO.NE.0 ) THEN
!$$$         CALL XERBLA( 'SSTEQR', -INFO )
!$$$         RETURN
!$$$      END IF
!
!    *** New starting with version 2.5 ***
!
      icompz = 2
!    *************************************
!
!     quick return if possible
!
      if ( n .eq. 0 ) return
!
      if ( n .eq. 1 ) then
         if (icompz.eq.2)  z(1) = one
         return
      endif
!
!     determine the unit roundoff and over/underflow thresholds.
!
      eps = qplamch (
     $      'E' )
      eps2 = eps**2
      safmin = qplamch (
     $         'S' )
      safmax = one / safmin
      ssfmax = sqrt(safmax) / three
      ssfmin = sqrt(safmin) / eps2
!
!     compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
!$$      if (icompz.eq.2)
!$$$     $   call qplaset('full', n, n, zero, one, z, ldz)
!
!     *** New starting with version 2.5 ***
!
      if ( icompz .eq. 2 ) then
         do j = 1, n-1
            z(j) = zero
         enddo 
         z(n) = one
      endif
!     *************************************
!
      nmaxit = n*maxit
      jtot = 0
!
!     determine where the matrix splits and choose ql or qr iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      l1 = 1
      nm1 = n - 1
!
   10 continue
      if ( l1 .gt. n ) goto 160
      if ( l1 .gt. 1 ) e(l1-1) = zero
      if ( l1 .le. nm1 ) then
         do m = l1, nm1
            tst = abs(e(m))
            if ( tst .eq. zero ) goto 30
            if ( tst .le. ( sqrt(abs(d(m))) *
     $                      sqrt(abs(d(m+1))) ) * eps ) then
               e(m) = zero
               goto 30
            endif
         enddo 
      endif
      m = n
!
   30 continue
      l = l1
      lsv = l
      lend = m
      lendsv = lend
      l1 = m + 1
      if ( lend .eq. l ) goto 10
!
!     scale submatrix in rows and columns l to lend
!
      anorm = qplanst (
     $        'I', lend-l+1, d(l), e(l) )
      iscale = 0
      if ( anorm .eq. zero ) goto 10
      if ( anorm .gt. ssfmax ) then
         iscale = 1
         call qplascl (
     $        'G', 0_wip, 0_wip, anorm, ssfmax, lend-l+1, 1_wip, 
     $         d(l), n, info )
         call qplascl (
     $        'G', 0_wip, 0_wip, anorm, ssfmax, lend-l, 1_wip, 
     $        e(l), n, info )
      else if (anorm.lt.ssfmin) then
         iscale = 2
         call qplascl (
     $        'G', 0_wip, 0_wip, anorm, ssfmin, lend-l+1, 1_wip, 
     $        d(l), n, info )
         call qplascl (
     $        'G', 0_wip, 0_wip, anorm, ssfmin, lend-l, 1_wip, 
     $        e(l), n, info )
      endif
!
!     choose between ql and qr iteration
!
      if (abs(d(lend)).lt.abs(d(l))) then
         lend = lsv
         l = lendsv
      endif
!
      if ( lend .gt. l ) then
!
!        ql iteration
!
!        look for small subdiagonal element.
!
   40    continue
         if ( l .ne. lend ) then
            lendm1 = lend - 1
            do m = l, lendm1
               tst = abs(e(m))**2
               if ( tst .le. ( eps2 * abs(d(m)) ) * 
     $                       abs(d(m+1)) + safmin ) goto 60
            enddo 
         endif
!
         m = lend
!
   60    continue
         if ( m .lt. lend ) e(m) = zero
         p = d(l)
         if ( m .eq. l ) goto 80
!
!        if remaining matrix is 2-by-2, use qplae2 or qplaev2
!        to compute its eigensystem.
!
         if ( m .eq. l+1 ) then
            if ( icompz .gt. 0 ) then
               call qplaev2 (
     $              d(l), e(l), d(l+1), rt1, rt2, c, s )
               work(l) = c
               work(n-1+l) = s
!$$$               call qplasr('R', 'V', 'B', n, 2, work(l),
!$$$     $                     work(n-1+l), z(1, l), ldz)
!
!              *** New starting with version 2.5 ***
!
               tst    = z(l+1)
               z(l+1) = c*tst - s*z(l)
               z(l)   = s*tst + c*z(l)
!              *************************************
            else
               call qplae2 (
     $              d(l), e(l), d(l+1), rt1, rt2 )
            endif
            d(l) = rt1
            d(l+1) = rt2
            e(l) = zero
            l = l + 2
            if ( l .le. lend ) goto 40
            goto 140
         endif
!
         if ( jtot .eq. nmaxit ) goto 140
         jtot = jtot + 1
!
!        form shift.
!
         g = ( d(l+1) - p ) / ( two*e(l) )
         r = qplapy2 (
     $       g, one )
         g = d(m) - p + (e(l) / (g+sign(r, g)))
!
         s = one
         c = one
         p = zero
!
!        inner loop
!
         mm1 = m - 1
         do i = mm1, l, -1
            f = s*e(i)
            b = c*e(i)
            call qplartg (
     $           g, f, c, s, r)
            if (i.ne.m-1)
     $         e(i+1) = r
            g = d(i+1) - p
            r = (d(i)-g)*s + two*c*b
            p = s*r
            d(i+1) = g + p
            g = c*r - b
!
!           if eigenvectors are desired, then save rotations.
!
            if ( icompz .gt. 0 ) then
               work(i) = c
               work(n-1+i) = -s
            endif
!
         enddo 
!
!        if eigenvectors are desired, then apply saved rotations.
!
         if ( icompz .gt. 0 ) then
            mm = m - l + 1
!$$$            call qplasr('R', 'V', 'B', n, mm, work(l), work(n-1+l),
!$$$     $                  z(1, l), ldz)
!
!             *** New starting with version 2.5 ***
!
              call qplasr (
     $             'R', 'V', 'B', 1_wip, mm, work(l), 
     $             work(n-1+l), z(l), 1_wip )
!             *************************************                             
         endif
!
         d(l) = d(l) - p
         e(l) = g
         goto 40
!
!        eigenvalue found.
!
   80    continue
         d(l) = p
!
         l = l + 1
         if ( l .le. lend ) goto 40
         goto 140
!
      else
!
!        qr iteration
!
!        look for small superdiagonal element.
!
   90    continue
         if ( l .ne. lend ) then
            lendp1 = lend + 1
            do m = l, lendp1, -1
               tst = abs(e(m-1))**2
               if ( tst .le. ( eps2 * abs(d(m)) ) *
     $                         abs(d(m-1)) + safmin ) 
     $            goto 110
            enddo 
         endif
!
         m = lend
!
  110    continue
         if ( m .gt. lend ) e(m-1) = zero
         p = d(l)
         if ( m .eq. l ) goto 130
!
!        if remaining matrix is 2-by-2, use qplae2 or qplaev2
!        to compute its eigensystem.
!
         if ( m .eq. l-1 ) then
            if ( icompz .gt. 0 ) then
               call qplaev2 (
     $              d(l-1), e(l-1), d(l), rt1, rt2, c, s)
!$$$               work(m) = c
!$$$               work(n-1+m) = s
!$$$               call qplasr('R', 'V', 'F', n, 2, work(m),
!$$$     $                     work(n-1+m), z(1, l-1), ldz)
!
!               *** New starting with version 2.5 ***
!
                tst    = z(l)
                z(l)   = c*tst - s*z(l-1)
                z(l-1) = s*tst + c*z(l-1)
!               ************************************* 
            else
               call qplae2 (
     $              d(l-1), e(l-1), d(l), rt1, rt2)
            endif
            d(l-1) = rt1
            d(l) = rt2
            e(l-1) = zero
            l = l - 2
            if ( l .ge. lend ) goto 90
            goto 140
         endif
!
         if ( jtot .eq. nmaxit ) goto 140
         jtot = jtot + 1
!
!        form shift.
!
         g = ( d(l-1) - p ) / ( two * e(l-1) )
         r = qplapy2 (
     $       g, one )
         g = d(m) - p + ( e(l-1) /( g + sign(r,g) ) )
!
         s = one
         c = one
         p = zero
!
!        inner loop
!
         lm1 = l - 1
         do i = m, lm1
            f = s*e(i)
            b = c*e(i)
            call qplartg (
     $           g, f, c, s, r)
            if ( i .ne. m ) e(i-1) = r
            g = d(i) - p
            r = ( d(i+1) - g )*s + two*c*b
            p = s*r
            d(i) = g + p
            g = c*r - b
!
!           if eigenvectors are desired, then save rotations.
!
            if ( icompz .gt. 0 ) then
               work(i) = c
               work(n-1+i) = s
            endif
         enddo 
!
!        if eigenvectors are desired, then apply saved rotations.
!
         if ( icompz .gt. 0 ) then
            mm = l - m + 1
!$$$            call qplasr('R', 'V', 'F', n, mm, work(m), work(n-1+m),
!$$$     $                  z(1, m), ldz)
!
!           *** New starting with version 2.5 ***
!
            call qplasr (
     $           'R', 'V', 'F', 1_wip, mm, work(m), work(n-1+m),
     $           z(m), 1_wip )
!           *************************************                             
         endif
!
         d(l) = d(l) - p
         e(lm1) = g
         goto 90
!
!        eigenvalue found.
!
  130    continue
         d(l) = p
!
         l = l - 1
         if ( l .ge. lend ) goto 90
         goto 140
!
      endif
!
!     undo scaling if necessary
!
  140 continue
      if ( iscale .eq. 1 ) then
         call qplascl (
     $        'G', 0_wip, 0_wip, ssfmax, anorm, lendsv-lsv+1, 1_wip,
     $        d(lsv), n, info )
         call qplascl(
     $        'G', 0_wip, 0_wip, ssfmax, anorm, lendsv-lsv, 1_wip, 
     $        e(lsv), n, info )
      else if ( iscale .eq. 2 ) then
         call qplascl (
     $        'G', 0_wip, 0_wip, ssfmin, anorm, lendsv-lsv+1, 1_wip,
     $        d(lsv), n, info )
         call qplascl (
     $        'G', 0_wip, 0_wip, ssfmin, anorm, lendsv-lsv, 1_wip, 
     $        e(lsv), n, info )
      endif
!
!     check for no convergence to an eigenvalue after a total
!     of n*maxit iterations.
!
      if ( jtot .lt. nmaxit ) goto 10
      do i = 1, n-1
         if ( e(i) .ne. zero ) info = info + 1
      enddo 
      goto 190
!
!     order eigenvalues and eigenvectors.
!
  160 continue
      if ( icompz .eq. 0 ) then
!
!        use quick sort
!
         call qplasrt (
     $        'I', n, d, info )
!
      else
!
!        use selection sort to minimize swaps of eigenvectors
!
         do ii = 2, n
            i = ii - 1
            k = i
            p = d(i)
            do j = ii, n
               if ( d(j) .lt. p ) then
                  k = j
                  p = d(j)
               endif
            enddo 
            if ( k .ne. i ) then
               d(k) = d(i)
               d(i) = p
!$$$               call qpswap(n, z(1, i), 1, z(1, k), 1)
!           *** New starting with version 2.5 ***
!
               p    = z(k)
               z(k) = z(i)
               z(i) = p
!           *************************************
            endif
         enddo 
      endif
!
  190 continue
      return
      end subroutine qpstqrb
!
!***********************************************************************
!-----------------------------------------------------------------------
!***********************************************************************
!
!@@@  LAPACK::
!
!-----------------------------------------------------------------------
!
      subroutine qplae2 ( 
     $           a, b, c, rt1, rt2 )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      real(wrp)          
     $   a, b, c, rt1, rt2
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, and RT2
!  is the eigenvalue of smaller absolute value.
!
!  Arguments
!  =========
!
!  A       (input) REAL(wrp)       
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) REAL(wrp)       
!          The (1,2) and (2,1) elements of the 2-by-2 matrix.
!
!  C       (input) REAL(wrp)       
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) REAL(wrp)       
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) REAL(wrp)       
!          The eigenvalue of smaller absolute value.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     DEPENDENCES:
      intrinsic          
     $   abs, sqrt
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   one, two, zero, half
      parameter ( 
     $   one  = 1.0_wrp, two  = 2.0_wrp, 
     $   zero = 0.0_wrp, half = 0.5_wrp )
      real(wrp)          
     $   ab, acmn, acmx, adf, df, rt, sm, tb
!
!     EXECUTABLE STATEMENTS:
!     Compute the eigenvalues
!
      sm = a+c
      df = a-c
      adf = abs( df )
      tb = b+b
      ab = abs( tb )
      if ( abs(a).gt.abs(c) ) then
         acmx = a
         acmn = c
      else
         acmx = c
         acmn = a
      endif
      if ( adf.gt.ab ) then
         rt = adf*sqrt( one + ( ab / adf )**2 )
      else if ( adf.lt.ab ) then
         rt = ab*sqrt( one + ( adf / ab )**2 )
      else
!
!        Includes case AB=ADF=0
!
         rt = ab*sqrt( two )
      endif
      if ( sm.lt.zero ) then
         rt1 = half*( sm-rt )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
      else if ( sm.gt.zero ) then
         rt1 = half*( sm + rt )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
      else
!
!        Includes case RT1 = RT2 = 0
!
         rt1 = half*rt
         rt2 =-half*rt
      endif
      return
      end subroutine qplae2
!
!***********************************************************************
!
      subroutine qplaev2 ( 
     $           a, b, c, rt1, rt2, cs1, sn1 )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      real(wrp)          
     $   a, b, c, cs1, rt1, rt2, sn1
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!  eigenvector for RT1, giving the decomposition
!
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!
!  Arguments
!  =========
!
!  A       (input) REAL(wrp)       
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) REAL(wrp)       
!          The (1,2) element and the conjugate of the (2,1) element of
!          the 2-by-2 matrix.
!
!  C       (input) REAL(wrp)       
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) REAL(wrp)       
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) REAL(wrp)       
!          The eigenvalue of smaller absolute value.
!
!  CS1     (output) REAL(wrp)       
!  SN1     (output) REAL(wrp)       
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     DEPENDENCES:
      intrinsic          
     $   abs, sqrt
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   one, two, zero, half
      parameter ( 
     $   one  = 1.0_wrp, two  = 2.0_wrp, 
     $   zero = 0.0_wrp, half = 0.5_wrp )
      integer
     $   sgn1, sgn2
      real(wrp)          
     $   ab, acmn, acmx, acs, adf, cs, ct, df, 
     $   rt, sm, tb, tn
!
!     EXECUTABLE STATEMENTS:
!     Compute the eigenvalues
!
      sm = a + c
      df = a - c
      adf = abs( df )
      tb = b + b
      ab = abs( tb )
      if ( abs(a) .gt. abs(c) ) then
         acmx = a
         acmn = c
      else
         acmx = c
         acmn = a
      endif
      if ( adf .gt. ab ) then
         rt = adf*sqrt( one + ( ab / adf )**2 )
      else if ( adf .lt. ab ) then
         rt = ab*sqrt( one + ( adf / ab )**2 )
      else
!
!        Includes case AB=ADF=0
!
         rt = ab*sqrt( two )
      endif
      if ( sm .lt. zero ) then
         rt1 = half*( sm - rt )
         sgn1 = -1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
      else if ( sm .gt. zero ) then
         rt1 = half*( sm + rt )
         sgn1 = 1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
      else
!
!        Includes case RT1 = RT2 = 0
!
         rt1 =  half*rt
         rt2 = -half*rt
         sgn1 = 1
      endif
!
!     Compute the eigenvector
!
      if ( df .ge. zero ) then
         cs = df + rt
         sgn2 = 1
      else
         cs = df - rt
         sgn2 = -1
      endif
      acs = abs( cs )
      if ( acs .gt. ab ) then
         ct  = -tb / cs
         sn1 = one / sqrt( one + ct*ct )
         cs1 = ct*sn1
      else
         if ( ab .eq. zero ) then
            cs1 = one
            sn1 = zero
         else
            tn  = -cs / tb
            cs1 = one / sqrt( one + tn*tn )
            sn1 = tn*cs1
         endif
      endif
      if ( sgn1 .eq. sgn2 ) then
         tn  =  cs1
         cs1 = -sn1
         sn1 =  tn
      endif
      return
      end subroutine qplaev2
!
!***********************************************************************
!
      function qplamch ( 
     $         cmach )
!
      implicit none
      include 'epcode_inc_qp.h'
!
      real(wrp)          
     $   qplamch
!
!     ARGUMENTS:
      character          
     $   cmach
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAMCH determines REAL(wrp)        machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by WFLAMCH:
!          = 'E' or 'E',   WFLAMCH := eps       = epsilon(.)
!          = 'S' or 's ,   WFLAMCH := sfmin
!          = 'B' or 'B',   WFLAMCH := base
!          = 'P' or 'P',   WFLAMCH := eps*base
!          = 'N' or 'N',   WFLAMCH := t
!          = 'R' or 'R',   WFLAMCH := rnd
!          = 'M' or 'M',   WFLAMCH := emin
!          = 'U' or 'U',   WFLAMCH := rmin      = tiny(.)
!          = 'L' or 'L',   WFLAMCH := emax
!          = 'O' or 'O',   WFLAMCH := rmax      = huge(.)
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold-base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold -(base**emax)*(1-eps)
!
! =====================================================================
!     
!     DEPENDENCES: 
      logical            
     $   wlsame
      external           
     $   wlsame
      external           
     $   qplamc2
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      logical            
     $   first, lrnd
      integer
     $   beta, imax, imin, it
      real(wrp)          
     $   base, emax, emin, eps, prec, rmach, 
     $   rmax, rmin, rnd, sfmin, small, t
!
!     SAVE:
      save               
     $   first, eps, sfmin, base, t, rnd, 
     $   emin, rmin, emax, rmax, prec
!     
!     DATA:
      data               
     $   first / .true. /
!
!     EXECUTABLE STATEMENTS:
!
      if ( first ) then
         first = .false.
         call qplamc2 ( 
     $        beta, it, lrnd, eps, imin, rmin, imax, rmax )
         base = beta
         t = it
         if ( lrnd ) then
            rnd = one
            eps = ( base**( 1-it ) ) / 2
         else
            rnd = zero
            eps = base**( 1-it )
         endif
         prec = eps*base
         emin = imin
         emax = imax
         sfmin = rmin
         small = one / rmax
         if ( small.ge.sfmin ) then
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            sfmin = small*( one+eps )
         endif
      endif
!
      if (      wlsame ( cmach, 'E' ) ) then
         rmach = eps
      else if ( wlsame ( cmach, 'S' ) ) then
         rmach = sfmin
      else if ( wlsame ( cmach, 'B' ) ) then
         rmach = base
      else if ( wlsame ( cmach, 'P' ) ) then
         rmach = prec
      else if ( wlsame ( cmach, 'N' ) ) then
         rmach = t
      else if ( wlsame ( cmach, 'R' ) ) then
         rmach = rnd
      else if ( wlsame ( cmach, 'M' ) ) then
         rmach = emin
      else if ( wlsame ( cmach, 'U' ) ) then
         rmach = rmin
      else if ( wlsame ( cmach, 'L' ) ) then
         rmach = emax
      else if ( wlsame ( cmach, 'O' ) ) then
         rmach = rmax
      endif
!
      qplamch = rmach
      return
!
!     End of WFLAMCH
!
      end function qplamch 
!
!***********************************************************************
!
      subroutine qplamc1 ( 
     $           beta, t, rnd, ieee1 )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      logical            
     $   ieee1, rnd
      integer
     $   beta, t
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
! =====================================================================
!
!     DEPENDENCES:
      real(wrp)          
     $   qplamc3
      external           
     $   qplamc3
!
!     LOCAL VARIABLES:
      logical            
     $   first, lieee1, lrnd
      integer
     $   lbeta, lt
      real(wrp)          
     $   a, b, c, f, one, qtr, savec, t1, t2
!
!     SAVE:
      save               
     $   first, lieee1, lbeta, lrnd, lt
!     
!     DATA:
      data               
     $   first /.true./
!
!     EXECUTABLE STATEMENTS:
!
      if ( first ) then
         first = .false.
         one = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  WFLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a+1.0 ) = a.
!
         a = 1
         c = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   10    continue
         if ( c .eq. one ) then
            a = 2*a
            c = qplamc3 ( 
     $          a, one )
            c = qplamc3 ( 
     $          c, -a )
            goto 10
         endif
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a+b ) .gt. a.
!
         b = 1
         c = qplamc3 ( 
     $       a, b )
!
!+       WHILE( C.EQ.A )LOOP
   20    continue
         if ( c.eq.a ) then
            b = 2*b
            c = qplamc3 ( 
     $          a, b )
            goto 20
         endif
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t+1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta-1 ).
!
         qtr = one / 4
         savec = c
         c = qplamc3 ( 
     $       c, -a )
         lbeta = c+qtr
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         b = lbeta
         f = qplamc3 ( 
     $       b/2, -b/100 )
         c = qplamc3 ( 
     $       f, a )
         if ( c.eq.a ) then
            lrnd = .true.
         else
            lrnd = .false.
         endif
         f = qplamc3 ( 
     $       b/2, b/100 )
         c = qplamc3 ( 
     $       f, a )
         if ( ( lrnd ) .and. ( c.eq.a ) )
     $      lrnd = .false.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         t1 = qplamc3 ( 
     $        b/2, a )
         t2 = qplamc3 ( 
     $        b/2, savec )
         lieee1 = ( t1.eq.a ) .and. ( t2.gt.savec ) .and. lrnd
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t+1.0 ) = 1.0.
!
         lt = 0
         a = 1
         c = 1
!
!+       WHILE( C.EQ.ONE )LOOP
   30    continue
         if ( c.eq.one ) then
            lt = lt+1
            a = a*lbeta
            c = qplamc3 ( 
     $          a, one )
            c = qplamc3 ( 
     $          c, -a )
            goto 30
         endif
!+       END WHILE
!
      endif
!
      beta = lbeta
      t = lt
      rnd = lrnd
      ieee1 = lieee1
      return
!
!     End of WFLAMC1
!
      end subroutine qplamc1
!
!***********************************************************************
!
      subroutine qplamc2 ( 
     $           beta, t, rnd, eps, emin, rmin, emax, rmax )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      logical            
     $   rnd
      integer
     $   beta, emax, emin, t
      real(wrp)          
     $   eps, rmax, rmin
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) REAL(wrp)       
!          The smallest positive number such that
!
!             fl( 1.0-EPS ) .LT. 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) REAL(wrp)       
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN-1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) REAL(wrp)       
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1-EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     DEPENDENCES:
      real(wrp)          
     $   qplamc3
      external           
     $   qplamc3
      external           
     $   qplamc1, 
     $   qplamc4, 
     $   qplamc5
      intrinsic          
     $   abs, max, min
!
!     LOCAL VARIABLES:
      logical            
     $   first, ieee, iwarn, lieee1, lrnd
      integer
     $   gnmin, gpmin, i, lbeta, lemax, lemin,
     $   lt, ngnmin, ngpmin
      real(wrp)          
     $   a, b, c, half, leps, lrmax, lrmin, one,
     $   rbase, sixth, small, third, two, zero
!     
!     SAVE:
      save               
     $   first, iwarn, lbeta, lemax, lemin, 
     $   leps, lrmax, lrmin, lt
!
!     DATA:
      data               
     $   first /.true./, 
     $   iwarn /.false./
!
!     EXECUTABLE STATEMENTS:
!
      if ( first ) then
         first = .false.
         zero = 0
         one = 1
         two = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  WFLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        WFLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         call qplamc1 ( 
     $        lbeta, lt, lrnd, lieee1 )
!
!        Start to find EPS.
!
         b = lbeta
         a = b**( -lt )
         leps = a
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         b = two / 3
         half = one / 2
         sixth = qplamc3 ( 
     $           b, -half )
         third = qplamc3 ( 
     $           sixth, sixth )
         b = qplamc3 ( 
     $       third, -half )
         b = qplamc3 ( 
     $       b, sixth )
         b = abs(b)
         if ( b.lt.leps ) b = leps
!
         leps = 1
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    continue
         if ( ( leps.gt.b ) .and. ( b.gt.zero ) ) then
            leps = b
            c = qplamc3 ( 
     $          half*leps, ( two**5 )*( leps**2 ) )
            c = qplamc3 ( 
     $          half, -c )
            b = qplamc3 ( 
     $          half, c )
            c = qplamc3 ( 
     $          half, -b )
            b = qplamc3 ( 
     $          half, c )
            goto 10
         endif
!+       END WHILE
!
         if ( a.lt.leps ) leps = a
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A =+or-1, and+or-(1+BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         rbase = one / lbeta
         small = one
         do i = 1, 3
            small = qplamc3 ( 
     $              small*rbase, zero )
         enddo 
         a = qplamc3 ( 
     $       one, small )
         call qplamc4 ( 
     $        ngpmin, one, lbeta )
         call qplamc4 ( 
     $        ngnmin, -one, lbeta )
         call qplamc4 ( 
     $        gpmin, a, lbeta )
         call qplamc4 ( 
     $        gnmin, -a, lbeta )
         ieee = .false.
!
         if ( ( ngpmin.eq.ngnmin ) .and. 
     $        ( gpmin.eq.gnmin )        ) then
            if ( ngpmin.eq.gpmin ) then
               lemin = ngpmin
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            else if ( ( gpmin-ngpmin ).eq.3 ) then
               lemin = ngpmin-1+lt
               ieee = .true.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            else
               lemin = min( ngpmin, gpmin )
!            ( A guess; no known machine )
               iwarn = .true.
            endif
!
         else if ( ( ngpmin.eq.gpmin ) .and. 
     $             ( ngnmin.eq.gnmin )     ) then
            if ( abs( ngpmin-ngnmin ).eq.1 ) then
               lemin = max( ngpmin, ngnmin )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            else
               lemin = min( ngpmin, ngnmin )
!            ( A guess; no known machine )
               iwarn = .true.
            endif
!
         else if ( ( abs( ngpmin-ngnmin ).eq.1 ) .and.
     $             ( gpmin.eq.gnmin )                 ) then
            if ( ( gpmin-min( ngpmin, ngnmin ) ).eq.3 ) then
               lemin = max( ngpmin, ngnmin )-1+lt
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            else
               lemin = min( ngpmin, ngnmin )
!            ( A guess; no known machine )
               iwarn = .true.
            endif
!
         else
            lemin = min( ngpmin, ngnmin, gpmin, gnmin )
!         ( A guess; no known machine )
            iwarn = .true.
         endif
!**
! Comment out this if block if EMIN is ok
         if ( iwarn ) then
            first = .true.
            write( 6, fmt = 9999 ) lemin
         endif
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine WFLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
         ieee = ieee .or. lieee1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN-1 ),  but some machines underflow during
!        this computation.
!
         lrmin = 1
         do i = 1, 1-lemin
            lrmin = qplamc3 ( 
     $              lrmin*rbase, zero )
         enddo 
!
!        Finally, call WFLAMC5 to compute EMAX and RMAX.
!
         call qplamc5 ( 
     $        lbeta, lt, lemin, ieee, lemax, lrmax )
      endif
!
      beta = lbeta
      t = lt
      rnd = lrnd
      eps = leps
      emin = lemin
      rmin = lrmin
      emax = lemax
      rmax = lrmax
!
      return
 9999 format( / / ' warning. the value emin may be incorrect:-',
     $      '  emin = ', i8, /
     $      ' if, after inspection, the value emin looks',
     $      ' acceptable please comment out ',
     $      / ' the if block as marked within the code of routine',
     $      ' qplamc2,', / ' otherwise supply emin explicitly.', / )
      end subroutine qplamc2
!
!***********************************************************************
!
      function qplamc3 ( 
     $         a, b )
!
      implicit none
      include 'epcode_inc_qp.h'
!
      real(wrp)          
     $   qplamc3
!
!     ARGUMENTS:
      real(wrp)          
     $   a, b
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) REAL(wrp)       
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      qplamc3 = a+b
!
      return
      end function qplamc3
!
!***********************************************************************
!
      subroutine qplamc4 ( 
     $           emin, start, base )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer
     $   base, emin
      real(wrp)          
     $   start
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAMC4 is a service routine for WFLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) REAL(wrp)       
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!     
!     DEPENDENCES:
      real(wrp)          
     $   qplamc3
      external           
     $   qplamc3
!
!     LOCAL VARIABLES:
      integer
     $   i
      real(wrp)          
     $   a, b1, b2, c1, c2, d1, d2, 
     $   one, rbase, zero
!
!     EXECUTABLE STATEMENTS:
!
      a = start
      one = 1
      rbase = one / base
      zero = 0
      emin = 1
      b1 = qplamc3 ( a*rbase, zero )
      c1 = a
      c2 = a
      d1 = a
      d2 = a
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 continue
      if ( ( c1.eq.a ) .and. 
     $     ( c2.eq.a ) .and. 
     $     ( d1.eq.a ) .and. 
     $     ( d2.eq.a )      ) then
         emin = emin-1
         a = b1
         b1 = qplamc3 ( 
     $        a / base, zero )
         c1 = qplamc3 ( 
     $        b1*base, zero )
         d1 = zero
         do i = 1, base
            d1 = d1+b1
         enddo 
         b2 = qplamc3 ( 
     $        a*rbase, zero )
         c2 = qplamc3 ( 
     $        b2 / rbase, zero )
         d2 = zero
         do i = 1, base
            d2 = d2+b2
         enddo 
         goto 10
      endif
!+    END WHILE
!
      return
      end subroutine qplamc4 
!
!***********************************************************************
!
      subroutine qplamc5 ( 
     $           beta, p, emin, ieee, emax, rmax )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      logical            
     $   ieee
      integer
     $   beta, emax, emin, p
      real(wrp)          
     $   rmax
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX+abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) REAL(wrp)       
!          The largest machine floating-point number.
!
! =====================================================================
!     
!     DEPENDENCES:
      real(wrp)          
     $   qplamc3
      external           
     $   qplamc3
      intrinsic          
     $   mod
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   zero, one
      parameter ( 
     $   zero = 0.0_wrp, 
     $   one  = 1.0_wrp )
      integer(wip)       
     $   exbits, expsum, i, lexp, nbits, try, uexp
      real(wrp)          
     $   oldy, recbas, y, z
!
!     EXECUTABLE STATEMENTS:
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX+abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      lexp = 1
      exbits = 1
   10 continue
      try = lexp*2
      if ( try.le.( -emin ) ) then
         lexp = try
         exbits = exbits+1
         goto 10
      endif
      if ( lexp.eq.-emin ) then
         uexp = lexp
      else
         uexp = try
         exbits = exbits+1
      endif
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      if ( ( uexp+emin ).gt.( -lexp-emin ) ) then
         expsum = 2*lexp
      else
         expsum = 2*uexp
      endif
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX-EMIN+1 .
!
      emax = expsum+emin-1
      nbits = 1+exbits+p
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
      if ( ( mod( nbits, 2 ).eq.1 ) .and. ( beta.eq.2 ) ) then
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         emax = emax-1
      endif
!
      if ( ieee ) then
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         emax = emax-1
      endif
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0-BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0-BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      recbas = one / beta
      z = beta-one
      y = zero
      do i = 1, p
         z = z*recbas
         if ( y.lt.one ) oldy = y
         y = qplamc3 ( 
     $       y, z )
      enddo 
      if ( y.ge.one ) y = oldy
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      do i = 1, emax
         y = qplamc3 ( 
     $       y*beta, zero )
      enddo 
!
      rmax = y
!
      return
      end subroutine qplamc5 
!
!***********************************************************************
!
      function qplanst ( 
     $         norm, n, d, e )
!
      implicit none
      include 'epcode_inc_qp.h'
!
      real(wrp)          
     $   qplanst
!
!     ARGUMENTS:
      character          
     $   norm
      integer(wip)       
     $   n
      real(wrp)          
     $   d(*), e(*)
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLANST  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric tridiagonal matrix A.
!
!  Description
!  ===========
!
!  WFLANST returns the value
!
!    WFLANST = ( max(abs(A(i,j))), NORM = 'M' or 'M'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'O'
!              (
!              ( normI(A),         NORM = 'I' or 'I'
!              (
!              ( normF(A),         NORM = 'F', 'F', 'E' or 'E'
!  where  
!     norm1  denotes the  one norm of a matrix (maximum column sum),
!     normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!     normF  denotes the  Frobenius norm of a matrix (square root of sum of
!            squares).  
!  Note that  
!     max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in WFLANST as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, WFLANST is
!          set to zero.
!
!  D       (input) REAL(wrp)        array, dimension (N)
!          The diagonal elements of A.
!
!  E       (input) REAL(wrp)        array, dimension (N-1)
!          The (n-1) sub-diagonal or super-diagonal elements of A.
!
!  =====================================================================
!
!     DEPENDENCES:
      logical            
     $   wlsame
      external           
     $   wlsame
      external           
     $   qplassq
      intrinsic          
     $   abs, max, sqrt
!
!     LOCAL VARIABLES:
      real(wrp)  
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      integer(wip)       
     $   i
      real(wrp)          
     $   anorm, rscale, rsum 
!
!     EXECUTABLE STATEMENTS:
!
      if ( n.le.0 ) then
         anorm = zero
      else if ( wlsame ( norm, 'M' ) ) then
!
!        Find max(abs(A(i,j))).
!
         anorm = abs( d(n) )
         do i = 1, n-1
            anorm = max( anorm, abs( d(i) ) )
            anorm = max( anorm, abs( e(i) ) )
         enddo 
      else if ( norm.eq.'1'         .or.
     $         wlsame ( norm, 'O' ) .or. 
     $         wlsame ( norm, 'I' )     ) then
!
!        Find norm1(A).
!
         if ( n.eq.1 ) then
            anorm = abs( d(1) )
         else
            anorm = max( abs(d(1))  +abs(e(1)),
     $                   abs(e(n-1))+abs(d(n)) )
            do i = 2, n-1
               anorm = max( anorm, abs(d(i))+abs(e(i)) +
     $                             abs(e(i-1)) )
            enddo 
         endif
      else if ( ( wlsame ( norm, 'F' ) ) .or. 
     $          ( wlsame ( norm, 'E' ) )     ) then
!
!        Find normF(A).
!
         rscale = zero
         rsum = one
         if ( n.gt.1 ) then
            call qplassq ( 
     $           n-1, e, 1_wip, rscale, rsum )
            rsum = 2*rsum 
         endif
         call qplassq ( 
     $        n, d, 1_wip, rscale, rsum )
         anorm = rscale*sqrt( rsum )
      endif
!
      qplanst = anorm
      return
      end function qplanst 
!
!***********************************************************************
!
      subroutine qplassq ( 
     $           n, x, incx, rscale, sumsq )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:   
      integer(wip) 
     $   n, incx
      real(wrp)  
     $   x(*), 
     $   rscale, sumsq
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x(1)**2 +...+ x(n)**2+( scale**2 )*sumsq,
!
!  where  x(i) = X( 1+( i-1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x(i) ) ).
!
!  scale and sumsq must be supplied in RSCALE and SUMSQ and
!  scl and smsq are overwritten on RSCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) REAL(wrp)        array, dimension (N)
!          The vector for which a scaled sum of squares is computed.
!             x(i)  = X( 1+( i-1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  RSCALE  (input/output) REAL(wrp)       
!          On entry, the value  scale  in the equation above.
!          On exit, RSCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) REAL(wrp)       
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
! =====================================================================
!
!     DEPENDENCES: 
      intrinsic   
     $   abs
!
!     LOCAL VARIABLES: 
      real(wrp)   
     $   zero
      parameter ( 
     $   zero = 0.0_wrp )
      integer(wip) 
     $   ix
      real(wrp) 
     $   absxi
!
!     EXECUTABLE STATEMENTS:
!
      if ( n .gt. 0_wip ) then
         do ix = 1, 1+(n-1)*incx, incx
            if ( x(ix) .ne. zero ) then
               absxi = abs( x(ix) )
               if ( rscale .lt. absxi ) then
                  sumsq = 1 + sumsq*( rscale / absxi )**2
                  rscale = absxi
               else
                  sumsq = sumsq + ( absxi / rscale )**2
               endif
            endif
         enddo 
      endif
      return
      end subroutine qplassq 
!
!***********************************************************************
!
      function qplapy2 ( 
     $         x, y )
!
      implicit none
      include 'epcode_inc_qp.h'
!
      real(wrp)   
     $   qplapy2
!
!     ARGUMENTS:
      real(wrp)
     $   x, y
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) REAL(wrp)       
!  Y       (input) REAL(wrp)       
!          X and Y specify the values x and y.
!
!  =====================================================================
!     
!     DEPENDENCES: 
      intrinsic          
     $   abs, max, min, sqrt
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   zero, one
      parameter ( 
     $   zero = 0.0_wrp, 
     $   one  = 1.0_wrp )
      real(wrp)          
     $   w, xabs, yabs, z
!
!     EXECUTABLE STATEMENTS:
!
      xabs = abs(x)
      yabs = abs(y)
      w = max( xabs, yabs )
      z = min( xabs, yabs )
      if ( z.eq.zero ) then
         qplapy2 = w
      else
         qplapy2 = w*sqrt( one + ( z / w )**2 )
      endif
      return
      end function qplapy2  
!
!
!***********************************************************************
!
      subroutine qplartg ( 
     $           f, g, cs, sn, r )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      real(wrp) 
     $   f, g, cs, sn, r
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLARTG generate a plane rotation so that
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2+SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a slower, more accurate version of the BLAS1 routine DROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in DBDSQR when
!        there are zeros on the diagonal).
!
!  If F exceeds G in magnitude, CS will be positive.
!
!  Arguments
!  =========
!
!  F       (input) REAL(wrp)       
!          The first component of vector to be rotated.
!
!  G       (input) REAL(wrp)       
!          The second component of vector to be rotated.
!
!  CS      (output) REAL(wrp)       
!          The cosine of the rotation.
!
!  SN      (output) REAL(wrp)       
!          The sine of the rotation.
!
!  R       (output) REAL(wrp)       
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!    
!     DEPENDENCES:
      real(wrp)   
     $   qplamch
      external 
     $   qplamch
      intrinsic  
     $   abs, int, log, max, sqrt
!
!     LOCAL VARIABLES:
      real(wrp)  
     $   zero, one, two
      parameter ( 
     $   zero = 0.0_wrp, 
     $   one  = 1.0_wrp, 
     $   two  = 2.0_wrp )
      logical       
     $   first
      integer(wip) 
     $   icount, i
      real(wrp)
     $   eps, f1, g1, safmin, safmn2, safmx2, rscale 
!
!     SAVE:
      save  
     $   first, safmx2, safmin, safmn2
!
!     DATA:
      data 
     $   first / .true. /
!
!     EXECUTABLE STATEMENTS:
!
      if ( first ) then
         first = .false.
         safmin = qplamch ('S')
         eps = qplamch ('E')
         safmn2 = qplamch ('B')**int( 
     $            log(safmin/eps) / log( qplamch ('B') ) / two )
         safmx2 = one / safmn2
      endif
      if ( g.eq.zero ) then
         cs = one
         sn = zero
         r = f
      else if ( f.eq.zero ) then
         cs = zero
         sn = one
         r = g
      else
         f1 = f
         g1 = g
         rscale = max( abs( f1 ), abs( g1 ) )
         if ( rscale.ge.safmx2 ) then
            icount = 0
   10       continue
            icount = icount+1
            f1 = f1*safmn2
            g1 = g1*safmn2
            rscale = max( abs( f1 ), abs( g1 ) )
            if ( rscale.ge.safmx2 ) goto 10
            r = sqrt( f1**2 + g1**2 )
            cs = f1 / r
            sn = g1 / r
            do i = 1, icount 
               r = r*safmx2
            enddo 
         else if ( rscale.le.safmn2 ) then
            icount = 0
   30       continue
            icount = icount+1
            f1 = f1*safmx2
            g1 = g1*safmx2
            rscale = max( abs( f1 ), abs( g1 ) )
            if ( rscale.le.safmn2 ) goto 30
            r = sqrt( f1**2 + g1**2 )
            cs = f1 / r
            sn = g1 / r
            do i = 1, icount 
               r = r*safmn2
            enddo 
         else
            r  = sqrt( f1**2 + g1**2 )
            cs = f1 / r
            sn = g1 / r
         endif
         if ( abs(f).gt.abs(g) .and. cs.lt.zero ) then
            cs = -cs
            sn = -sn
            r  = -r
         endif
      endif
      return
      end subroutine qplartg
!
!***********************************************************************
!
      subroutine qplascl ( 
     $           ctype, kl, ku, cfrom, cto, m, n, a, lda, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character  
     $   ctype 
      integer(wip) 
     $   kl, ku, lda, m, n
      real(wrp) 
     $   cfrom, cto
      real(wrp)  
     $   a(lda,*)
      integer(wip) 
     $   info
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLASCL multiplies the M by N real matrix A by the real scalar
!  CTO/CFROM.  This is done without over/underflow as long as the final
!  result CTO*A(I,J)/CFROM does not over/underflow. CTYPE specifies that
!  A may be full, upper triangular, lower triangular, upper Hessenberg,
!  or banded.
!
!  Arguments
!  =========
!
!  CTYPE    (input) CHARACTER*1
!          CTYPE indices the storage type of the input matrix.
!          = 'G':  A is a full matrix.
!          = 'L':  A is a lower triangular matrix.
!          = 'U':  A is an upper triangular matrix.
!          = 'H':  A is an upper Hessenberg matrix.
!          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the lower
!                  half stored.
!          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the upper
!                  half stored.
!          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!                  bandwidth KU.
!
!  KL      (input) INTEGER
!          The lower bandwidth of A.  Referenced only if CTYPE = 'B',
!          'Q' or 'Z'.
!
!  KU      (input) INTEGER
!          The upper bandwidth of A.  Referenced only if CTYPE = 'B',
!          'Q' or 'Z'.
!
!  CFROM   (input) REAL(wrp)       
!  CTO     (input) REAL(wrp)       
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,M)
!          The matrix to be multiplied by CTO/CFROM.  See CTYPE for the
!          storage type.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  INFO    (output) INTEGER
!          0 -successful exit
!          <0-if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     DEPENDENCES:
      logical    
     $   wlsame
      real(wrp)    
     $   qplamch
      external  
     $   wlsame, 
     $   qplamch
      intrinsic   
     $   abs, max, min
      external  
     $   wxerbla
!
!     LOCAL VARIABLES: 
      real(wrp)  
     $   zero, one
      parameter ( 
     $   zero = 0.0_wrp, 
     $   one  = 1.0_wrp )
      logical     
     $   done
      integer(wip) 
     $   i, itype, j, k1, k2, k3, k4
      real(wrp) 
     $   bignum, cfrom1, cfromc, cto1, 
     $   ctoc, mul, smlnum
!
!     EXECUTABLE STATEMENTS:
!     Test the input arguments
!
      info = 0
!
      if (      wlsame ( ctype, 'G' ) ) then
         itype = 0
      else if ( wlsame ( ctype, 'L' ) ) then
         itype = 1
      else if ( wlsame ( ctype, 'U' ) ) then
         itype = 2
      else if ( wlsame ( ctype, 'H' ) ) then
         itype = 3
      else if ( wlsame ( ctype, 'B' ) ) then
         itype = 4
      else if ( wlsame ( ctype, 'Q' ) ) then
         itype = 5
      else if ( wlsame ( ctype, 'Z' ) ) then
         itype = 6
      else
         itype = -1
      endif
!
      if ( itype.eq.-1 ) then
         info = -1
      else if ( cfrom.eq.zero ) then
         info = -4
      else if ( m.lt.0 ) then
         info = -6
      else if ( n.lt.0 .or. ( itype.eq.4 .and. n.ne.m ) .or.
     $         ( itype.eq.5 .and. n.ne.m ) ) then
         info = -7
      else if ( itype.le.3 .and. lda.lt.max( 1_wip, m ) ) then
         info = -9
      else if ( itype.ge.4 ) then
         if ( kl.lt.0 .or. kl.gt.max( m-1, 0_wip ) ) then
            info = -2
         else if ( ku.lt.0 .or. ku.gt.max( n-1, 0_wip ) .or.
     $            ((itype.eq.4 .or. itype.eq.5) .and. kl.ne.ku) )
     $             then
            info = -3
         else if ( ( itype.eq.4 .and. lda.lt.kl+1 ) .or.
     $             ( itype.eq.5 .and. lda.lt.ku+1 ) .or.
     $             ( itype.eq.6 .and. lda.lt.2*kl+ku+1 ) ) then
            info = -9
         endif
      endif
!
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wflascl', -info )
         return
      endif
!
!     Quick return if possible
!
      if ( n.eq.0 .or. m.eq.0 ) return
!
!     Get machine parameters
!
      smlnum = qplamch ( 'S' )
      bignum = one / smlnum
!
      cfromc = cfrom
      ctoc = cto
!
   10 continue
      cfrom1 = cfromc*smlnum
      cto1 = ctoc / bignum
      if ( abs( cfrom1 ).gt.abs( ctoc ) .and. ctoc.ne.zero ) then
         mul = smlnum
         done = .false.
         cfromc = cfrom1
      else if ( abs( cto1 ).gt.abs( cfromc ) ) then
         mul = bignum
         done = .false.
         ctoc = cto1
      else
         mul = ctoc / cfromc
         done = .true.
      endif
!
      if ( itype.eq.0 ) then
!
!        Full matrix
!
         do j = 1, n
            do i = 1, m
               a(i,j) = a(i,j)*mul
            enddo 
         enddo 
!
      else if ( itype.eq.1 ) then
!
!        Lower triangular matrix
!
         do j = 1, n
            do i = j, m
               a(i,j) = a(i,j)*mul
            enddo 
         enddo 
!
      else if ( itype.eq.2 ) then
!
!        Upper triangular matrix
!
         do j = 1, n
            do i = 1, min(j,m)
               a(i,j) = a(i,j)*mul
            enddo 
         enddo 
!
      else if ( itype.eq.3 ) then
!
!        Upper Hessenberg matrix
!
         do j = 1, n
            do i = 1, min(j+1,m)
               a(i,j) = a(i,j)*mul
            enddo 
         enddo 
!
      else if ( itype.eq.4 ) then
!
!        Lower half of a symmetric band matrix
!
         k3 = kl+1
         k4 = n+1
         do j = 1, n
            do i = 1, min(k3,k4-j)
               a(i,j) = a(i,j)*mul
            enddo 
         enddo 
!
      else if ( itype.eq.5 ) then
!
!        Upper half of a symmetric band matrix
!
         k1 = ku+2
         k3 = ku+1
         do j = 1, n
            do i = max(k1-j,1_wip), k3
               a(i,j) = a(i,j)*mul
            enddo 
         enddo 
!
      else if ( itype.eq.6 ) then
!
!        Band matrix
!
         k1 = kl + ku + 2
         k2 = kl + 1
         k3 = 2*kl + ku + 1
         k4 = kl + ku + 1 + m
         do j = 1, n
            do i = max(k1-j,k2), min(k3,k4-j)
               a(i,j) = a(i,j)*mul
            enddo 
         enddo 
!
      endif
!
      if ( .not.done ) goto 10
!
      return
      end subroutine qplascl 
!
!***********************************************************************
!
      subroutine qplaset ( 
     $           uplo, m, n, alpha, beta, a, lda )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character     
     $   uplo
      integer(wip)  
     $   m, n, lda 
      real(wrp)  
     $   alpha, beta,
     $   a(lda,*)
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) REAL(wrp)       
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) REAL(wrp)       
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================
!
!     DEPENDENCES:
      logical     
     $   wlsame
      external  
     $   wlsame
      intrinsic   
     $   min
!
!     LOCAL VARIABLES:
      integer(wip)    
     $   i, j
!
!     EXECUTABLE STATEMENTS:
!
      if ( wlsame ( uplo, 'U' ) ) then
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         do j = 2, n
            do i = 1, min(j-1,m)
               a(i,j) = alpha
            enddo 
         enddo 
!
      else if ( wlsame ( uplo, 'L' ) ) then
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         do j = 1, min(m,n)
            do i = j+1, m
               a(i,j) = alpha
            enddo 
         enddo 
!
      else
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         do j = 1, n
            do i = 1, m
               a(i,j) = alpha
            enddo 
         enddo 
      endif
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      do i = 1, min(m,n)
         a(i,i) = beta
      enddo 
!
      return
      end subroutine qplaset
!
!***********************************************************************
!
      subroutine qplasr ( 
     $           side, pivot, cdirect, m, n, c, s, a, lda )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character  
     $   side, pivot, cdirect 
      integer(wip)  
     $   m, n, lda 
      real(wrp)    
     $   a(lda,*), c(*), s(*)
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLASR   performs the transformation
!
!     A := P*A,   when SIDE = 'L' or 'L'  (  Left-hand side )
!
!     A := A*P',  when SIDE = 'R' or 'R'  ( Right-hand side )
!
!  where A is an m by n real matrix and P is an orthogonal matrix,
!  consisting of a sequence of plane rotations determined by the
!  parameters PIVOT and CDIRECT as follows ( z = m when SIDE = 'L' or 'L'
!  and z = n when SIDE = 'R' or 'R' ):
!
!  When  CDIRECT = 'F' or 'F'  ( Forward sequence ) then
!
!     P = P( z-1 )*...*P(2)*P(1),
!
!  and when CDIRECT = 'B' or 'B'  ( Backward sequence ) then
!
!     P = P(1)*P(2)*...*P( z-1 ),
!
!  where  P(k) is a plane rotation matrix for the following planes:
!
!     when  PIVOT = 'V' or 'V'  ( Variable pivot ),
!        the plane (k, k+1)
!
!     when  PIVOT = 'T' or 'T'  ( Top pivot ),
!        the plane (1, k+1)
!
!     when  PIVOT = 'B' or 'B'  ( Bottom pivot ),
!        the plane (k, z)
!
!  c(k) and s(k)  must contain the  cosine and sine that define the
!  matrix  P(k).  The two by two plane rotation part of the matrix
!  P(k), R(k), is assumed to be of the form
!
!     R(k) = (  c(k)  s(k) ).
!              ( -s(k)  c(k) )
!
!  This version vectorises across rows of the array A when SIDE = 'L'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P'
!
!  CDIRECT (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P( z-1 )*...*P(2)*P(1)
!          = 'B':  Backward, P = P(1)*P(2)*...*P( z-1 )
!
!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
!  C, S    (input) REAL(wrp)        arrays, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          c(k) and s(k) contain the cosine and sine that define the
!          matrix P(k).  The two by two plane rotation part of the
!          matrix P(k), R(k), is assumed to be of the form
!          R(k) = (  c(k)  s(k) ).
!                   ( -s(k)  c(k) )
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          The m by n matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P' if SIDE = 'L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     DEPENDENCES:
      logical  
     $   wlsame
      external 
     $   wlsame
      external   
     $   wxerbla
      intrinsic   
     $   max
!
!     LOCAL VARIABLES:
      real(wrp) 
     $   one, zero
      parameter (
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      integer(wip)  
     $   i, j, info
      real(wrp)  
     $   ctemp, stemp, temp
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters
!
      info = 0
      if ( .not.( wlsame ( side, 'L' ) .or. 
     $            wlsame ( side, 'R' )     ) ) then
         info = 1
      else if ( .not.( wlsame ( pivot, 'V' ) .or. 
     $                 wlsame ( pivot, 'T' ) .or. 
     $                 wlsame ( pivot, 'B' )     ) ) then
         info = 2
      else if ( .not.( wlsame ( cdirect, 'F') .or. 
     $                 wlsame ( cdirect, 'B')     ) ) then
         info = 3
      else if ( m.lt.0 ) then
         info = 4
      else if ( n.lt.0 ) then
         info = 5
      else if ( lda.lt.max( 1_wip, m ) ) then
         info = 9
      endif
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wflasr ', info )
         return
      endif
!
!     Quick return if possible
!
      if ( ( m.eq.0 ) .or. ( n.eq.0 ) ) return
      if ( wlsame ( side, 'L' ) ) then
!
!        Form  P * A
!
         if ( wlsame ( pivot, 'V' ) ) then
            if ( wlsame ( cdirect, 'F' ) ) then
               do j = 1, m-1
                  ctemp = c(j)
                  stemp = s(j)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, n
                        temp = a(j+1,i)
                        a(j+1,i) = ctemp*temp - stemp*a(j,i)
                        a(  j,i) = stemp*temp + ctemp*a(j,i)
                     enddo 
                  endif
               enddo 
            else if ( wlsame ( cdirect, 'B' ) ) then
               do j = m-1, 1, -1
                  ctemp = c(j)
                  stemp = s(j)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, n
                        temp = a(j+1,i)
                        a(j+1,i) = ctemp*temp - stemp*a(j,i)
                        a(  j,i) = stemp*temp + ctemp*a(j,i)
                     enddo 
                  endif
               enddo 
            endif
         else if ( wlsame ( pivot, 'T' ) ) then
            if ( wlsame ( cdirect, 'F' ) ) then
               do j = 2, m
                  ctemp = c(j-1)
                  stemp = s(j-1)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, n
                        temp = a(j,i)
                        a(j,i) = ctemp*temp - stemp*a(1,i)
                        a(1,i) = stemp*temp + ctemp*a(1,i)
                     enddo 
                  endif
               enddo 
            else if ( wlsame ( cdirect, 'B' ) ) then
               do j = m, 2, -1
                  ctemp = c( j-1 )
                  stemp = s( j-1 )
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, n
                        temp = a(j,i)
                        a(j,i) = ctemp*temp - stemp*a(1,i)
                        a(1,i) = stemp*temp + ctemp*a(1,i)
                     enddo 
                  endif
               enddo 
            endif
         else if ( wlsame ( pivot, 'B' ) ) then
            if ( wlsame ( cdirect, 'F' ) ) then
               do j = 1, m-1
                  ctemp = c(j)
                  stemp = s(j)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, n
                        temp = a(j,i)
                        a(j,i) = stemp*a(m,i) + ctemp*temp
                        a(m,i) = ctemp*a(m,i) - stemp*temp
                     enddo 
                  endif
               enddo 
            else if ( wlsame ( cdirect, 'B' ) ) then
               do j = m-1, 1, -1
                  ctemp = c(j)
                  stemp = s(j)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, n
                        temp = a(j,i)
                        a(j,i) = stemp*a(m,i) + ctemp*temp
                        a(m,i) = ctemp*a(m,i) - stemp*temp
                     enddo 
                  endif
               enddo 
            endif
         endif
      else if ( wlsame ( side, 'R' ) ) then
!
!        Form A * P'
!
         if ( wlsame ( pivot, 'V' ) ) then
            if ( wlsame ( cdirect, 'F' ) ) then
               do j = 1, n-1
                  ctemp = c(j)
                  stemp = s(j)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, m
                        temp = a(i,j+1)
                        a(i,j+1) = ctemp*temp - stemp*a(i,j)
                        a(i,  j) = stemp*temp + ctemp*a(i,j)
                     enddo 
                  endif
               enddo
            else if ( wlsame ( cdirect, 'B' ) ) then
               do j = n-1, 1, -1
                  ctemp = c(j)
                  stemp = s(j)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, m
                        temp = a(i,j+1)
                        a(i,j+1) = ctemp*temp - stemp*a(i,j)
                        a(i,  j) = stemp*temp + ctemp*a(i,j)
                     enddo 
                  endif
               enddo 
            endif
         else if ( wlsame ( pivot, 'T' ) ) then
            if ( wlsame ( cdirect, 'F' ) ) then
               do j = 2, n
                  ctemp = c(j-1)
                  stemp = s(j-1)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, m
                        temp = a(i,j)
                        a(i,j) = ctemp*temp - stemp*a(i,1)
                        a(i,1) = stemp*temp + ctemp*a(i,1)
                     enddo 
                  endif
               enddo 
            else if ( wlsame ( cdirect, 'B' ) ) then
               do j = n, 2, -1
                  ctemp = c(j-1)
                  stemp = s(j-1)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, m
                        temp = a(i,j)
                        a(i,j) = ctemp*temp - stemp*a(i,1)
                        a(i,1) = stemp*temp + ctemp*a(i,1)
                     enddo 
                  endif
               enddo 
            endif
         else if ( wlsame ( pivot, 'B' ) ) then
            if ( wlsame ( cdirect, 'F' ) ) then
               do j = 1, n-1
                  ctemp = c(j)
                  stemp = s(j)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, m
                        temp = a(i,j)
                        a(i,j) = stemp*a(i,n) + ctemp*temp
                        a(i,n) = ctemp*a(i,n) - stemp*temp
                     enddo 
                  endif
               enddo 
            else if ( wlsame ( cdirect, 'B' ) ) then
               do j = n-1, 1, -1
                  ctemp = c(j)
                  stemp = s(j)
                  if ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) ) then
                     do i = 1, m
                        temp = a(i,j)
                        a(i,j) = stemp*a(i,n) + ctemp*temp
                        a(i,n) = ctemp*a(i,n) - stemp*temp
                     enddo 
                  endif
               enddo 
            endif
         endif
      endif
!
      return
      end subroutine qplasr 
!
!***********************************************************************
!
!     Sorting arrays using stack
!
      subroutine qplasrt ( 
     $           id, n, d, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      character     
     $   id
      integer(wip) 
     $   n
      real(wrp)  
     $   d(*)
      integer(wip) 
     $   info
! 
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments
!  =========
!
!  ID      (input) CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) INTEGER
!          The length of the array D.
!
!  D       (input/output) REAL(wrp)        array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     DEPENDENCES: 
      logical   
     $   wlsame
      external  
     $   wlsame
      external   
     $   wxerbla
!
!     LOCAL VARIABLES: 
      integer
     $   iselect 
      parameter ( 
     $   iselect = 20 )
      integer(wip) 
     $   dir, endd, i, j, start, stkpnt 
      real(wrp)   
     $   d1, d2, d3, dmnmx, tmp
!
!!    STACK: (to limits N to about 2**32)
!!    integer(wip) 
!!   $   stack ( 2, 32 )
!
!     STACK: (it limits N to about 2**62 < huge(1_wip) = 2**63-1 ???)
      integer(wip) 
     $   stack ( 2, 62 )
!
!     EXECUTABLE STATEMENTS:
!     Test the input paramters.
!
      info = 0
      dir = -1
      if (      wlsame ( id, 'D' ) ) then
         dir = 0
      else if ( wlsame ( id, 'I' ) ) then
         dir = 1
      endif
      if ( dir .eq. -1 ) then
         info = -1
      else if ( n .lt. 0 ) then
         info = -2
      endif
      if ( info .ne. 0 ) then
         call wxerbla ( 
     $        'wflasrt', -info )
         return
      endif
!
!     Quick return if possible
!
      if ( n .le. 1_wip ) return
!
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = n
!
   10 continue
!
      start  = stack(1,stkpnt)
      endd   = stack(2,stkpnt)
      stkpnt = stkpnt - 1
!
      if ( endd-start .le. iselect .and. endd-start .gt. 0 ) then
!
!        Do Insertion sort on D( START:ENDD )
!
         if ( dir .eq. 0 ) then
!
!           Sort into decreasing order
!
            do i = start+1, endd
               do j = i, start+1, -1
                  if ( d(j) .gt. d(j-1) ) then
                     dmnmx  = d(j)
                     d(j)   = d(j-1)
                     d(j-1) = dmnmx
                  else
                     exit 
                  endif
               enddo 
            enddo 
!
         else
!
!           Sort into increasing order
!
            do i = start+1, endd
               do j = i, start+1, -1
                  if ( d(j) .lt. d(j-1) ) then
                     dmnmx  = d(j)
                     d(j)   = d(j-1)
                     d(j-1) = dmnmx
                  else
                     exit 
                  endif
               enddo 
            enddo 
!
         endif
!
      else if ( endd-start .gt. iselect ) then
!
!        Do Quick sort 
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
         d1 = d(start)
         d2 = d(endd)
         i  = ( start + endd ) / 2
         d3 = d(i)
         if ( d1 .lt. d2 ) then
            if ( d3 .lt. d1 ) then
               dmnmx = d1
            else if ( d3 .lt. d2 ) then
               dmnmx = d3
            else
               dmnmx = d2
            endif
         else
            if ( d3 .lt. d2 ) then
               dmnmx = d2
            else if ( d3 .lt. d1 ) then
               dmnmx = d3
            else
               dmnmx = d1
            endif
         endif
!
         if ( dir .eq. 0 ) then
!
!           Sort into decreasing order
!
            i = start - 1
            j = endd  + 1
   60       continue
   70       continue
            j = j-1
            if ( d(j) .lt. dmnmx ) goto 70
   80       continue
            i = i+1
            if ( d(i) .gt. dmnmx ) goto 80
            if ( i .lt. j ) then
               tmp  = d(i)
               d(i) = d(j)
               d(j) = tmp
               goto 60
            endif
            if ( j-start .gt. endd-j-1 ) then
               stkpnt = stkpnt+1
               stack(1,stkpnt) = start
               stack(2,stkpnt) = j
               stkpnt = stkpnt+1
               stack(1,stkpnt) = j+1
               stack(2,stkpnt) = endd
            else
               stkpnt = stkpnt+1
               stack(1,stkpnt) = j+1
               stack(2,stkpnt) = endd
               stkpnt = stkpnt+1
               stack(1,stkpnt) = start
               stack(2,stkpnt) = j
            endif
         else
!
!           Sort into increasing order
!
            i = start - 1
            j = endd  + 1
   90       continue
  100       continue
            j = j-1
            if ( d(j) .gt. dmnmx ) goto 100
  110       continue
            i = i+1
            if ( d(i) .lt. dmnmx ) goto 110
            if ( i .lt. j ) then
               tmp  = d(i)
               d(i) = d(j)
               d(j) = tmp
               goto 90
            endif
            if ( j-start .gt. endd-j-1 ) then
               stkpnt = stkpnt+1
               stack(1,stkpnt) = start
               stack(2,stkpnt) = j
               stkpnt = stkpnt+1
               stack(1,stkpnt) = j+1
               stack(2,stkpnt) = endd
            else
               stkpnt = stkpnt+1
               stack(1,stkpnt) = j+1
               stack(2,stkpnt) = endd
               stkpnt = stkpnt+1
               stack(1,stkpnt) = start
               stack(2,stkpnt) = j
            endif
         endif
      endif
      if ( stkpnt .gt. 0 ) goto 10
      return
      end subroutine qplasrt
!
!***********************************************************************
!
      subroutine qpsteqr ( 
     $           compz, n, d, e, z, ldz, work, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      character   
     $   compz
      integer(wip) 
     $   n, ldz 
      real(wrp)   
     $   d(*), e(*), work(*), z(ldz,*)
      integer(wip) 
     $   info
! 
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band symmetric matrix can also be found
!  if QPSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
!  tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  symmetric matrix.  On entry, Z must contain the
!                  orthogonal matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) REAL(wrp)        array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) REAL(wrp)        array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) REAL(wrp)        array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,N).
!
!  WORK    (workspace) REAL(wrp)        array, dimension (max(1,2*N-2))
!          If COMPZ = 'N', then WORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm has failed to find all the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero; on exit, D
!                and E contain the elements of a symmetric tridiagonal
!                matrix which is orthogonally similar to the original
!                matrix.
!
!  =====================================================================
!
!     DEPENDENCES: 
      logical 
     $   wlsame
      real(wrp) 
     $   qplamch, 
     $   qplanst, 
     $   qplapy2
      external
     $   wlsame, 
     $   qplamch, 
     $   qplanst, 
     $   qplapy2
      external
     $   qplae2,  
     $   qplaev2, 
     $   qplartg, 
     $   qplascl, 
     $   qplaset, 
     $   qplasr,
     $   qplasrt, 
     $   qpswap, 
     $   wxerbla
      intrinsic
     $   abs, max, sign, sqrt
!
!     LOCAL VARIABLES: 
      real(wrp)  
     $   zero, one, two, three
      parameter ( 
     $   zero  = 0.0_wrp, 
     $   one   = 1.0_wrp, 
     $   two   = 2.0_wrp,  
     $   three = 3.0_wrp )
      integer
     $   maxit
      parameter ( 
     $   maxit = 30 )
      integer
     $   icompz
      integer(wip) 
     $   i, ii, iscale, j, jtot, k, l, 
     $   l1, lend, lendm1, lendp1, lendsv, 
     $   lm1, lsv, m, mm, mm1, nm1, nmaxit
      real(wrp)  
     $   anorm, b, c, eps, eps2, f, g, p, r, 
     $   rt1, rt2, s, safmax, safmin, 
     $   ssfmax, ssfmin, tst
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      info = 0
!
      if (      wlsame ( compz, 'N' ) ) then
         icompz = 0
      else if ( wlsame ( compz, 'V' ) ) then
         icompz = 1
      else if ( wlsame ( compz, 'I' ) ) then
         icompz = 2
      else
         icompz = -1
      endif
      if ( icompz .lt. 0  ) then
         info = -1
      else if ( n .lt. 0_wip  ) then
         info = -2
      else if ( ( ldz .lt. 1_wip  )     .or. 
     $          ( icompz .gt. 0  .and. 
     $            ldz .lt. max(1_wip,n) )    ) then
         info = -6
      endif
      if ( info .ne. 0 ) then
         call wxerbla ( 
     $        'wfsteqr', -info )
         return
      endif
!
!     Quick return if possible
!
      if ( n .eq. 0_wip  )
     $   return
!
      if ( n.eq.1_wip  ) then
         if ( icompz.eq.2  )
     $      z(1,1) = one
         return
      endif
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      eps  = qplamch ( 'E' )
      eps2 = eps**2
      safmin = qplamch ( 'S' )
      safmax = one / safmin
      ssfmax = sqrt( safmax ) / three
      ssfmin = sqrt( safmin ) / eps2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      if ( icompz .eq. 2 )
     $   call qplaset ( 
     $        'Full', n, n, zero, one, z, ldz )
!
      nmaxit = n*maxit
      jtot = 0_wip 
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      l1 = 1_wip 
      nm1 = n - 1
!
   10 continue
      if ( l1 .gt. n ) goto 160
      if ( l1 .gt. 1_wip ) e(l1-1) = zero
      if ( l1 .le. nm1 ) then
         do m = l1, nm1
            tst = abs( e(m) )
            if ( tst .eq. zero ) goto 30
            if ( tst .le. ( sqrt( abs(d(m  )) )*
     $                      sqrt( abs(d(m+1)) ) )*eps ) then
               e(m) = zero
               goto 30
            endif
         enddo 
      endif
      m = n
!
   30 continue
      l = l1
      lsv = l
      lend = m
      lendsv = lend
      l1 = m + 1
      if ( lend .eq. l ) goto 10
!
!     Scale submatrix in rows and columns L to LEND
!
      anorm = qplanst ( 
     $        'I', lend-l+1, d(l), e(l) )
      iscale = 0
      if ( anorm .eq. zero ) goto 10
      if ( anorm .gt. ssfmax ) then
         iscale = 1_wip 
         call qplascl ( 
     $        'G', 0_wip, 0_wip, anorm, ssfmax, lend-l+1, 
     $        1_wip, d(l), n, info )
         call qplascl ( 
     $        'G', 0_wip, 0_wip, anorm, ssfmax, lend-l, 
     $        1_wip, e(l), n, info )
      else if ( anorm .lt. ssfmin ) then
         iscale = 2_wip 
         call qplascl ( 
     $        'G', 0_wip, 0_wip, anorm, ssfmin, lend-l+1, 
     $        1_wip, d(l), n, info )
         call qplascl ( 
     $        'G', 0_wip, 0_wip, anorm, ssfmin, lend-l, 
     $        1_wip, e(l), n, info )
      endif
!
!     Choose between QL and QR iteration
!
      if ( abs( d(lend) ) .lt. abs( d(l) ) ) then
         lend = lsv
         l = lendsv
      endif
!
      if ( lend .gt. l ) then
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    continue
         if ( l .ne. lend ) then
            lendm1 = lend-1
            do m = l, lendm1
               tst = abs( e(m) )**2
               if ( tst .le. ( eps2*abs(d(m)) ) *
     $                       abs(d(m+1)) + safmin ) goto 60
            enddo 
         endif
!
         m = lend
!
   60    continue
         if ( m .lt. lend ) e(m) = zero
         p = d(l)
         if ( m .eq. l ) goto 80
!
!        If remaining matrix is 2-by-2, use WFLAE2 or SLAEV2
!        to compute its eigensystem.
!
         if ( m .eq. l+1 ) then
            if ( icompz .gt. 0 ) then
               call qplaev2 ( 
     $              d(l), e(l), d(l+1), rt1, rt2, c, s )
               work(l)     = c
               work(n-1+l) = s
               call qplasr ( 
     $              'R', 'V', 'B', n, 2_wip, work(l),
     $              work(n-1+l), z(1,l), ldz )
            else
               call qplae2 ( 
     $              d(l), e(l), d(l+1), rt1, rt2 )
            endif
            d(l)   = rt1
            d(l+1) = rt2
            e(l)   = zero
            l = l + 2
            if ( l .le. lend ) goto 40
            goto 140
         endif
!
         if ( jtot .eq. nmaxit ) goto 140
         jtot = jtot + 1
!
!        Form shift.
!
         g = ( d(l+1) - p ) / ( two*e(l) )
         r = qplapy2 ( 
     $       g, one )
         g = d(m) - p + ( e(l) / ( g+sign(r, g) ) )
!
         s = one
         c = one
         p = zero
!
!        Inner loop
!
         mm1 = m - 1
         do i = mm1, l, -1
            f = s*e(i)
            b = c*e(i)
            call qplartg ( 
     $           g, f, c, s, r )
            if ( i .ne. m-1 ) e(i+1) = r
            g = d(i+1) - p
            r = ( d(i) - g )*s + two*c*b
            p = s*r
            d(i+1) = g + p
            g = c*r - b
!
!           If eigenvectors are desired, then save rotations.
!
            if ( icompz .gt. 0 ) then
               work(i)     = c
               work(n-1+i) = -s
            endif
         enddo 
!
!        If eigenvectors are desired, then apply saved rotations.
!
         if ( icompz .gt. 0 ) then
            mm = m-l+1
            call qplasr ( 
     $           'R', 'V', 'B', n, mm, work(l), work(n-1+l),
     $           z(1,l), ldz )
         endif
!
         d(l) = d(l) - p
         e(l) = g
         goto 40
!
!        Eigenvalue found.
!
   80    continue
         d(l) = p
!
         l = l + 1
         if ( l .le. lend ) goto 40
         goto 140
!
      else
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    continue
         if ( l .ne. lend ) then
            lendp1 = lend + 1
            do m = l, lendp1, -1
               tst = abs( e(m-1) )**2
               if ( tst .le. ( eps2*abs(d(m)) )*
     $             abs(d(m-1)) + safmin ) goto 110
            enddo 
         endif
!
         m = lend
!
  110    continue
         if ( m .gt. lend ) e(m-1) = zero
         p = d(l)
         if ( m .eq. l ) goto 130
!
!        If remaining matrix is 2-by-2, use WFLAE2 or SLAEV2
!        to compute its eigensystem.
!
         if ( m .eq. l-1 ) then
            if ( icompz .gt. 0 ) then
               call qplaev2 ( 
     $              d(l-1), e(l-1), d(l), rt1, rt2, c, s )
               work(m)     = c
               work(n-1+m) = s
               call qplasr( 
     $              'R', 'V', 'F', n, 2_wip, work(m),
     $              work(n-1+m), z(1,l-1), ldz )
            else
               call qplae2 ( 
     $              d(l-1), e(l-1), d(l), rt1, rt2 )
            endif
            d(l-1) = rt1
            d(l)   = rt2
            e(l-1) = zero
            l = l - 2
            if ( l .ge. lend ) goto 90
            goto 140
         endif
!
         if ( jtot .eq. nmaxit ) goto 140
         jtot = jtot + 1
!
!        Form shift.
!
         g = ( d(l-1) - p ) / ( two*e(l-1) )
         r = qplapy2 ( 
     $       g, one )
         g = d(m) - p + ( e(l-1) / ( g + sign(r, g) ) )
!
         s = one
         c = one
         p = zero
!
!        Inner loop
!
         lm1 = l - 1
         do i = m, lm1
            f = s*e(i)
            b = c*e(i)
            call qplartg ( 
     $           g, f, c, s, r )
            if ( i .ne. m ) e(i-1) = r
            g = d(i) - p
            r = ( d(i+1) - g )*s + two*c*b
            p = s*r
            d(i) = g + p
            g = c*r - b
!
!           If eigenvectors are desired, then save rotations.
!
            if ( icompz .gt. 0 ) then
               work(i)     = c
               work(n-1+i) = s
            endif
         enddo 
!
!        If eigenvectors are desired, then apply saved rotations.
!
         if ( icompz .gt. 0 ) then
            mm = l-m+1
            call qplasr ( 
     $           'R', 'V', 'F', n, mm, work(m), work(n-1+m),
     $           z(1,m), ldz )
         endif
!
         d(l) = d(l) - p
         e(lm1) = g
         goto 90
!
!        Eigenvalue found.
!
  130    continue
         d(l) = p
!
         l = l - 1
         if ( l .ge. lend ) goto 90
         goto 140
!
      endif
!
!     Undo scaling if necessary
!
  140 continue
      if ( iscale .eq. 1 ) then
         call qplascl ( 
     $        'G', 0_wip, 0_wip, ssfmax, anorm, lendsv-lsv+1, 
     $        1_wip, d(lsv), n, info )
         call qplascl ( 
     $        'G', 0_wip, 0_wip, ssfmax, anorm, lendsv-lsv, 
     $        1_wip, e(lsv), n, info )
      else if ( iscale .eq. 2 ) then
         call qplascl ( 
     $        'G', 0_wip, 0_wip, ssfmin, anorm, lendsv-lsv+1, 
     $        1_wip, d(lsv), n, info )
         call qplascl ( 
     $        'G', 0_wip, 0_wip, ssfmin, anorm, lendsv-lsv, 
     $        1_wip, e(lsv), n, info )
      endif
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      if ( jtot .lt. nmaxit ) goto 10
      do i = 1, n-1
         if ( e(i) .ne. zero ) info = info + 1
      enddo 
      goto 190
!
!     Order eigenvalues and eigenvectors.
!
  160 continue
      if ( icompz .eq. 0 ) then
!
!        Use Quick Sort
!
         call qplasrt ( 
     $        'I', n, d, info )
!
      else
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         do ii = 2, n
            i = ii-1
            k = i
            p = d(i)
            do j = ii, n
               if ( d(j) .lt. p ) then
                  k = j
                  p = d(j)
               endif
            enddo 
            if ( k .ne. i ) then
               d(k) = d(i)
               d(i) = p
               call qpswap ( 
     $              n, z(1,i), 1_wip, z(1,k), 1_wip )
            endif
         enddo 
      endif
!
  190 continue
      return
      end subroutine qpsteqr 
!
!***********************************************************************
!
!***  BLAS routines required are moved to dev_blas.f temporally.
!
!     We already have all the BLAS required by both LAPACK and ARPACK.
!
!     Here we add more LAPACK routines, which are required by ARPACK.
!     They are as follows:
!
!     ?  dlacpy  LAPACK matrix copy routine.
!     ?  dlarnv  LAPACK routine for generating a random vector.
!     ?  dgeqr2  LAPACK routine that computes the QR factorization of a matrix.
!     ?  dorm2r  LAPACK routine that applies an orthogonal matrix in factored form.
!
!***********************************************************************
!
      subroutine qplacpy ( 
     $           uplo, m, n, a, lda, b, ldb )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      character 
     $   uplo
      integer(wip)
     $   m, n, lda, ldb 
      real(wrp)  
     $   a(lda,*), b(ldb,*)
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  DLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================
!
!     DEPENDENCES: 
      logical
     $   wlsame
      external
     $   wlsame
      intrinsic
     $   min
!
!     LOCAL VARIABLES: 
      integer(wip)
     $   i, j
!    
!     EXECUTABLE STATEMENTS:
!
      if ( wlsame ( uplo, 'U' ) ) then
         do j = 1, n
            do i = 1, min(j,m)
               b(i,j) = a(i,j)
            enddo 
         enddo 
      else if ( wlsame ( uplo, 'L' ) ) then
         do j = 1, n
            do i = j, m
               b(i,j) = a(i,j)
            enddo 
         enddo 
      else
         do j = 1, n
            do i = 1, m
               b(i,j) = a(i,j)
            enddo 
         enddo 
      endif
      return
      end subroutine qplacpy  
!
!***********************************************************************
!
      subroutine qplarnv ( 
     $           idist, iseed, n, x )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer
     $   idist, 
     $   iseed(4)
      integer(wip)
     $   n
      real(wrp) 
     $   x(*)
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!     
!  Purpose
!  =======
!
!  DLARNV returns a vector of n random real numbers from a uniform or
!  normal distribution.
!
!  Arguments
!  =========
!
!  IDIST   (input) INTEGER
!          Specifies the distribution of the random numbers:
!          = 1:  uniform (0,1)
!          = 2:  uniform (-1,1)
!          = 3:  normal (0,1)
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated.
!
!  X       (output) DOUBLE PRECISION array, dimension (N)
!          The generated random numbers.
!
!  Further Details
!  ===============
!
!  This routine calls the auxiliary routine DLARUV to generate random
!  real numbers from a uniform (0,1) distribution, in batches of up to
!  128 using vectorisable code. The Box-Muller method is used to
!  transform numbers from a uniform to a normal distribution.
!
!  =====================================================================
!    
!     DEPENDENCES: 
      intrinsic
     $   cos, log, min, sqrt
      external 
     $   qplaruv
!
!     LOCAL VARIABLES: 
      real(wrp)  
     $   one, two
      parameter ( 
     $   one = 1.0_wrp, 
     $   two = 2.0_wrp )
      integer
     $   lv
      parameter ( 
     $   lv = 128 )
      real(wrp)
     $   u(lv)
      real(wrp)   
     $   twopi
      parameter ( 
     $   twopi = 6.2831853071795864769252867663_wrp )
      integer(wip)
     $   i, il, iv, il2
!     
!     EXECUTABLE STATEMENTS:
!
      do iv = 1, n, lv / 2
         il = min( lv / 2, n-iv+1 )
         if ( idist .eq. 3 ) then
            il2 = 2*il
         else
            il2 = il
         endif
!
!        Call DLARUV to generate IL2 numbers from a uniform (0,1)
!        distribution (IL2 <= LV)
!
         call qplaruv ( 
     $        iseed, il2, u )
!
         if ( idist .eq. 1 ) then
!
!           Copy generated numbers
!
            do i = 1, il
               x(iv+i-1) = u(i)
            enddo 
         else if ( idist .eq. 2 ) then
!
!           Convert generated numbers to uniform (-1,1) distribution
!
            do i = 1, il
               x(iv+i-1) = two*u(i) - one
            enddo 
         else if ( idist .eq. 3 ) then
!
!           Convert generated numbers to normal (0,1) distribution
!
            do i = 1, il
               x(iv+i-1) = sqrt( -two*log( u(2*i-1) ) )*
     $                     cos( twopi*u(2*i) )
            enddo 
         endif
      enddo 
      return
      end subroutine qplarnv 
!
!***********************************************************************
!
      subroutine qplaruv ( 
     $           iseed, n, x )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer
     $   iseed(4)
      integer(wip)
     $   n
      real(wrp) 
     $   x(n)
!
!  -- lapack auxiliary routine (version 3.0) --
!     univ. of tennessee, univ. of california berkeley, nag ltd.,
!     courant institute, argonne national lab, and rice university
!     october 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  purpose
!  =======
!
!  qplaruv returns a vector of n random real numbers from a uniform (0,1)
!  distribution (n <= 128).
!
!  this is an auxiliary routine called by dlarnv and zlarnv.
!
!  arguments
!  =========
!
!  iseed   (input/output) integer array, dimension (4)
!          on entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and iseed(4) must be
!          odd.
!          on exit, the seed is updated.
!
!  n       (input) integer
!          the number of random numbers to be generated. n <= 128.
!
!  x       (output) real(wrp) array, dimension (n)
!          the generated random numbers.
!
!  further details
!  ===============
!
!  this routine uses a multiplicative congruential method with modulus
!  2**48 and multiplier 33952834046453 (see g.s.fishman,
!  'multiplicative congruential random number generators with modulus
!  2**b: an exhaustive analysis for b = 32 and a partial analysis for
!  b = 48', math. comp. 189, pp 331-344, 1990).
!
!  48-bit integers are stored in 4 integer array elements with 12 bits
!  per element. hence the routine is portable across machines with
!  integers of 32 bits or more.
!
!  =====================================================================
!    
!     DEPENDENCES:  
      intrinsic
     $   min, mod, real 
!
!     LOCAL VARIABLES: 
      real(wrp)
     $   one
      parameter ( 
     $   one = 1.0_wrp )
      integer
     $   lv, ipw2
      parameter ( 
     $   lv = 128, 
     $   ipw2 = 4096 ) 
      real(wrp)
     $   r
      parameter ( 
     $   r = one / ipw2 )
      integer  
     $   mm( lv, 4 )
      integer
     $   i, i1, i2, i3, i4, it1, it2, it3, it4, j
!     
!     DATA:
      data   
     $   ( mm(  1, j ), j = 1, 4 ) /  494,  322, 2508, 2549 /,
     $   ( mm(  2, j ), j = 1, 4 ) / 2637,  789, 3754, 1145 /,
     $   ( mm(  3, j ), j = 1, 4 ) /  255, 1440, 1766, 2253 /,
     $   ( mm(  4, j ), j = 1, 4 ) / 2008,  752, 3572,  305 /,
     $   ( mm(  5, j ), j = 1, 4 ) / 1253, 2859, 2893, 3301 /,
     $   ( mm(  6, j ), j = 1, 4 ) / 3344,  123,  307, 1065 /,
     $   ( mm(  7, j ), j = 1, 4 ) / 4084, 1848, 1297, 3133 /,
     $   ( mm(  8, j ), j = 1, 4 ) / 1739,  643, 3966, 2913 /,
     $   ( mm(  9, j ), j = 1, 4 ) / 3143, 2405,  758, 3285 /,
     $   ( mm( 10, j ), j = 1, 4 ) / 3468, 2638, 2598, 1241 /,
     $   ( mm( 11, j ), j = 1, 4 ) /  688, 2344, 3406, 1197 /,
     $   ( mm( 12, j ), j = 1, 4 ) / 1657,   46, 2922, 3729 /,
     $   ( mm( 13, j ), j = 1, 4 ) / 1238, 3814, 1038, 2501 /,
     $   ( mm( 14, j ), j = 1, 4 ) / 3166,  913, 2934, 1673 /,
     $   ( mm( 15, j ), j = 1, 4 ) / 1292, 3649, 2091,  541 /,
     $   ( mm( 16, j ), j = 1, 4 ) / 3422,  339, 2451, 2753 /,
     $   ( mm( 17, j ), j = 1, 4 ) / 1270, 3808, 1580,  949 /,
     $   ( mm( 18, j ), j = 1, 4 ) / 2016,  822, 1958, 2361 /,
     $   ( mm( 19, j ), j = 1, 4 ) /  154, 2832, 2055, 1165 /,
     $   ( mm( 20, j ), j = 1, 4 ) / 2862, 3078, 1507, 4081 /,
     $   ( mm( 21, j ), j = 1, 4 ) /  697, 3633, 1078, 2725 /,
     $   ( mm( 22, j ), j = 1, 4 ) / 1706, 2970, 3273, 3305 /,
     $   ( mm( 23, j ), j = 1, 4 ) /  491,  637,   17, 3069 /,
     $   ( mm( 24, j ), j = 1, 4 ) /  931, 2249,  854, 3617 /,
     $   ( mm( 25, j ), j = 1, 4 ) / 1444, 2081, 2916, 3733 /,
     $   ( mm( 26, j ), j = 1, 4 ) /  444, 4019, 3971,  409 /,
     $   ( mm( 27, j ), j = 1, 4 ) / 3577, 1478, 2889, 2157 /,
     $   ( mm( 28, j ), j = 1, 4 ) / 3944,  242, 3831, 1361 /,
     $   ( mm( 29, j ), j = 1, 4 ) / 2184,  481, 2621, 3973 /,
     $   ( mm( 30, j ), j = 1, 4 ) / 1661, 2075, 1541, 1865 /,
     $   ( mm( 31, j ), j = 1, 4 ) / 3482, 4058,  893, 2525 /,
     $   ( mm( 32, j ), j = 1, 4 ) /  657,  622,  736, 1409 /,
     $   ( mm( 33, j ), j = 1, 4 ) / 3023, 3376, 3992, 3445 /,
     $   ( mm( 34, j ), j = 1, 4 ) / 3618,  812,  787, 3577 /,
     $   ( mm( 35, j ), j = 1, 4 ) / 1267,  234, 2125,   77 /,
     $   ( mm( 36, j ), j = 1, 4 ) / 1828,  641, 2364, 3761 /,
     $   ( mm( 37, j ), j = 1, 4 ) /  164, 4005, 2460, 2149 /,
     $   ( mm( 38, j ), j = 1, 4 ) / 3798, 1122,  257, 1449 /,
     $   ( mm( 39, j ), j = 1, 4 ) / 3087, 3135, 1574, 3005 /,
     $   ( mm( 40, j ), j = 1, 4 ) / 2400, 2640, 3912,  225 /,
     $   ( mm( 41, j ), j = 1, 4 ) / 2870, 2302, 1216,   85 /,
     $   ( mm( 42, j ), j = 1, 4 ) / 3876,   40, 3248, 3673 /,
     $   ( mm( 43, j ), j = 1, 4 ) / 1905, 1832, 3401, 3117 /,
     $   ( mm( 44, j ), j = 1, 4 ) / 1593, 2247, 2124, 3089 /,
     $   ( mm( 45, j ), j = 1, 4 ) / 1797, 2034, 2762, 1349 /,
     $   ( mm( 46, j ), j = 1, 4 ) / 1234, 2637,  149, 2057 /,
     $   ( mm( 47, j ), j = 1, 4 ) / 3460, 1287, 2245,  413 /,
     $   ( mm( 48, j ), j = 1, 4 ) / 328, 1691,  166,    65 /,
     $   ( mm( 49, j ), j = 1, 4 ) / 2861, 496,  466,  1845 /,
     $   ( mm( 50, j ), j = 1, 4 ) / 1950, 1597, 4018,  697 /,
     $   ( mm( 51, j ), j = 1, 4 ) /  617, 2394, 1399, 3085 /,
     $   ( mm( 52, j ), j = 1, 4 ) / 2070, 2584,  190, 3441 /,
     $   ( mm( 53, j ), j = 1, 4 ) / 3331, 1843, 2879, 1573 /,
     $   ( mm( 54, j ), j = 1, 4 ) /  769,  336,  153, 3689 /,
     $   ( mm( 55, j ), j = 1, 4 ) / 1558, 1472, 2320, 2941 /,
     $   ( mm( 56, j ), j = 1, 4 ) / 2412, 2407,   18,  929 /,
     $   ( mm( 57, j ), j = 1, 4 ) / 2800,  433,  712,  533 /,
     $   ( mm( 58, j ), j = 1, 4 ) /  189, 2096, 2159, 2841 /,
     $   ( mm( 59, j ), j = 1, 4 ) /  287, 1761, 2318, 4077 /,
     $   ( mm( 60, j ), j = 1, 4 ) / 2045, 2810, 2091,  721 /,
     $   ( mm( 61, j ), j = 1, 4 ) / 1227,  566, 3443, 2821 /,
     $   ( mm( 62, j ), j = 1, 4 ) / 2838,  442, 1510, 2249 /,
     $   ( mm( 63, j ), j = 1, 4 ) /  209,   41,  449, 2397 /,
     $   ( mm( 64, j ), j = 1, 4 ) / 2770, 1238, 1956, 2817 /,
     $   ( mm( 65, j ), j = 1, 4 ) / 3654, 1086, 2201,  245 /,
     $   ( mm( 66, j ), j = 1, 4 ) / 3993,  603, 3137, 1913 /,
     $   ( mm( 67, j ), j = 1, 4 ) /  192,  840, 3399, 1997 /,
     $   ( mm( 68, j ), j = 1, 4 ) / 2253, 3168, 1321, 3121 /,
     $   ( mm( 69, j ), j = 1, 4 ) / 3491, 1499, 2271,  997 /,
     $   ( mm( 70, j ), j = 1, 4 ) / 2889, 1084, 3667, 1833 /,
     $   ( mm( 71, j ), j = 1, 4 ) / 2857, 3438, 2703, 2877 /,
     $   ( mm( 72, j ), j = 1, 4 ) / 2094, 2408,  629, 1633 /,
     $   ( mm( 73, j ), j = 1, 4 ) / 1818, 1589, 2365,  981 /,
     $   ( mm( 74, j ), j = 1, 4 ) /  688, 2391, 2431, 2009 /,
     $   ( mm( 75, j ), j = 1, 4 ) / 1407,  288, 1113,  941 /,
     $   ( mm( 76, j ), j = 1, 4 ) /  634,   26, 3922, 2449 /,
     $   ( mm( 77, j ), j = 1, 4 ) / 3231,  512, 2554,  197 /,
     $   ( mm( 78, j ), j = 1, 4 ) /  815, 1456,  184, 2441 /,
     $   ( mm( 79, j ), j = 1, 4 ) / 3524,  171, 2099,  285 /,
     $   ( mm( 80, j ), j = 1, 4 ) / 1914, 1677, 3228, 1473 /,
     $   ( mm( 81, j ), j = 1, 4 ) /  516, 2657, 4012, 2741 /,
     $   ( mm( 82, j ), j = 1, 4 ) /  164, 2270, 1921, 3129 /,
     $   ( mm( 83, j ), j = 1, 4 ) /  303, 2587, 3452,  909 /,
     $   ( mm( 84, j ), j = 1, 4 ) / 2144, 2961, 3901, 2801 /,
     $   ( mm( 85, j ), j = 1, 4 ) / 3480, 1970,  572,  421 /,
     $   ( mm( 86, j ), j = 1, 4 ) /  119, 1817, 3309, 4073 /,
     $   ( mm( 87, j ), j = 1, 4 ) / 3357,  676, 3171, 2813 /,
     $   ( mm( 88, j ), j = 1, 4 ) /  837, 1410,  817, 2337 /,
     $   ( mm( 89, j ), j = 1, 4 ) / 2826, 3723, 3039, 1429 /,
     $   ( mm( 90, j ), j = 1, 4 ) / 2332, 2803, 1696, 1177 /,
     $   ( mm( 91, j ), j = 1, 4 ) / 2089, 3185, 1256, 1901 /,
     $   ( mm( 92, j ), j = 1, 4 ) / 3780,  184, 3715,   81 /,
     $   ( mm( 93, j ), j = 1, 4 ) / 1700,  663, 2077, 1669 /,
     $   ( mm( 94, j ), j = 1, 4 ) / 3712,  499, 3019, 2633 /,
     $   ( mm( 95, j ), j = 1, 4 ) /  150, 3784, 1497, 2269 /,
     $   ( mm( 96, j ), j = 1, 4 ) / 2000, 1631, 1101,  129 /,
     $   ( mm( 97, j ), j = 1, 4 ) / 3375, 1925,  717, 1141 /,
     $   ( mm( 98, j ), j = 1, 4 ) / 1621, 3912,   51,  249 /,
     $   ( mm( 99, j ), j = 1, 4 ) / 3090, 1398,  981, 3917 /,
     $   ( mm(100, j ), j = 1, 4 ) / 3765, 1349, 1978, 2481 /,
     $   ( mm(101, j ), j = 1, 4 ) / 1149, 1441, 1813, 3941 /,
     $   ( mm(102, j ), j = 1, 4 ) / 3146, 2224, 3881, 2217 /,
     $   ( mm(103, j ), j = 1, 4 ) /   33, 2411,   76, 2749 /,
     $   ( mm(104, j ), j = 1, 4 ) / 3082, 1907, 3846, 3041 /,
     $   ( mm(105, j ), j = 1, 4 ) / 2741, 3192, 3694, 1877 /,
     $   ( mm(106, j ), j = 1, 4 ) /  359, 2786, 1682,  345 /,
     $   ( mm(107, j ), j = 1, 4 ) / 3316,  382,  124, 2861 /,
     $   ( mm(108, j ), j = 1, 4 ) / 1749,   37, 1660, 1809 /,
     $   ( mm(109, j ), j = 1, 4 ) /  185,  759, 3997, 3141 /,
     $   ( mm(110, j ), j = 1, 4 ) / 2784, 2948,  479, 2825 /,
     $   ( mm(111, j ), j = 1, 4 ) / 2202, 1862, 1141,  157 /,
     $   ( mm(112, j ), j = 1, 4 ) / 2199, 3802,  886, 2881 /,
     $   ( mm(113, j ), j = 1, 4 ) / 1364, 2423, 3514, 3637 /,
     $   ( mm(114, j ), j = 1, 4 ) / 1244, 2051, 1301, 1465 /,
     $   ( mm(115, j ), j = 1, 4 ) / 2020, 2295, 3604, 2829 /,
     $   ( mm(116, j ), j = 1, 4 ) / 3160, 1332, 1888, 2161 /,
     $   ( mm(117, j ), j = 1, 4 ) / 2785, 1832, 1836, 3365 /,
     $   ( mm(118, j ), j = 1, 4 ) / 2772, 2405, 1990,  361 /,
     $   ( mm(119, j ), j = 1, 4 ) / 1217, 3638, 2058, 2685 /,
     $   ( mm(120, j ), j = 1, 4 ) / 1822, 3661,  692, 3745 /,
     $   ( mm(121, j ), j = 1, 4 ) / 1245,  327, 1194, 2325 /,
     $   ( mm(122, j ), j = 1, 4 ) / 2252, 3660,   20, 3609 /,
     $   ( mm(123, j ), j = 1, 4 ) / 3904,  716, 3285, 3821 /,
     $   ( mm(124, j ), j = 1, 4 ) / 2774, 1842, 2046, 3537 /,
     $   ( mm(125, j ), j = 1, 4 ) /  997, 3987, 2107,  517 /,
     $   ( mm(126, j ), j = 1, 4 ) / 2573, 1368, 3508, 3017 /,
     $   ( mm(127, j ), j = 1, 4 ) / 1148, 1848, 3525, 2141 /,
     $   ( mm(128, j ), j = 1, 4 ) /  545, 2366, 3801, 1537 /
!     
!     EXECUTABLE STATEMENTS:
!
      i1 = iseed(1)
      i2 = iseed(2)
      i3 = iseed(3)
      i4 = iseed(4)
!
      do i = 1, min(n,lv)
!
!        multiply the seed by i-th power of the multiplier modulo 2**48
!
         it4 = i4*mm(i,4)
         it3 = it4 / ipw2
         it4 = it4 - ipw2*it3
         it3 = it3 + i3*mm(i,4) + i4*mm(i,3)
         it2 = it3 / ipw2
         it3 = it3 - ipw2*it2
         it2 = it2 + i2*mm(i,4) + i3*mm(i,3) + i4*mm(i,2)
         it1 = it2 / ipw2
         it2 = it2 - ipw2*it1
         it1 = it1 + i1*mm(i,4) + i2*mm(i,3) + i3*mm(i,2) +
     $               i4*mm(i,1)
         it1 = mod( it1, ipw2 )
!
!        convert 48-bit integer to a real number in the interval (0,1)
!
         x(i) = r*( real(it1, kind=wrp) +
     $               r*( real(it2, kind=wrp) +
     $                  r*( real(it3, kind=wrp) +
     $                     r*real(it4, kind=wrp) ) ) )
!
      enddo 
!
!     return final value of seed
!
      iseed(1) = it1
      iseed(2) = it2
      iseed(3) = it3
      iseed(4) = it4
      return
      end subroutine qplaruv
!
!***********************************************************************
!
      subroutine qpgeqr2 ( 
     $           m, n, a, lda, tau, work, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)
     $   m, n, lda 
      real(wrp)  
     $   a(lda,*), tau(*), work(*)
      integer(wip)
     $   info
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  DGEQR2 computes a QR factorization of a real m by n matrix A:
!  A = Q * R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(m,n) by n upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I-tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!  =====================================================================
!     
!     DEPENDENCES: 
      external
     $   qplarf, 
     $   qplarfg, 
     $   wxerbla
      intrinsic
     $   max, min
!
!     LOCAL VARIABLES: 
      real(wrp) 
     $   one
      parameter ( 
     $   one = 1.0_wrp )
      integer(wip) 
     $   i, k
      real(wrp) 
     $   aii
!    
!     EXECUTABLE STATEMENTS:
!     Test the input arguments
!
      info = 0
      if ( m .lt. 0 ) then
         info = -1
      else if ( n .lt. 0 ) then
         info = -2
      else if ( lda .lt. max(1,m) ) then
         info = -4
      endif
      if ( info .ne. 0 ) then
         call wxerbla ( 
     $        'wfgeqr2', -info )
         return
      endif
!
      k = min(m,n)
!
      do i = 1, k
!
!        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
         call qplarfg (
     $        m-i+1, a(i,i), a(min(i+1,m),i), 1_wip,
     $        tau(i) )
         if ( i .lt. n ) then
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
            aii = a(i,i)
            a(i,i) = one
            call qplarf (
     $           'Left', m-i+1, n-i, a(i,i), 1_wip, tau(i),
     $           a(i,i+1), lda, work )
            a(i,i) = aii
         endif
      enddo 
      return
      end subroutine qpgeqr2 
!
!***********************************************************************
!
      subroutine qplarf ( 
     $           side, m, n, v, incv, tau, c, ldc, work )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character          
     $   side
      integer(wip)       
     $   m, n, incv, ldc
      real(wrp)          
     $   tau
      real(wrp) 
     $   c(ldc,*), v(*), work(*)
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I-tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) REAL(wrp)        array, dimension
!                     (1+(M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1+(N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) REAL(wrp)       
!          The value tau in the representation of H.
!
!  C       (input/output) REAL(wrp)        array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) REAL(wrp)        array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     DEPENDENCES:
      external
     $   qpgemv, 
     $   qpger
      logical 
     $   wlsame
      external
     $   wlsame
!
!     LOCAL VARIABLES:
      real(wrp)  
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!
!     EXECUTABLE STATEMENTS:
!
      if ( wlsame ( side, 'L' ) ) then
!
!        Form  H * C
!
         if ( tau.ne.zero ) then
!
!           w := C' * v
!
            call qpgemv( 
     $           'Transpose', 
     $           m, n, one, c, ldc, v, incv, zero, work, 1_wip )
!
!           C := C-v * w'
!
            call qpger ( 
     $           m, n, -tau, v, incv, work, 1_wip, c, ldc )
         endif
      else
!
!        Form  C * H
!
         if ( tau.ne.zero ) then
!
!           w := C * v
!
            call qpgemv ( 
     $           'No transpose', 
     $           m, n, one, c, ldc, v, incv, zero, work, 1_wip )
!
!           C := C-w * v'
!
            call qpger ( 
     $           m, n, -tau, work, 1_wip, v, incv, c, ldc )
         endif
      endif
      return
      end subroutine qplarf
!
!***********************************************************************
!
      subroutine qplarfg ( 
     $           n, alpha, x, incx, tau )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer(wip)       
     $   n, incx
      real(wrp)   
     $   alpha, tau, 
     $   x(*)
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I-tau * (1) * ( 1 v' ) ,
!                      (v)
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) REAL(wrp)       
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) REAL(wrp)        array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) REAL(wrp)       
!          The value tau.
!
!  =====================================================================
!
!     DEPENDENCES:
      real(wrp)  
     $   qplamch, 
     $   qplapy2, 
     $   qpnrm2
      external 
     $   qplamch, 
     $   qplapy2, 
     $   qpnrm2
      intrinsic  
     $   abs, sign
      external 
     $   qpscal
!
!     LOCAL VARIABLES:
      real(wrp) 
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      integer(wip)  
     $   j, knt
      real(wrp)   
     $   beta, rsafmn, safmin, xnorm
!
!     EXECUTABLE STATEMENTS:
!
      if ( n.le.1 ) then
         tau = zero
         return
      endif
!
      xnorm = qpnrm2 ( 
     $        n-1, x, incx )
!
      if ( xnorm.eq.zero ) then
!
!        H  =  I
!
         tau = zero
      else
!
!        general case
!
         beta = -sign( qplapy2 ( 
     $                 alpha, xnorm ), alpha )
         safmin = qplamch ('S') / qplamch ('E')
         if ( abs( beta ).lt.safmin ) then
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            rsafmn = one / safmin
            knt = 0
   10       continue
            knt = knt+1
            call qpscal ( 
     $           n-1, rsafmn, x, incx )
            beta  = beta *rsafmn
            alpha = alpha*rsafmn
            if ( abs( beta ).lt.safmin )
     $         goto 10
!
!           New BETA is at most 1, at least SAFMIN
!
            xnorm = qpnrm2 ( 
     $              n-1, x, incx )
            beta = -sign( qplapy2 ( 
     $                    alpha, xnorm ), alpha )
            tau = ( beta-alpha ) / beta
            call qpscal ( 
     $           n-1, one / ( alpha-beta ), x, incx )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
            alpha = beta
            do j = 1, knt
               alpha = alpha*safmin
            enddo
         else
            tau = ( beta-alpha ) / beta
            call qpscal ( 
     $           n-1, one / ( alpha-beta ), x, incx )
            alpha = beta
         endif
      endif
!
      return
      end subroutine qplarfg 
!
!***********************************************************************
!
      subroutine qporm2r ( 
     $           side, trans, m, n, k, a, lda, tau, c, ldc,
     $           work, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      character
     $   side, trans
      integer(wip)
     $   m, n, k, lda, ldc 
      real(wrp) 
     $   a(lda,*), c(ldc,*), tau(*), work(*)
      integer(wip)
     $   info
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  DORM2R overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!     
!     DEPENDENCES: 
      logical
     $   wlsame
      external 
     $   wlsame
      external
     $   qplarf, 
     $   wxerbla
      intrinsic
     $   max
!
!     LOCAL VARIABLES: 
      real(wrp)  
     $   one
      parameter ( 
     $   one = 1.0_wrp )
      logical 
     $   left, notran
      integer(wip)    
     $   i, i1, i2, i3, ic, jc, mi, ni, nq
      real(wrp) 
     $   aii
!     
!     EXECUTABLE STATEMENTS:
!     Test the input arguments
!
      info   = 0
      left   = wlsame ( side, 'L' )
      notran = wlsame ( trans, 'N' )
!
!     NQ is the order of Q
!
      if ( left ) then
         nq = m
      else
         nq = n
      endif
      if ( .not.left .and. .not.wlsame ( side, 'R' ) ) then
         info = -1
      else if ( .not.notran .and. .not.wlsame ( trans, 'T' ) ) then
         info = -2
      else if ( m .lt. 0 ) then
         info = -3
      else if ( n .lt. 0 ) then
         info = -4
      else if ( k .lt. 0 .or. k.gt.nq ) then
         info = -5
      else if ( lda .lt. max(1,nq) ) then
         info = -7
      else if ( ldc .lt. max(1,m) ) then
         info = -10
      endif
      if ( info .ne. 0 ) then
         call wxerbla (
     $        'wform2r', -info )
         return
      endif
!
!     Quick return if possible
!
      if ( m.eq.0 .or. n.eq.0 .or. k.eq.0 ) return
!
      if ( ( left .and. .not.notran ) .or. 
     $     ( .not.left .and. notran )    )  then
         i1 = 1
         i2 = k
         i3 = 1
      else
         i1 = k
         i2 = 1
         i3 = -1
      endif
!
      if ( left ) then
         ni = n
         jc = 1
      else
         mi = m
         ic = 1
      endif
!
      do i = i1, i2, i3
         if ( left ) then
!
!           H(i) is applied to C(i:m,1:n)
!
            mi = m-i+1
            ic = i
         else
!
!           H(i) is applied to C(1:m,i:n)
!
            ni = n-i+1
            jc = i
         endif
!
!        Apply H(i)
!
         aii = a(i,i)
         a(i,i) = one
         call qplarf (
     $        side, mi, ni, a(i,i), 1_wip, tau(i), c(ic,jc),
     $        ldc, work )
         a(i,i) = aii
      enddo 
      return
      end subroutine qporm2r
!
!***********************************************************************
!-----------------------------------------------------------------------
!***********************************************************************
!
!@@@  BLAS::
! 
!-----------------------------------------------------------------------
!
      subroutine qpaxpy ( 
     $           n, da, dx, incx, dy, incy )
!      
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)  
     $   n, incx, incy
      real(wrp)  
     $   da, dx(*), dy(*)
!
!     DAXPY constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     Jack Dongarra, linpack, 3/11/78.
!     +  modified 12/3/93, array(1) declarations changed to array(*)
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     LOCAL VARIABLES: 
      integer(wip)  
     $   i, ix, iy, m, mp1
!
!     EXECUTABLE STATEMENTS:
!
      if ( n .le. 0_wip ) return
      if ( da .eq. 0.0_wrp ) return
      if ( incx .eq. 1_wip .and. incy .eq. 1_wip ) goto 20
!
!        code for unequal increments or equal increments
!        not equal to 1
!
      ix = 1_wip
      iy = 1_wip
      if ( incx .lt. 0_wip ) ix = (-n+1)*incx + 1
      if ( incy .lt. 0_wip ) iy = (-n+1)*incy + 1
      do i = 1,n
         dy(iy) = dy(iy) + da*dx(ix)
         ix = ix + incx
         iy = iy + incy
      enddo 
      return
!
!        code for both increments equal to 1
!        clean-up loop
!
   20 continue 
      m = mod( n, 4_wip )
      if ( m .eq. 0_wip ) goto 40
      do i = 1,m
         dy(i) = dy(i) + da*dx(i)
      enddo 
      if ( n .lt. 4_wip ) return
   40 continue 
      mp1 = m + 1
      do i = mp1, n, 4
         dy(i  ) = dy(i  ) + da*dx(i  )
         dy(i+1) = dy(i+1) + da*dx(i+1)
         dy(i+2) = dy(i+2) + da*dx(i+2)
         dy(i+3) = dy(i+3) + da*dx(i+3)
      enddo 
      return
      end subroutine qpaxpy  
! 
!***********************************************************************
!
      subroutine qpcopy (
     $           n, dx, incx, dy, incy )
!      
      implicit none 
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer(wip)  
     $   n, incx, incy
      real(wrp)     
     $   dx(*), dy(*)
!
!     DCOPY copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to one.
!     Jack Dongarra, linpack, 3/11/78.
!     +  modified 12/3/93, array(1) declarations changed to array(*)
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     LOCAL VARIABLES:
      integer(wip)  
     $   i, ix, iy, m, mp1
!
!     EXECUTABLE STATEMENTS:
!
      if ( n .le. 0_wip ) return
      if ( incx .eq. 1_wip .and. incy .eq. 1_wip ) goto 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if ( incx .lt. 0_wip ) ix = (-n+1)*incx + 1
      if ( incy .lt. 0_wip ) iy = (-n+1)*incy + 1
      do i = 1,n
         dy(iy) = dx(ix)
         ix = ix + incx
         iy = iy + incy
      enddo 
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 continue 
      m = mod( n, 7_wip )
      if ( m .eq. 0_wip ) goto 40
      do i = 1,m
         dy(i) = dx(i)
      enddo 
      if ( n .lt. 7_wip ) return
   40 continue 
      mp1 = m + 1
      do i = mp1, n, 7
         dy(i  ) = dx(i  )
         dy(i+1) = dx(i+1)
         dy(i+2) = dx(i+2)
         dy(i+3) = dx(i+3)
         dy(i+4) = dx(i+4)
         dy(i+5) = dx(i+5)
         dy(i+6) = dx(i+6)
      enddo 
      return
      end subroutine qpcopy 
!
!***********************************************************************
!
      function qpdot ( 
     $         n, dx, incx, dy, incy )
!
      implicit none
      include 'epcode_inc_qp.h'
!
      real(wrp) 
     $   qpdot
!
!     ARGUMENTS:
      real(wrp) 
     $   dx(*), dy(*)
      integer(wip) 
     $   n, incx, incy
!
!     DDOT forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     Jack Dongarra, linpack, 3/11/78.
!     +  modified 12/3/93, array(1) declarations changed to array(*)
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     LOCAL VARIABLES:
      real(wrp) 
     $   dtemp
      integer(wip) 
     $   i, ix, iy, m, mp1
!
!     EXECUTABLE STATEMENTS:
!
      qpdot = 0.0_wrp 
      dtemp = 0.0_wrp 
      if ( n.le.0_wip ) return
      if ( incx.eq.1_wip .and. incy.eq.1_wip ) goto 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if ( incx.lt.0_wip ) ix = (-n+1)*incx + 1
      if ( incy.lt.0_wip ) iy = (-n+1)*incy + 1
      do i = 1,n
         dtemp = dtemp + dx(ix)*dy(iy)
         ix = ix + incx
         iy = iy + incy
      enddo 
      qpdot = dtemp
      return
!
!        code for both increments equal to 1
!
!        clean-up loop
!
   20 continue 
      m = mod( n, 5_wip )
      if ( m .eq. 0_wip ) goto 40
      do i = 1,m
         dtemp = dtemp + dx(i)*dy(i)
      enddo 
      if ( n .lt. 5_wip ) goto 60
   40 continue 
      mp1 = m + 1
      do i = mp1, n, 5
         dtemp = dtemp + 
     $           dx(i  )*dy(i  ) + 
     $           dx(i+1)*dy(i+1) +
     $           dx(i+2)*dy(i+2) + 
     $           dx(i+3)*dy(i+3) + 
     $           dx(i+4)*dy(i+4)
      enddo 
   60 continue 
      qpdot = dtemp
      return
      end function qpdot
!
!***********************************************************************
!
      function qpnrm2 ( 
     $         n, x, incx )
!
      implicit none
      include 'epcode_inc_qp.h'
!
      real(wrp)   
     $   qpnrm2 
!
!     ARGUMENTS:
      integer(wip)  
     $   incx, n
      real(wrp) 
     $   x(*)
!
!     DNRM2 returns the euclidean norm of a vector via the function
!     name, so that
!        DNRM2 := sqrt( x'*x )
!  -- This version qpitten on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to WFLASSQ.
!     Sven Hammarling, Nag Ltd.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      intrinsic 
     $   abs, sqrt
!
!     LOCAL VARIABLES: 
      real(wrp)
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      integer(wip)
     $   ix
      real(wrp) 
     $   absxi, norm, rscale, ssq
!     
!     EXECUTABLE STATEMENTS:
!
      if ( n.lt.1_wip .or. incx.lt.1_wip ) then
         norm  = zero
      else if ( n.eq.1_wip ) then
         norm  = abs ( x(1) )
      else
         rscale = zero
         ssq    = one
!
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL WFLASSQ( N, X, INCX, SCALE, SSQ )
!
         do ix = 1, 1+(n-1)*incx, incx
            if ( x(ix) .ne. zero ) then
               absxi = abs ( x(ix) )
               if ( rscale .lt. absxi ) then
                  ssq    = one + ssq*( rscale/absxi )**2
                  rscale = absxi
               else
                  ssq    = ssq +     ( absxi/rscale )**2
               endif
            endif
         enddo 
         norm  = rscale * sqrt ( ssq )
      endif
      qpnrm2 = norm
      return
      end function qpnrm2 
!
!***********************************************************************
!
      subroutine qpscal ( 
     $           n, da, dx, incx )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      integer(wip) 
     $   n, incx
      real(wrp) 
     $   da, dx(*)
!
!     DSCAL scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     Jack Dongarra, linpack, 3/11/78.
!     +  modified 3/93 to return if incx .le. 0.
!     +  modified 12/3/93, array(1) declarations changed to array(*)
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     LOCAL VARIABLES:
      integer(wip) 
     $   i, m, mp1, nincx
!
!     EXECUTABLE STATEMENTS:
!
      if ( n .le. 0_wip .or. incx .le. 0_wip ) return
      if ( incx .eq. 1_wip ) goto 20
!
!     code for increment not equal to 1
!
      nincx = n*incx
      do i = 1, nincx, incx
         dx(i) = da*dx(i)
      enddo 
      return
!
!     code for increment equal to 1
!     clean-up loop
!
   20 continue 
      m = mod( n, 5_wip )
      if ( m .eq. 0_wip ) goto 40
      do i = 1,m
         dx(i) = da*dx(i)
      enddo 
      if ( n .lt. 5_wip ) return
   40 continue 
      mp1 = m + 1
      do i = mp1, n, 5
         dx(i  ) = da*dx(i  )
         dx(i+1) = da*dx(i+1)
         dx(i+2) = da*dx(i+2)
         dx(i+3) = da*dx(i+3)
         dx(i+4) = da*dx(i+4)
      enddo 
      return
      end subroutine qpscal 
!
!***********************************************************************
!
      subroutine qpswap (
     $           n, dx, incx, dy, incy )
!      
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      integer(wip) 
     $   n, incx, incy
      real(wrp) 
     $   dx(*), dy(*)
!
!     DSWAP interchanges two vectors.
!     uses unrolled loops for increments equal one.
!     Jack Dongarra, linpack, 3/11/78.
!     +  modified 12/3/93, array(1) declarations changed to array(*)
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     LOCAL VARIABLES:
      real(wrp) 
     $   dtemp
      integer(wip) 
     $   i, ix, iy, m, mp1
!
!     EXECUTABLE STATEMENTS:
!
      if ( n.le.0_wip ) return
      if ( incx.eq.1_wip .and. incy.eq.1_wip ) goto 20
!
!     code for unequal increments or equal increments
!     not equal to 1
!
      ix = 1
      iy = 1
      if ( incx.lt.0_wip ) ix = (-n+1)*incx + 1
      if ( incy.lt.0_wip ) iy = (-n+1)*incy + 1
      do i = 1,n
         dtemp  = dx(ix)
         dx(ix) = dy(iy)
         dy(iy) = dtemp
         ix = ix + incx
         iy = iy + incy
      enddo 
      return
!
!     code for both increments equal to 1
!     clean-up loop
!
   20 continue 
      m = mod ( n, 3_wip )
      if ( m .eq. 0_wip ) goto 40
      do i = 1,m
         dtemp = dx(i)
         dx(i) = dy(i)
         dy(i) = dtemp
      enddo 
      if ( n .lt. 3_wip ) return
   40 continue 
      mp1 = m + 1
      do i = mp1, n, 3
         dtemp   = dx(i)
         dx(i)   = dy(i)
         dy(i)   = dtemp
         dtemp   = dx(i+1)
         dx(i+1) = dy(i+1)
         dy(i+1) = dtemp
         dtemp   = dx(i+2)
         dx(i+2) = dy(i+2)
         dy(i+2) = dtemp
      enddo 
      return
      end subroutine qpswap
!
!***********************************************************************
!
!     WFGEMV is called many times by other routines 
!
      subroutine qpgemv ( 
     $           trans, m, n, alpha, a, lda, x, incx,
     $           beta, y, incy )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      character(len=1) 
     $   trans
      integer(wip) 
     $   m, n, lda, incx, incy 
      real(wrp)   
     $   alpha, beta,
     $   a(lda,*), x(*), y(*)
!
!  Purpose
!  =======
!
!  WFGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(wrp)       .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(wrp)        array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - REAL(wrp)        array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - REAL(wrp)       .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL(wrp)        array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      logical 
     $   wlsame
      external
     $   wlsame
      external 
     $   wxerbla
      intrinsic
     $   max
!
!     LOCAL VARIABLES:   
      real(wrp)  
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      real(wrp)
     $   temp
      integer(wip) 
     $   i, ix, iy, j, jx, jy, 
     $   kx, ky, lenx, leny, 
     $   info
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      info = 0
      if ( .not. wlsame ( trans, 'N' ) .and.
     $     .not. wlsame ( trans, 'T' ) .and.
     $     .not. wlsame ( trans, 'C' )  ) then
         info = 1
      else if ( m .lt. 0_wip ) then
         info = 2
      else if ( n .lt. 0_wip ) then
         info = 3
      else if ( lda .lt. max( 1_wip, m ) ) then
         info = 6
      else if ( incx .eq. 0_wip ) then
         info = 8
      else if ( incy .eq. 0_wip ) then
         info = 11
      endif
      if ( info .ne. 0 ) then
         call wxerbla ( 
     $        'wfgemv ', info )
         return
      endif
!
!     Quick return if possible.
!
      if ( ( m .eq. 0_wip ) .or. ( n .eq. 0_wip ) .or.
     $    ( ( alpha .eq. zero ) .and. ( beta .eq. one ) ) )
     $   return
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      if ( wlsame ( trans, 'N' ) ) then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      endif
      if ( incx .gt. 0_wip ) then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      endif
      if ( incy .gt. 0_wip ) then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      endif
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      if ( beta.ne.one ) then
         if ( incy.eq.1_wip ) then
            if ( beta.eq.zero ) then
               do i = 1, leny
                  y(i) = zero
               enddo 
            else
               do i = 1, leny
                  y(i) = beta*y(i)
               enddo 
            endif
         else
            iy = ky
            if ( beta .eq. zero ) then
               do i = 1, leny
                  y(iy) = zero
                  iy    = iy + incy
               enddo 
            else
               do i = 1, leny
                  y(iy) = beta*y(iy)
                  iy    = iy + incy
               enddo 
            endif
         endif
      endif
      if ( alpha.eq.zero ) return
      if ( wlsame ( trans, 'N' ) ) then
!
!        Form  y := alpha*A*x + y.
!
         jx = kx
         if ( incy.eq.1_wip ) then
            do j = 1, n
               if ( x(jx) .ne. zero ) then
                  temp = alpha*x(jx)
                  do i = 1, m
                     y(i) = y(i) + temp*a(i,j)
                  enddo 
               endif
               jx = jx + incx
            enddo 
         else
            do j = 1, n
               if ( x(jx) .ne. zero ) then
                  temp = alpha*x(jx)
                  iy   = ky
                  do i = 1, m
                     y(iy) = y(iy) + temp*a(i,j)
                     iy    = iy + incy
                  enddo 
               endif
               jx = jx + incx
            enddo 
         endif
      else
!
!        Form  y := alpha*A'*x + y.
!
         jy = ky
         if ( incx.eq.1_wip ) then
            do j = 1, n
               temp = zero
               do i = 1, m
                  temp = temp + a(i,j)*x(i)
               enddo 
               y(jy) = y(jy) + alpha*temp
               jy    = jy + incy
            enddo 
         else
            do j = 1, n
               temp = zero
               ix   = kx
               do i = 1, m
                  temp = temp + a(i,j)*x(ix)
                  ix   = ix + incx
               enddo 
               y(jy) = y(jy) + alpha*temp
               jy    = jy + incy
            enddo 
         endif
      endif
      return
      end subroutine qpgemv 
!
!***********************************************************************
!
      subroutine qpger ( 
     $           m, n, alpha, x, incx, y, incy, a, lda )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      real(wrp)  
     $   alpha
      integer(wip)  
     $   m, n, incx, incy, lda
      real(wrp) 
     $   a(lda,*), x(*), y(*)
!
!  Purpose
!  =======
!
!  WFGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(wrp)       .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL(wrp)        array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL(wrp)        array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL(wrp)        array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      external  
     $   wxerbla
      intrinsic  
     $   max
!
!     LOCAL VARIABLES:
      real(wrp)   
     $   zero
      parameter ( 
     $   zero = 0.0_wrp )
      real(wrp) 
     $   temp
      integer(wip) 
     $   i, ix, j, jy, kx, 
     $   info
! 
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      info = 0
      if      ( m.lt.0_wip ) then
         info = 1
      else if ( n.lt.0_wip ) then
         info = 2
      else if ( incx.eq.0_wip ) then
         info = 5
      else if ( incy.eq.0_wip ) then
         info = 7
      else if ( lda.lt.max( 1_wip, m ) ) then
         info = 9
      endif
      if ( info.ne.0 ) then
         call wxerbla (
     $        'wfger  ', info )
         return
      endif
!
!     Quick return if possible.
!
      if ( ( m.eq.0_wip ).or.( n.eq.0_wip ).or.( alpha.eq.zero ) )
     $   return
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      if ( incy.gt.0_wip ) then
         jy = 1
      else
         jy = 1 - ( n - 1 )*incy
      endif
      if ( incx.eq.1_wip ) then
         do j = 1, n
            if ( y(jy) .ne. zero ) then
               temp = alpha*y(jy)
               do i = 1, m
                  a(i,j) = a(i,j) + x(i)*temp
               enddo 
            endif
            jy = jy + incy
         enddo 
      else
         if ( incx.gt.0_wip ) then
            kx = 1_wip
         else
            kx = 1 - ( m - 1 )*incx
         endif
         do j = 1, n
            if ( y(jy) .ne. zero ) then
               temp = alpha*y(jy)
               ix   = kx
               do i = 1, m
                  a(i,j) = a(i,j) + x(ix)*temp
                  ix     = ix + incx
               enddo 
            endif
            jy = jy + incy
         enddo 
      endif
      return
      end subroutine qpger 
!
!***********************************************************************
!
      subroutine qpsyr2 ( 
     $           uplo, n, alpha, x, incx, y, incy, a, lda )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      character(len=1) 
     $   uplo  
      real(wrp) 
     $   alpha
      integer(wip) 
     $   incx, incy, lda, n
      real(wrp)  
     $   a(lda,*), x(*), y(*)
!
!  Purpose
!  =======
!
!  WFSYR2  performs the symmetric rank 2 operation
!
!     A := alpha*x*y' + alpha*y*x' + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an n
!  by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(wrp)       .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - REAL(wrp)        array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - REAL(wrp)        array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - REAL(wrp)        array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      logical  
     $   wlsame
      external 
     $   wlsame
      external  
     $   wxerbla
      intrinsic
     $   max
!
!     LOCAL VARIABLES:
      real(wrp)  
     $   zero
      parameter ( 
     $   zero = 0.0_wrp )
      real(wrp) 
     $   temp1, temp2
      integer(wip) 
     $   i, ix, iy, j, jx, jy, kx, ky,
     $   info
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      info = 0
      if ( .not. wlsame ( uplo, 'U' ) .and.
     $     .not. wlsame ( uplo, 'L' )      ) then
         info = 1
      else if ( n.lt.0_wip ) then
         info = 2
      else if ( incx.eq.0_wip ) then
         info = 5
      else if ( incy.eq.0_wip ) then
         info = 7
      else if ( lda.lt.max( 1_wip, n ) ) then
         info = 9
      endif
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wfsyr2 ', info )
         return
      endif
!
!     Quick return if possible.
!
      if ( ( n.eq.0_wip ).or.( alpha.eq.zero ) ) return
!
!     Set up the start points in X and Y if the increments are not both
!     unity.
!
      if ( ( incx.ne.1_wip ).or.( incy.ne.1_wip ) ) then
         if ( incx.gt.0_wip ) then
            kx = 1_wip
         else
            kx = 1 - ( n - 1 )*incx
         endif
         if ( incy.gt.0_wip ) then
            ky = 1_wip
         else
            ky = 1 - ( n - 1 )*incy
         endif
         jx = kx
         jy = ky
      endif
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      if ( wlsame ( uplo, 'U' ) ) then
!
!        Form  A  when A is stored in the upper triangle.
!
         if ( ( incx.eq.1_wip ).and.( incy.eq.1_wip ) ) then
            do j = 1, n
               if ( ( x(j) .ne. zero ) .or. 
     $              ( y(j) .ne. zero )    ) then
                  temp1 = alpha*y(j)
                  temp2 = alpha*x(j)
                  do i = 1, j
                     a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
                  enddo 
               endif
            enddo 
         else
            do j = 1, n
               if ( ( x(jx) .ne. zero ) .or. 
     $              ( y(jy) .ne. zero )    ) then
                  temp1 = alpha*y(jy)
                  temp2 = alpha*x(jx)
                  ix    = kx
                  iy    = ky
                  do i = 1, j
                     a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                     ix     = ix + incx
                     iy     = iy + incy
                  enddo 
               endif
               jx = jx + incx
               jy = jy + incy
            enddo 
         endif
      else
!
!        Form  A  when A is stored in the lower triangle.
!
         if ( ( incx.eq.1_wip ).and.( incy.eq.1_wip ) ) then
            do j = 1, n
               if ( ( x(j) .ne. zero ) .or.
     $              ( y(j) .ne. zero )    ) then
                  temp1 = alpha*y(j)
                  temp2 = alpha*x(j)
                  do i = j, n
                     a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
                  enddo 
               endif
            enddo 
         else
            do j = 1, n
               if ( ( x(jx) .ne. zero ).or.
     $              ( y(jy) .ne. zero )   ) then
                  temp1 = alpha*y(jy)
                  temp2 = alpha*x(jx)
                  ix    = jx
                  iy    = jy
                  do i = j, n
                     a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                     ix     = ix + incx
                     iy     = iy + incy
                  enddo 
               endif
               jx = jx + incx
               jy = jy + incy
            enddo 
         endif
      endif
      return
      end subroutine qpsyr2 
!
!***********************************************************************
!
      subroutine qpgemm ( 
     $           transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $           beta, c, ldc )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      character(len=1)  
     $   transa, transb
      integer(wip)
     $   m, n, k, lda, ldb, ldc
      real(wrp)
     $   alpha, beta
      real(wrp)
     $   a(lda,*), b(ldb,*), c(ldc,*)
!
!  Purpose
!  =======
!
!  WFGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(wrp)       .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(wrp)        array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL(wrp)        array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - REAL(wrp)       .
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - REAL(wrp)        array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      logical 
     $   wlsame
      external 
     $   wlsame
      external  
     $   wxerbla
      intrinsic 
     $   max
!
!     LOCAL VARIABLES:
      logical  
     $   nota, notb
      integer(wip) 
     $   i, j, l, ncola, nrowa, nrowb
      integer(wip)
     $   info
      real(wrp)  
     $   temp
      real(wrp)
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!     
!     EXECUTABLE STATEMENTS:
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
      nota = wlsame ( transa, 'N' )
      notb = wlsame ( transb, 'N' )
      if ( nota ) then
         nrowa = m
         ncola = k
      else
         nrowa = k
         ncola = m
      endif
      if ( notb ) then
         nrowb = k
      else
         nrowb = n
      endif
!
!     Test the input parameters.
!
      info = 0
      if ( ( .not. nota                   ) .and.
     $     ( .not. wlsame ( transa, 'C' ) ) .and.
     $     ( .not. wlsame ( transa, 'T' ) )      ) then
         info = 1
      else if (( .not. notb                   ) .and.
     $         ( .not. wlsame ( transb, 'C' ) ) .and.
     $         ( .not. wlsame ( transb, 'T' ) )      ) then
         info = 2
      else if ( m  .lt.0_wip               ) then
         info = 3
      else if ( n  .lt.0_wip               ) then
         info = 4 
      else if ( k  .lt.0_wip               ) then
         info = 5
      else if ( lda.lt.max( 1_wip, nrowa ) ) then
         info = 8
      else if ( ldb.lt.max( 1_wip, nrowb ) ) then
         info = 10
      else if ( ldc.lt.max( 1_wip, m     ) ) then
         info = 13
      endif
      if ( info.ne.0_wip ) then
         call wxerbla ( 
     $        'wfgemm ', info )
         return
      endif
!
!     Quick return if possible.
!
      if ( ( m.eq.0_wip ).or.( n.eq.0_wip ).or.
     $     ( ( ( alpha.eq.zero ).or.( k.eq.0_wip ) ).and.
     $       ( beta.eq.one ) ) )
     $   return
!
!     And if  alpha.eq.zero.
!
      if ( alpha.eq.zero ) then
         if ( beta.eq.zero ) then
            do j = 1, n
               do i = 1, m
                  c(i,j) = zero
               enddo 
            enddo 
         else
            do j = 1, n
               do i = 1, m
                  c(i,j) = beta*c(i,j)
               enddo 
            enddo 
         endif
         return
      endif
!
!     Start the operations.
!
      if ( notb ) then
         if ( nota ) then
!
!           Form  C := alpha*A*B + beta*C.
!
            do j = 1, n
               if ( beta .eq. zero ) then
                  do i = 1, m
                     c(i,j) = zero
                  enddo 
               else if ( beta .ne. one ) then
                  do i = 1, m
                     c(i,j) = beta*c(i,j)
                  enddo 
               endif
               do l = 1, k
                  if ( b(l,j) .ne. zero ) then
                     temp = alpha*b(l,j)
                     do i = 1, m
                        c(i,j) = c(i,j) + temp*a(i,l)
                     enddo 
                  endif
               enddo 
            enddo 
         else
!
!           Form  C := alpha*A'*B + beta*C
!
            do j = 1, n
               do i = 1, m
                  temp = zero
                  do l = 1, k
                     temp = temp + a(l,i)*b(l,j)
                  enddo 
                  if ( beta .eq. zero ) then
                     c(i,j) = alpha*temp
                  else
                     c(i,j) = alpha*temp + beta*c(i,j)
                  endif
               enddo 
            enddo 
         endif
      else
         if ( nota ) then
!
!           Form  C := alpha*A*B' + beta*C
!
            do j = 1, n
               if ( beta .eq. zero ) then
                  do i = 1, m
                     c(i,j) = zero
                  enddo 
               else if ( beta .ne. one ) then
                  do i = 1, m
                     c(i,j) = beta*c(i,j)
                  enddo 
               endif
               do l = 1, k
                  if ( b(j,l) .ne. zero ) then
                     temp = alpha*b(j,l)
                     do i = 1, m
                        c(i,j) = c(i,j) + temp*a(i,l)
                     enddo 
                  endif
               enddo 
            enddo 
         else
!
!           Form  C := alpha*A'*B' + beta*C
!
            do j = 1, n
               do i = 1, m
                  temp = zero
                  do l = 1, k
                     temp = temp + a(l,i)*b(j,l)
                  enddo 
                  if ( beta .eq. zero ) then
                     c(i,j) = alpha*temp
                  else
                     c(i,j) = alpha*temp + beta*c(i,j)
                  endif
               enddo 
            enddo 
         endif
      endif
      return
      end subroutine qpgemm  
!
!***********************************************************************
!
      subroutine qpsyr2k ( 
     $           uplo, trans, n, k, alpha, a, lda, b, ldb,
     $           beta, c, ldc )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      character(len=1) 
     $   uplo, trans
      integer(wip) 
     $   n, k, lda, ldb, ldc
      real(wrp) 
     $   alpha, beta,
     $   a(lda,*), b(ldb,*), c(ldc,*)
!     
!
!  Purpose
!  =======
!
!  WFSYR2K  performs one of the symmetric rank 2k operations
!
!     C := alpha*A*B' + alpha*B*A' + beta*C,
!
!  or
!
!     C := alpha*A'*B + alpha*B'*A + beta*C,
!
!  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
!  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
!  matrices in the second case.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of the  array  C  is to be  referenced  as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry,  TRANS  specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
!                                        beta*C.
!
!              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
!                                        beta*C.
!
!              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
!                                        beta*C.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N specifies the order of the matrix C.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!           of  columns  of the  matrices  A and B,  and on  entry  with
!           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
!           of rows of the matrices  A and B.  K must be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(wrp)       .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(wrp)        array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by n  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!           be at least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - REAL(wrp)        array of DIMENSION ( LDB, kb ), where kb is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  k by n  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDB must be at least  max( 1, n ), otherwise  LDB must
!           be at least  max( 1, k ).
!           Unchanged on exit.
!
!  BETA   - REAL(wrp)       .
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.
!
!  C      - REAL(wrp)        array of DIMENSION ( LDC, n ).
!           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!           upper triangular part of the array C must contain the upper
!           triangular part  of the  symmetric matrix  and the strictly
!           lower triangular part of C is not referenced.  On exit, the
!           upper triangular part of the array  C is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!           lower triangular part of the array C must contain the lower
!           triangular part  of the  symmetric matrix  and the strictly
!           upper triangular part of C is not referenced.  On exit, the
!           lower triangular part of the array  C is overwritten by the
!           lower triangular part of the updated matrix.
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, n ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      logical
     $   wlsame
      external
     $   wlsame
      external
     $   wxerbla
      intrinsic
     $   max
!
!     LOCAL VARIABLES:
      logical 
     $   upper
      integer(wip)  
     $   i, j, l, nrowa,
     $   info
      real(wrp) 
     $   temp1, temp2
      real(wrp) 
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!    
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      if ( wlsame ( trans, 'N' ) ) then
         nrowa = n
      else
         nrowa = k
      endif
      upper = wlsame ( uplo, 'U' )
!
      info = 0
      if ( ( .not. upper                 ) .and.
     $     ( .not. wlsame ( uplo, 'L' ) )     ) then
         info = 1
      else if (( .not. wlsame ( trans, 'N' ) ) .and.
     $         ( .not. wlsame ( trans, 'T' ) ) .and.
     $         ( .not. wlsame ( trans, 'C' ) )      ) then
         info = 2
      else if ( n  .lt.0_wip               ) then
         info = 3
      else if ( k  .lt.0_wip               ) then
         info = 4
      else if ( lda.lt.max( 1_wip, nrowa ) ) then
         info = 7
      else if ( ldb.lt.max( 1_wip, nrowa ) ) then
         info = 9
      else if ( ldc.lt.max( 1_wip, n     ) ) then
         info = 12
      endif
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wfsyr2k', info )
         return
      endif
!
!     Quick return if possible.
!
      if ( ( n.eq.0_wip ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0_wip ) ).and.
     $      ( beta.eq.one ) ) )
     $   return
!
!     And when  alpha.eq.zero.
!
      if ( alpha .eq. zero ) then
         if ( upper ) then
            if ( beta .eq. zero ) then
               do j = 1, n
                  do i = 1, j
                     c(i,j) = zero
                  enddo 
               enddo 
            else
               do j = 1, n
                  do i = 1, j
                     c(i,j) = beta*c(i,j)
                  enddo 
               enddo 
            endif
         else
            if ( beta .eq. zero ) then
               do j = 1, n
                  do i = j, n
                     c(i,j) = zero
                  enddo 
               enddo 
            else
               do j = 1, n
                  do i = j, n
                     c(i,j) = beta*c(i,j)
                  enddo 
               enddo 
            endif
         endif
         return
      endif
!
!     Start the operations.
!
      if ( wlsame ( trans, 'N' ) ) then
!
!        Form  C := alpha*A*B' + alpha*B*A' + C.
!
         if ( upper ) then
            do j = 1, n
               if ( beta .eq. zero ) then
                  do i = 1, j
                     c(i,j) = zero
                  enddo 
               else if ( beta .ne. one ) then
                  do i = 1, j
                     c(i,j) = beta*c(i,j)
                  enddo 
               endif
               do l = 1, k
                  if ( ( a(j,l) .ne. zero ) .or.
     $                 ( b(j,l) .ne. zero )    ) then
                     temp1 = alpha*b(j,l)
                     temp2 = alpha*a(j,l)
                     do i = 1, j
                        c(i,j) = c(i,j)       +
     $                           a(i,l)*temp1 + 
     $                           b(i,l)*temp2
                     enddo 
                  endif
               enddo 
            enddo 
         else
            do j = 1, n
               if ( beta .eq. zero ) then
                  do i = j, n
                     c(i,j) = zero
                  enddo 
               else if ( beta .ne. one ) then
                  do i = j, n
                     c(i,j) = beta*c(i,j)
                  enddo 
               endif
               do l = 1, k
                  if ( ( a(j,l) .ne. zero ).or.
     $                 ( b(j,l) .ne. zero )     ) then
                     temp1 = alpha*b(j,l)
                     temp2 = alpha*a(j,l)
                     do i = j, n
                        c(i,j) = c(i,j)       +
     $                           a(i,l)*temp1 + 
     $                           b(i,l)*temp2
                     enddo 
                  endif
               enddo 
            enddo 
         endif
      else
!
!        Form  C := alpha*A'*B + alpha*B'*A + C.
!
         if ( upper ) then
            do j = 1, n
               do i = 1, j
                  temp1 = zero
                  temp2 = zero
                  do l = 1, k
                     temp1 = temp1 + a(l,i)*b(l,j)
                     temp2 = temp2 + b(l,i)*a(l,j)
                  enddo 
                  if ( beta.eq.zero ) then
                     c(i,j) = alpha*temp1 + alpha*temp2
                  else
                     c(i,j) = beta *c(i,j) +
     $                        alpha*temp1 + alpha*temp2
                  endif
               enddo 
            enddo 
         else
            do j = 1, n
               do i = j, n
                  temp1 = zero
                  temp2 = zero
                  do l = 1, k
                     temp1 = temp1 + a(l,i)*b(l,j)
                     temp2 = temp2 + b(l,i)*a(l,j)
                  enddo 
                  if ( beta.eq.zero ) then
                     c(i,j) = alpha*temp1 + alpha*temp2
                  else
                     c(i,j) = beta *c(i,j) +
     $                        alpha*temp1 + alpha*temp2
                  endif
               enddo 
            enddo 
         endif
      endif
      return
      end subroutine qpsyr2k 
!
!***********************************************************************
!
      subroutine qptrmm ( 
     $           side, uplo, transa, diag, m, n, alpha, 
     $           a, lda, b, ldb )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      character(len=1)   
     $   side, uplo, transa, diag
      integer(wip)     
     $   m, n, lda, ldb
      real(wrp) 
     $   alpha,
     $   a(lda,*), b(ldb,*)
!
!  Purpose
!  =======
!
!  WFTRMM  performs one of the matrix-matrix operations
!
!     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  Parameters
!  ==========
!
!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!
!           Unchanged on exit.
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(wrp)       .
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - REAL(wrp)        array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - REAL(wrp)        array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      logical 
     $   wlsame
      external 
     $   wlsame
      external 
     $   wxerbla
      intrinsic
     $   max
!
!     LOCAL VARIABLES:
      logical 
     $   lside, nounit, upper
      integer(wip) 
     $   i, j, k, nrowa
      integer(wip)
     $   info
      real(wrp) 
     $   temp
      real(wrp)  
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
!     
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      lside  = wlsame ( side, 'L' )
      if ( lside ) then
         nrowa = m
      else
         nrowa = n
      endif
      nounit = wlsame ( diag, 'N' )
      upper  = wlsame ( uplo, 'U' )
!
      info   = 0
      if ( ( .not. lside                ) .and.
     $     ( .not. wlsame ( side, 'R' ) )   ) then
         info = 1
      else if (( .not. upper                ) .and.
     $         ( .not. wlsame ( uplo, 'L' ) )      ) then
         info = 2
      else if (( .not. wlsame ( transa, 'N' ) ) .and.
     $         ( .not. wlsame ( transa, 'T' ) ) .and.
     $         ( .not. wlsame ( transa, 'C' ) )      ) then
         info = 3
      else if (( .not.wlsame ( diag, 'U' ) ) .and.
     $         ( .not.wlsame ( diag, 'N' ) )      ) then
         info = 4
      else if ( m  .lt.0_wip               ) then
         info = 5
      else if ( n  .lt.0_wip               ) then
         info = 6
      else if ( lda.lt.max( 1_wip, nrowa ) ) then
         info = 9
      else if ( ldb.lt.max( 1_wip, m     ) ) then
         info = 11
      endif
      if ( info .ne. 0 ) then
         call wxerbla ( 
     $        'wftrmm ', info )
         return
      endif
!
!     Quick return if possible.
!
      if ( n .eq. 0_wip ) return
!
!     And when  alpha.eq.zero.
!
      if ( alpha .eq. zero ) then
         do j = 1, n
            do i = 1, m
               b(i,j) = zero
            enddo 
         enddo 
         return
      endif
!
!     Start the operations.
!
      if ( lside ) then
         if ( wlsame ( transa, 'N' ) ) then
!
!           Form  B := alpha*A*B.
!
            if ( upper ) then
               do j = 1, n
                  do k = 1, m
                     if ( b(k,j) .ne. zero ) then
                        temp = alpha*b(k,j)
                        do i = 1, k - 1
                           b(i,j) = b(i,j) + temp*a(i,k)
                        enddo 
                        if ( nounit ) temp = temp*a(k,k)
                        b(k,j) = temp
                     endif
                  enddo 
               enddo 
            else
               do j = 1, n
                  do k = m, 1, -1
                     if ( b(k,j) .ne. zero ) then
                        temp   = alpha*b(k,j)
                        b(k,j) = temp
                        if ( nounit ) b(k,j) = b(k,j)*a(k,k)
                        do i = k + 1, m
                           b(i,j) = b(i,j) + temp*a(i,k)
                        enddo 
                     endif
                  enddo 
               enddo 
            endif
         else
!
!           Form  B := alpha*A'*B.
!
            if ( upper ) then
               do j = 1, n
                  do i = m, 1, -1
                     temp = b(i,j)
                     if ( nounit ) temp = temp*a(i,i)
                     do k = 1, i - 1
                        temp = temp + a(k,i)*b(k,j)
                     enddo 
                     b(i,j) = alpha*temp
                  enddo 
               enddo 
            else
               do j = 1, n
                  do i = 1, m
                     temp = b(i,j)
                     if ( nounit ) temp = temp*a(i,i)
                     do k = i + 1, m
                        temp = temp + a(k,i)*b(k,j)
                     enddo 
                     b(i,j) = alpha*temp
                  enddo 
               enddo 
            endif
         endif
      else
         if ( wlsame ( transa, 'N' ) ) then
!
!           Form  B := alpha*B*A.
!
            if ( upper ) then
               do j = n, 1, -1
                  temp = alpha
                  if ( nounit ) temp = temp*a(j,j)
                  do i = 1, m
                     b(i,j) = temp*b(i,j)
                  enddo 
                  do k = 1, j - 1
                     if ( a(k,j) .ne. zero ) then
                        temp = alpha*a(k,j)
                        do i = 1, m
                           b(i,j) = b(i,j) + temp*b(i,k)
                        enddo 
                     endif
                  enddo 
               enddo 
            else
               do j = 1, n
                  temp = alpha
                  if ( nounit ) temp = temp*a(j,j)
                  do i = 1, m
                     b(i,j) = temp*b(i,j)
                  enddo 
                  do k = j + 1, n
                     if ( a(k,j) .ne. zero ) then
                        temp = alpha*a(k,j)
                        do i = 1, m
                           b(i,j) = b(i,j) + temp*b(i,k)
                        enddo 
                     endif
                  enddo 
               enddo 
            endif
         else
!
!           Form  B := alpha*B*A'.
!
            if ( upper ) then
               do k = 1, n
                  do j = 1, k - 1
                     if ( a(j,k) .ne. zero ) then
                        temp = alpha*a(j,k)
                        do i = 1, m
                           b(i,j) = b(i,j) + temp*b(i,k)
                        enddo 
                     endif
                  enddo 
                  temp = alpha
                  if ( nounit ) temp = temp*a(k,k)
                  if ( temp .ne. one ) then
                     do i = 1, m
                        b(i,k) = temp*b(i,k)
                     enddo 
                  endif
               enddo 
            else
               do k = n, 1, -1
                  do j = k + 1, n
                     if ( a(j,k) .ne. zero ) then
                        temp = alpha*a(j,k)
                        do i = 1, m
                           b(i,j) = b(i,j) + temp*b(i,k)
                        enddo 
                     endif
                  enddo 
                  temp = alpha
                  if ( nounit ) temp = temp*a(k,k)
                  if ( temp .ne. one ) then
                     do i = 1, m
                        b(i,k) = temp*b(i,k)
                     enddo 
                  endif
               enddo 
            endif
         endif
      endif
      return
      end subroutine qptrmm
!
!-----------------------------------------------------------------------
!***********************************************************************
!
!     LAPACK DSYEV::
!
!     Remove the rest when you call only ARPACK (e.g., for NS=0).
!
!     For the cases NS>0, we must keep this part since we need to 
!     call DSYEV for moderate LJJ.
!
!-----------------------------------------------------------------------
!
!@@begin of DSYEV
!
!-----------------------------------------------------------------------
!
      subroutine qpsyev ( 
     $           jobz, uplo, n, a, lda, w, work, lwork, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character          
     $   jobz, uplo
      integer(wip)       
     $   n, lda, lwork
      real(wrp) 
     $   a(lda,*), w(*), work(*)
      integer(wip)  
     $   info
! 
!  -- LAPACK driver routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
!  real SYMMETRIC matrix A.
!
!  Arguments
!  =========
!
!  JOBZ    (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  W       (output) REAL(wrp)        array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
!  WORK    (workspace/output) REAL(wrp)        array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The length of the array WORK.  LWORK >= max(1,3*N-1).
!          For optimal efficiency, LWORK >= (NB+2)*N,
!          where NB is the blocksize for QPSYTRD returned by WILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by WXERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.
!
!  =====================================================================
!     
!     DEPENDENCES:
      logical            
     $   wlsame
      integer(wip)       
     $   wilaenv
      real(wrp)          
     $   qplamch, 
     $   qplansy
      external           
     $   wlsame, 
     $   wilaenv, 
     $   qplamch, 
     $   qplansy
      external           
     $   qplascl, 
     $   qporgtr, 
     $   qpscal, 
     $   qpsteqr, 
     $   qpsterf, 
     $   qpsytrd, 
     $   wxerbla
      intrinsic
     $   max, sqrt
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   zero, one
      parameter ( 
     $   zero = 0.0_wrp, 
     $   one  = 1.0_wrp )
      logical            
     $   lower, lquery, wantz
      integer(wip)       
     $   imax, inde, indtau, indwrk, iscale,
     $   llwork, lopt, lwkopt, nb, iinfo
      real(wrp)          
     $   anrm, bignum, eps, rmax, rmin, safmin, 
     $   sigma, smlnum
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      wantz = wlsame ( jobz, 'V' )
      lower = wlsame ( uplo, 'L' )
      lquery = ( lwork.eq.-1 )
!
      info = 0
      if ( .not.( wantz .or. wlsame ( jobz, 'N' ) ) ) then
         info = -1
      else if ( .not.( lower .or. wlsame ( uplo, 'U' ) ) ) then
         info = -2
      else if ( n.lt.0 ) then
         info = -3
      else if ( lda.lt.max( 1_wip, n ) ) then
         info = -5
      else if ( lwork.lt.max( 1_wip, 3*n-1 ) .and. .not.lquery ) then
         info = -8
      endif
!
      if ( info.eq.0 ) then
         nb = wilaenv ( 
     $        1, 'WRSYTRD', uplo, 
     $        n, -1_wip, -1_wip, -1_wip )
         lwkopt = max( 1_wip, (nb+2)*n )
         work(1) = lwkopt
      endif
!
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wfsyev ', -info )
         return
      else if ( lquery ) then
!
!        To obtain only the length of work(:), work(1) = lwkopt,
!        and exit.
         return
!
      endif
!
!     Quick return if possible
!
      if ( n.eq.0 ) then
         work(1) = 1
         return
      endif
!
      if ( n.eq.1 ) then
         w(1) = a(1,1)
         work(1) = 3
         if ( wantz )
     $      a(1,1) = one
         return
      endif
!
!     Get machine constants.
!
      safmin = qplamch ( 'safe minimum' )
      eps = qplamch ( 'precision' )
      smlnum = safmin / eps
      bignum = one / smlnum
      rmin = sqrt( smlnum )
      rmax = sqrt( bignum )
!
!     Scale matrix to allowable range, if necessary.
!
      anrm = qplansy ( 
     $       'M', uplo, n, a, lda, work )
      iscale = 0
      if ( anrm.gt.zero .and. anrm.lt.rmin ) then
         iscale = 1
         sigma = rmin / anrm
      else if ( anrm.gt.rmax ) then
         iscale = 1
         sigma = rmax / anrm
      endif
      if ( iscale.eq.1 )
     $   call qplascl ( 
     $        uplo, 0_wip, 0_wip, one, sigma, 
     $        n, n, a, lda, info )
!
!     Call QPSYTRD to reduce symmetric matrix to tridiagonal form.
!
      inde = 1
      indtau = inde + n
      indwrk = indtau + n
      llwork = lwork - indwrk + 1
      call qpsytrd ( 
     $     uplo, n, a, lda, w, work(inde), work(indtau),
     $     work(indwrk), llwork, iinfo )
      lopt = 2*n + work(indwrk)
!
!     For eigenvalues only, call WFSTERF.  For eigenvectors, first call
!     WFORGTR to generate the orthogonal matrix, then call WFSTEQR.
!
      if ( .not.wantz ) then
         call qpsterf ( 
     $        n, w, work(inde), info )
      else
         call qporgtr ( 
     $        uplo, n, a, lda, work(indtau), work(indwrk),
     $        llwork, iinfo )
         call qpsteqr ( 
     $        jobz, n, w, work(inde), a, lda, work(indtau),
     $        info )
      endif
!
!     If matrix was scaled, then rescale eigenvalues appropriately.
!
      if ( iscale.eq.1 ) then
         if ( info.eq.0 ) then
            imax = n
         else
            imax = info-1
         endif
         call qpscal ( 
     $        imax, one / sigma, w, 1_wip )
      endif
!
!     Set WORK(1) to optimal workspace size.
!
      work(1) = lwkopt
!
      return
      end subroutine qpsyev
!
!***********************************************************************
!
      subroutine qpsytd2 ( 
     $           uplo, n, a, lda, d, e, tau, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character          
     $   uplo
      integer(wip)       
     $   n, lda
      real(wrp)          
     $   a(lda,*), d(*), e(*), tau(*)
      integer(wip)       
     $   info
! 
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
!  form T by an orthogonal similarity transformation: Q' * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          qpitten by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) REAL(wrp)        array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) REAL(wrp)        array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) REAL(wrp)        array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I-tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I-tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!     
!     DEPENDENCES:
      external           
     $   qpaxpy, 
     $   qplarfg, 
     $   qpsymv, 
     $   qpsyr2, 
     $   wxerbla
      logical            
     $   wlsame
      real(wrp)          
     $   qpdot
      external           
     $   wlsame, 
     $   qpdot
      intrinsic          
     $   max, min
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   one, zero, half
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp,
     $   half = 0.5_wrp  )
      logical            
     $   upper
      integer(wip)       
     $   i
      real(wrp)          
     $   alpha, taui
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters
!
      info = 0
      upper = wlsame ( uplo, 'U' )
      if ( .not.upper .and. .not.wlsame ( uplo, 'L' ) ) then
         info = -1
      else if ( n.lt.0_wip ) then
         info = -2
      else if ( lda.lt.max( 1_wip, n ) ) then
         info = -4
      endif
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wfsytd2', -info )
         return
      endif
!
!     Quick return if possible
!
      if ( n.le.0_wip )
     $   return
!
      if ( upper ) then
!
!        Reduce the upper triangle of A
!
         do i = n-1, 1, -1
!
!           Generate elementary reflector H(i) = I-tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
            call qplarfg ( 
     $           i, a(i,i+1), a(1,i+1), 1_wip, taui )
            e(i) = a(i,i+1)
!
            if ( taui.ne.zero ) then
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               a(i,i+1) = one
!
!              Compute  x := tau * A * v  storing x in TAU(1:i)
!
               call qpsymv ( 
     $              uplo, i, taui, a, lda, a(1,i+1), 
     $              1_wip, zero, tau, 1_wip )
!
!              Compute  w := x-1/2 * tau * (x'*v) * v
!
               alpha = -half * taui * 
     $                  qpdot( 
     $                  i, tau, 1_wip, a(1,i+1), 1_wip )
!     
               call qpaxpy ( 
     $              i, alpha, a(1,i+1), 1_wip, tau, 1_wip )
!
!              Apply the transformation as a rank-2 update:
!                 A := A-v * w'-w * v'
!
               call qpsyr2 ( 
     $              uplo, i, -one, a(1,i+1), 
     $              1_wip, tau, 1_wip, a, lda )
!
               a(i,i+1) = e(i)
            endif
            d(i+1) = a(i+1,i+1)
            tau(i) = taui
         enddo 
         d(1) = a(1,1)
      else
!
!        Reduce the lower triangle of A
!
         do i = 1, n-1
!
!           Generate elementary reflector H(i) = I-tau * v * v'
!           to annihilate A(i+2:n,i)
!
            call qplarfg ( 
     $           n-i, a(i+1,i), a(min(i+2,n),i), 
     $           1_wip, taui )
            e(i) = a(i+1,i)
!
            if ( taui.ne.zero ) then
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               a(i+1,i) = one
!
!              Compute  x := tau * A * v  storing y in TAU(i:n-1)
!
               call qpsymv ( 
     $              uplo, n-i, taui, a(i+1,i+1), lda,
     $              a(i+1,i), 1_wip, zero, tau(i), 1_wip )
!
!              Compute  w := x-1/2 * tau * (x'*v) * v
!
               alpha = -half * taui *
     $                  qpdot ( 
     $                  n-i, tau(i), 1_wip, 
     $                  a(i+1,i), 1_wip )
!     
               call qpaxpy ( 
     $              n-i, alpha, a(i+1,i), 1_wip, 
     $              tau(i), 1_wip )
!
!              Apply the transformation as a rank-2 update:
!                 A := A-v * w'-w * v'
!
               call qpsyr2 ( 
     $              uplo, n-i, -one, a(i+1,i), 1_wip, 
     $              tau(i), 1_wip, a(i+1,i+1), lda )
!
               a(i+1,i) = e(i)
            endif
            d(i) = a(i,i)
            tau(i) = taui
         enddo 
         d(n) = a(n,n)
      endif
      return
      end subroutine qpsytd2 
!
!***********************************************************************
!
      subroutine qpsytrd ( 
     $           uplo, n, a, lda, d, e, tau, work, lwork, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character          
     $   uplo
      integer(wip)       
     $   n, lda, lwork
      real(wrp)          
     $   a(lda,*), d(*), e(*), tau(*), work(*)
      integer(wip)       
     $   info
! 
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  DSYTRD reduces a real symmetric matrix A to real symmetric
!  tridiagonal form T by an orthogonal similarity transformation:
!  Q**T * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          qpitten by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) REAL(wrp)        array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) REAL(wrp)        array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) REAL(wrp)        array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) REAL(wrp)        array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= 1.
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by WXERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I-tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I-tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     DEPENDECNES:
      external           
     $   qplatrd, 
     $   qpsyr2k, 
     $   qpsytd2, 
     $   wxerbla
      intrinsic          
     $   max
      logical            
     $   wlsame
      integer(wip)       
     $   wilaenv
      external           
     $   wlsame, 
     $   wilaenv
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   one
      parameter ( 
     $   one = 1.0_wrp )
      logical            
     $   lquery, upper
      integer(wip)       
     $   i, iws, j, kk, ldwork, 
     $   lwkopt, nb, nbmin, nx, 
     $   iinfo
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters
!
      info = 0
      upper = wlsame ( uplo, 'U' )
      lquery = ( lwork.eq.-1 )
      if ( .not.upper .and. .not.wlsame ( uplo, 'L' ) ) then
         info = -1
      else if ( n.lt.0_wip ) then
         info = -2
      else if ( lda.lt.max( 1_wip, n ) ) then
         info = -4
      else if ( lwork.lt.1_wip .and. .not.lquery ) then
         info = -9
      endif
!
      if ( info.eq.0 ) then
!
!        Determine the block size.
!
         nb = wilaenv ( 
     $        1, 'WRSYTRD', uplo, 
     $        n, -1_wip, -1_wip, -1_wip )
         lwkopt = n*nb
         work(1) = lwkopt
      endif
!
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wfsytrd', -info )
         return
      else if ( lquery ) then
         return
      endif
!
!     Quick return if possible
!
      if ( n.eq.0_wip ) then
         work(1) = 1
         return
      endif
!
      nx  = n
      iws = 1
      if ( nb.gt.1_wip .and. nb.lt.n ) then
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
         nx = max ( nb, 
     $              wilaenv ( 
     $              3, 'WRSYTRD', uplo, 
     $              n, -1_wip, -1_wip, -1_wip ) )
         if ( nx.lt.n ) then
!
!           Determine if workspace is large enough for blocked code.
!
            ldwork = n
            iws = ldwork*nb
            if ( lwork.lt.iws ) then
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code by setting NX = N.
!
               nb = max( lwork / ldwork, 1_wip )
               nbmin = wilaenv ( 
     $                 2, 'WRSYTRD', uplo, 
     $                 n, -1_wip, -1_wip, -1_wip )
               if ( nb .lt. nbmin )
     $            nx = n
            endif
         else
            nx = n
         endif
      else
         nb = 1
      endif
!
      if ( upper ) then
!
!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.
!
         kk = n - ( (n-nx+nb-1)/nb )*nb
         do i = n-nb+1, kk+1, -nb
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            call qplatrd ( 
     $           uplo, i+nb-1, nb, a, lda, e, tau, work,
     $           ldwork )
!
!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A-V*W'-W*V'
!
            call qpsyr2k ( 
     $           uplo, 'No transpose', i-1, nb, -one, a(1,i),
     $           lda, work, ldwork, one, a, lda )
!
!           Copy superdiagonal elements back into A, and diagonal
!           elements into D
!
            do j = i, i+nb-1
               a(j-1,j) = e(j-1)
               d(j) = a(j,j)
            enddo 
         enddo 
!
!        Use unblocked code to reduce the last or only block
!
         call qpsytd2 ( 
     $        uplo, kk, a, lda, d, e, tau, iinfo )
      else
!
!        Reduce the lower triangle of A
!
         do i = 1, n-nx, nb
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            call qplatrd ( 
     $           uplo, n-i+1, nb, a(i,i), lda, e(i),
     $           tau(i), work, ldwork )
!
!           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
!           an update of the form:  A := A-V*W'-W*V'
!
            call qpsyr2k ( 
     $           uplo, 'No transpose', n-i-nb+1, nb, -one,
     $           a(i+nb,i), lda, work(nb+1), ldwork, one,
     $           a(i+nb,i+nb), lda )
!
!           Copy subdiagonal elements back into A, and diagonal
!           elements into D
!
            do j = i, i+nb-1
               a(j+1,j) = e(j)
               d(j) = a(j,j)
            enddo 
         enddo 
!
!        Use unblocked code to reduce the last or only block
!
         call qpsytd2 ( 
     $        uplo, n-i+1, a(i,i), lda, d(i), e(i),
     $        tau(i), iinfo )
      endif
!
      work(1) = lwkopt
      return
      end subroutine qpsytrd  
!
!***********************************************************************
!
      subroutine qplarfb ( 
     $           side, trans, cdirect, storev, m, n, k, v, ldv,
     $           t, ldt, c, ldc, work, ldwork )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character 
     $   side, trans, cdirect, storev
      integer(wip) 
     $   k, ldc, ldt, ldv, ldwork, m, n
      real(wrp) 
     $   c(ldc,*), t(ldt,*), v(ldv,*),
     $   work(ldwork,*)
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)
!
!  CDIRECT (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) REAL(wrp)        array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) REAL(wrp)        array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) REAL(wrp)        array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) REAL(wrp)        array, dimension (LDWORK,K)
!
!  LDWORK  (input) INTEGER
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
!  =====================================================================
!
!     DEPENDENCES:
      logical  
     $   wlsame
      external 
     $   wlsame
      external 
     $   qpcopy, 
     $   qpgemm, 
     $   qptrmm
!
!     LOCAL VARIABLES:
      real(wrp)  
     $   one
      parameter (
     $   one = 1.0_wrp )
      character   
     $   transt
      integer(wip) 
     $   i, j
!
!     EXECUTABLE STATEMENTS:
!     Quick return if possible
!
      if ( m.le.0 .or. n.le.0 ) return
!
      if ( wlsame ( trans, 'N' ) ) then
         transt = 'T'
      else
         transt = 'N'
      endif
!
      if ( wlsame ( storev, 'C' ) ) then
!
         if ( wlsame ( cdirect, 'F' ) ) then
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is Unit Lower triangular.
!
            if ( wlsame ( side, 'L' ) ) then
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1+C2'*V2)  (stored in WORK)
!
!              W := C1'
!
               do j = 1, k
                  call qpcopy ( 
     $                 n, c(j,1), ldc, work(1,j), 1_wip )
               enddo 
!
!              W := W * V1
!
               call qptrmm ( 
     $              'Right', 'Lower', 'No transpose', 'Unit', 
     $              n, k, one, v, ldv, work, ldwork )
               if ( m.gt.k ) then
!
!                 W := W+C2'*V2
!
                  call qpgemm ( 
     $                 'Transpose', 'No transpose', 
     $                 n, k, m-k, one, c(k+1,1), ldc, 
     $                 v(k+1,1), ldv, one, work, ldwork )
               endif
!
!              W := W * T'  or  W * T
!
               call qptrmm ( 
     $              'Right', 'Upper', transt, 'Non-unit', 
     $              n, k, one, t, ldt, work, ldwork )
!
!              C := C-V * W'
!
               if ( m.gt.k ) then
!
!                 C2 := C2-V2 * W'
!
                  call qpgemm ( 
     $                 'No transpose', 'Transpose', 
     $                 m-k, n, k, -one, v(k+1,1), ldv, 
     $                 work, ldwork, one, c(k+1,1), ldc )
               endif
!
!              W := W * V1'
!
               call qptrmm ( 
     $              'Right', 'Lower', 'Transpose', 'Unit', 
     $              n, k, one, v, ldv, work, ldwork )
!
!              C1 := C1-W'
!
               do j = 1, k
                  do i = 1, n
                     c(j,i) = c(j,i) - work(i,j)
                  enddo 
               enddo 
!
            else if ( wlsame ( side, 'R' ) ) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1+C2*V2)  (stored in WORK)
!
!              W := C1
!
               do j = 1, k
                  call qpcopy ( 
     $                 m, c(1,j), 1_wip, work(1,j), 1_wip )
               enddo 
!
!              W := W * V1
!
               call qptrmm ( 
     $              'Right', 'Lower', 'No transpose', 'Unit', 
     $              m, k, one, v, ldv, work, ldwork )
               if ( n.gt.k ) then
!
!                 W := W+C2 * V2
!
                  call qpgemm ( 
     $                 'No transpose', 'No transpose', 
     $                 m, k, n-k, one, c(1,k+1), ldc, 
     $                 v(k+1,1), ldv, one, work, ldwork )
               endif
!
!              W := W * T  or  W * T'
!
               call qptrmm ( 
     $              'Right', 'Upper', trans, 'Non-unit',
     $               m, k, one, t, ldt, work, ldwork )
!
!              C := C-W * V'
!
               if ( n.gt.k ) then
!
!                 C2 := C2-W * V2'
!
                  call qpgemm ( 
     $                 'No transpose', 'Transpose', 
     $                 m, n-k, k, -one, work, ldwork, 
     $                 v(k+1,1), ldv, one, c(1,k+1), ldc )
               endif
!
!              W := W * V1'
!
               call qptrmm ( 
     $              'Right', 'Lower', 'Transpose', 'Unit',
     $               m, k, one, v, ldv, work, ldwork )
!
!              C1 := C1-W
!
               do j = 1, k
                  do i = 1, m
                     c(i,j) = c(i,j) - work(i,j)
                  enddo 
               enddo 
            endif
!
         else
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is Unit upper triangular.
!
            if ( wlsame ( side, 'L' ) ) then
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V  =  (C1'*V1+C2'*V2)  (stored in WORK)
!
!              W := C2'
!
               do j = 1, k
                  call qpcopy ( 
     $                 n, c(m-k+j,1), ldc, work(1,j), 1_wip )
               enddo 
!
!              W := W * V2
!
               call qptrmm ( 
     $              'Right', 'Upper', 'No transpose', 'Unit',
     $              n, k, one, v(m-k+1,1), ldv, work, ldwork )
               if ( m.gt.k ) then
!
!                 W := W+C1'*V1
!
                  call qpgemm ( 
     $                 'Transpose', 'No transpose', n, k, m-k,
     $                 one, c, ldc, v, ldv, one, work, ldwork )
               endif
!
!              W := W * T'  or  W * T
!
               call qptrmm ( 
     $              'Right', 'Lower', transt, 'Non-unit',
     $              n, k, one, t, ldt, work, ldwork )
!
!              C := C-V * W'
!
               if ( m.gt.k ) then
!
!                 C1 := C1-V1 * W'
!
                  call qpgemm ( 
     $                 'No transpose', 'Transpose',
     $                 m-k, n, k, -one, v, ldv, work, ldwork, 
     $                 one, c, ldc )
               endif
!
!              W := W * V2'
!
               call qptrmm ( 
     $              'Right', 'Upper', 'Transpose', 'Unit',
     $              n, k, one, v(m-k+1,1), ldv, work, ldwork )
!
!              C2 := C2-W'
!
               do j = 1, k
                  do i = 1, n
                     c(m-k+j,i) = c(m-k+j,i) - work(i,j)
                  enddo 
               enddo 
!
            else if ( wlsame ( side, 'R' ) ) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V  =  (C1*V1+C2*V2)  (stored in WORK)
!
!              W := C2
!
               do j = 1, k
                  call qpcopy ( 
     $                 m, c(1,n-k+j), 1_wip, work(1,j), 1_wip )
               enddo 
!
!              W := W * V2
!
               call qptrmm ( 
     $              'Right', 'Upper', 'No transpose', 'Unit',
     $              m, k, one, v(n-k+1,1), ldv, work, ldwork )
               if ( n.gt.k ) then
!
!                 W := W+C1 * V1
!
                  call qpgemm ( 
     $                 'No transpose', 'No transpose', m, k, n-k,
     $                 one, c, ldc, v, ldv, one, work, ldwork )
               endif
!
!              W := W * T  or  W * T'
!
               call qptrmm ( 
     $              'Right', 'Lower', trans, 'Non-unit', m, k,
     $              one, t, ldt, work, ldwork )
!
!              C := C-W * V'
!
               if ( n.gt.k ) then
!
!                 C1 := C1-W * V1'
!
                  call qpgemm ( 
     $                 'No transpose', 'Transpose', m, n-k, k,
     $                 -one, work, ldwork, v, ldv, one, c, ldc )
               endif
!
!              W := W * V2'
!
               call qptrmm ( 
     $              'Right', 'Upper', 'Transpose', 'Unit', m, k,
     $              one, v(n-k+1,1), ldv, work, ldwork )
!
!              C2 := C2-W
!
               do j = 1, k
                  do i = 1, m
                     c(i,n-k+j) = c(i,n-k+j) - work(i,j)
                  enddo 
               enddo 
            endif
         endif
!
      else if ( wlsame ( storev, 'R' ) ) then
!
         if ( wlsame ( cdirect, 'F' ) ) then
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is Unit upper triangular.
!
            if ( wlsame ( side, 'L' ) ) then
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1'+C2'*V2') (stored in WORK)
!
!              W := C1'
!
               do j = 1, k
                  call qpcopy ( 
     $                 n, c(j,1), ldc, work(1,j), 1_wip )
               enddo 
!
!              W := W * V1'
!
               call qptrmm ( 
     $              'Right', 'Upper', 'Transpose', 'Unit',
     $               n, k, one, v, ldv, work, ldwork )
               if ( m.gt.k ) then
!
!                 W := W+C2'*V2'
!
                  call qpgemm ( 
     $                 'Transpose', 'Transpose', n, k, m-k, one,
     $                 c(k+1,1), ldc, v(1,k+1), ldv, one,
     $                 work, ldwork )
               endif
!
!              W := W * T'  or  W * T
!
               call qptrmm ( 
     $              'Right', 'Upper', transt, 'Non-unit',
     $              n, k, one, t, ldt, work, ldwork )
!
!              C := C-V' * W'
!
               if ( m.gt.k ) then
!
!                 C2 := C2-V2' * W'
!
                  call qpgemm ( 
     $                 'Transpose', 'Transpose', m-k, n, k, -one,
     $                 v(1,k+1), ldv, work, ldwork, one,
     $                 c(k+1,1), ldc )
               endif
!
!              W := W * V1
!
               call qptrmm ( 
     $              'Right', 'Upper', 'No transpose', 'Unit',
     $              n, k, one, v, ldv, work, ldwork )
!
!              C1 := C1-W'
!
               do j = 1, k
                  do i = 1, n
                     c(j,i) = c(j,i) - work(i,j)
                  enddo 
               enddo 
!
            else if ( wlsame ( side, 'R' ) ) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1'+C2*V2')  (stored in WORK)
!
!              W := C1
!
               do j = 1, k
                  call qpcopy ( 
     $                 m, c(1,j), 1_wip, work(1,j), 1_wip )
               enddo 
!
!              W := W * V1'
!
               call qptrmm ( 
     $              'Right', 'Upper', 'Transpose', 'Unit',
     $               m, k, one, v, ldv, work, ldwork )
               if ( n.gt.k ) then
!
!                 W := W+C2 * V2'
!
                  call qpgemm ( 
     $                 'No transpose', 'Transpose', m, k, n-k,
     $                 one, c(1,k+1), ldc, v(1,k+1), ldv,
     $                 one, work, ldwork )
               endif
!
!              W := W * T  or  W * T'
!
               call qptrmm ( 
     $              'Right', 'Upper', trans, 'Non-unit',
     $              m, k, one, t, ldt, work, ldwork )
!
!              C := C-W * V
!
               if ( n.gt.k ) then
!
!                 C2 := C2-W * V2
!
                  call qpgemm ( 
     $                 'No transpose', 'No transpose', m, n-k, k,
     $                 -one, work, ldwork, v(1,k+1), ldv, one,
     $                 c(1,k+1), ldc )
               endif
!
!              W := W * V1
!
               call qptrmm ( 
     $              'Right', 'Upper', 'No transpose', 'Unit',
     $              m, k, one, v, ldv, work, ldwork )
!
!              C1 := C1-W
!
               do j = 1, k
                  do i = 1, m
                     c(i,j) = c(i,j) - work(i,j)
                  enddo 
               enddo 
!
            endif
!
         else
!
!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is Unit Lower triangular.
!
            if ( wlsame ( side, 'L' ) ) then
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W := C' * V'  =  (C1'*V1'+C2'*V2') (stored in WORK)
!
!              W := C2'
!
               do j = 1, k
                  call qpcopy ( 
     $                 n, c(m-k+j,1), ldc, work(1,j), 1_wip )
               enddo 
!
!              W := W * V2'
!
               call qptrmm ( 
     $              'Right', 'Lower', 'Transpose', 'Unit', 
     $              n, k, one, v(1,m-k+1), ldv, work, ldwork )
               if ( m.gt.k ) then
!
!                 W := W+C1'*V1'
!
                  call qpgemm ( 
     $                 'Transpose', 'Transpose', n, k, m-k, one,
     $                 c, ldc, v, ldv, one, work, ldwork )
               endif
!
!              W := W * T'  or  W * T
!
               call qptrmm ( 
     $              'Right', 'Lower', transt, 'Non-unit',
     $              n, k, one, t, ldt, work, ldwork )
!
!              C := C-V' * W'
!
               if ( m.gt.k ) then
!
!                 C1 := C1-V1' * W'
!
                  call qpgemm( 
     $                 'Transpose', 'Transpose', m-k, n, k, -one,
     $                 v, ldv, work, ldwork, one, c, ldc )
               endif
!
!              W := W * V2
!
               call qptrmm ( 
     $              'Right', 'Lower', 'No transpose', 'Unit',
     $              n, k, one, v(1,m-k+1), ldv, work, ldwork )
!
!              C2 := C2-W'
!
               do j = 1, k
                  do i = 1, n
                     c(m-k+j,i) = c(m-k+j,i) - work(i,j)
                  enddo 
               enddo 
!
            else if ( wlsame ( side, 'R' ) ) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W := C * V'  =  (C1*V1'+C2*V2')  (stored in WORK)
!
!              W := C2
!
               do j = 1, k
                  call qpcopy ( 
     $                 m, c(1,n-k+j), 1_wip, work(1,j), 1_wip )
               enddo 
!
!              W := W * V2'
!
               call qptrmm ( 
     $              'Right', 'Lower', 'Transpose', 'Unit',
     $               m, k, one, v(1,n-k+1), ldv, work, ldwork )
               if ( n.gt.k ) then
!
!                 W := W+C1 * V1'
!
                  call qpgemm ( 
     $                 'No transpose', 'Transpose', m, k, n-k,
     $                 one, c, ldc, v, ldv, one, work, ldwork )
               endif
!
!              W := W * T  or  W * T'
!
               call qptrmm ( 
     $              'Right', 'Lower', trans, 'Non-unit',
     $              m, k, one, t, ldt, work, ldwork )
!
!              C := C-W * V
!
               if ( n.gt.k ) then
!
!                 C1 := C1-W * V1
!
                  call qpgemm ( 
     $                 'No transpose', 'No transpose', m, n-k, k,
     $                 -one, work, ldwork, v, ldv, one, c, ldc )
               endif
!
!              W := W * V2
!
               call qptrmm ( 
     $              'Right', 'Lower', 'No transpose', 'Unit', 
     $              m, k, one, v(1,n-k+1), ldv, work, ldwork )
!
!              C1 := C1-W
!
               do j = 1, k
                  do i = 1, m
                     c(i,n-k+j) = c(i,n-k+j) - work(i,j)
                  enddo 
               enddo 
!
            endif
!
         endif
      endif
!
      return
      end subroutine qplarfb 
!
!***********************************************************************
!
      subroutine qplarft ( 
     $           cdirect, storev, n, k, v, ldv, tau, t, ldt )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS:
      character   
     $   cdirect, storev
      integer(wip) 
     $   n, k, ldv, ldt
      real(wrp)   
     $   t(ldt,*), tau(*), v(ldv,*)
!
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If CDIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If CDIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I-V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I-V' * T * V
!
!  Arguments
!  =========
!
!  CDIRECT (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.
!
!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) REAL(wrp)        array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) REAL(wrp)        array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) REAL(wrp)        array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If CDIRECT = 'F', T is upper triangular; if CDIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!  ===============
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  CDIRECT = 'F' and STOREV = 'C':         CDIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  CDIRECT = 'B' and STOREV = 'C':         CDIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
!  =====================================================================
!
!     DEPENDENCES:
      external  
     $   qpgemv, 
     $   qptrmv
      logical  
     $   wlsame
      external 
     $   wlsame
!
!     LOCAL VARIABLES:
      real(wrp)     
     $   one, zero
      parameter  ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      integer(wip) 
     $   i, j
      real(wrp) 
     $   vii
!
!     EXECUTABLE STATEMENTS:
!     Quick return if possible
!
      if ( n.eq.0 ) return
!
      if ( wlsame ( cdirect, 'F' ) ) then
         do i = 1, k
            if ( tau(i).eq.zero ) then
!
!              H(i)  =  I
!
               do j = 1, i
                  t(j,i) = zero
               enddo 
            else
!
!              general case
!
               vii = v(i,i)
               v(i,i) = one
               if ( wlsame ( storev, 'C' ) ) then
!
!                 T(1:i-1,i) :=-tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
                  call qpgemv ( 
     $                 'Transpose', n-i+1, i-1, -tau(i),
     $                 v(i,1), ldv, v(i,i), 1_wip, zero,
     $                 t(1,i), 1_wip )
               else
!
!                 T(1:i-1,i) :=-tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
                  call qpgemv ( 
     $                 'No transpose', i-1, n-i+1, -tau(i),
     $                 v(1,i), ldv, v(i,i), ldv, zero,
     $                 t(1,i), 1_wip )
               endif
               v(i,i) = vii
!
!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
!
               call qptrmv ( 
     $              'Upper', 'No transpose', 'Non-unit',
     $              i-1, t, ldt, t(1,i), 1_wip )
               t(i,i) = tau(i)
            endif
         enddo 
      else
         do i = k, 1, -1
            if ( tau(i).eq.zero ) then
!
!              H(i)  =  I
!
               do j = i, k
                  t(j,i) = zero
               enddo 
            else
!
!              general case
!
               if ( i.lt.k ) then
                  if ( wlsame ( storev, 'C' ) ) then
                     vii = v(n-k+i,i)
                     v(n-k+i,i) = one
!
!                    T(i+1:k,i) :=
!                           -tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                     call qpgemv ( 
     $                    'Transpose', n-k+i, k-i, -tau(i),
     $                    v(1,i+1), ldv, v(1,i), 1_wip, zero,
     $                    t(i+1,i), 1_wip )
                     v(n-k+i,i) = vii
                  else
                     vii = v(i,n-k+i)
                     v(i,n-k+i) = one
!
!                    T(i+1:k,i) :=
!                           -tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                     call qpgemv ( 
     $                    'No transpose', k-i, n-k+i, -tau(i),
     $                    v(i+1,1), ldv, v(i,1), ldv, zero,
     $                    t(i+1,i), 1_wip )
                     v(i,n-k+i) = vii
                  endif
!
!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
!
                  call qptrmv ( 
     $                 'Lower', 'No transpose', 'Non-unit', k-i,
     $                 t(i+1,i+1), ldt, t(i+1,i), 1_wip )
               endif
               t(i,i) = tau(i)
            endif
         enddo 
      endif
      return
      end subroutine qplarft
!
!***********************************************************************
!
      subroutine qplatrd ( 
     $           uplo, n, nb, a, lda, e, tau, w, ldw )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      character   
     $   uplo
      integer(wip)  
     $   n, nb, lda, ldw 
      real(wrp) 
     $   a(lda,*), e(*), tau(*), w(ldw,*)
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLATRD reduces NB rows and columns of a real symmetric matrix A to
!  symmetric tridiagonal form by an orthogonal similarity
!  transformation Q' * A * Q, and returns the matrices V and W which are
!  needed to apply the transformation to the unreduced part of A.
!
!  If UPLO = 'U', WFLATRD reduces the last NB rows and columns of a
!  matrix, of which the upper triangle is supplied;
!  if UPLO = 'L', WFLATRD reduces the first NB rows and columns of a
!  matrix, of which the lower triangle is supplied.
!
!  This is an auxiliary routine called by QPSYTRD.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U': Upper triangular
!          = 'L': Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.
!
!  NB      (input) INTEGER
!          The number of rows and columns to be reduced.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit:
!          if UPLO = 'U', the last NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements above the diagonal
!            with the array TAU, represent the orthogonal matrix Q as a
!            product of elementary reflectors;
!          if UPLO = 'L', the first NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements below the diagonal
!            with the array TAU, represent the  orthogonal matrix Q as a
!            product of elementary reflectors.
!          See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= (1,N).
!
!  E       (output) REAL(wrp)        array, dimension (N-1)
!          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
!          elements of the last NB columns of the reduced matrix;
!          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
!          the first NB columns of the reduced matrix.
!
!  TAU     (output) REAL(wrp)        array, dimension (N-1)
!          The scalar factors of the elementary reflectors, stored in
!          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
!          See Further Details.
!
!  W       (output) REAL(wrp)        array, dimension (LDW,NB)
!          The n-by-nb matrix W required to update the unreduced part
!          of A.
!
!  LDW     (input) INTEGER
!          The leading dimension of the array W. LDW >= max(1,N).
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n) H(n-1) . . . H(n-nb+1).
!
!  Each H(i) has the form
!
!     H(i) = I-tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
!  and tau in TAU(i-1).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(nb).
!
!  Each H(i) has the form
!
!     H(i) = I-tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!  and tau in TAU(i).
!
!  The elements of the vectors v together form the n-by-nb matrix V
!  which is needed, with W, to apply the transformation to the unreduced
!  part of the matrix, using a symmetric rank-2k update of the form:
!  A := A-V*W'-W*V'.
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5 and nb = 2:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  a   a   a   v4  v5 )              (  d                  )
!    (      a   a   v4  v5 )              (  1   d              )
!    (          a   1   v5 )              (  v1  1   a          )
!    (              d   1  )              (  v1  v2  a   a      )
!    (                  d  )              (  v1  v2  a   a   a  )
!
!  where d denotes a diagonal element of the reduced matrix, a denotes
!  an element of the original matrix that is unchanged, and vi denotes
!  an element of the vector defining H(i).
!
!  =====================================================================
!
!     DEPENDENCES: 
      external
     $   qpaxpy, 
     $   qpgemv, 
     $   qplarfg, 
     $   qpscal, 
     $   qpsymv
      logical  
     $   wlsame
      real(wrp)
     $   qpdot
      external  
     $   wlsame, 
     $   qpdot
      intrinsic 
     $   min
!
!     LOCAL VARIABLES: 
      real(wrp) 
     $   zero, one, half
      parameter ( 
     $   zero = 0.0_wrp, 
     $   one  = 1.0_wrp, 
     $   half = 0.5_wrp )
      integer(wip)
     $   i, iw
      real(wrp)  
     $   alpha
!
!     EXECUTABLE STATEMENTS:
!     Quick return if possible
!
      if ( n .le. 0_wip ) return
!
      if ( wlsame ( uplo, 'U' ) ) then
!
!        Reduce last NB columns of upper triangle
!
         do i = n, n-nb+1, -1
            iw = i-n+nb
            if ( i .lt. n ) then
!
!              Update A(1:i,i)
!
               call qpgemv ( 
     $              'No transpose', i, n-i, -one, a(1,i+1),
     $              lda, w(i,iw+1), ldw, one, a(1,i), 1_wip )
               call qpgemv ( 
     $              'No transpose', i, n-i, -one, w(1,iw+1),
     $              ldw, a(i,i+1), lda, one, a(1,i), 1_wip )
            endif
            if ( i .gt. 1_wip ) then
!
!              Generate elementary reflector H(i) to annihilate
!              A(1:i-2,i)
!
               call qplarfg ( 
     $              i-1, a(i-1,i), a(1,i), 1_wip, tau(i-1) )
               e(i-1)   = a(i-1,i)
               a(i-1,i) = one
!
!              Compute W(1:i-1,i)
!
               call qpsymv ( 
     $              'Upper', i-1, one, a, lda, a(1,i), 1_wip,
     $              zero, w(1,iw), 1_wip )
               if ( i.lt.n ) then
                  call qpgemv ( 
     $                 'Transpose', i-1, n-i, one, w(1,iw+1),
     $                 ldw, a(1,i), 1_wip, zero, 
     $                 w(i+1,iw), 1_wip )
                  call qpgemv ( 
     $                 'No transpose', i-1, n-i, -one,
     $                 a(1,i+1), lda, w(i+1,iw), 1_wip, one,
     $                 w(1,iw), 1_wip )
                  call qpgemv ( 
     $                 'Transpose', i-1, n-i, one, a(1,i+1),
     $                 lda, a(1,i), 1_wip, zero, 
     $                 w(i+1,iw), 1_wip )
                  call qpgemv ( 
     $                 'No transpose', i-1, n-i, -one,
     $                 w(1,iw+1), ldw, w(i+1,iw), 1_wip, one,
     $                 w(1,iw), 1_wip )
               endif
               call qpscal ( 
     $              i-1, tau(i-1), w(1,iw), 1_wip )
               alpha = -half * tau(i-1) * 
     $                  qpdot ( 
     $                  i-1, w(1,iw), 1_wip, a(1,i), 1_wip )
               call qpaxpy ( 
     $              i-1, alpha, a(1,i), 1_wip, w(1,iw), 1_wip )
            endif
!
         enddo 
      else
!
!        Reduce first NB columns of lower triangle
!
         do i = 1, nb
!
!           Update A(i:n,i)
!
            call qpgemv ( 
     $           'No transpose', n-i+1, i-1, -one, a(i,1),
     $           lda, w(i,1), ldw, one, a(i,i), 1_wip )
            call qpgemv ( 
     $           'No transpose', n-i+1, i-1, -one, w(i,1),
     $           ldw, a(i,1), lda, one, a(i,i), 1_wip )
            if ( i .lt. n ) then
!
!              Generate elementary reflector H(i) to annihilate
!              A(i+2:n,i)
!
               call qplarfg ( 
     $              n-i, a(i+1,i), a(min(i+2,n),i), 1_wip,
     $              tau(i) )
               e(i)     = a(i+1,i)
               a(i+1,i) = one
!
!              Compute W(i+1:n,i)
!
               call qpsymv ( 
     $              'Lower', n-i, one, a(i+1,i+1), lda,
     $              a(i+1,i), 1_wip, zero, w(i+1,i), 1_wip )
               call qpgemv ( 
     $              'Transpose', n-i, i-1, one, w(i+1,1), ldw,
     $              a(i+1,i), 1_wip, zero, w(1,i), 1_wip )
               call qpgemv ( 
     $              'No transpose', n-i, i-1, -one, a(i+1,1),
     $              lda, w(1,i), 1_wip, one, w(i+1,i), 1_wip )
               call qpgemv ( 
     $              'Transpose', n-i, i-1, one, a(i+1,1), lda,
     $              a(i+1,i), 1_wip, zero, w(1,i), 1_wip )
               call qpgemv ( 
     $              'No transpose', n-i, i-1, -one, w(i+1,1),
     $              ldw, w(1,i), 1_wip, one, w(i+1,i), 1_wip )
               call qpscal ( 
     $              n-i, tau(i), w(i+1,i), 1_wip )
               alpha = -half * tau(i) * 
     $                  qpdot( 
     $                  n-i, w(i+1,i), 1_wip, a(i+1,i), 1_wip )
               call qpaxpy ( 
     $              n-i, alpha, a(i+1,i), 1_wip, 
     $              w(i+1,i), 1_wip )
            endif
!
         enddo 
      endif
!
      return
      end subroutine qplatrd
!
!***********************************************************************
!

!
      subroutine qporg2l ( 
     $           m, n, k, a, lda, tau, work, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip) 
     $   m, n, k, lda 
      real(wrp)   
     $   a(lda,*), tau(*), work(*)
      integer(wip) 
     $   info
! 
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFORG2L generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the last n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the m by n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) REAL(wrp)        array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  WORK    (workspace) REAL(wrp)        array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     DEPENDENCES: 
      external 
     $   qplarf, 
     $   qpscal, 
     $   wxerbla
      intrinsic 
     $   max
!
!     LOCAL VARIABLES: 
      real(wrp) 
     $   one, zero
      parameter ( 
     $   one = 1.0_wrp, zero = 0.0_wrp )
      integer(wip)  
     $   i, ii, j, l
!
!     EXECUTABLE STATEMENTS:
!     Test the input arguments
!
      info = 0
      if ( m .lt. 0_wip ) then
         info = -1
      else if ( n .lt. 0_wip .or. n .gt. m ) then
         info = -2
      else if ( k .lt. 0_wip .or. k .gt. n ) then
         info = -3
      else if ( lda .lt. max(1_wip, m) ) then
         info = -5
      endif
      if ( info .ne. 0 ) then
         call wxerbla ( 
     $        'wforg2l', -info )
         return
      endif
!
!     Quick return if possible
!
      if ( n .le. 0_wip ) return
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
      do j = 1, n-k
         do l = 1, m
            a(l,j) = zero
         enddo 
         a(m-n+j,j) = one
      enddo 
!
      do i = 1, k
         ii = n-k+i
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
         a(m-n+ii,ii) = one
         call qplarf ( 
     $        'Left', m-n+ii, ii-1, a(1,ii), 1_wip, tau(i), a,
     $        lda, work )
         call qpscal ( 
     $        m-n+ii-1, -tau(i), a(1,ii), 1_wip )
         a(m-n+ii,ii) = one - tau(i)
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
         do l = m-n+ii+1, m
            a(l,ii) = zero
         enddo 
      enddo 
      return
      end subroutine qporg2l 
!
!***********************************************************************
!
      subroutine qporg2r ( 
     $           m, n, k, a, lda, tau, work, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)  
     $   m, n, k, lda 
      real(wrp)  
     $   a(lda,*), tau(*), work(*)
      integer(wip)  
     $   info
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFORG2R generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) REAL(wrp) array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) REAL(wrp) array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace) REAL(wrp) array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     DEPENDENCES: 
      external
     $   qplarf, 
     $   qpscal, 
     $   wxerbla
      intrinsic 
     $   max
!
!     LOCAL VARIABLES: 
      real(wrp)  
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      integer(wip)  
     $   i, j, l
!
!     EXECUTABLE STATEMENTS:
!     Test the input arguments
!
      info = 0
      if ( m .lt. 0_wip ) then
         info = -1
      else if ( n .lt. 0_wip .or. n .gt. m ) then
         info = -2
      else if ( k .lt. 0_wip .or. k .gt. n ) then
         info = -3
      else if ( lda .lt. max(1_wip, m) ) then
         info = -5
      endif
      if ( info .ne. 0 ) then
         call wxerbla ( 
     $        'wforg2r', -info )
         return
      endif
!
!     Quick return if possible
!
      if ( n .le. 0_wip ) return
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      do j = k+1, n
         do l = 1, m
            a(l,j) = zero
         enddo 
         a(j,j) = one
      enddo 
!
      do i = k, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         if ( i .lt. n ) then
            a(i,i) = one
            call qplarf ( 
     $           'Left', m-i+1, n-i, a(i,i), 1_wip, tau(i),
     $           a(i,i+1), lda, work )
         endif
         if ( i .lt. m )
     $      call qpscal ( 
     $           m-i, -tau(i), a(i+1,i), 1_wip )
         a(i,i) = one - tau(i)
!
!        Set A(1:i-1,i) to zero
!
         do l = 1, i-1
            a(l,i) = zero
         enddo 
      enddo 
      return
      end subroutine qporg2r  
!
!***********************************************************************
!
      subroutine qporgql ( 
     $           m, n, k, a, lda, tau, work, lwork, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip) 
     $   m, n, k, lda, lwork 
      real(wrp)   
     $   a(lda,*), tau(*), work(*)
      integer(wip) 
     $   info
! 
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFORGQL generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the last N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) REAL(wrp)        array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  WORK    (workspace/output) REAL(wrp)        array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by WXERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     DEPENDENCES: 
      external   
     $   qplarfb, 
     $   qplarft, 
     $   qporg2l, 
     $   wxerbla
      intrinsic  
     $   max, min
      integer(wip) 
     $   wilaenv
      external  
     $   wilaenv
!
!     LOCAL VARIABLES: 
      real(wrp)      
     $   zero
      parameter ( 
     $   zero = 0.0_wrp )
      logical      
     $   lquery
      integer(wip) 
     $   i, ib, iws, j, kk, l, ldwork, lwkopt, 
     $   nb, nbmin, nx, iinfo
!
!     EXECUTABLE STATEMENTS:
!     Test the input arguments
!
      info = 0
      nb = wilaenv ( 
     $     1, 'WFORGQL', ' ', m, n, k, -1_wip )
      lwkopt = max( 1_wip, n )*nb
      work(1) = lwkopt
      lquery = ( lwork .eq. -1 )
      if ( m .lt. 0 ) then
         info = -1
      else if ( n.lt.0 .or. n.gt.m ) then
         info = -2
      else if ( k.lt.0 .or. k.gt.n ) then
         info = -3
      else if ( lda.lt.max( 1_wip, m ) ) then
         info = -5
      else if ( lwork.lt.max( 1_wip, n ) .and. .not.lquery ) then
         info = -8
      endif
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wforgql', -info )
         return
      else if ( lquery ) then
         return
      endif
!
!     Quick return if possible
!
      if ( n.le.0 ) then
         work(1) = 1
         return
      endif
!
      nbmin = 2
      nx = 0
      iws = n
      if ( nb.gt.1 .and. nb.lt.k ) then
!
!        Determine when to cross over from blocked to unblocked code.
!
         nx = max( 0_wip, 
     $             wilaenv ( 
     $             3, 'WFORGQL', ' ', m, n, k, -1_wip ) )
         if ( nx.lt.k ) then
!
!           Determine if workspace is large enough for blocked code.
!
            ldwork = n
            iws = ldwork*nb
            if ( lwork.lt.iws ) then
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               nb = lwork / ldwork
               nbmin = max( 2_wip, 
     $                      wilaenv ( 
     $                      2, 'WFORGQL', ' ', m, n, k, -1_wip ) )
            endif
         endif
      endif
!
      if ( nb.ge.nbmin .and. nb.lt.k .and. nx.lt.k ) then
!
!        Use blocked code after the first block.
!        The last kk columns are handled by the block method.
!
         kk = min( k, ( ( k-nx+nb-1 ) / nb )*nb )
!
!        Set A(m-kk+1:m,1:n-kk) to zero.
!
         do j = 1, n-kk
            do i = m-kk+1, m
               a(i,j) = zero
            enddo 
         enddo 
      else
         kk = 0
      endif
!
!     Use unblocked code for the first or only block.
!
      call qporg2l ( 
     $     m-kk, n-kk, k-kk, a, lda, tau, work, iinfo )
!
      if ( kk.gt.0 ) then
!
!        Use blocked code
!
         do i = k-kk+1, k, nb
            ib = min( nb, k-i+1 )
            if ( n-k+i.gt.1 ) then
!
!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)
!
               call qplarft ( 
     $              'Backward', 'Columnwise', m-k+i+ib-1, ib,
     $              a(1,n-k+i), lda, tau(i), work, ldwork )
!
!              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
!
               call qplarfb ( 
     $              'Left', 'No transpose', 'Backward',
     $              'Columnwise', m-k+i+ib-1, n-k+i-1, ib,
     $              a(1,n-k+i), lda, work, ldwork, a, lda,
     $              work(ib+1), ldwork )
            endif
!
!           Apply H to rows 1:m-k+i+ib-1 of current block
!
            call qporg2l ( 
     $           m-k+i+ib-1, ib, ib, a(1,n-k+i), lda,
     $           tau(i), work, iinfo )
!
!           Set rows m-k+i+ib:m of current block to zero
!
            do j = n-k+i, n-k+i+ib-1
               do l = m-k+i+ib, m
                  a(l,j) = zero
               enddo 
            enddo 
         enddo 
      endif
!
      work(1) = iws
      return
      end subroutine qporgql
!
!***********************************************************************
!
      subroutine qporgqr ( 
     $           m, n, k, a, lda, tau, work, lwork, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)   
     $   m, n, k, lda, lwork 
      real(wrp) 
     $   a(lda,*), tau(*), work(*)
      integer(wip)   
     $   info
!
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFORGQR generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the first N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) REAL(wrp)        array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace/output) REAL(wrp)        array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by WXERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
!  =====================================================================
!
!     DEPENDENCES: 
      external
     $   qplarfb, 
     $   qplarft, 
     $   qporg2r, 
     $   wxerbla
      intrinsic   
     $   max, min
      integer(wip) 
     $   wilaenv
      external  
     $   wilaenv
!
!     LOCAL VARIABLES: 
      real(wrp) 
     $   zero
      parameter ( 
     $   zero = 0.0_wrp )
      logical    
     $   lquery
      integer(wip) 
     $   i, ib, iws, j, ki, kk, l, ldwork, lwkopt, 
     $   nb, nbmin, nx, iinfo
!
!     EXECUTABLE STATEMENTS:
!     Test the input arguments
!
      info = 0
      nb = wilaenv ( 
     $     1, 'WFORGQR', ' ', m, n, k, -1_wip )
      lwkopt = max(1_wip, n)*nb
      work(1) = lwkopt
      lquery = ( lwork.eq.-1 )
      if ( m.lt.0 ) then
         info = -1
      else if ( n.lt.0 .or. n.gt.m ) then
         info = -2
      else if ( k.lt.0 .or. k.gt.n ) then
         info = -3
      else if ( lda .lt. max(1_wip, m) ) then
         info = -5
      else if ( lwork .lt. max(1_wip, n) .and. .not.lquery ) then
         info = -8
      endif
      if ( info .ne. 0 ) then
         call wxerbla ( 
     $        'wforgqr', -info )
         return
      else if ( lquery ) then
         return
      endif
!
!     Quick return if possible
!
      if ( n .le. 0 ) then
         work(1) = 1
         return
      endif
!
      nbmin = 2
      nx = 0
      iws = n
      if ( nb.gt.1 .and. nb.lt.k ) then
!
!        Determine when to cross over from blocked to unblocked code.
!
         nx = max( 0_wip, 
     $             wilaenv ( 
     $             3, 'WFORGQR', ' ', m, n, k, -1_wip ) )
         if ( nx.lt.k ) then
!
!           Determine if workspace is large enough for blocked code.
!
            ldwork = n
            iws = ldwork*nb
            if ( lwork .lt. iws ) then
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
               nb = lwork / ldwork
               nbmin = max( 2_wip, 
     $                      wilaenv ( 
     $                      2, 'WFORGQR', ' ', m, n, k, -1_wip ) )
            endif
         endif
      endif
!
      if ( nb.ge.nbmin .and. nb.lt.k .and. nx.lt.k ) then
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
         ki = ( ( k-nx-1 ) / nb )*nb
         kk = min( k, ki+nb )
!
!        Set A(1:kk,kk+1:n) to zero.
!
         do j = kk+1, n
            do i = 1, kk
               a(i,j) = zero
            enddo 
         enddo 
      else
         kk = 0
      endif
!
!     Use unblocked code for the last or only block.
!
      if ( kk.lt.n )
     $   call qporg2r( 
     $        m-kk, n-kk, k-kk, a(kk+1,kk+1), lda,
     $        tau(kk+1), work, iinfo )
!
      if ( kk.gt.0 ) then
!
!        Use blocked code
!
         do i = ki+1, 1, -nb
            ib = min( nb, k-i+1 )
            if ( i+ib.le.n ) then
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
               call qplarft ( 
     $              'Forward', 'Columnwise', m-i+1, ib,
     $              a(i,i), lda, tau(i), work, ldwork )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
               call qplarfb ( 
     $              'Left', 'No transpose', 'Forward',
     $              'Columnwise', m-i+1, n-i-ib+1, ib,
     $              a(i,i), lda, work, ldwork, a(i,i+ib),
     $              lda, work(ib+1), ldwork )
            endif
!
!           Apply H to rows i:m of current block
!
            call qporg2r ( 
     $           m-i+1, ib, ib, a(i,i), lda, tau(i), work,
     $           iinfo )
!
!           Set rows 1:i-1 of current block to zero
!
            do j = i, i+ib-1
               do l = 1, i-1
                  a(l,j) = zero
               enddo 
            enddo 
         enddo 
      endif
!
      work(1) = iws
      return
      end subroutine qporgqr
!
!***********************************************************************
!
      subroutine qporgtr ( 
     $           uplo, n, a, lda, tau, work, lwork, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      character 
     $   uplo
      integer(wip) 
     $   n, lda, lwork 
      real(wrp)   
     $   a(lda,*), tau(*), work(*)
      integer(wip) 
     $   info
! 
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFORGTR generates a real orthogonal matrix Q which is defined as the
!  product of n-1 elementary reflectors of order N, as returned by
!  QPSYTRD:
!
!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangle of A contains elementary reflectors
!                 from QPSYTRD;
!          = 'L': Lower triangle of A contains elementary reflectors
!                 from QPSYTRD.
!
!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.
!
!  A       (input/output) REAL(wrp)        array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by QPSYTRD.
!          On exit, the N-by-N orthogonal matrix Q.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).
!
!  TAU     (input) REAL(wrp)        array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by QPSYTRD.
!
!  WORK    (workspace/output) REAL(wrp)        array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK. LWORK >= max(1,N-1).
!          For optimum performance LWORK >= (N-1)*NB, where NB is
!          the optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by WXERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     DEPENDENCES: 
      logical 
     $   wlsame
      integer(wip) 
     $   wilaenv
      external  
     $   wlsame, 
     $   wilaenv
      external 
     $   qporgql, 
     $   qporgqr, 
     $   wxerbla
      intrinsic
     $   max
!
!     LOCAL VARIABLES: 
      real(wrp) 
     $   zero, one
      parameter ( 
     $   zero = 0.0_wrp, 
     $   one  = 1.0_wrp )
      logical   
     $   lquery, upper
      integer(wip) 
     $   i, j, lwkopt, nb, iinfo
!
!     EXECUTABLE STATEMENTS:
!     Test the input arguments
!
      info = 0
      lquery = ( lwork.eq.-1 )
      upper = wlsame ( uplo, 'U' )
      if ( .not.upper .and. .not.wlsame ( uplo, 'L' ) ) then
         info = -1
      else if ( n.lt.0 ) then
         info = -2
      else if ( lda.lt.max( 1_wip, n ) ) then
         info = -4
      else if ( lwork.lt.max( 1_wip, n-1 ) .and. 
     $         .not.lquery ) then
         info = -7
      endif
!
      if ( info.eq.0 ) then
         if ( upper ) then
            nb = wilaenv ( 
     $           1, 'WFORGQL', ' ', n-1, n-1, n-1, -1_wip )
         else
            nb = wilaenv ( 
     $           1, 'WFORGQR', ' ', n-1, n-1, n-1, -1_wip )
         endif
         lwkopt = max( 1_wip, n-1 )*nb
         work(1) = lwkopt
      endif
!
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'WFORGTR', -info )
         return
      else if ( lquery ) then
         return
      endif
!
!     Quick return if possible
!
      if ( n.eq.0 ) then
         work(1) = 1
         return
      endif
!
      if ( upper ) then
!
!        Q was determined by a call to QPSYTRD with UPLO = 'U'
!
!        Shift the vectors which define the elementary reflectors one
!        column to the left, and set the last row and column of Q to
!        those of the unit matrix
!
         do j = 1, n-1
            do i = 1, j-1
               a(i,j) = a(i,j+1)
            enddo 
            a(n,j) = zero
         enddo 
         do i = 1, n-1
            a(i,n) = zero
         enddo 
         a(n,n) = one
!
!        Generate Q(1:n-1,1:n-1)
!
         call qporgql ( 
     $        n-1, n-1, n-1, a, lda, tau, work, lwork, iinfo )
!
      else
!
!        Q was determined by a call to QPSYTRD with UPLO = 'L'.
!
!        Shift the vectors which define the elementary reflectors one
!        column to the right, and set the first row and column of Q to
!        those of the unit matrix
!
         do j = n, 2, -1
            a(1,j) = zero
            do i = j+1, n
               a(i,j) = a(i,j-1)
            enddo 
         enddo 
         a(1,1) = one
         do i = 2, n
            a(i,1) = zero
         enddo 
         if ( n.gt.1 ) then
!
!           Generate Q(2:n,2:n)
!
            call qporgqr ( 
     $           n-1, n-1, n-1, a(2, 2), lda, tau, work,
     $           lwork, iinfo )
         endif
      endif
      work(1) = lwkopt
      return
      end subroutine qporgtr
!
!***********************************************************************
!
      subroutine qpsterf ( 
     $           n, d, e, info )
!
      implicit none
      include 'epcode_inc_qp.h'
!
!     ARGUMENTS: 
      integer(wip)  
     $   n
      real(wrp) 
     $   d(*), e(*)
      integer(wip)  
     $   info
! 
!  -- LAPACK routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     June 30, 1999
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) REAL(wrp)        array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) REAL(wrp)        array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm failed to find all of the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero.
!
!  =====================================================================
!
!     DEPENDENCES: 
      real(wrp)
     $   qplamch, 
     $   qplanst, 
     $   qplapy2
      external 
     $   qplamch, 
     $   qplanst, 
     $   qplapy2
      external 
     $   qplae2, 
     $   qplascl, 
     $   qplasrt, 
     $   wxerbla
      intrinsic 
     $   abs, sign, sqrt
!
!     LOCAL VARIABLES: 
      real(wrp) 
     $   zero, one, two, three
      parameter ( 
     $   zero  = 0.0_wrp, 
     $   one   = 1.0_wrp, 
     $   two   = 2.0_wrp,
     $   three = 3.0_wrp )
      integer
     $   maxit
      parameter ( 
     $   maxit = 30 )
      integer(wip)  
     $   i, iscale, jtot, l, l1, lend, 
     $   lendsv, lsv, m, nmaxit
      real(wrp)
     $   alpha, anorm, bb, c, eps, 
     $   eps2, rgamma, oldc, oldgam, 
     $   p, r, rt1, rt2, rte, s, 
     $   safmax, safmin, sigma, 
     $   ssfmax, ssfmin
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      info = 0
!
!     Quick return if possible
!
      if ( n .lt. 0 ) then
         info = -1
         call wxerbla ( 
     $        'wfsterf', -info )
         return
      endif
      if ( n .le. 1 ) return
!
!     Determine the unit roundoff for this environment.
!
      eps = qplamch ( 'E' )
      eps2 = eps**2
      safmin = qplamch ( 'S' )
      safmax = one / safmin
      ssfmax = sqrt( safmax ) / three
      ssfmin = sqrt( safmin ) / eps2
!
!     Compute the eigenvalues of the tridiagonal matrix.
!
      nmaxit = n*maxit
      sigma = zero
      jtot = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      l1 = 1
!
   10 continue
      if ( l1 .gt. n ) goto 170
      if ( l1 .gt. 1 ) e(l1-1) = zero
      do m = l1, n-1
         if ( abs(e(m)) .le. ( sqrt(abs(d(m  )))*
     $                         sqrt(abs(d(m+1))) )*eps ) then
            e(m) = zero
            goto 30
         endif
      enddo 
      m = n
!
   30 continue
      l = l1
      lsv = l
      lend = m
      lendsv = lend
      l1 = m + 1
      if ( lend .eq. l ) goto 10
!
!     Scale submatrix in rows and columns L to LEND
!
      anorm = qplanst ( 
     $        'I', lend-l+1, d(l), e(l) )
      iscale = 0
      if ( anorm .gt. ssfmax ) then
         iscale = 1
         call qplascl ( 
     $        'G', 0_wip, 0_wip, anorm, ssfmax, lend-l+1, 
     $        1_wip, d(l), n, info )
         call qplascl ( 
     $        'G', 0_wip, 0_wip, anorm, ssfmax, lend-l, 
     $        1_wip, e(l), n, info )
      else if ( anorm .lt. ssfmin ) then
         iscale = 2
         call qplascl ( 
     $        'G', 0_wip, 0_wip, anorm, ssfmin, lend-l+1, 
     $        1_wip, d(l), n, info )
         call qplascl ( 
     $        'G', 0_wip, 0_wip, anorm, ssfmin, lend-l, 
     $        1_wip, e(l), n, info )
      endif
!
      do i = l, lend-1
         e(i) = e(i)**2
      enddo 
!
!     Choose between QL and QR iteration
!
      if ( abs( d(lend) ) .lt. abs( d(l) ) ) then
         lend = lsv
         l = lendsv
      endif
!
      if ( lend .ge. l ) then
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   50    continue
         if ( l .ne. lend ) then
            do m = l, lend-1
               if ( abs(e(m)) .le. eps2 * abs(d(m)*d(m+1)) )
     $            goto 70
            enddo 
         endif
         m = lend
!
   70    continue
         if ( m .lt. lend ) e(m) = zero
         p = d(l)
         if ( m .eq. l ) goto 90
!
!        If remaining matrix is 2 by 2, use WFLAE2 to compute its
!        eigenvalues.
!
         if ( m .eq. l+1 ) then
            rte = sqrt( e(l) )
            call qplae2 ( 
     $           d(l), rte, d(l+1), rt1, rt2 )
            d(l)   = rt1
            d(l+1) = rt2
            e(l)   = zero
            l = l + 2
            if ( l .le. lend ) goto 50
            goto 150
         endif
!
         if ( jtot .eq. nmaxit ) goto 150
         jtot = jtot + 1
!
!        Form shift.
!
         rte = sqrt( e(l) )
         sigma = ( d(l+1) - p ) / ( two*rte )
         r = qplapy2 ( 
     $       sigma, one )
         sigma = p - ( rte / ( sigma + sign(r,sigma) ) )
!
         c = one
         s = zero
         rgamma = d(m) - sigma
         p = rgamma * rgamma 
!
!        Inner loop
!
         do i = m-1, l, -1
            bb = e(i)
            r = p + bb
            if ( i .ne. m-1 ) e(i+1) = s*r
            oldc = c
            c = p / r
            s = bb / r
            oldgam = rgamma 
            alpha = d(i)
            rgamma = c*( alpha - sigma ) - s*oldgam
            d(i+1) = oldgam + ( alpha - rgamma )
            if ( c .ne. zero ) then
               p = ( rgamma * rgamma ) / c
            else
               p = oldc*bb
            endif
         enddo 
!
         e(l) = s*p
         d(l) = sigma + rgamma 
         goto 50
!
!        Eigenvalue found.
!
   90    continue
         d(l) = p
!
         l = l + 1
         if ( l .le. lend ) goto 50
         goto 150
!
      else
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
  100    continue
         do m = l, lend+1, -1
            if ( abs(e(m-1)) .le. eps2 * abs(d(m)*d(m-1)) )
     $         goto 120
         enddo 
         m = lend
!
  120    continue
         if ( m .gt. lend ) e(m-1) = zero
         p = d(l)
         if ( m .eq. l ) goto 140
!
!        If remaining matrix is 2 by 2, use WFLAE2 to compute its
!        eigenvalues.
!
         if ( m .eq. l-1 ) then
            rte = sqrt( e(l-1) )
            call qplae2 ( 
     $           d(l), rte, d(l-1), rt1, rt2 )
            d(l)   = rt1
            d(l-1) = rt2
            e(l-1) = zero
            l = l - 2
            if ( l .ge. lend ) goto 100
            goto 150
         endif
!
         if ( jtot .eq. nmaxit ) goto 150
         jtot = jtot + 1
!
!        Form shift.
!
         rte = sqrt( e(l-1) )
         sigma = ( d(l-1) - p ) / ( two*rte )
         r = qplapy2 ( 
     $       sigma, one )
         sigma = p - ( rte / ( sigma + sign(r,sigma) ) )
!
         c = one
         s = zero
         rgamma = d(m) - sigma
         p = rgamma * rgamma 
!
!        Inner loop
!
         do i = m, l-1
            bb = e(i)
            r = p + bb
            if ( i .ne. m ) e(i-1) = s*r
            oldc = c
            c = p / r
            s = bb / r
            oldgam = rgamma
            alpha = d(i+1)
            rgamma = c*( alpha - sigma ) - s*oldgam
            d(i) = oldgam + ( alpha - rgamma )
            if ( c .ne. zero ) then
               p = ( rgamma*rgamma ) / c
            else
               p = oldc*bb
            endif
         enddo 
!
         e(l-1) = s*p
         d(l)   = sigma+rgamma
         goto 100
!
!        Eigenvalue found.
!
  140    continue
         d(l) = p
!
         l = l - 1
         if ( l .ge. lend ) goto 100
         goto 150
!
      endif
!
!     Undo scaling if necessary
!
  150 continue
      if ( iscale .eq. 1 )
     $   call qplascl ( 
     $        'G', 0_wip, 0_wip, ssfmax, anorm, lendsv-lsv+1, 1_wip,
     $        d(lsv), n, info )
      if ( iscale.eq.2 )
     $   call qplascl ( 
     $        'G', 0_wip, 0_wip, ssfmin, anorm, lendsv-lsv+1, 1_wip,
     $        d(lsv), n, info )
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      if ( jtot .lt. nmaxit ) goto 10
      do i = 1, n-1
         if ( e(i) .ne. zero ) info = info + 1
      enddo 
      goto 180
!
!     Sort eigenvalues in increasing order.
!
  170 continue
      call qplasrt ( 
     $     'I', n, d, info )
!
  180 continue
      return
      end subroutine qpsterf 
!
!***********************************************************************
!
      subroutine qpsymv ( 
     $           uplo, n, alpha, a, lda, x, incx,
     $           beta, y, incy )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      character(len=1) 
     $   uplo
      real(wrp) 
     $   alpha, beta
      integer(wip) 
     $   incx, incy, lda, n
      real(wrp) 
     $   a(lda,*), x(*), y(*)
!
!  Purpose
!  =======
!
!  WFSYMV  performs the matrix-vector  operation
!
!     y := alpha*A*x + beta*y,
!
!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL(wrp)       .
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - REAL(wrp)        array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - REAL(wrp)        array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - REAL(wrp)       .
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - REAL(wrp)        array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      logical
     $   wlsame
      external 
     $   wlsame
      external
     $   wxerbla
      intrinsic
     $   max
!
!     LOCAL VARIABLES:
      real(wrp)  
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      real(wrp)
     $   temp1, temp2
      integer(wip) 
     $   i, ix, iy, j, jx, jy, kx, ky,
     $   info
!    
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      info = 0
      if ( .not. wlsame ( uplo, 'U' ) .and.
     $     .not. wlsame ( uplo, 'L' )      ) then
         info = 1
      else if ( n.lt.0_wip ) then
         info = 2
      else if ( lda.lt.max( 1_wip, n ) ) then
         info = 5
      else if ( incx.eq.0_wip ) then
         info = 7
      else if ( incy.eq.0_wip ) then
         info = 10
      endif
      if ( info.ne.0 ) then
         call wxerbla ( 
     $        'wfsymv ', info )
         return
      endif
!
!     Quick return if possible.
!
      if ( ( n.eq.0_wip ) .or. 
     $   ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
!
!     Set up the start points in  X  and  Y.
!
      if ( incx.gt.0_wip ) then
         kx = 1_wip
      else
         kx = 1 - ( n - 1 )*incx
      endif
      if ( incy.gt.0_wip ) then
         ky = 1_wip
      else
         ky = 1 - ( n - 1 )*incy
      endif
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
!     First form  y := beta*y.
!
      if ( beta.ne.one ) then
         if ( incy.eq.1_wip ) then
            if ( beta.eq.zero ) then
               do i = 1, n
                  y(i) = zero
               enddo 
            else
               do i = 1, n
                  y(i) = beta*y(i)
               enddo 
            endif
         else
            iy = ky
            if ( beta.eq.zero ) then
               do i = 1, n
                  y(iy) = zero
                  iy    = iy + incy
               enddo 
            else
               do i = 1, n
                  y(iy) = beta*y(iy)
                  iy    = iy + incy
               enddo 
            endif
         endif
      endif
      if ( alpha .eq. zero ) return
      if ( wlsame ( uplo, 'U' ) ) then
!
!        Form  y  when A is stored in upper triangle.
!
         if ( ( incx.eq.1_wip ).and.( incy.eq.1_wip ) ) then
            do j = 1, n
               temp1 = alpha*x(j)
               temp2 = zero
               do i = 1, j - 1
                  y(i)  = y(i)  + temp1*a(i,j)
                  temp2 = temp2 + a(i,j)*x(i)
               enddo 
               y(j) = y(j) + temp1*a(j,j) + alpha*temp2
            enddo 
         else
            jx = kx
            jy = ky
            do j = 1, n
               temp1 = alpha*x(jx)
               temp2 = zero
               ix    = kx
               iy    = ky
               do i = 1, j - 1
                  y(iy) = y(iy) + temp1*a(i,j)
                  temp2 = temp2 + a(i,j)*x(ix)
                  ix    = ix + incx
                  iy    = iy + incy
               enddo 
               y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
               jx    = jx + incx
               jy    = jy + incy
            enddo 
         endif
      else
!
!        Form  y  when A is stored in lower triangle.
!
         if ( ( incx.eq.1_wip ).and.( incy.eq.1_wip ) ) then
            do j = 1, n
               temp1 = alpha*x(j)
               temp2 = zero
               y(j)  = y(j) + temp1*a(j,j)
               do i = j + 1, n
                  y(i)  = y(i)  + temp1*a(i,j)
                  temp2 = temp2 + a(i,j)*x(i)
               enddo 
               y(j) = y(j) + alpha*temp2
            enddo 
         else
            jx = kx
            jy = ky
            do j = 1, n
               temp1 = alpha*x(jx)
               temp2 = zero
               y(jy) = y(jy) + temp1*a(j,j)
               ix    = jx
               iy    = jy
               do i = j + 1, n
                  ix    = ix + incx
                  iy    = iy + incy
                  y(iy) = y(iy) + temp1*a(i,j)
                  temp2 = temp2 + a(i,j)*x(ix)
               enddo 
               y(jy) = y(jy) + alpha*temp2
               jx    = jx + incx
               jy    = jy + incy
            enddo 
         endif
      endif
      return
      end subroutine qpsymv
!
!***********************************************************************
!
      subroutine qptrmv ( 
     $           uplo, trans, diag, n, a, lda, x, incx )
!
      implicit none
      include 'epcode_inc_qp.h'
! 
!     ARGUMENTS:
      character(len=1) 
     $   uplo, trans, diag 
      integer(wip) 
     $   n, lda, incx
      real(wrp) 
     $   a(lda,*), x(*)
!
!  Purpose
!  =======
!
!  WFTRMV  performs one of the matrix-vector operations
!
!     x := A*x,   or   x := A'*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - REAL(wrp)        array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - REAL(wrp)        array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!     Tran Quoc Viet, a portable blas with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!     DEPENDENCES:
      logical 
     $   wlsame
      external 
     $   wlsame
      external 
     $   wxerbla
      intrinsic 
     $   max
!
!     LOCAL VARIABLES:
      real(wrp) 
     $   zero
      parameter ( 
     $   zero = 0.0_wrp )
      real(wrp)   
     $   temp
      integer(wip) 
     $   i, ix, j, jx, kx
      integer(wip)
     $   info
      logical
     $   nounit
!
!     EXECUTABLE STATEMENTS:
!     Test the input parameters.
!
      info = 0
      if ( .not. wlsame ( uplo, 'U' ) .and.
     $     .not. wlsame ( uplo, 'L' )      ) then
         info = 1
      else if ( .not. wlsame ( trans, 'N' ) .and.
     $          .not. wlsame ( trans, 'T' ) .and.
     $          .not. wlsame ( trans, 'C' )      ) then
         info = 2
      else if ( .not. wlsame ( diag, 'U' ) .and.
     $          .not. wlsame ( diag, 'N' )      ) then
         info = 3
      else if ( n.lt.0_wip ) then
         info = 4
      else if ( lda.lt.max( 1_wip, n ) ) then
         info = 6
      else if ( incx.eq.0_wip ) then
         info = 8
      endif
      if ( info.ne.0_wip ) then
         call wxerbla ( 
     $        'wftrmv ', info )
         return
      endif
!
!     Quick return if possible.
!
      if ( n.eq.0_wip ) return
!
      nounit = wlsame ( diag, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      if ( incx.le.0_wip ) then
         kx = 1 - ( n - 1 )*incx
      else if ( incx.ne.1 ) then
         kx = 1
      endif
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      if ( wlsame ( trans, 'N' ) ) then
!
!        Form  x := A*x.
!
         if ( wlsame ( uplo, 'U' ) ) then
            if ( incx.eq.1_wip ) then
               do j = 1, n
                  if ( x(j) .ne. zero ) then
                     temp = x(j)
                     do i = 1, j - 1
                        x(i) = x(i) + temp*a(i,j)
                     enddo 
                     if ( nounit ) x(j) = x(j)*a(j,j)
                  endif
               enddo 
            else
               jx = kx
               do j = 1, n
                  if ( x(jx) .ne. zero ) then
                     temp = x(jx)
                     ix   = kx
                     do i = 1, j - 1
                        x(ix) = x(ix) + temp*a(i,j)
                        ix    = ix + incx
                     enddo
                     if ( nounit ) x(jx) = x(jx)*a(j,j)
                  endif
                  jx = jx + incx
               enddo 
            endif
         else
            if ( incx .eq. 1_wip ) then
               do j = n, 1, -1
                  if ( x(j) .ne. zero ) then
                     temp = x(j)
                     do i = n, j + 1, -1
                        x(i) = x(i) + temp*a(i,j)
                     enddo 
                     if ( nounit ) x(j) = x(j)*a(j,j)
                  endif
               enddo 
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do j = n, 1, -1
                  if ( x(jx) .ne. zero ) then
                     temp = x(jx)
                     ix   = kx
                     do i = n, j + 1, -1
                        x(ix) = x(ix) + temp*a(i,j)
                        ix    = ix - incx
                     enddo 
                     if ( nounit ) x(jx) = x(jx)*a(j,j)
                  endif
                  jx = jx - incx
               enddo 
            endif
         endif
      else
!
!        Form  x := A'*x.
!
         if ( wlsame ( uplo, 'U' ) ) then
            if ( incx.eq.1_wip ) then
               do j = n, 1, -1
                  temp = x(j)
                  if ( nounit ) temp = temp*a(j,j)
                  do i = j - 1, 1, -1
                     temp = temp + a(i,j)*x(i)
                  enddo 
                  x(j) = temp
               enddo 
            else
               jx = kx + ( n - 1 )*incx
               do j = n, 1, -1
                  temp = x(jx)
                  ix   = jx
                  if ( nounit ) temp = temp*a(j,j)
                  do i = j - 1, 1, -1
                     ix   = ix - incx
                     temp = temp + a(i,j)*x(ix)
                  enddo 
                  x(jx) = temp
                  jx    = jx - incx
               enddo 
            endif
         else
            if ( incx.eq.1_wip ) then
               do j = 1, n
                  temp = x(j)
                  if ( nounit ) temp = temp*a(j,j)
                  do i = j + 1, n
                     temp = temp + a(i,j)*x(i)
                  enddo 
                  x(j) = temp
               enddo 
            else
               jx = kx
               do j = 1, n
                  temp = x(jx)
                  ix   = jx
                  if ( nounit ) temp = temp*a(j,j)
                  do i = j + 1, n
                     ix   = ix + incx
                     temp = temp + a(i,j)*x(ix)
                  enddo 
                  x(jx) = temp
                  jx    = jx + incx
               enddo 
            endif
         endif
      endif
      return
      end subroutine qptrmv
!
!***********************************************************************
!
      function qplansy ( 
     $         norm, uplo, n, a, lda, work )
!
      implicit none
      include 'epcode_inc_qp.h'
!
      real(wrp)          
     $   qplansy
!
!     ARGUMENTS:
      character          
     $   norm, uplo
      integer(wip)       
     $   lda, n
      real(wrp)          
     $   a(lda,*), work(*)
! 
!  -- LAPACK auxiliary routine (version 3.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!     Tran Quoc Viet, a portable lapack with flexible precision on demand.
!     +  2024/03/30, customizing precision for arguments.
!
!  Purpose
!  =======
!
!  WFLANSY  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric matrix A.
!
!  Description
!  ===========
!
!  WFLANSY returns the value
!
!    WFLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'M'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'O'
!              (
!              ( normI(A),         NORM = 'I' or 'I'
!              (
!              ( normF(A),         NORM = 'F', 'F', 'E' or 'E'
!
!  where  
!     norm1  denotes the  one norm of a matrix (maximum column sum),
!     normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!     normF  denotes the  Frobenius norm of a matrix (square root of sum of
!            squares).  
!  Note that  
!     max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in WFLANSY as described
!          above.
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is to be referenced.
!          = 'U':  Upper triangular part of A is referenced
!          = 'L':  Lower triangular part of A is referenced
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, WFLANSY is
!          set to zero.
!
!  A       (input) REAL(wrp)        array, dimension (LDA,N)
!          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) REAL(wrp)        array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
!          WORK is not referenced.
!
! =====================================================================
!
!     DEPENDENCES:
      external           
     $   qplassq
      logical            
     $   wlsame
      external           
     $   wlsame
      intrinsic          
     $   abs, max, sqrt
!
!     LOCAL VARIABLES:
      real(wrp)          
     $   one, zero
      parameter ( 
     $   one  = 1.0_wrp, 
     $   zero = 0.0_wrp )
      integer(wip)   
     $   i, j
      real(wrp)          
     $   absa, rscale, rsum, rval 
!
!     EXECUTABLE STATEMENTS:
!
      if ( n.eq.0 ) then
         rval = zero
      else if ( wlsame ( norm, 'M' ) ) then
!
!        Find max(abs(A(i,j))).
!
         rval = zero
         if ( wlsame ( uplo, 'U' ) ) then
            do j = 1, n
               do i = 1, j
                  rval = max( rval, abs( a(i,j) ) )
               enddo 
            enddo 
         else
            do j = 1, n
               do i = j, n
                  rval = max( rval, abs( a(i,j) ) )
               enddo 
            enddo 
         endif
      else if ( ( wlsame ( norm, 'I' ) ) .or. 
     $          ( wlsame ( norm, 'O' ) ) .or.
     $          ( norm.eq.'1' )            ) then
!
!        Find normI(A) ( = norm1(A), since A is symmetric).
!
         rval = zero
         if ( wlsame ( uplo, 'U' ) ) then
            do j = 1, n
               rsum = zero
               do i = 1, j-1
                  absa = abs( a(i,j) )
                  rsum = rsum+absa
                  work(i) = work(i)+absa
               enddo 
               work(j) = rsum+abs( a(j,j) )
            enddo 
            do i = 1, n
               rval = max( rval, work(i) )
            enddo 
         else
            do i = 1, n
               work(i) = zero
            enddo 
            do j = 1, n
               rsum = work(j)+abs( a(j,j) )
               do i = j+1, n
                  absa = abs( a(i,j) )
                  rsum = rsum+absa
                  work(i) = work(i)+absa
               enddo 
               rval = max( rval, rsum )
            enddo 
         endif
      else if ( ( wlsame ( norm, 'F' ) ) .or. 
     $          ( wlsame ( norm, 'E' ) )     ) then
!
!        Find normF(A).
!
         rscale = zero
         rsum = one
         if ( wlsame ( uplo, 'U' ) ) then
            do j = 2, n
               call qplassq ( 
     $              j-1, a(1,j), 1_wip, rscale, rsum )
            enddo 
         else
            do j = 1, n-1
               call qplassq ( 
     $              n-j, a(j+1,j), 1_wip, rscale, rsum )
            enddo 
         endif
         rsum = 2*rsum 
         call qplassq ( 
     $        n, a, lda+1, rscale, rsum )
         rval = rscale*sqrt( rsum )
      endif
!
      qplansy = rval 
      return
      end function qplansy 
!
!-----------------------------------------------------------------------
!
!@@end of DSYEV
!
!-----------------------------------------------------------------------
!
