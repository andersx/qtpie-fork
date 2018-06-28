
!>
!!    \param N : an integer specifying the size of the problem
!!    \param A : a real, positive definite NxN matrix
!!    \param b : a real vector with N elements
!!    \param x : (Output) solution to matrix equation
!<
      subroutine solver(N, A, b, x)
      implicit none
      integer, intent(in) :: N
      double precision, dimension(N,N), intent( in) :: A
      double precision, dimension(N)  , intent( in) :: b
      double precision, dimension(N)  , intent(inout) :: x

      double precision, external :: dnrm2
      external :: dcg
      double precision, external :: ConditionNumber

      integer i

      if (abs(x(1)) .lt. 1.0d-8) x = 0.0d0
!     Check that b is not zero vector, else return 0
      if (dnrm2(N,b,1) .gt. 1.0d-8) then
!     Use conjugate gradients routine
        call dcg(N, A, b, x)
!     Use LAPACK SVD-based solver
!      call lapack_svdsolver(N, A, b, x)
      else
        print *, 'WARNING, b = 0'
        x = 0.0d0
      end if
!      print *, ConditionNumber(N, A)
!      print *, "Condition number = ", ConditionNumber(N, A)
!      print *, "Norm of b = ", dnrm2(N,b,1)

!        print *, "b = "
!      do i=1,N
!        print *, b(i)
!        end do       
      end subroutine solver

      subroutine lapack_svdsolver(N, AA, b, x)
      implicit none

      integer, intent(in) :: N
      double precision, dimension(N,N) :: A
      double precision, dimension(N,N), intent( in) :: AA
      double precision, dimension(N)  , intent( in) :: b
      double precision, dimension(N)  , intent(inout) :: x

!     Used for LAPACK solver
      double precision, dimension(N) :: S !< Matrix of singular values
      integer :: Rank, stat, WorkSize
      double precision, dimension(:), allocatable :: WORK
      double precision, parameter :: RCond = 1.0d-8
      external :: dgelss

!     Avoids bug where LAPACK could overwrite matrix
      A = AA

      x = b
C     First find optimal workspace size
      allocate(WORK(1))
      call dgelss(N, N, 1, A, N, x, N, S, RCond, Rank, WORK, 
     &     -1, stat)
      WorkSize = WORK(1)
      deallocate(WORK)

      allocate(WORK(WorkSize))
      call dgelss(N, N, 1, A, N, x, N, S, RCond, Rank, WORK,
     &     WorkSize, stat)
      deallocate(WORK)
      end subroutine

!>
!! Use LAPACK routine to calculate condition number of the NxN matrix A
!<
      double precision function ConditionNumber(N, AA)
      implicit none

      integer, intent(in) :: N
      double precision, dimension(N,N), intent( in) :: AA
      double precision, dimension(N,N) :: A

      double precision, dimension(N) :: c

!     Used for LAPACK solver
      double precision, dimension(N) :: S !< Matrix of singular values
      integer :: Rank, stat, WorkSize
      double precision, dimension(:), allocatable :: WORK
      double precision, parameter :: RCond = 1.0d-8
      external :: dgelss

      integer :: i
      c = 0.0d0

c     Known bug: can mess around with the matrix A for some reason
c     This is a workaround
      A = AA
C     First find optimal workspace size
      allocate(WORK(1))
      call dgelss(N, N, 1, A, N, c, N, S, RCond, Rank, WORK, 
     &     -1, stat)
      WorkSize = WORK(1)
      deallocate(WORK)

      allocate(WORK(WorkSize))
      call dgelss(N, N, 1, A, N, c, N, S, RCond, Rank, WORK,
     &     WorkSize, stat)
      deallocate(WORK)

      if (Rank .eq. N) then
         ConditionNumber = S(1)/S(N)
      else
         print *, "Matrix found to be singular"
         ConditionNumber = S(1)/S(Rank)
      end if

      print *, "Singular values:"
      do i =1, N
        print *, S(i)
      end do

      end function ConditionNumber

!>
!!    Double precision conjugate gradient solver with Jacobi preconditioner
!!
!!     Solves the matrix problem Ax = b for x
!!    Implemented from Golub and van Loan's stuff
!!    \author Jiahao Chen
!!    \date   2008-01-28
!!    \param N : an integer specifying the size of the problem
!!    \param A : a real, positive definite NxN matrix
!!    \param b : a real vector with N elements
!!    \param x : (Output) solution to matrix equation
!!    On input, contains initial guess
!<
      subroutine dcg(N, A, b, x)
      implicit none
      integer, intent(in) :: N
      double precision, dimension(N,N), intent( in) :: A
      double precision, dimension(N)  , intent( in) :: b
      double precision, dimension(N)  , intent(inout) :: x

      integer :: max_k = 100000 !< maximum number of iterations
      double precision, parameter :: tol = 1.0d-7 !< Convergence tolerance

      integer :: k !< Iteration loop counter
      external :: Precondition

c     Residual vector, p, q, z
      double precision, dimension(N) :: r, p, q, z
      double precision :: alpha, norm, critical_norm, gamma, gamma0

c     BLAS routines
      double precision, external :: ddot, dnrm2
      external :: dcopy, daxpy, dgemv

c      logical, parameter :: Verbose = .True.
      logical, parameter :: Verbose = .False.

*     Termination criterion norm
      critical_norm = tol * dnrm2(N,b,1)

*     Calculate initial guess x from diagonal part P(A) x = b 
!     The secret code to want an initial guess calculated is to pass an
!     initial guess with the first entry equal to floating-point zero.
!     If not, we'll just use the pre-specified initial guess that's already in x

      if (abs(x(1)) .lt. 1.0d0-10) call Precondition(N,A,b,x)
      
*     Calculate residual r = b - Ax
c     r = b (Copy b into r)
      call dcopy(N,b,1,r,1)
c     r = r - Ax
      call dgemv('N',N,N,-1.0d0,A,N,x,1,1.0d0,r,1)

c     Calculate norm
      norm = dnrm2(N,r,1)
      if (Verbose) then
         print *, "Iteration",0,":",norm, norm/critical_norm
      end if

      if (norm.lt.critical_norm) goto 1

      do k=1,max_k
*        Generate preconditioned z from P(A) z = r
         call Precondition(N,A,r,z)
c         z = r
*        Propagate old vectors
         if (k.ne.1) gamma0 = gamma
c        gamma  = r . z
         gamma  = ddot(N,r,1,z,1)

         if (k.ne.1) then
c           p = z + gamma/gamma0 * p
c           With BLAS, first overwrite z,then copy result from z to p
c           z = z + gamma/gamma0 * p
            call daxpy(N,gamma/gamma0,p,1,z,1)
         end if
c        p = z
         call dcopy(N,z,1,p,1)
*        Form matrix-vector product
c        q = A p
         call dgemv('N',N,N,1.0d0,A,N,p,1,0.0d0,q,1)

*        Calculate step size
c        alpha = gamma / p.q
         alpha = gamma / ddot(N,p,1,q,1)

*        Propagate by step size
c        x = x + alpha * p
         call daxpy(N, alpha,p,1,x,1)
c        r = r - alpha * q
         call daxpy(N,-alpha,q,1,r,1)

*        Calculate new norm of residual
         norm = dnrm2(N,r,1)

*        If requested, print convergence information
         if (Verbose) then
            print *, "Iteration",k,":",norm, norm/critical_norm
         end if
*        Check termination criterion

c        Done if || r || < tol || b ||
         if (norm.lt.critical_norm) goto 1
      end do

c     Oops, reached maximum iterations without convergence
      print *, "dcg: maximum iterations reached."
      print *, "ERROR: Solution is not converged."
      stop
c     Finally, return the answer
 1    if (Verbose) then
         print *, "dcg: solution found with residual", norm
      end if
      end subroutine dcg

!>
!! Calculates the solution x of the approximate preconditioned problem 
!! \f[
!! \mathbf{P}(\mathbf{A}) \vec{x} = \vec{b}
!! \f]
!<
      subroutine Precondition(N,A,b,x)
      implicit none
      integer, intent(in) :: N
      double precision, dimension(N,N), intent( in) :: A
      double precision, dimension(N)  , intent( in) :: b
      double precision, dimension(N)  , intent(out) :: x

      integer :: i
c      double precision :: ReciprocalSumOfDiagonals, MatrixElement

c     Use NO preconditioning
c      x = b

c     Use Jacobi preconditioning
      do i=1,N
         if (abs(A(i,i)) .gt. 1.0d-10) then
         x(i) = b(i) / A(i,i)
         else
         print *, "Error: divide by zero!"
         print *, "Error in column", i,":", A(i,i)
         stop
        end if
      end do

!     Use Gauss-Siedel preconditioning
!      do i=N,1,-1
!         x(i) = b(i)
!         do j=i+1,N
!           x(i) = x(i) - x(j) * A(j,i)
!         end do
!         x(i) = x(i) / A(i,i)
!      end do

*     Add in exact solution for last column and last row
!     ( 0 v ) ( x ) = (   y v   )
!     ( v w ) ( y ) = ( x.v + wy) 
!      ReciprocalSumOfDiagonals = 0.0d0
!
!      do i = 1,N-1
!         ReciprocalSumOfDiagonals = ReciprocalSumOfDiagonals
!     &        + 1.0d0/A(i,i)
!      end do
!
!      x(N) = - b(N) / ReciprocalSumOfDiagonals
!
!      do i = 1,N-1
!         MatrixElement = 1.0d0/(ReciprocalSumOfDiagonals * A(i,i))
!         x(i) = x(i) + b(N) * MatrixElement
!         x(N) = x(N) + b(i) * MatrixElement
!      end do

*     Calculate initial guess from approximate inverse
*     W is the inverse of the preconditioning matrix
*     W = approximate inverse of A
c      call ApproximateInverse(A,W,N)
*     Form matrix-vector product x = W b
c      call dgemv('N',N,N,1.0d0,W,N,b,1,0.0d0,x,1)

      end subroutine
!>
!!    Calculates an approximate inverse to a matrix of the form
!!    \f[
!!    \mathbf{M}=\left(\begin{array}{cc}\mathbf{J} & 1\\
!!                   1 & 0\end{array}\right)
!!    \f]
!!
!!    The inverse is calculated by approximating J by its diagonal, in which
!!    case an exact inverse can be constructed.
!<
      subroutine ApproximateInverse(M, W, N)
      implicit none
      integer, intent(in) :: N !< Size of matrix
      double precision, dimension(N,N), intent(in) :: M !< Matrix to invert
      double precision, dimension(N,N), intent(out) :: W !< Approximate inverse matrix

      integer :: i!, j
      double precision :: ReciprocalSum
      double precision :: MatrixElement!, OffDiagonalMatrixElement
 
      W = 0.0d0

      ReciprocalSum = 0.0d0
      do i = 1,N-1
         ReciprocalSum = ReciprocalSum + 1.0d0/M(i,i)
      end do

      W(N,N) = -1.0d0/ReciprocalSum
      do i = 1,N-1
         MatrixElement = 1.0d0/(ReciprocalSum * M(i,i))
         W(i,N) = MatrixElement
         W(N,i) = MatrixElement

c This is an approximation to the approximate problem, replacing that which follows

         W(i,i) = 1.0d0/M(i,i)

c This code computes the exact solution to the approximation, but exhibits
c slower convergence when used as a preconditioner. Go figure.
c         do j = 1, i-1
c            OffDiagonalMatrixElement = -MatrixElement/M(j,j)
c            W(i,j) = OffDiagonalMatrixElement
c            W(j,i) = OffDiagonalMatrixElement
c         end do
c
c         W(i,i) = 0.0d0
c         do j = 1,N-1
c            W(i,i) = W(i,i) - W(i,j)
c         end do
      end do
      end subroutine
