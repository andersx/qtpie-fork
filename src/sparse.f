      module SparseMatrix
      implicit none
!>
!! Data type for compressed sparse row matrix format
!! \note NO range bounds checking 
!<
        type CSRMatrix
          integer :: RowDim, ColDim
          integer, dimension(:), allocatable :: RowStart
          integer, dimension(:), allocatable :: ColIdx
          double precision, dimension(:), allocatable :: Value
        end type CSRMatrix

      contains

!>
!!    Initializes a CSR matrix
!<
      subroutine CSRNew(A, RowDim, ColDim, MaxNumVals)
      implicit none
      integer, intent(in) :: RowDim, ColDim, MaxNumVals
      type(CSRMatrix), intent(inout) :: A

      integer :: stat

      call CSRDelete(A)

      A%RowDim = RowDim
      A%ColDim = ColDim

      allocate(A%RowStart(RowDim+1), A%ColIdx(MaxNumVals),
     &         A%Value(MaxNumVals), STAT=stat)

      if (stat .ne. 0) then
         print *, "CSRNew: Error allocating new sparse matrix"
         print *, "Error code =",stat
         stop
      end if

      A%RowStart = 0

      end subroutine CSRNew

!>
!!    Deletes a CSR matrix
!<
      subroutine CSRDelete(A)
      implicit none
      type(CSRMatrix), intent(inout) :: A
      integer :: stat

      if (allocated(A%Value)) then
         deallocate(A%RowStart, A%ColIdx, A%Value, STAT=stat)

         if (stat .ne. 0) then
            print *, "CSRDelete: Error deallocating sparse matrix"
            print *, "Error code =",stat
            stop
         end if
      end if

      end subroutine CSRDelete

!>
!!    Prints a CSR matrix
!<
      subroutine CSRPrint(A, rowidx)
      implicit none
      type(CSRMatrix), intent(inout) :: A
      integer, optional :: rowidx

      integer :: i, low, upp

      if (present(rowidx)) then
         low = rowidx
         upp = rowidx
      else
         low = 1
         upp = A%RowDim
      end if
      do i = low, upp
         print *, "Row", i, "runs from",A%RowStart(i),
     &        "to",A%RowStart(i+1)-1
         print *, " "
         print *, "Column indices"
         print *, "--------------"
         print *, A%ColIdx(A%RowStart(i):A%RowStart(i+1)-1)
         print *, "Matrix elements"
         print *, "---------------"
         print *, A%Value(A%RowStart(i):A%RowStart(i+1)-1)
      end do
      
      end subroutine CSRPrint
!>
!! Does sparse matrix-vector multiplies
!<
      subroutine MatrixVectorMultiply(N, A, x, y)
      implicit none
      integer, intent(in) :: N
      type(CSRMatrix), intent(in) :: A
      double precision, dimension(N), intent(in) :: x
      double precision, dimension(N), intent(out) :: y

      integer :: i, j
      double precision :: ThisElement

      do i = 1, N
         ThisElement = 0.0d0
         do j = A%RowStart(i), A%RowStart(i+1) - 1
            ThisElement = ThisElement + A%Value(j) * x(A%ColIdx(j))
         end do
         y(i) = ThisElement
      end do
      
      end subroutine MatrixVectorMultiply

!>
!! Sums over the entire row of a CSR.
!<
      double precision function SumRow(A, RowIdx) 
      implicit none
      type(CSRMatrix), intent(in) :: A
      integer, intent(in) :: RowIdx !<Index to find sum of

      integer :: i

      SumRow = 0.0d0
      do i = A%RowStart(RowIdx), A%RowStart(RowIdx+1) - 1
         SumRow = SumRow + A%Value(i)
      end do

      end function SumRow

      end module SparseMatrix
