!> Reads XYZ file
!!
!! loads an external file containing a XYZ geometry
!! \param fileName: name of the XYZ geometry file
!! \return A Molecule data structure
!<
      function loadXYZ(fileName)
        use Parameters
        implicit none
        character (len=*), intent(in) :: fileName
        type(Molecule) :: loadXYZ
        character (len=2) :: AtomSymbol
        integer :: j, l, stat
C       file handle
        integer :: fXYZ = 101 

        logical :: isParameterized

        open(unit=fXYZ, status="old", action="read", iostat=stat,
     &       file=fileName)
        if (stat.ne.0) then 
          print *,"Problem loading geometry file ", fileName
          stop
        end if
C       First line says how many atoms there are
        read (unit=fXYZ, fmt=*) loadXYZ%NumAtoms
C       Second line may contain a comment, skip it
        read (unit=fXYZ, fmt=*)

C       Allocate memory for atoms and coordinate data
        call NewAtoms(loadXYZ%Atoms, loadXYZ%NumAtoms)

        do j=1,loadXYZ%NumAtoms
          read (unit=fXYZ, fmt=*) AtomSymbol, loadXYZ%Atoms(j)%Position 
C         Convert units from Angstroms to atomic units (Bohr)
          loadXYZ%Atoms(j)%Position = loadXYZ%Atoms(j)%Position
     &         * Angstrom
C         look up AtomSymbol to assign parameters
          isParameterized = .False.
          do l=1,numParameterizedAtoms
            if (AtomSymbol.eq.ParameterizedAtoms(l)%Symbol) then
              loadXYZ%Atoms(j)%Element = ParameterizedAtoms(l)
              isParameterized = .True.
            end if
          end do

C         assign basis set
          if (isParameterized) then
C           Assign a Gaussian basis
            call AssignsGTOBasis(loadXYZ%Atoms(j))
            call AssignFittedGTOBasis(loadXYZ%Atoms(j))
C           Replace with this line to assign STO
c            call AssignsSTOBasis(loadXYZ%Atoms(j))
          else
             print *, "Error: Unknown element type: ", AtomSymbol
             stop
          end if
        end do
c       By default, assign zero total charge
        ! loadXYZ%TotalCharge = 1.01d0
        loadXYZ%TotalCharge = 0.0d0
        close(fXYZ)        
      end function loadXYZ

!>
!! Read in parameters from an external file
!<
      subroutine UpdateParameters(filename, Mol)
        use Parameters
        implicit none
        character (len=*), intent(in) :: fileName
        type(AtomData), dimension(:), allocatable :: ParameterSet
        integer :: i, j, N, stat

        type(Molecule), intent(inout), optional :: Mol
        logical :: isParameterized
        double precision, dimension(:), allocatable ::
     &     CustomGaussianExponent
        double precision, external:: sGTOFromHardness

C       file handle
        integer :: fPar = 1002

        open(unit=fPar, status="old", action="read", iostat=stat,
     &       file=fileName)
        if (stat.ne.0) then 
          print *,"Problem loading parameter file ", fileName
          stop
        end if

C       First line says how many Parameters there are
        read (unit=fPar, fmt=*) N

c       Allocate
        if (allocated(ParameterSet)) deallocate(ParameterSet)
        allocate(ParameterSet(N))
        
        call NewVector(CustomGaussianExponent,N)

        do i=1,N
          read (unit=fPar, fmt=*) ParameterSet(i)%Symbol,
     &      ParameterSet(i)%Z, ParameterSet(i)%FormalCharge,
     &      ParameterSet(i)%Electronegativity, ParameterSet(i)%Hardness,
     &      CustomGaussianExponent(i)
c     Assume units of electron volts are specified
          ParameterSet(i)%Electronegativity =
     ,         ParameterSet(i)%Electronegativity * eV
          ParameterSet(i)%Hardness = ParameterSet(i)%Hardness * eV
c     If exponent specified is zero, then calculate it automatically from hardness relation
          if (abs(CustomGaussianExponent(i)) .lt. 1.0d-16) then
            CustomGaussianExponent(i) = 
     &          sGTOFromHardness(ParameterSet(i)%Hardness)
            print *, "Automatically generated Gaussian exponent"
            print *, ParameterSet(i)%Symbol,
     &         ParameterSet(i)%Z, ParameterSet(i)%FormalCharge,
     &         CustomGaussianExponent(i)

          end if
        end do

        close(fPar)

c     If Molecule is specified, update its parameters
        if (present(Mol)) then
           do i=1,Mol%NumAtoms

C     look up AtomSymbol to assign parameters
              isParameterized = .False.
              do j=1,N
                 if (Mol%Atoms(i)%Element%Symbol
     &                .eq.ParameterSet(j)%Symbol) then
                    Mol%Atoms(i)%Element = ParameterSet(j)
                    Mol%Atoms(i)%Basis%zeta = CustomGaussianExponent(j)
                    isParameterized = .True.
                 end if
              end do
              
C     assign basis set
c              if (isParameterized) then
c                 call AssignsGTOBasis(Mol%Atoms(i))
c              else
c                 print *, "Warning, could not parameterize atom", i
c              end if
           end do
        end if
      end subroutine UpdateParameters

!>
!! dumps QTPIE calculation results
!! \param Mol: molecule data structure
!! \param fileName Name of the log file to write or append to
!<
      subroutine WriteLog(Mol, fileName)
        use Parameters
        character (len=*), intent(in) :: fileName
        type(Molecule), intent(in) :: Mol
C       file handle
        integer :: fXYZ = 101, stat

        open(unit=fXYZ, status="new", action="write", iostat=stat,
     &       file=fileName)
        if (stat.ne.0) then
c	gfortran's code is 17, pgf90's is 208, ifort's is 10
          if ((stat.eq.208) .or. (stat.eq.17) .or. (stat.eq.10)) then
C           File already exists
            open(unit=fXYZ, status="old", action="write", iostat=stat,
     &           position="append", file=fileName)
          else
            print *,"Problem opening log file"
            print *,"Status code = ", stat
            stop
          end if
        end if

C       write out charges
        write (unit=fXYZ, fmt=1) Mol%Energy, Mol%Atoms(:)%Charge
 1     format(99999f10.5)

        close(fXYZ)
      end subroutine writelog

!>
!! dumps molecular geometry from QTPIE in XYZ formal
!! \param Mol: molecule data structure
!! \param fileName: Name of geometry file to write or append to 
!<
      subroutine WriteXYZ(Mol, fileName)
        use Parameters
        character (len=*), intent(in) :: fileName
        type(Molecule), intent(in) :: Mol
C       file handle
        integer :: fXYZ = 102, stat, j

        open(unit=fXYZ, status="new", action="write", iostat=stat,
     &       file=fileName)
        if (stat.ne.0) then
          if ((stat.eq.208) .or. (stat .eq. 17) .or. (stat.eq.10)) then
C           File already exists
            open(unit=fXYZ, status="old", action="write", iostat=stat,
     &           position="append", file=fileName)
          else
            print *,"Problem writing geometry file ",fileName
            print *,"Status code = ", stat
            stop
          end if
        end if
C       write file
        write (unit=fXYZ, fmt=*) Mol%NumAtoms
        write (unit=fXYZ, fmt=*) "Written by QTPIE : WriteXYZ()"
        do j=1,Mol%NumAtoms
          write (unit=fXYZ, fmt=2) Mol%Atoms(j)%Element%Symbol, 
     &            Mol%Atoms(j)%Position/Angstrom
        end do
        close(fXYZ)
 
 2     format(a2,3f15.10)
 
       end subroutine WriteXYZ

