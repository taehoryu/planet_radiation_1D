! DO NOT EDIT THIS FILE!!!
!
! This file is automatically generated by write_probin.py at
! compile-time.
!
! To add a runtime parameter, do so by editting the appropriate _parameters
! file.

! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module extern_probin_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  private

  logical, allocatable, public :: eos_input_is_constant
  !$acc declare create(eos_input_is_constant)
  real (kind=rt), allocatable, public :: eos_gamma
  !$acc declare create(eos_gamma)
  logical, allocatable, public :: eos_assume_neutral
  !$acc declare create(eos_assume_neutral)
  real (kind=rt), allocatable, public :: small_x
  !$acc declare create(small_x)
  logical, allocatable, public :: use_tables
  !$acc declare create(use_tables)
  logical, allocatable, public :: use_c12ag_deboer17
  !$acc declare create(use_c12ag_deboer17)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: eos_input_is_constant
  attributes(managed) :: eos_gamma
  attributes(managed) :: eos_assume_neutral
  attributes(managed) :: small_x
  attributes(managed) :: use_tables
  attributes(managed) :: use_c12ag_deboer17
#endif

end module extern_probin_module

subroutine runtime_init(name, namlen)

  use extern_probin_module
  use amrex_error_module, only: amrex_error

  implicit none

  integer :: namlen
  integer :: name(namlen)
  
  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin


  namelist /extern/ eos_input_is_constant
  namelist /extern/ eos_gamma
  namelist /extern/ eos_assume_neutral
  namelist /extern/ small_x
  namelist /extern/ use_tables
  namelist /extern/ use_c12ag_deboer17

  allocate(eos_input_is_constant)
  allocate(eos_gamma)
  allocate(eos_assume_neutral)
  allocate(small_x)
  allocate(use_tables)
  allocate(use_c12ag_deboer17)

  eos_input_is_constant = .true.
  eos_gamma = 5.d0/3.d0
  eos_assume_neutral = .true.
  small_x = 1.d-30
  use_tables = .false.
  use_c12ag_deboer17 = .false.


  ! create the filename
  if (namlen > maxlen) then
     call amrex_error("probin file name too long")
  endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do


  ! read in the namelist
  un = 9
  open (unit=un, file=probin(1:namlen), form='formatted', status='old')
  read (unit=un, nml=extern, iostat=status)

  if (status < 0) then
     ! the namelist does not exist, so we just go with the defaults
     continue

  else if (status > 0) then
     ! some problem in the namelist
     call amrex_error("ERROR: problem in the extern namelist")
  endif

  close (unit=un)

  !$acc update &
  !$acc device(eos_input_is_constant, eos_gamma, eos_assume_neutral) &
  !$acc device(small_x, use_tables, use_c12ag_deboer17)

end subroutine runtime_init

subroutine runtime_pretty_print(name, namlen) bind(C, name="runtime_pretty_print")

  use amrex_constants_module
  use extern_probin_module
  use amrex_error_module, only: amrex_error

  implicit none

  integer :: unit, i

  integer :: namlen
  integer :: name(namlen)

  logical :: ltest

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  if (namlen > maxlen) then
     call amrex_error("probin file name too long")
  endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! open the job info file in Fortran
  open(newunit=unit, file=probin(1:namlen), status="old", position="append")

100 format (1x, a3, 2x, a32, 1x, "=", 1x, a)
101 format (1x, a3, 2x, a32, 1x, "=", 1x, i10)
102 format (1x, a3, 2x, a32, 1x, "=", 1x, g20.10)
103 format (1x, a3, 2x, a32, 1x, "=", 1x, l)

  ltest = eos_input_is_constant .eqv. .true.
  write (unit,103) merge("   ", "[*]", ltest), &
 "eos_input_is_constant", eos_input_is_constant

  ltest = eos_gamma == 5.d0/3.d0
  write (unit,102) merge("   ", "[*]", ltest), &
 "eos_gamma", eos_gamma

  ltest = eos_assume_neutral .eqv. .true.
  write (unit,103) merge("   ", "[*]", ltest), &
 "eos_assume_neutral", eos_assume_neutral

  ltest = small_x == 1.d-30
  write (unit,102) merge("   ", "[*]", ltest), &
 "small_x", small_x

  ltest = use_tables .eqv. .false.
  write (unit,103) merge("   ", "[*]", ltest), &
 "use_tables", use_tables

  ltest = use_c12ag_deboer17 .eqv. .false.
  write (unit,103) merge("   ", "[*]", ltest), &
 "use_c12ag_deboer17", use_c12ag_deboer17

  close(unit=unit)

end subroutine runtime_pretty_print
