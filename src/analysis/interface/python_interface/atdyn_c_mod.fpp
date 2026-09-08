!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Module  atdyn_c_mod
!! @brief  C-callable wrapper for atdyn MD engine
!! @authors Claude Code
!
!  (c) Copyright 2024 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module atdyn_c_mod
  use, intrinsic :: iso_c_binding
  use at_setup_atdyn_mod
  use at_setup_mpi_mod
  use at_dynamics_mod
  use at_minimize_mod
  use at_control_mod
  use at_energy_mod
  use at_energy_pme_mod
  use at_output_mod
  use at_dynvars_str_mod
  use at_dynamics_str_mod
  use at_minimize_str_mod
  use at_ensemble_str_mod
  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_enefunc_str_mod
  use at_energy_str_mod
  use at_constraints_str_mod
  use at_output_str_mod
  use at_remd_str_mod
  use at_rpath_str_mod
  use at_vibration_str_mod
  use at_morph_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use timers_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  use error_mod

  implicit none
  private

  public :: atdyn_md_c
  public :: atdyn_min_c
  public :: deallocate_atdyn_results_c
  public :: reset_atdyn_state_c

  ! Module-level pointers for results (to be deallocated later)
  real(wp), pointer, save :: energies_ptr(:,:) => null()
  real(wp), pointer, save :: final_coords_ptr(:,:) => null()

  ! Context shared between a bind(C) entry point and its guarded body (see
  ! run_guarded in error_mod). The entry point copies the control text in,
  ! runs the body under the library-mode error guard and copies the outputs
  ! back out. A GENESIS error_msg inside the engine therefore reaches Python
  ! as a status/message pair instead of terminating the process.
  type :: t_atdyn_ctx
    character(kind=c_char), allocatable :: ctrl_text(:)   ! input control text
    integer       :: ctrl_len       = 0
    integer       :: nframes        = 0                   ! outputs
    integer       :: nterms         = 0
    integer       :: natom          = 0
    integer       :: converged      = 0
    real(wp)      :: final_gradient = 0.0_wp
    type(s_error) :: err
  end type t_atdyn_ctx

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    atdyn_md_c
  !> @brief        Run MD simulation with control text (string)
  !! @authors      Claude Code
  !
  !  On failure status/msg carry the GENESIS error and every result output is
  !  null/zero. The engine data of a failed run is not released (the unwind
  !  skips the engine's own cleanup), so callers that need a clean process
  !  after a failure should use the *_isolated Python variants.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine atdyn_md_c(ctrl_text, ctrl_len, &
                        result_energies, result_nframes, result_nterms, &
                        result_final_coords, result_natom, &
                        status, msg, msglen) &
        bind(C, name="atdyn_md_c")

    implicit none

    ! Input - control text
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len

    ! Output - energies
    type(c_ptr), intent(out) :: result_energies
    integer(c_int), intent(out) :: result_nframes
    integer(c_int), intent(out) :: result_nterms

    ! Output - final coordinates
    type(c_ptr), intent(out) :: result_final_coords
    integer(c_int), intent(out) :: result_natom

    ! Status
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_atdyn_ctx), target :: c

    result_energies     = c_null_ptr
    result_nframes      = 0
    result_nterms       = 0
    result_final_coords = c_null_ptr
    result_natom        = 0

    call pack_ctrl_text(c, ctrl_text, ctrl_len)

    call run_guarded(atdyn_md_body, c, c%err, status, msg, msglen)
    if (status /= 0) return

    result_energies     = c_loc(energies_ptr)
    result_nframes      = c%nframes
    result_nterms       = c%nterms
    result_final_coords = c_loc(final_coords_ptr)
    result_natom        = c%natom

  end subroutine atdyn_md_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    atdyn_md_body
  !> @brief        Guarded body of atdyn_md_c
  !! @authors      Claude Code
  !! @param[in]    ctx : C pointer to the caller's t_atdyn_ctx
  !
  !  Module procedure by design: run_guarded takes c_funloc() of it, and an
  !  internal procedure would need a trampoline (see error_mod).
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine atdyn_md_body(ctx) bind(C, name="atdyn_md_body")

    implicit none

    ! Arguments
    type(c_ptr), value :: ctx

    ! Local variables
    type(t_atdyn_ctx), pointer :: c
    type(s_ctrl_data)      :: ctrl_data
    type(s_molecule)       :: molecule
    type(s_enefunc)        :: enefunc
    type(s_dynvars)        :: dynvars
    type(s_pairlist)       :: pairlist
    type(s_boundary)       :: boundary
    type(s_constraints)    :: constraints
    type(s_ensemble)       :: ensemble
    type(s_dynamics)       :: dynamics
    type(s_output)         :: output
    type(s_remd)           :: remd
    type(s_rpath)          :: rpath

    integer                :: nframes, nterms, natom
    integer                :: i
#ifdef OMP
    integer                :: omp_get_max_threads
#endif

    call c_f_pointer(ctx, c)

    if (c%ctrl_len <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     'atdyn_md_c: control text is empty')
      return
    end if

    ! Reset timers for clean state in library mode
    call reset_timers()

    ! Initialize MPI variables for non-MPI build
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.
    my_world_rank = 0
    nproc_world   = 1

    ! Initialize nthread (must be done before any setup calls)
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    ! Initialize timers
    call timer(TimerTotal, TimerOn)

    ! [Step1] Read control from string
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters'
    write(MsgOut,'(A)') ' '

    call control_md_from_string(c%ctrl_text, c%ctrl_len, ctrl_data)

    ! [Step2] Setup MPI (no-op for non-MPI)
    write(MsgOut,'(A)') '[STEP2] Setup MPI'
    write(MsgOut,'(A)') ' '
    call setup_mpi_md(ctrl_data%ene_info)

    ! [Step3] Setup simulation
    write(MsgOut,'(A)') '[STEP3] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup_atdyn_md(ctrl_data, output, molecule, enefunc, pairlist, &
                        dynvars, dynamics, constraints, ensemble, boundary)

    natom = molecule%num_atoms

    ! [Step4] Compute initial energy
    write(MsgOut,'(A)') '[STEP4] Compute Single Point Energy for Molecules'
    write(MsgOut,'(A)') ' '

    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        enefunc%nonb_limiter,  &
                        dynvars%coord,     &
                        dynvars%trans,     &
                        dynvars%coord_pbc, &
                        dynvars%energy,    &
                        dynvars%temporary, &
                        dynvars%force,     &
                        dynvars%force_omp, &
                        dynvars%virial,    &
                        dynvars%virial_extern, &
                        constraints)

    call output_energy(dynvars%step, enefunc, dynvars%energy)

    ! [Step5] Run MD
    write(MsgOut,'(A)') '[STEP5] Perform Molecular Dynamics Simulation'
    write(MsgOut,'(A)') ' '

    call timer(TimerDynamics, TimerOn)
    call run_md(output, molecule, enefunc, dynvars, dynamics, &
                pairlist, boundary, constraints, ensemble)
    call timer(TimerDynamics, TimerOff)

    ! Collect results
    nframes = 1  ! Just final frame for now
    nterms = 8   ! total, bond, angle, urey_bradley, dihedral, improper, electrostatic, vdw

    ! Deallocate previous results if any
    if (associated(energies_ptr)) deallocate(energies_ptr)
    if (associated(final_coords_ptr)) deallocate(final_coords_ptr)

    allocate(energies_ptr(nterms, nframes))
    energies_ptr(1, 1) = dynvars%energy%total
    energies_ptr(2, 1) = dynvars%energy%bond
    energies_ptr(3, 1) = dynvars%energy%angle
    energies_ptr(4, 1) = dynvars%energy%urey_bradley
    energies_ptr(5, 1) = dynvars%energy%dihedral
    energies_ptr(6, 1) = dynvars%energy%improper
    energies_ptr(7, 1) = dynvars%energy%electrostatic
    energies_ptr(8, 1) = dynvars%energy%van_der_waals

    allocate(final_coords_ptr(3, natom))
    do i = 1, natom
      final_coords_ptr(1, i) = dynvars%coord(1, i)
      final_coords_ptr(2, i) = dynvars%coord(2, i)
      final_coords_ptr(3, i) = dynvars%coord(3, i)
    end do

    c%nframes = nframes
    c%nterms  = nterms
    c%natom   = natom

    ! [Step6] Deallocate simulation data (but keep results)
    write(MsgOut,'(A)') ' '
    write(MsgOut,'(A)') '[STEP6] Deallocate Arrays'
    write(MsgOut,'(A)') ' '

    call dealloc_pme
    call dealloc_constraints_all(constraints)
    call dealloc_boundary_all(boundary)
    call dealloc_pairlist_all(pairlist)
    call dealloc_energy_all(dynvars%energy)
    call dealloc_dynvars_all(dynvars)
    call dealloc_enefunc_all(enefunc)
    call dealloc_molecules_all(molecule)

    call timer(TimerTotal, TimerOff)
    call output_time

  end subroutine atdyn_md_body

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    atdyn_min_c
  !> @brief        Run energy minimization with control text (string)
  !! @authors      Claude Code
  !
  !  Error reporting and cleanup follow atdyn_md_c.
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine atdyn_min_c(ctrl_text, ctrl_len, &
                         result_energies, result_nsteps, result_nterms, &
                         result_final_coords, result_natom, &
                         result_converged, result_final_gradient, &
                         status, msg, msglen) &
        bind(C, name="atdyn_min_c")

    implicit none

    ! Input - control text
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer(c_int), value :: ctrl_len

    ! Output - energies
    type(c_ptr), intent(out) :: result_energies
    integer(c_int), intent(out) :: result_nsteps
    integer(c_int), intent(out) :: result_nterms

    ! Output - final coordinates
    type(c_ptr), intent(out) :: result_final_coords
    integer(c_int), intent(out) :: result_natom

    ! Output - convergence info
    integer(c_int), intent(out) :: result_converged
    real(c_double), intent(out) :: result_final_gradient

    ! Status
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_atdyn_ctx), target :: c

    result_energies       = c_null_ptr
    result_nsteps         = 0
    result_nterms         = 0
    result_final_coords   = c_null_ptr
    result_natom          = 0
    result_converged      = 0
    result_final_gradient = 0.0_c_double

    call pack_ctrl_text(c, ctrl_text, ctrl_len)

    call run_guarded(atdyn_min_body, c, c%err, status, msg, msglen)
    if (status /= 0) return

    result_energies       = c_loc(energies_ptr)
    result_nsteps         = c%nframes
    result_nterms         = c%nterms
    result_final_coords   = c_loc(final_coords_ptr)
    result_natom          = c%natom
    result_converged      = c%converged
    result_final_gradient = c%final_gradient

  end subroutine atdyn_min_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    atdyn_min_body
  !> @brief        Guarded body of atdyn_min_c
  !! @authors      Claude Code
  !! @param[in]    ctx : C pointer to the caller's t_atdyn_ctx
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine atdyn_min_body(ctx) bind(C, name="atdyn_min_body")

    implicit none

    ! Arguments
    type(c_ptr), value :: ctx

    ! Local variables
    type(t_atdyn_ctx), pointer :: c
    type(s_ctrl_data)      :: ctrl_data
    type(s_molecule)       :: molecule
    type(s_enefunc)        :: enefunc
    type(s_dynvars)        :: dynvars
    type(s_pairlist)       :: pairlist
    type(s_boundary)       :: boundary
    type(s_minimize)       :: minimize
    type(s_output)         :: output

    integer                :: nsteps, nterms, natom
    integer                :: i
#ifdef OMP
    integer                :: omp_get_max_threads
#endif

    call c_f_pointer(ctx, c)

    if (c%ctrl_len <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     'atdyn_min_c: control text is empty')
      return
    end if

    ! Reset timers for clean state in library mode
    call reset_timers()

    ! Initialize MPI variables
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.
    my_world_rank = 0
    nproc_world   = 1

    ! Initialize nthread (must be done before any setup calls)
#ifdef OMP
    nthread = omp_get_max_threads()
#else
    nthread = 1
#endif

    ! Initialize timers
    call timer(TimerTotal, TimerOn)

    ! [Step1] Read control from string
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Minimization'
    write(MsgOut,'(A)') ' '

    call control_min_from_string(c%ctrl_text, c%ctrl_len, ctrl_data)

    ! [Step2] Setup MPI
    write(MsgOut,'(A)') '[STEP2] Setup MPI'
    write(MsgOut,'(A)') ' '
    call setup_mpi_md(ctrl_data%ene_info)

    ! [Step3] Setup simulation
    write(MsgOut,'(A)') '[STEP3] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup_atdyn_min(ctrl_data, output, molecule, enefunc, pairlist, &
                         dynvars, minimize, boundary)

    natom = molecule%num_atoms

    ! [Step4] Compute initial energy
    write(MsgOut,'(A)') '[STEP4] Compute Single Point Energy for Molecules'
    write(MsgOut,'(A)') ' '

    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                        enefunc%nonb_limiter,  &
                        dynvars%coord,     &
                        dynvars%trans,     &
                        dynvars%coord_pbc, &
                        dynvars%energy,    &
                        dynvars%temporary, &
                        dynvars%force,     &
                        dynvars%force_omp, &
                        dynvars%virial,    &
                        dynvars%virial_extern)

    call output_energy(dynvars%step, enefunc, dynvars%energy)

    ! [Step5] Run minimization
    write(MsgOut,'(A)') '[STEP5] Perform Energy Minimization'
    write(MsgOut,'(A)') ' '

    call timer(TimerDynamics, TimerOn)
    call run_min(output, molecule, enefunc, dynvars, minimize, &
                 pairlist, boundary)
    call timer(TimerDynamics, TimerOff)

    ! Collect results
    nsteps = 1
    nterms = 8

    ! Deallocate previous results if any
    if (associated(energies_ptr)) deallocate(energies_ptr)
    if (associated(final_coords_ptr)) deallocate(final_coords_ptr)

    allocate(energies_ptr(nterms, nsteps))
    energies_ptr(1, 1) = dynvars%energy%total
    energies_ptr(2, 1) = dynvars%energy%bond
    energies_ptr(3, 1) = dynvars%energy%angle
    energies_ptr(4, 1) = dynvars%energy%urey_bradley
    energies_ptr(5, 1) = dynvars%energy%dihedral
    energies_ptr(6, 1) = dynvars%energy%improper
    energies_ptr(7, 1) = dynvars%energy%electrostatic
    energies_ptr(8, 1) = dynvars%energy%van_der_waals

    allocate(final_coords_ptr(3, natom))
    do i = 1, natom
      final_coords_ptr(1, i) = dynvars%coord(1, i)
      final_coords_ptr(2, i) = dynvars%coord(2, i)
      final_coords_ptr(3, i) = dynvars%coord(3, i)
    end do

    c%nframes        = nsteps
    c%nterms         = nterms
    c%natom          = natom
    c%converged      = 0  ! TODO: get actual convergence status
    c%final_gradient = dynvars%rms_gradient

    ! [Step6] Deallocate simulation data
    write(MsgOut,'(A)') ' '
    write(MsgOut,'(A)') '[STEP6] Deallocate Arrays'
    write(MsgOut,'(A)') ' '

    call dealloc_pme
    call dealloc_boundary_all(boundary)
    call dealloc_pairlist_all(pairlist)
    call dealloc_energy_all(dynvars%energy)
    call dealloc_dynvars_all(dynvars)
    call dealloc_enefunc_all(enefunc)
    call dealloc_molecules_all(molecule)

    call timer(TimerTotal, TimerOff)
    call output_time

  end subroutine atdyn_min_body

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    pack_ctrl_text
  !> @brief        Copy the C control text into the context
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine pack_ctrl_text(c, ctrl_text, ctrl_len)

    implicit none

    type(t_atdyn_ctx),      intent(inout) :: c
    character(kind=c_char), intent(in)    :: ctrl_text(*)
    integer(c_int),         intent(in)    :: ctrl_len

    integer :: i

    c%ctrl_len = max(0, int(ctrl_len))
    allocate(c%ctrl_text(max(c%ctrl_len, 1)))
    c%ctrl_text(1) = c_null_char
    do i = 1, c%ctrl_len
      c%ctrl_text(i) = ctrl_text(i)
    end do

  end subroutine pack_ctrl_text

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_atdyn_results_c
  !> @brief        Deallocate result arrays
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine deallocate_atdyn_results_c() bind(C, name="deallocate_atdyn_results_c")
    implicit none

    if (associated(energies_ptr)) then
      deallocate(energies_ptr)
      nullify(energies_ptr)
    end if
    if (associated(final_coords_ptr)) then
      deallocate(final_coords_ptr)
      nullify(final_coords_ptr)
    end if

  end subroutine deallocate_atdyn_results_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    reset_atdyn_state_c
  !> @brief        Reset atdyn global state for multiple sequential runs
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine reset_atdyn_state_c() bind(C, name="reset_atdyn_state_c")
    implicit none

    ! Reset timer state
    call reset_timers()

    ! Deallocate any previous results
    if (associated(energies_ptr)) then
      deallocate(energies_ptr)
      nullify(energies_ptr)
    end if
    if (associated(final_coords_ptr)) then
      deallocate(final_coords_ptr)
      nullify(final_coords_ptr)
    end if

  end subroutine reset_atdyn_state_c

end module atdyn_c_mod
