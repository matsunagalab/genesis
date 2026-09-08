!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  ma_main
!! @brief   MBAR analysis tool
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module mbar_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use mbar_analyze_mod

  use mbar_control_mod
  use mbar_option_str_mod
  use output_str_mod
  use input_str_mod
  use molecules_str_mod
  use fileio_control_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  !> State shared between mbar_analysis_c and its guarded body
  !! mbar_analysis_body (see run_guarded in error_mod).
  type :: t_mbar_analysis_ctx
    ! inputs
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    logical     :: return_weights = .false.
    ! outputs (fene_f / weights_f are handed to Python on success)
    real(wp), pointer :: fene_f(:,:) => null()
    real(wp), pointer :: weights_f(:,:) => null()
    integer :: n_replica = 0
    integer :: n_blocks = 0
    integer :: n_weight_replica = 0
    integer :: n_weight_step = 0
    ! error state
    type(s_error) :: err
  end type t_mbar_analysis_ctx

  private :: t_mbar_analysis_ctx

 contains
  subroutine mbar_analysis_c(ctrl_text, ctrl_len, return_weights, result_fene, &
                             n_replica, n_blocks, result_weights, &
                             n_weight_replica, n_weight_step, &
                             status, msg, msglen) &
        bind(C, name="mbar_analysis_c")
    use conv_f_c_util
    implicit none
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    integer(c_int), value :: return_weights
    type(c_ptr), intent(out) :: result_fene
    integer(c_int), intent(out) :: n_replica
    integer(c_int), intent(out) :: n_blocks
    type(c_ptr), intent(out) :: result_weights
    integer(c_int), intent(out) :: n_weight_replica
    integer(c_int), intent(out) :: n_weight_step
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    type(t_mbar_analysis_ctx), target :: c

    result_fene    = c_null_ptr
    result_weights = c_null_ptr
    c%ctrl_text_ptr  = c_loc(ctrl_text)
    c%ctrl_len       = ctrl_len
    c%return_weights = (return_weights /= 0)

    ! Run under the library-mode error guard (see run_guarded in error_mod)
    ! so a fatal error_msg becomes a catchable error rather than exit(1).
    call run_guarded(mbar_analysis_body, c, c%err, status, msg, msglen)

    n_replica        = c%n_replica
    n_blocks         = c%n_blocks
    n_weight_replica = c%n_weight_replica
    n_weight_step    = c%n_weight_step
    if (error_has(c%err)) return

    if (associated(c%fene_f))    result_fene    = c_loc(c%fene_f)
    if (associated(c%weights_f)) result_weights = c_loc(c%weights_f)
  end subroutine mbar_analysis_c

  !> Guarded body of mbar_analysis_c (see run_guarded in error_mod).
  subroutine mbar_analysis_body(ctx) bind(C, name="mbar_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_mbar_analysis_ctx), pointer :: c
    character(kind=c_char), pointer :: ctrl_text(:)

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%ctrl_text_ptr, ctrl_text, [c%ctrl_len])

    call mbar_analysis_main( &
        ctrl_text, c%ctrl_len, c%return_weights, &
        c%fene_f, c%n_replica, c%n_blocks, c%weights_f, &
        c%n_weight_replica, c%n_weight_step, c%err)
  end subroutine mbar_analysis_body

  subroutine mbar_analysis_main( &
          ctrl_text, ctrl_len, return_weights, result_fene, &
          n_replica, n_blocks, result_weights, &
          n_weight_replica, n_weight_step, err)
    implicit none
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    logical,                intent(in) :: return_weights
    real(wp), pointer, intent(out) :: result_fene(:,:)
    integer,           intent(out) :: n_replica
    integer,           intent(out) :: n_blocks
    real(wp), pointer, intent(out) :: result_weights(:,:)
    integer,           intent(out) :: n_weight_replica
    integer,           intent(out) :: n_weight_step
    type(s_error),     intent(inout) :: err

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_molecule)       :: molecule
    type(s_option)         :: option
    type(s_input)          :: input
    type(s_output)         :: output
    type(s_mbar_result)    :: res
    real(wp)               :: out_unit
    real(wp), allocatable  :: f_k2d(:,:)
    integer                :: i, j, nbrella, nrep_x, nrep_y

    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    nullify(result_fene)
    nullify(result_weights)
    n_replica = 0
    n_blocks = 0
    n_weight_replica = 0
    n_weight_step = 0

    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '
    call control_from_string(ctrl_text, ctrl_len, ctrl_data)

    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '
    call setup(ctrl_data, molecule, option, input, output)

    if (return_weights .and. option%nblocks /= 1) then
      call error_set(err, ERROR_BLOCK_NOT_SUPP, &
                     'Analyze> in-memory weights require nblocks = 1.')
      return
    end if

    if (.not. option%check_only) then

      write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
      write(MsgOut,'(A)') ' '
      call analyze_mbar_unified(molecule, input, option, res)

      ! free energies in the Python layout: (n_blocks, n_replica) in Fortran
      ! order, i.e. (n_replica, n_blocks) for NumPy; the values the CLI
      ! writes to fenefile
      out_unit = mbar_output_unit(option)
      nbrella  = size(res%f_k(1)%v)
      if (option%dimension == 1) then
        n_replica = option%num_replicas
        n_blocks  = option%nblocks
        allocate(result_fene(n_blocks, n_replica))
        do i = 1, nbrella
          result_fene(:, i) = [(res%f_k(j)%v(i)*out_unit, j=1,option%nblocks)]
        end do
      else
        if (option%nblocks > 1) then
          call error_set(err, ERROR_BLOCK_NOT_SUPP, &
                         'Output_MBar> n-block > 1 is not supported in 2D')
          return
        end if
        nrep_x = option%rest_nreplica(option%rest_func_no(1, 1))
        nrep_y = option%rest_nreplica(option%rest_func_no(2, 1))
        n_replica = nrep_x
        n_blocks  = nrep_y
        allocate(f_k2d(nrep_y, nrep_x))
        f_k2d = reshape(res%f_k(1)%v, (/nrep_y, nrep_x/))
        allocate(result_fene(n_blocks, n_replica))
        do j = 1, nrep_y
          result_fene(:, j) = [(f_k2d(j,i)*out_unit, i=1,nrep_x)]
        end do
        deallocate(f_k2d)
      end if

      ! MBAR weights of every sample (nstep, nbrella)
      if (return_weights) then
        if (.not. allocated(res%weight_k)) then
          call error_set(err, ERROR_NOT_SUPPORTED, &
                         'Analyze> weights are unavailable for this MBAR input type.')
          return
        end if
        n_weight_step    = size(res%weight_k, 1)
        n_weight_replica = size(res%weight_k, 2)
        allocate(result_weights(n_weight_step, n_weight_replica))
        result_weights(:,:) = res%weight_k(:,:)
      end if

    end if

    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '
    call dealloc_option(option)
    call dealloc_molecules_all(molecule)

  end subroutine mbar_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in MBAR_ANALYSIS
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, option, input, output)
    use mbar_control_mod
    use mbar_option_mod
    use mbar_option_str_mod
    use output_mod
    use input_mod
    use output_str_mod
    use input_str_mod
    use select_mod
    use molecules_mod
    use molecules_str_mod
    use fileio_grocrd_mod
    use fileio_grotop_mod
    use fileio_ambcrd_mod
    use fileio_prmtop_mod
    use fileio_psf_mod
    use fileio_pdb_mod
    use constants_mod
    implicit none


    ! formal arguments
    type(s_ctrl_data),       intent(in)    :: ctrl_data
    type(s_molecule),        intent(inout) :: molecule
    type(s_option),          intent(inout) :: option
    type(s_input),           intent(inout) :: input
    type(s_output),          intent(inout) :: output


    ! setup input
    !
    call setup_input(ctrl_data%inp_info, input)


    ! setup selection
    !
    call setup_selection(ctrl_data%sel_info, molecule)


    ! setup option
    !
    call setup_option(ctrl_data%opt_info, ctrl_data%sel_info, &
                      molecule, option)


    ! setup output
    !
    call setup_output(ctrl_data%out_info, output)


    return

  end subroutine setup

end module mbar_c_mod
