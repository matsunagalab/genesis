!--------1---------2---------3---------4---------5---------6---------7---------8
! 
!> Program  wa_main
!! @brief   WHAM analysis tool
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module wham_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use wham_impl_mod

  use wa_control_mod
  use wa_option_str_mod
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

  !> State shared between wa_analysis_c and its guarded body wa_analysis_body
  !! (see run_guarded in error_mod).
  type :: t_wa_analysis_ctx
    ! inputs
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    ! outputs (pmf_f is handed to Python on success)
    real(wp), pointer :: pmf_f(:,:) => null()
    integer :: n_bins = 0
    integer :: n_bin_x = 0
    ! error state
    type(s_error) :: err
  end type t_wa_analysis_ctx

  private :: t_wa_analysis_ctx
  private :: wa_analysis_body

contains
  subroutine wa_analysis_c(ctrl_text, ctrl_len, result_pmf, n_bins, n_bin_x,    &
                           status, msg, msglen) &
        bind(C, name="wa_analysis_c")
    use conv_f_c_util
    implicit none
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int), value :: ctrl_len
    type(c_ptr), intent(out)    :: result_pmf
    integer(c_int), intent(out) :: n_bins
    integer(c_int), intent(out) :: n_bin_x
    integer(c_int),          intent(out) :: status
    character(kind=c_char),  intent(out) :: msg(*)
    integer(c_int),          value       :: msglen

    type(t_wa_analysis_ctx), target :: c

    result_pmf = c_null_ptr
    c%ctrl_text_ptr = c_loc(ctrl_text)
    c%ctrl_len      = ctrl_len

    ! Run the analysis under the library-mode error guard so that a fatal
    ! error_msg (e.g. a missing cvfile) is turned into a catchable error
    ! instead of aborting the host process (see run_guarded in error_mod).
    call run_guarded(wa_analysis_body, c, c%err, status, msg, msglen)

    n_bins  = c%n_bins
    n_bin_x = c%n_bin_x
    if (error_has(c%err)) return

    if (associated(c%pmf_f)) result_pmf = c_loc(c%pmf_f)
  end subroutine wa_analysis_c

  !> Guarded body of wa_analysis_c (see run_guarded in error_mod).
  subroutine wa_analysis_body(ctx) bind(C, name="wa_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_wa_analysis_ctx), pointer :: c
    character(kind=c_char), pointer :: ctrl_text(:)

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%ctrl_text_ptr, ctrl_text, [c%ctrl_len])

    call wa_analysis_main( &
        ctrl_text, c%ctrl_len, c%pmf_f, c%n_bins, c%n_bin_x, c%err)
  end subroutine wa_analysis_body

  subroutine wa_analysis_main( &
          ctrl_text, ctrl_len, result_pmf, n_bins, n_bin_x, err)
    implicit none
    character(kind=c_char), intent(in) :: ctrl_text(*)
    integer,                intent(in) :: ctrl_len
    real(wp), pointer, intent(out) :: result_pmf(:,:)
    integer,           intent(out) :: n_bins
    integer,           intent(out) :: n_bin_x
    type(s_error),                   intent(inout) :: err


    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_option)         :: option
    type(s_input)          :: input
    type(s_output)         :: output
    type(s_molecule)       :: molecule


    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.


    ! [Step1] Read control file
    !
    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '

    call control_from_string(ctrl_text, ctrl_len, ctrl_data)


    ! [Step2] Set relevant variables and structures 
    !
    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '

    call setup(ctrl_data, molecule, option, input, output)


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    ! call analyze(molecule, s_trajes_c, ana_period, input, output, option)
    call analyze(molecule, input, output, option, result_pmf, n_bins, n_bin_x, err)
    if (error_has(err)) return


    ! [Step4] Deallocate memory
    !
    write(MsgOut,'(A)') '[STEP4] Deallocate memory'
    write(MsgOut,'(A)') ' '

    call dealloc_option(option)
    call dealloc_molecules_all(molecule)
end subroutine wa_analysis_main

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup
  !> @brief        setup variables and structures in WHAM_ANALYSIS
  !! @authors      NT
  !! @param[in]    ctrl_data  : information of control parameters
  !! @param[inout] molecule   : molecule information
  !! @param[inout] option     : option information
  !! @param[inout] input      : input information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup(ctrl_data, molecule, option, input, output)
    use wa_control_mod
    use wa_option_mod
    use wa_option_str_mod
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

end module wham_c_mod

