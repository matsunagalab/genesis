!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Module   pmf_c_mod
!! @brief   C-callable wrapper for PMF analysis (Python interface)
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module pmf_c_mod
  use, intrinsic :: iso_c_binding
  use pmf_impl_mod

  use pm_control_mod
  use pm_setup_mod
  use pm_option_str_mod
  use output_str_mod
  use input_str_mod
  use error_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  !> State shared between pmf_analysis_c and its guarded body
  !! pmf_analysis_body (see run_guarded in error_mod).
  type :: t_pmf_analysis_ctx
    ! inputs
    type(c_ptr) :: ctrl_text_ptr = c_null_ptr
    integer     :: ctrl_len = 0
    ! outputs (pmf_f is handed to Python on success)
    real(wp), pointer :: pmf_f(:,:) => null()
    integer :: n_out1 = 0
    integer :: n_out2 = 0
    ! error state
    type(s_error) :: err
  end type t_pmf_analysis_ctx

  private :: t_pmf_analysis_ctx
  private :: pmf_analysis_body

contains

  subroutine pmf_analysis_c(ctrl_text, ctrl_len, result_pmf, n_out1, n_out2, &
                            status, msg, msglen) &
        bind(C, name="pmf_analysis_c")
    implicit none
    character(kind=c_char), intent(in), target :: ctrl_text(*)
    integer(c_int),         value         :: ctrl_len
    type(c_ptr),            intent(out)   :: result_pmf
    integer(c_int),         intent(out)   :: n_out1
    integer(c_int),         intent(out)   :: n_out2
    integer(c_int),         intent(out)   :: status
    character(kind=c_char), intent(out)   :: msg(*)
    integer(c_int),         value         :: msglen

    type(t_pmf_analysis_ctx), target :: c

    result_pmf = c_null_ptr
    c%ctrl_text_ptr = c_loc(ctrl_text)
    c%ctrl_len      = ctrl_len

    ! Run the analysis under the library-mode error guard so that a fatal
    ! error_msg (e.g. a missing cvfile) is turned into a catchable error
    ! instead of aborting the host process (see run_guarded in error_mod).
    call run_guarded(pmf_analysis_body, c, c%err, status, msg, msglen)

    n_out1 = c%n_out1
    n_out2 = c%n_out2
    if (error_has(c%err)) return

    if (associated(c%pmf_f)) result_pmf = c_loc(c%pmf_f)
  end subroutine pmf_analysis_c

  !> Guarded body of pmf_analysis_c (see run_guarded in error_mod).
  subroutine pmf_analysis_body(ctx) bind(C, name="pmf_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_pmf_analysis_ctx), pointer :: c
    character(kind=c_char), pointer :: ctrl_text(:)

    call c_f_pointer(ctx, c)
    call c_f_pointer(c%ctrl_text_ptr, ctrl_text, [c%ctrl_len])

    call pmf_analysis_main( &
        ctrl_text, c%ctrl_len, c%pmf_f, c%n_out1, c%n_out2, c%err)
  end subroutine pmf_analysis_body

  subroutine pmf_analysis_main( &
          ctrl_text, ctrl_len, result_pmf, n_out1, n_out2, err)
    implicit none
    character(kind=c_char), intent(in)    :: ctrl_text(*)
    integer,                intent(in)    :: ctrl_len
    real(wp), pointer,      intent(out)   :: result_pmf(:,:)
    integer,                intent(out)   :: n_out1
    integer,                intent(out)   :: n_out2
    type(s_error),          intent(inout) :: err

    ! local variables
    type(s_ctrl_data)      :: ctrl_data
    type(s_option)         :: option
    type(s_input)          :: input
    type(s_output)         :: output


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

    call setup(ctrl_data, option, input, output)


    ! [Step3] Analyze trajectory
    !
    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '

    call analyze(input, output, option, result_pmf, n_out1, n_out2, err)
    if (error_has(err)) return

  end subroutine pmf_analysis_main

end module pmf_c_mod
