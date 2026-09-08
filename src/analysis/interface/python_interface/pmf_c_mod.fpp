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
  use pm_analyze_mod

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
    real(wp), allocatable  :: pmf1d(:,:), pmf2d(:,:)
    integer                :: i, j

    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    nullify(result_pmf)
    n_out1 = 0
    n_out2 = 0

    write(MsgOut,'(A)') '[STEP1] Read Control Parameters for Analysis'
    write(MsgOut,'(A)') ' '
    call control_from_string(ctrl_text, ctrl_len, ctrl_data)

    write(MsgOut,'(A)') '[STEP2] Set Relevant Variables and Structures'
    write(MsgOut,'(A)') ' '
    call setup(ctrl_data, option, input, output)

    if (option%check_only) return

    write(MsgOut,'(A)') '[STEP3] Analysis trajectory files'
    write(MsgOut,'(A)') ' '
    call analyze_pmf_unified(input, option, pmf1d, pmf2d)

    ! Python layout: 1-D -> rows (center, standard PMF, Gaussian PMF) per bin;
    ! 2-D -> (nbin_y, nbin_x) so that NumPy sees (nbin_x, nbin_y)
    if (option%dimension == 1) then
      n_out1 = size(pmf1d, 2)      ! number of bins (Python rows)
      n_out2 = 3
      allocate(result_pmf(n_out2, n_out1))
      do i = 1, n_out1
        result_pmf(1, i) = option%center(1, i)
        result_pmf(2, i) = pmf1d(1, i)
        result_pmf(3, i) = pmf1d(2, i)
      end do
    else if (option%dimension == 2) then
      n_out1 = size(pmf2d, 1)      ! nbin_x (Python rows)
      n_out2 = size(pmf2d, 2)      ! nbin_y (Python cols)
      allocate(result_pmf(n_out2, n_out1))
      do i = 1, n_out1
        do j = 1, n_out2
          result_pmf(j, i) = pmf2d(i, j)
        end do
      end do
    else
      call error_set(err, ERROR_DIM_NOT_SUPP, &
                     'pmf_analysis> dimension must be 1 or 2')
    end if

  end subroutine pmf_analysis_main

end module pmf_c_mod
