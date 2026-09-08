!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   diffusion_c_mod
!> @brief   Python (bind(C)) entry point of diffusion_analysis
!! @authors Claude Code
!
!  The fit itself lives in da_analyze_mod (shared with the CLI program).
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module diffusion_c_mod
  use, intrinsic :: iso_c_binding
  use constants_mod
  use da_analyze_mod
  use error_mod
  use messages_mod
  implicit none

  public :: diffusion_analysis_c

  !> State shared between diffusion_analysis_c and its guarded body.
  type :: t_diffusion_ctx
    ! inputs (copied from the bind(C) arguments)
    type(c_ptr)    :: msd_ptr = c_null_ptr
    integer        :: ndata = 0
    integer        :: ncols = 0
    real(c_double) :: time_step = 1.0_c_double
    real(c_double) :: distance_unit = 1.0_c_double
    integer        :: ndofs = 0
    integer        :: start_step = 0
    integer        :: stop_step = 0
    type(c_ptr)    :: out_data_ptr = c_null_ptr
    integer        :: out_data_size = 0
    type(c_ptr)    :: diff_coeff_ptr = c_null_ptr
    integer        :: n_sets = 0
    ! error state
    type(s_error)  :: err
  end type t_diffusion_ctx

  private :: t_diffusion_ctx

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    diffusion_analysis_c
  !> @brief        Diffusion analysis (least-squares fit of MSD data)
  !! @authors      Claude Code
  !! @param[in]    msd_ptr        : (ncols, ndata) data: time and MSD sets
  !! @param[in]    time_step      : time unit of column 1 (ps)
  !! @param[in]    distance_unit  : length unit of the MSD columns (angstrom)
  !! @param[in]    ndofs          : degrees of freedom of every MSD set
  !! @param[in]    start_step     : first data index of the fit range
  !! @param[in]    stop_step      : last data index of the fit range
  !! @param[in]    out_data_ptr   : pre-allocated (2*n_sets+1, ndata) results
  !! @param[in]    diff_coeff_ptr : pre-allocated (n_sets) coefficients (cm^2/s)
  !! @param[out]   status, msg    : error status and message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine diffusion_analysis_c(msd_ptr, ndata, ncols, &
                                  time_step, distance_unit, ndofs, &
                                  start_step, stop_step, &
                                  out_data_ptr, out_data_size, &
                                  diff_coeff_ptr, n_sets, &
                                  status, msg, msglen) &
        bind(C, name="diffusion_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: msd_ptr
    integer(c_int), value :: ndata
    integer(c_int), value :: ncols
    real(c_double), value :: time_step
    real(c_double), value :: distance_unit
    integer(c_int), value :: ndofs
    integer(c_int), value :: start_step
    integer(c_int), value :: stop_step
    type(c_ptr), value :: out_data_ptr
    integer(c_int), value :: out_data_size
    type(c_ptr), value :: diff_coeff_ptr
    integer(c_int), value :: n_sets
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(t_diffusion_ctx), target :: c

    c%msd_ptr        = msd_ptr
    c%ndata          = ndata
    c%ncols          = ncols
    c%time_step      = time_step
    c%distance_unit  = distance_unit
    c%ndofs          = ndofs
    c%start_step     = start_step
    c%stop_step      = stop_step
    c%out_data_ptr   = out_data_ptr
    c%out_data_size  = out_data_size
    c%diff_coeff_ptr = diff_coeff_ptr
    c%n_sets         = n_sets

    call run_guarded(diffusion_analysis_body, c, c%err, status, msg, msglen)

  end subroutine diffusion_analysis_c

  !> Guarded body of diffusion_analysis_c (see run_guarded in error_mod).
  subroutine diffusion_analysis_body(ctx) bind(C, name="diffusion_analysis_body")
    implicit none
    type(c_ptr), value :: ctx

    type(t_diffusion_ctx), pointer :: c
    real(wp), pointer :: msd_f(:,:)
    real(wp), pointer :: out_data_f(:,:)
    real(wp), pointer :: diff_coeff_f(:)
    integer, allocatable :: ndofs(:)
    type(s_fitting_result), allocatable :: fittings(:)
    integer :: iset

    call c_f_pointer(ctx, c)

    ! Validate inputs
    if (.not. c_associated(c%msd_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: msd_ptr is null")
      return
    end if
    if (c%ndata <= 0) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: ndata must be positive")
      return
    end if
    if (c%ncols < 2) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: ncols must be >= 2")
      return
    end if
    if (c%n_sets /= c%ncols - 1) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: n_sets must equal ncols - 1")
      return
    end if
    if (.not. c_associated(c%out_data_ptr) .or. &
        c%out_data_size < (2 * c%n_sets + 1) * c%ndata) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: out_data buffer is null or too small")
      return
    end if
    if (.not. c_associated(c%diff_coeff_ptr)) then
      call error_set(c%err, ERROR_INVALID_PARAM, &
                     "diffusion_analysis_c: diff_coeff_ptr is null")
      return
    end if

    ! Zero-copy views of the NumPy buffers
    call c_f_pointer(c%msd_ptr, msd_f, [c%ncols, c%ndata])
    call c_f_pointer(c%out_data_ptr, out_data_f, [2 * c%n_sets + 1, c%ndata])
    call c_f_pointer(c%diff_coeff_ptr, diff_coeff_f, [c%n_sets])

    allocate(ndofs(c%n_sets), fittings(c%n_sets))
    ndofs(:) = c%ndofs

    write(MsgOut,'(A)') '[STEP1] Diffusion Analysis'
    write(MsgOut,'(A)') ' '

    call analyze_diffusion_unified(msd_f, real(c%time_step, wp), &
                                   real(c%distance_unit, wp), ndofs, &
                                   c%start_step, c%stop_step, &
                                   out_data_f, fittings)

    do iset = 1, c%n_sets
      diff_coeff_f(iset) = fittings(iset)%coeff(2) * 1e-4_wp
    end do

  end subroutine diffusion_analysis_body

end module diffusion_c_mod
