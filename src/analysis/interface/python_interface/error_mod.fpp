!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> @brief   utilities for error handling
!! @authors CK
!
!--------1---------2---------3---------4---------5---------6---------7---------8
module error_mod
  use iso_c_binding, only: c_int, c_char, c_null_char, c_funptr, c_ptr, &
                           c_funloc, c_loc
  implicit none
  private
  public :: s_error, error_init, error_clear, error_set, error_has, &
            fi_msg_len, error_to_c, error_finish_to_c, &
            fi_error_guard_run_ctx, error_from_pending, &
            guard_body_iface, run_guarded

  ! Library-mode error guard (implemented in lib.a / fileio_data_.c).
  !
  ! fi_error_guard_run_ctx(body, ctx) calls body(ctx) under a setjmp guard and
  ! returns a nonzero value if the body aborted via error_msg (see
  ! messages_mod). The bind(C) wrapper then turns the pending error into an
  ! s_error (error_from_pending) and reports it to the C caller, instead of
  ! the process being killed by exit(1).
  !
  ! Wrappers normally call the run_guarded helper below instead of using
  ! fi_error_guard_run_ctx directly. Usage pattern (see rmsd_analysis_lazy_c
  ! for a complete example):
  !   * put every input, output and resource the body touches into a
  !     derived-type variable with the TARGET attribute in the wrapper;
  !   * make the body a bind(C) MODULE procedure taking type(c_ptr), value,
  !     and recover the context with c_f_pointer;
  !   * report failures only through the s_error in the context, and let the
  !     wrapper release the context's resources after the guard returns.
  ! Never pass c_funloc() of an INTERNAL procedure: host association forces
  ! gfortran to build a trampoline on the stack, which makes the shared
  ! library require an executable stack; dlopen() refuses that on
  ! glibc >= 2.41, and it crashes outright on non-executable stacks.
  interface
    function fi_error_guard_run_ctx(body, ctx) &
        bind(C, name="fi_error_guard_run_ctx") result(rc)
      import :: c_int, c_funptr, c_ptr
      type(c_funptr), value :: body
      type(c_ptr),    value :: ctx
      integer(c_int)        :: rc
    end function fi_error_guard_run_ctx
  end interface

  ! Interface of a guarded body (see run_guarded): a bind(C) MODULE procedure
  ! that receives the wrapper's context as an opaque pointer.
  abstract interface
    subroutine guard_body_iface(ctx) bind(C)
      import :: c_ptr
      type(c_ptr), value :: ctx
    end subroutine guard_body_iface
  end interface

  ! Legacy error code (kept for backward compatibility)
  integer, public, parameter :: ERROR_CODE = 101

  ! Memory errors (100-199)
  integer, public, parameter :: ERROR_ALLOC         = 101
  integer, public, parameter :: ERROR_DEALLOC       = 102
  integer, public, parameter :: ERROR_NOT_ALLOCATED = 103

  ! File errors (200-299)
  integer, public, parameter :: ERROR_FILE_NOT_FOUND = 201
  integer, public, parameter :: ERROR_FILE_FORMAT    = 202
  integer, public, parameter :: ERROR_FILE_READ      = 203

  ! Validation errors (300-399)
  integer, public, parameter :: ERROR_INVALID_PARAM  = 301
  integer, public, parameter :: ERROR_MISSING_PARAM  = 302
  integer, public, parameter :: ERROR_DIMENSION      = 303
  integer, public, parameter :: ERROR_ATOM_COUNT     = 304
  integer, public, parameter :: ERROR_GRID_SIZE      = 305
  integer, public, parameter :: ERROR_BOND_INFO      = 306
  integer, public, parameter :: ERROR_SELECTION      = 307

  ! Data errors (400-499)
  integer, public, parameter :: ERROR_DATA_MISMATCH  = 401
  integer, public, parameter :: ERROR_NO_DATA        = 402
  integer, public, parameter :: ERROR_MASS_UNDEFINED = 403
  integer, public, parameter :: ERROR_PBC_BOX        = 404

  ! Not supported errors (500-599)
  integer, public, parameter :: ERROR_NOT_SUPPORTED  = 501
  integer, public, parameter :: ERROR_DIM_NOT_SUPP   = 502
  integer, public, parameter :: ERROR_FUNC_NOT_SUPP  = 503
  integer, public, parameter :: ERROR_BLOCK_NOT_SUPP = 504

  ! Internal errors (600-699)
  integer, public, parameter :: ERROR_INTERNAL       = 601
  integer, public, parameter :: ERROR_SYNTAX         = 602
  integer, public, parameter :: ERROR_GENERIC        = 699

  integer(c_int), parameter :: MSG_LEN_ = 2048

  type :: s_error
     integer :: code = 0
     character(len=:), allocatable :: msg
  end type s_error

contains
  integer(c_int) function fi_msg_len() bind(C, name="fi_msg_len")
    fi_msg_len = MSG_LEN_
  end function

  !> Populate err from the error that error_msg recorded before unwinding.
  !! Used on the failure path of run_guarded so the C caller receives
  !! the message and category code exactly as the aborting routine reported.
  subroutine error_from_pending(err)
    use messages_mod, only: fi_pending_msg, fi_pending_code
    type(s_error), intent(inout) :: err
    integer :: code

    code = fi_pending_code
    if (code == 0) code = ERROR_GENERIC   ! guarantee error_has(err) is true

    if (allocated(fi_pending_msg)) then
      call error_set(err, code, fi_pending_msg)
    else
      call error_set(err, code, 'GENESIS aborted')
    end if
  end subroutine error_from_pending

  subroutine error_init(e)
    type(s_error), intent(out) :: e
    e%code = 0
    if (allocated(e%msg)) deallocate(e%msg)
  end subroutine

  subroutine error_clear(e)
    type(s_error), intent(inout) :: e
    e%code = 0
    if (allocated(e%msg)) e%msg = ''
  end subroutine

  subroutine error_set(e, code, text)
    type(s_error), intent(inout) :: e
    integer, intent(in) :: code
    character(*), intent(in) :: text
    e%code = code
    e%msg  = text
  end subroutine

  logical function error_has(e)
    type(s_error), intent(in) :: e
    error_has = (e%code /= 0)
  end function

  subroutine error_to_c(err, status, msg, msglen)
    use iso_fortran_env, only: error_unit
    implicit none
    type(s_error),         intent(in)   :: err
    integer(c_int),       intent(inout) :: status
    character(kind=c_char), intent(inout) :: msg(*)
    integer(c_int),       value         :: msglen

    integer :: i, n
    character(len=:), allocatable       :: s

    status = err%code
    if (msglen > 0) msg(1) = c_null_char

    if (allocated(err%msg)) then
       s = trim(err%msg)
    else
       s = ''
    end if
    n = min(len_trim(s),msglen-1)
    do i = 1, n
      msg(i) = s(i:i)
    end do
    if (msglen > 0) msg(n+1) = c_null_char
    ! write(error_unit, '(A,I0,2A)') 'fi_error :', status, ' ',s
    ! flush(error_unit)
  end subroutine 

  !> Report either outcome of err to the C caller in one call.
  !! Success leaves status at 0 with an empty message, so a bind(C) wrapper
  !! needs no branch of its own unless it has cleanup to run on failure.
  subroutine error_finish_to_c(err, status, msg, msglen)
    implicit none
    type(s_error),          intent(in)    :: err
    integer(c_int),         intent(inout) :: status
    character(kind=c_char), intent(inout) :: msg(*)
    integer(c_int),         value         :: msglen

    if (error_has(err)) then
      call error_to_c(err, status, msg, msglen)
      return
    end if

    status = 0
    if (msglen > 0) msg(1) = c_null_char
  end subroutine

  !> Run body(ctx) under the library-mode error guard and report the outcome
  !! to the C caller in one call.
  !!
  !! ctx is the wrapper's context variable (any derived type; it needs the
  !! TARGET attribute). On a normal return err is whatever the body left in
  !! it; when the body aborted through error_msg, err receives the pending
  !! message/code. Either way status/msg are filled, so the wrapper only has
  !! to copy its outputs out of ctx and release the resources recorded there.
  subroutine run_guarded(body, ctx, err, status, msg, msglen)
    procedure(guard_body_iface)           :: body
    type(*), target,        intent(inout) :: ctx
    type(s_error),          intent(inout) :: err
    integer(c_int),         intent(inout) :: status
    character(kind=c_char), intent(inout) :: msg(*)
    integer(c_int),         value         :: msglen

    if (fi_error_guard_run_ctx(c_funloc(body), c_loc(ctx)) /= 0) then
      call error_from_pending(err)
    end if
    call error_finish_to_c(err, status, msg, msglen)
  end subroutine run_guarded

end module error_mod
