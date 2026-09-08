!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   result_sink_mod
!> @brief   Abstract result sink for unified analysis output
!! @authors Claude Code
!
!  (c) Copyright 2024 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module result_sink_mod

  use, intrinsic :: iso_c_binding
  use fileio_mod
  use messages_mod
  use constants_mod

  implicit none
  private

  ! Sink type enumeration
  integer, parameter, public :: SINK_FILE  = 1
  integer, parameter, public :: SINK_ARRAY = 2
  integer, parameter, public :: SINK_TEXT  = 3

  ! Abstract result sink type
  type, public :: s_result_sink
    integer :: sink_type = SINK_FILE

    ! File-based sink (CLI mode)
    integer            :: file_unit = -1
    logical            :: file_open = .false.

    ! Array-based sink (Python mode - zerocopy)
    real(wp), pointer  :: results(:) => null()
    integer            :: array_size = 0

    ! Row (vector) results: ncol values per frame, written as one file row
    ! or stored in an (ncol, nframe) array (init_sink_*_rows / write_result_row)
    integer            :: ncol = 0
    real(wp), pointer  :: results2d(:,:) => null()
    character(len=64)  :: row_format = '(i10,1x,100(f8.3,1x))'
    ! format of one scalar result line (index, value) of a file sink
    character(len=64)  :: value_format = '(i10,1x,f10.5)'

    ! Text results: lines are appended to the file or, for a text sink,
    ! accumulated in ``text`` (newline separated, see write_result_line)
    character(len=:), allocatable :: text
    integer            :: text_len = 0

    ! Common
    integer            :: current_index = 0
  end type s_result_sink

  ! Public subroutines
  public :: init_sink_file
  public :: init_sink_array
  public :: write_result
  public :: write_result_with_index
  public :: init_sink_file_rows
  public :: init_sink_array_rows
  public :: write_result_row
  public :: init_sink_text
  public :: write_result_line
  public :: sink_is_active
  public :: get_result_count
  public :: finalize_sink

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_sink_file
  !> @brief        Initialize file-based result sink (CLI mode)
  !! @authors      Claude Code
  !! @param[inout] sink     : result sink
  !! @param[in]    filename : output file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_sink_file(sink, filename, value_format)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    character(*),        intent(in)    :: filename
    character(*),        intent(in), optional :: value_format

    if (present(value_format)) sink%value_format = value_format

    sink%sink_type = SINK_FILE
    sink%current_index = 0

    if (len_trim(filename) > 0) then
      call open_file(sink%file_unit, filename, IOFileOutputNew)
      sink%file_open = .true.
    else
      sink%file_open = .false.
    end if

    return

  end subroutine init_sink_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_sink_array
  !> @brief        Initialize array-based result sink (Python zerocopy mode)
  !! @authors      Claude Code
  !! @param[inout] sink       : result sink
  !! @param[in]    results    : pre-allocated results array pointer
  !! @param[in]    array_size : size of results array
  !! @note         The results array must be pre-allocated by caller (Python).
  !!               Fortran writes directly to this array (zerocopy).
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_sink_array(sink, results, array_size)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    real(wp), target,    intent(in)    :: results(:)
    integer,             intent(in)    :: array_size

    sink%sink_type = SINK_ARRAY
    sink%results => results
    sink%array_size = array_size
    sink%current_index = 0

    return

  end subroutine init_sink_array

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_result
  !> @brief        Write a single result value (auto-increments index)
  !! @authors      Claude Code
  !! @param[inout] sink  : result sink
  !! @param[in]    value : result value to write
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_result(sink, value)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    real(wp),            intent(in)    :: value

    sink%current_index = sink%current_index + 1
    call write_result_with_index(sink, sink%current_index, value)

    return

  end subroutine write_result

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_result_with_index
  !> @brief        Write a result value with explicit index
  !! @authors      Claude Code
  !! @param[inout] sink  : result sink
  !! @param[in]    idx   : result index (1-indexed)
  !! @param[in]    value : result value to write
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_result_with_index(sink, idx, value)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    integer,             intent(in)    :: idx
    real(wp),            intent(in)    :: value

    select case(sink%sink_type)

    case(SINK_FILE)
      if (sink%file_open) then
        write(sink%file_unit, sink%value_format) idx, value
      end if

    case(SINK_ARRAY)
      if (idx >= 1 .and. idx <= sink%array_size) then
        sink%results(idx) = value
      else
        call error_msg('write_result_with_index> Index out of bounds')
      end if

    end select

    return

  end subroutine write_result_with_index

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_sink_file_rows
  !> @brief        Initialize a file sink that receives ncol values per frame
  !! @authors      Claude Code
  !! @param[inout] sink       : result sink
  !! @param[in]    filename   : output file name (empty: discard results)
  !! @param[in]    ncol       : number of values per frame
  !! @param[in]    row_format : optional Fortran format of one row
  !!                            (default '(i10,1x,100(f8.3,1x))')
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_sink_file_rows(sink, filename, ncol, row_format)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    character(*),        intent(in)    :: filename
    integer,             intent(in)    :: ncol
    character(*),        intent(in), optional :: row_format

    call init_sink_file(sink, filename)
    sink%ncol = ncol
    if (present(row_format)) sink%row_format = row_format

    return

  end subroutine init_sink_file_rows

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_sink_array_rows
  !> @brief        Initialize an array sink that receives ncol values per frame
  !! @authors      Claude Code
  !! @param[inout] sink    : result sink
  !! @param[in]    results : pre-allocated (ncol, nframe) array (zero-copy)
  !! @param[in]    ncol    : number of values per frame
  !! @param[in]    nframe  : number of frames the array can hold
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_sink_array_rows(sink, results, ncol, nframe)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    real(wp), target,    intent(in)    :: results(:,:)
    integer,             intent(in)    :: ncol
    integer,             intent(in)    :: nframe

    sink%sink_type = SINK_ARRAY
    sink%results2d => results
    sink%ncol = ncol
    sink%array_size = nframe
    sink%current_index = 0

    return

  end subroutine init_sink_array_rows

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_result_row
  !> @brief        Write the values of the next frame (index auto-incremented)
  !! @authors      Claude Code
  !! @param[inout] sink   : result sink
  !! @param[in]    values : ncol values of this frame
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_result_row(sink, values)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    real(wp),            intent(in)    :: values(:)

    sink%current_index = sink%current_index + 1

    select case(sink%sink_type)

    case(SINK_FILE)
      if (sink%file_open) then
        write(sink%file_unit, sink%row_format) sink%current_index, values
      end if

    case(SINK_ARRAY)
      if (sink%current_index >= 1 .and. sink%current_index <= sink%array_size) then
        sink%results2d(:, sink%current_index) = values
      else
        call error_msg('write_result_row> Index out of bounds')
      end if

    end select

    return

  end subroutine write_result_row

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_sink_text
  !> @brief        Initialize a sink that accumulates text lines in memory
  !!               (the Python interface returns sink%text(1:sink%text_len))
  !! @authors      Claude Code
  !! @param[inout] sink : result sink
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_sink_text(sink)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink

    sink%sink_type = SINK_TEXT
    if (allocated(sink%text)) deallocate(sink%text)
    allocate(character(len=1024) :: sink%text)
    sink%text_len = 0
    sink%current_index = 0

    return

  end subroutine init_sink_text

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_result_line
  !> @brief        Append one text line (trailing blanks removed)
  !! @authors      Claude Code
  !! @param[inout] sink : result sink (file or text; others ignore the line)
  !! @param[in]    line : the line
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_result_line(sink, line)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink
    character(*),        intent(in)    :: line

    ! local variables
    character(len=:), allocatable :: grown
    integer                       :: n, need

    sink%current_index = sink%current_index + 1

    select case(sink%sink_type)

    case(SINK_FILE)
      if (sink%file_open) then
        write(sink%file_unit, '(A)') trim(line)
      end if

    case(SINK_TEXT)
      n = len_trim(line)
      need = sink%text_len + n + 1
      if (.not. allocated(sink%text)) then
        allocate(character(len=max(1024, 2*need)) :: sink%text)
        sink%text_len = 0
      else if (need > len(sink%text)) then
        allocate(character(len=max(2*need, 2*len(sink%text))) :: grown)
        grown(1:sink%text_len) = sink%text(1:sink%text_len)
        call move_alloc(grown, sink%text)
      end if
      sink%text(sink%text_len+1:sink%text_len+n) = line(1:n)
      sink%text(sink%text_len+n+1:sink%text_len+n+1) = achar(10)
      sink%text_len = need

    end select

    return

  end subroutine write_result_line

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      sink_is_active
  !> @brief        .true. when the sink was initialized to receive results
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function sink_is_active(sink) result(active)

    ! formal arguments
    type(s_result_sink), intent(in) :: sink

    ! return value
    logical :: active

    active = sink%file_open .or. associated(sink%results) .or. &
             associated(sink%results2d) .or. sink%sink_type == SINK_TEXT

    return

  end function sink_is_active

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      get_result_count
  !> @brief        Get number of results written
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_result_count(sink) result(count)

    ! formal arguments
    type(s_result_sink), intent(in) :: sink
    integer :: count

    count = sink%current_index

    return

  end function get_result_count

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    finalize_sink
  !> @brief        Finalize and cleanup result sink
  !! @authors      Claude Code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine finalize_sink(sink)

    ! formal arguments
    type(s_result_sink), intent(inout) :: sink

    select case(sink%sink_type)

    case(SINK_FILE)
      if (sink%file_open) then
        call close_file(sink%file_unit)
        sink%file_open = .false.
      end if
      sink%file_unit = -1

    case(SINK_ARRAY)
      sink%results => null()
      sink%results2d => null()
      sink%array_size = 0

    case(SINK_TEXT)
      if (allocated(sink%text)) deallocate(sink%text)
      sink%text_len = 0

    end select

    sink%ncol = 0
    sink%current_index = 0

    return

  end subroutine finalize_sink

end module result_sink_mod
