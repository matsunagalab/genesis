module s_trajectories_c_mod
  use, intrinsic :: iso_c_binding
  use trajectory_str_mod
  use select_atoms_str_mod
  implicit none
  private

  ! multi frame s_trajectory
  type, public, bind(C) :: s_trajectories_c
    type(c_ptr) :: coords      ! double[nframe][natom][3]
    type(c_ptr) :: pbc_boxes   ! double[nframe][3][3]
    integer(c_int) :: nframe
    integer(c_int) :: natom
  end type s_trajectories_c

  public :: get_frame
  public :: set_frame
  public :: set_frame_filtered
  public :: deep_copy_s_trajectories_c
  public :: init_empty_s_trajectories_c
  public :: deallocate_s_trajectories_c
  public :: deallocate_s_trajectories_c_array
  ! bind(C) entry points called from Python must also be PUBLIC (gfortran 16.2.0
  ! hides PRIVATE bind(C) procedures, PR fortran/126872).
  public :: allocate_s_trajectories_c_array
  public :: join_s_trajectories_c

  ! Module-level save pointer for trajectory arrays (to prevent dangling pointer)
  ! Note: This is NOT thread-safe. GENESIS Python interface is single-threaded only.
  type(s_trajectories_c), pointer, save :: trajs_array_ptr(:) => null()

contains

  subroutine init_empty_s_trajectories_c(trajs, natom, nframe) &
      bind(C, name="init_empty_s_trajectories_c")
    implicit none
    type(s_trajectories_c), intent(inout) :: trajs
    integer(c_int), intent(in) :: natom
    integer(c_int), intent(in) :: nframe
    real(c_double), pointer :: coords(:,:,:), pbc_boxes(:,:,:)
    integer(c_int) :: i

    if (nframe > 0) then
      allocate(pbc_boxes(3, 3, nframe))
      pbc_boxes(:,:,:) = 0.0_C_double
      do i = 1, nframe
        pbc_boxes(1,1,i) = 1.0_C_double
        pbc_boxes(2,2,i) = 1.0_C_double
        pbc_boxes(3,3,i) = 1.0_C_double
      end do
      trajs%pbc_boxes = C_LOC(pbc_boxes)
      if (natom > 0) then
        allocate(coords(3, natom, nframe))
        trajs%coords = C_LOC(coords)
      else
        trajs%coords = c_null_ptr
      end if
    else
        trajs%pbc_boxes = c_null_ptr
      trajs%coords = c_null_ptr
    end if
    trajs%nframe = nframe
    trajs%natom = natom
  end subroutine init_empty_s_trajectories_c

  subroutine deallocate_s_trajectories_c(trajs) &
      bind(C, name="deallocate_s_trajectories_c")
    implicit none
    type(s_trajectories_c), intent(inout) :: trajs
    real(c_double), pointer :: coords(:,:,:), boxes(:,:,:)

    ! Check if coords pointer is valid before deallocating
    if (c_associated(trajs%coords) .and. trajs%natom > 0 .and. trajs%nframe > 0) then
      call C_F_POINTER(trajs%coords, coords, [3, trajs%natom, trajs%nframe])
      deallocate(coords)
      trajs%coords = c_null_ptr
    end if

    ! Check if pbc_boxes pointer is valid before deallocating
    if (c_associated(trajs%pbc_boxes) .and. trajs%nframe > 0) then
      call C_F_POINTER(trajs%pbc_boxes, boxes, [3, 3, trajs%nframe])
      deallocate(boxes)
      trajs%pbc_boxes = c_null_ptr
    end if

    trajs%natom = 0
    trajs%nframe = 0
  end subroutine deallocate_s_trajectories_c

  subroutine deep_copy_s_trajectories_c(src, dst) &
      bind(C, name="deep_copy_s_trajectories_c")
    implicit none
    type(s_trajectories_c), intent(in) :: src
    type(s_trajectories_c), intent(out) :: dst
    real(c_double), pointer :: src_coords(:,:,:)
    real(c_double), pointer :: src_boxes(:,:,:)
    real(c_double), pointer :: dst_coords(:,:,:)
    real(c_double), pointer :: dst_boxes(:,:,:)

    dst%natom = src%natom
    dst%nframe = src%nframe

    allocate(dst_coords(3, dst%natom, dst%nframe))
    call C_F_POINTER(src%coords, src_coords, [3, src%natom, src%nframe])
    dst_coords(:,:,:) = src_coords(:,:,:)
    dst%coords = c_loc(dst_coords)

    allocate(dst_boxes(3, 3, dst%nframe))
    call C_F_POINTER(src%pbc_boxes, src_boxes, [3, 3, src%nframe])
    dst_boxes(:,:,:) = src_boxes(:,:,:)
    dst%pbc_boxes = c_loc(dst_boxes)
  end subroutine deep_copy_s_trajectories_c

  ! Accessor for single frame
  subroutine get_frame(trajs_c, frame_idx, traj_fort)
    type(s_trajectories_c), intent(in) :: trajs_c
    integer, intent(in) :: frame_idx
    type(s_trajectory), intent(inout) :: traj_fort
    real(C_double), pointer :: coords(:,:,:), boxes(:,:,:)

    if (allocated(traj_fort%coord)) then
      if (.not. all(shape(traj_fort%coord) == [3, trajs_c%natom])) then
        deallocate(traj_fort%coord)
        allocate(traj_fort%coord(3,trajs_c%natom))
      end if
    else
      allocate(traj_fort%coord(3,trajs_c%natom))
    end if
    call C_F_POINTER(trajs_c%coords, coords, [3, trajs_c%natom, trajs_c%nframe])
    call C_F_POINTER(trajs_c%pbc_boxes, boxes, [3, 3, trajs_c%nframe])
    traj_fort%pbc_box = boxes(:, :, frame_idx)
    traj_fort%coord = coords(:, :, frame_idx)
  end subroutine

  subroutine set_frame(out_traj, src_traj, frame_idx)
    type(s_trajectories_c), intent(inout) :: out_traj
    type(s_trajectory), intent(in) :: src_traj
    integer, intent(in) :: frame_idx

    real(C_double), pointer :: coords(:,:,:), boxes(:,:,:)

    call C_F_POINTER(out_traj%coords, coords, [3, out_traj%natom, out_traj%nframe])
    call C_F_POINTER(out_traj%pbc_boxes, boxes, [3, 3, out_traj%nframe])
    coords(:, :, frame_idx) = src_traj%coord
    boxes(:, :, frame_idx) = src_traj%pbc_box
  end subroutine

  subroutine set_frame_filtered(out_traj, src_traj, frame_idx, selatoms)
    type(s_trajectories_c), intent(inout) :: out_traj
    type(s_trajectory), intent(in) :: src_traj
    integer, intent(in) :: frame_idx
    type(s_selatoms), intent(in) :: selatoms

    real(C_double), pointer :: coords(:,:,:), boxes(:,:,:)
    integer :: i, atom_idx

    call C_F_POINTER(out_traj%coords, coords, [3, out_traj%natom, out_traj%nframe])
    call C_F_POINTER(out_traj%pbc_boxes, boxes, [3, 3, out_traj%nframe])
    
    ! Copy only selected atoms
    do i = 1, size(selatoms%idx)
      atom_idx = selatoms%idx(i)
      coords(:, i, frame_idx) = src_traj%coord(:, atom_idx)
    end do
    
    ! Copy PBC box (unchanged)
    boxes(:, :, frame_idx) = src_traj%pbc_box
  end subroutine set_frame_filtered

  subroutine join_s_trajectories_c(trajs_array, len, joined) &
      bind(C, name="join_s_trajectories_c")
    implicit none
    type(c_ptr), intent(in) :: trajs_array ! array of s_trajectories_c
    integer(c_int), intent(in) :: len
    type(s_trajectories_c), intent(out) :: joined
    type(s_trajectories_c), pointer :: src(:)
    integer :: i, j
    integer :: nframe_sum
    integer :: natom
    real(C_double), pointer :: coords(:,:,:), boxes(:,:,:)
    real(C_double), pointer :: jcoords(:,:,:), jboxes(:,:,:)

    if (len > 0) then
        call C_F_POINTER(trajs_array, src, [len])
        natom = src(1)%natom
        nframe_sum = src(1)%nframe
        do i = 2, len
          nframe_sum = nframe_sum + src(i)%nframe
          if (natom /= src(i)%natom) then
            natom = 0
            exit
          end if
        end do
        if (natom > 0) then
            call init_empty_s_trajectories_c(joined, natom, nframe_sum)
            call C_F_POINTER(joined%coords, jcoords, [3, natom, joined%nframe])
            call C_F_POINTER(joined%pbc_boxes, jboxes, [3, 3, joined%nframe])
            j = 1
            do i = 1, len
              call C_F_POINTER(src(i)%coords, coords, [3, natom, src(i)%nframe])
              call C_F_POINTER(src(i)%pbc_boxes, boxes, [3, 3, src(i)%nframe])
              jcoords(:, :, j:j+src(i)%nframe-1) = coords(:,:,:)
              jboxes(:, :, j:j+src(i)%nframe-1) = boxes(:,:,:)
              j = j + src(i)%nframe
            end do
        end if
    else
        call init_empty_s_trajectories_c(joined, 0, 0)
    end if
  end subroutine join_s_trajectories_c

  subroutine allocate_s_trajectories_c_array(trajs_array, len) &
      bind(C, name="allocate_s_trajectories_c_array")
    implicit none
    integer(c_int), intent(in) :: len
    type(c_ptr), intent(out) :: trajs_array
    integer :: i

    ! Clean up any existing allocation first
    if (associated(trajs_array_ptr)) then
      do i = 1, size(trajs_array_ptr)
        call deallocate_s_trajectories_c(trajs_array_ptr(i))
      end do
      deallocate(trajs_array_ptr)
      nullify(trajs_array_ptr)
    end if

    ! Use module-level pointer to prevent dangling pointer
    allocate(trajs_array_ptr(len))
    trajs_array = c_loc(trajs_array_ptr)
  end subroutine allocate_s_trajectories_c_array

  subroutine deallocate_s_trajectories_c_array(trajs_array, len) &
      bind(C, name="deallocate_s_trajectories_c_array")
    implicit none
    type(c_ptr), intent(inout) :: trajs_array
    integer(c_int), intent(in) :: len
    integer :: i

    ! Use module-level pointer
    if (associated(trajs_array_ptr)) then
      do i = 1, size(trajs_array_ptr)
        call deallocate_s_trajectories_c(trajs_array_ptr(i))
      end do
      deallocate(trajs_array_ptr)
      nullify(trajs_array_ptr)
    end if
    trajs_array = c_null_ptr
  end subroutine deallocate_s_trajectories_c_array
end module s_trajectories_c_mod
