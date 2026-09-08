module s_molecule_c_mod
  use, intrinsic :: iso_c_binding
  use molecules_str_mod  ! Contains s_molecule definition
  implicit none
  private

  type, public, bind(C) :: s_molecule_c
    integer(c_int) :: num_deg_freedom

    integer(c_int) :: num_atoms
    integer(c_int) :: num_bonds
    integer(c_int) :: num_enm_bonds
    integer(c_int) :: num_angles
    integer(c_int) :: num_dihedrals
    integer(c_int) :: num_impropers
    integer(c_int) :: num_cmaps

    integer(c_int) :: num_residues
    integer(c_int) :: num_molecules
    integer(c_int) :: num_segments
    ! shift or not shift
    logical(c_bool) :: shift_origin
    logical(c_bool) :: special_hydrogen

    real(c_double)  :: total_charge
    ! atom
    type(c_ptr) :: atom_no ! c_int[num_atoms]
    type(c_ptr) :: segment_name ! c_char[num_atoms][4]
    type(c_ptr) :: segment_no ! c_int[num_atoms]
    type(c_ptr) :: residue_no ! c_int[num_atoms]
    type(c_ptr) :: residue_c_no ! c_int[num_atoms]
    type(c_ptr) :: residue_name ! c_char[num_atoms][6]
    type(c_ptr) :: atom_name ! c_char[num_atoms][4]
    type(c_ptr) :: atom_cls_name ! c_char[num_atoms][6]
    type(c_ptr) :: atom_cls_no ! c_int[num_atoms]
    type(c_ptr) :: charge ! c_double[num_atoms]
    type(c_ptr) :: mass ! c_double[num_atoms]
    type(c_ptr) :: inv_mass ! c_double[num_atoms]
    type(c_ptr) :: imove ! c_int[num_atoms]
    type(c_ptr) :: stokes_radius ! c_double[num_atoms]
    type(c_ptr) :: inv_stokes_radius ! c_double[num_atoms]

    type(c_ptr) :: chain_id ! c_char[num_atoms][1]
    type(c_ptr) :: atom_coord ! c_double[num_atoms][3]
    type(c_ptr) :: atom_occupancy ! c_double[num_atoms]
    type(c_ptr) :: atom_temp_factor ! c_double[num_atoms]
    type(c_ptr) :: atom_velocity ! c_double[num_atoms][3]
    type(c_ptr) :: light_atom_name ! c_bool[num_atoms]
    type(c_ptr) :: light_atom_mass ! c_bool[num_atoms]

    type(c_ptr) :: molecule_no ! c_int[num_atoms]

    ! bond
    type(c_ptr) :: bond_list ! c_int[num_bonds][2]
    ! enm
    type(c_ptr) :: enm_list ! c_int[num_enm_bonds][2]
    ! angle
    type(c_ptr) :: angl_list ! c_int[num_angles][3]
    ! dihedral
    type(c_ptr) :: dihe_list ! c_int[num_dihedrals][4]
    ! improper
    type(c_ptr) :: impr_list ! c_int[num_impropers][4]
    ! cmap
    type(c_ptr) :: cmap_list ! c_int[num_cmaps][8]
    ! molecule
    type(c_ptr) :: molecule_atom_no ! c_int[num_molecules]
    type(c_ptr) :: molecule_mass ! c_double[num_molecules]
    type(c_ptr) :: molecule_name ! c_char[num_molecules][10]
    ! ref coord
    type(c_ptr) :: atom_refcoord ! c_double[num_atoms][3]
    ! fit coord
    type(c_ptr) :: atom_fitcoord ! c_double[num_atoms][3]
    ! principal component mode
    integer(c_int) :: num_pc_modes
    type(c_ptr) :: pc_mode ! c_double[num_pc_modes]
    ! ! FEP
    integer(c_int)           :: fep_topology
    integer(c_int)           :: num_hbonds_singleA
    integer(c_int)           :: num_hbonds_singleB
    type(c_ptr)              :: num_atoms_fep     ! c_int[5]
    type(c_ptr)              :: num_bonds_fep     ! c_int[5]
    type(c_ptr)              :: num_angles_fep    ! c_int[5]
    type(c_ptr)              :: num_dihedrals_fep ! c_int[5]
    type(c_ptr)              :: num_impropers_fep ! c_int[5]
    type(c_ptr)              :: num_cmaps_fep     ! c_int[5]
    type(c_ptr) :: bond_list_fep ! c_int[5][nbnd_fep_max][2]
    type(c_ptr) :: angl_list_fep ! c_int[5][nangl_fep_max][3]
    type(c_ptr) :: dihe_list_fep ! c_int[5][ndihe_fep_max][4]
    type(c_ptr) :: impr_list_fep ! c_int[5][nimpr_fep_max][4]
    type(c_ptr) :: cmap_list_fep ! c_int[5][ncmap_fep_max][8]
    type(c_ptr) :: id_singleA ! c_int[size_id_singleA]
    type(c_ptr) :: id_singleB ! c_int[size_id_singleB]
    type(c_ptr) :: fepgrp ! c_int[size_fepgrp]
    type(c_ptr) :: fepgrp_bond ! c_int[5][5]
    type(c_ptr) :: fepgrp_angl ! c_int[5][5][5]
    type(c_ptr) :: fepgrp_dihe ! c_int[5][5][5][5]
    type(c_ptr) :: fepgrp_cmap ! c_int[5*5*5*5*5*5*5*5]

    ! for deallocate data
    integer :: nbnd_fep_max
    integer :: nangl_fep_max
    integer :: ndihe_fep_max
    integer :: nimpr_fep_max
    integer :: ncmap_fep_max
    integer :: size_id_singleA
    integer :: size_id_singleB
    integer :: size_fepgrp
  end type s_molecule_c

  public :: f2c_s_molecule
  public :: c2f_s_molecule
  public :: deallocate_s_molecule_c
  ! bind(C) entry points called from Python must also be PUBLIC: gfortran
  ! 16.2.0 gives PRIVATE bind(C) procedures hidden visibility (PR fortran/126872),
  ! so dlsym() cannot find them.
  public :: allocate_s_molecule_c

contains

  subroutine f2c_s_molecule(f_src, c_dst)
    use conv_f_c_util
    implicit none
    type(s_molecule), intent(in) :: f_src
    type(s_molecule_c), intent(out) :: c_dst
    integer, pointer :: ptr_int(:)
    integer :: i

    c_dst%num_deg_freedom = f_src%num_deg_freedom
    c_dst%num_atoms = f_src%num_atoms
    c_dst%num_bonds = f_src%num_bonds
    c_dst%num_enm_bonds = f_src%num_enm_bonds
    c_dst%num_angles = f_src%num_angles
    c_dst%num_dihedrals = f_src%num_dihedrals
    c_dst%num_impropers = f_src%num_impropers
    c_dst%num_cmaps     = f_src%num_cmaps

    c_dst%num_residues = f_src%num_residues
    c_dst%num_molecules = f_src%num_molecules
    c_dst%num_segments = f_src%num_segments

    c_dst%shift_origin = f_src%shift_origin
    c_dst%special_hydrogen = f_src%special_hydrogen

    c_dst%total_charge = f_src%total_charge

    c_dst%atom_no = f2c_int_array_nullcheck(f_src%atom_no)
    c_dst%segment_name  = f2c_string_array_nullcheck(f_src%segment_name)
    c_dst%segment_no    = f2c_int_array_nullcheck(f_src%segment_no)
    c_dst%residue_no    = f2c_int_array_nullcheck(f_src%residue_no)
    c_dst%residue_c_no  = f2c_int_array_nullcheck(f_src%residue_c_no)
    c_dst%residue_name  = f2c_string_array_nullcheck(f_src%residue_name)
    c_dst%atom_name     = f2c_string_array_nullcheck(f_src%atom_name)
    c_dst%atom_cls_name = f2c_string_array_nullcheck(f_src%atom_cls_name)
    c_dst%atom_cls_no   = f2c_int_array_nullcheck(f_src%atom_cls_no)
    c_dst%charge            = f2c_double_array_nullcheck(f_src%charge)
    c_dst%mass              = f2c_double_array_nullcheck(f_src%mass)
    c_dst%inv_mass          = f2c_double_array_nullcheck(f_src%inv_mass)
    c_dst%imove             = f2c_int_array_nullcheck(f_src%imove)
    c_dst%stokes_radius     = f2c_double_array_nullcheck(f_src%stokes_radius)
    c_dst%inv_stokes_radius = f2c_double_array_nullcheck(f_src%inv_stokes_radius)

    c_dst%chain_id = f2c_string_array_nullcheck(f_src%chain_id)
    c_dst%atom_coord = f2c_double_array_nullcheck(f_src%atom_coord)
    c_dst%atom_occupancy = f2c_double_array_nullcheck(f_src%atom_occupancy)
    c_dst%atom_temp_factor = f2c_double_array_nullcheck(f_src%atom_temp_factor)
    c_dst%atom_velocity = f2c_double_array_nullcheck(f_src%atom_velocity)
    c_dst%light_atom_name = f2c_bool_array_nullcheck(f_src%light_atom_name)
    c_dst%light_atom_mass = f2c_bool_array_nullcheck(f_src%light_atom_mass)

    c_dst%molecule_no = f2c_int_array_nullcheck(f_src%molecule_no)

    c_dst%bond_list = f2c_int_array_nullcheck(f_src%bond_list)
    c_dst%enm_list = f2c_int_array_nullcheck(f_src%enm_list)
    c_dst%angl_list = f2c_int_array_nullcheck(f_src%angl_list)
    c_dst%dihe_list = f2c_int_array_nullcheck(f_src%dihe_list)
    c_dst%impr_list = f2c_int_array_nullcheck(f_src%impr_list)
    c_dst%cmap_list = f2c_int_array_nullcheck(f_src%cmap_list)
    c_dst%molecule_atom_no = f2c_int_array_nullcheck(f_src%molecule_atom_no)
    c_dst%molecule_mass = f2c_double_array_nullcheck(f_src%molecule_mass)
    c_dst%molecule_name = f2c_string_array_nullcheck(f_src%molecule_name)
    c_dst%atom_refcoord = f2c_double_array_nullcheck(f_src%atom_refcoord)
    c_dst%atom_fitcoord = f2c_double_array_nullcheck(f_src%atom_fitcoord)
    c_dst%num_pc_modes = f_src%num_pc_modes
    c_dst%pc_mode = f2c_double_array_nullcheck(f_src%pc_mode)

    c_dst%fep_topology       = f_src%fep_topology
    c_dst%num_hbonds_singleA = f_src%num_hbonds_singleA
    c_dst%num_hbonds_singleB = f_src%num_hbonds_singleB

    c_dst%num_atoms_fep     = f2c_int_array(f_src%num_atoms_fep)
    c_dst%num_bonds_fep     = f2c_int_array(f_src%num_bonds_fep)
    c_dst%num_angles_fep    = f2c_int_array(f_src%num_angles_fep)
    c_dst%num_dihedrals_fep = f2c_int_array(f_src%num_dihedrals_fep)
    c_dst%num_impropers_fep = f2c_int_array(f_src%num_impropers_fep)
    c_dst%num_cmaps_fep     = f2c_int_array(f_src%num_cmaps_fep)

    if (allocated(f_src%bond_list_fep)) then
      c_dst%bond_list_fep = f2c_int_array(f_src%bond_list_fep)
      c_dst%nbnd_fep_max = size(f_src%bond_list_fep, dim=2)
    else
      c_dst%nbnd_fep_max = 0
    end if
    if (allocated(f_src%angl_list_fep)) then
      c_dst%angl_list_fep = f2c_int_array(f_src%angl_list_fep)
      c_dst%nangl_fep_max = size(f_src%angl_list_fep, dim=2)
    else
      c_dst%nangl_fep_max = 0
    end if
    if (allocated(f_src%dihe_list_fep)) then
      c_dst%dihe_list_fep= f2c_int_array(f_src%dihe_list_fep)
      c_dst%ndihe_fep_max = size(f_src%dihe_list_fep, dim=2)
    else
      c_dst%ndihe_fep_max = 0
    end if
    if (allocated(f_src%impr_list_fep)) then
      c_dst%impr_list_fep = f2c_int_array(f_src%impr_list_fep)
      c_dst%nimpr_fep_max = size(f_src%impr_list_fep, dim=2)
    else
      c_dst%nimpr_fep_max = 0
    end if
    if (allocated(f_src%cmap_list_fep)) then
      c_dst%cmap_list_fep = f2c_int_array(f_src%cmap_list_fep)
      c_dst%ncmap_fep_max = size(f_src%cmap_list_fep, dim=2)
    else
      c_dst%ncmap_fep_max = 0
    end if
    if (allocated(f_src%id_singleA)) then
      c_dst%id_singleA  = f2c_int_array(f_src%id_singleA)
      c_dst%size_id_singleA = size(f_src%id_singleA)
    else
      c_dst%size_id_singleA = 0
    end if
    if (allocated(f_src%id_singleB)) then
      c_dst%id_singleB  = f2c_int_array(f_src%id_singleB)
      c_dst%size_id_singleB = size(f_src%id_singleB)
    else
      c_dst%size_id_singleB = 0
    end if
    if (allocated(f_src%fepgrp)) then
      c_dst%fepgrp      = f2c_int_array(f_src%fepgrp)
      c_dst%size_fepgrp = size(f_src%fepgrp)
    else
      c_dst%size_fepgrp = 0
    end if
    c_dst%fepgrp_bond = f2c_int_array(f_src%fepgrp_bond)
    c_dst%fepgrp_angl = f2c_int_array(f_src%fepgrp_angl)
    c_dst%fepgrp_dihe = f2c_int_array(f_src%fepgrp_dihe)
    c_dst%fepgrp_cmap = f2c_int_array(f_src%fepgrp_cmap)
  end subroutine f2c_s_molecule

  subroutine c2f_s_molecule(c_src, f_dst)
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(in) :: c_src
    type(s_molecule), intent(out) :: f_dst
    integer(c_int), pointer :: work_int(:)

    f_dst%num_deg_freedom = c_src%num_deg_freedom
    f_dst%num_atoms = c_src%num_atoms
    f_dst%num_bonds = c_src%num_bonds
    f_dst%num_enm_bonds = c_src%num_enm_bonds
    f_dst%num_angles = c_src%num_angles
    f_dst%num_dihedrals = c_src%num_dihedrals
    f_dst%num_impropers = c_src%num_impropers
    f_dst%num_cmaps     = c_src%num_cmaps

    f_dst%num_residues = c_src%num_residues
    f_dst%num_molecules = c_src%num_molecules
    f_dst%num_segments = c_src%num_segments

    f_dst%shift_origin = c_src%shift_origin
    f_dst%special_hydrogen = c_src%special_hydrogen

    f_dst%total_charge = c_src%total_charge

    call c2f_int_array(f_dst%atom_no, c_src%atom_no, [c_src%num_atoms])

    call c2f_string_array(f_dst%segment_name, c_src%segment_name, c_src%num_atoms)
    call c2f_int_array(f_dst%segment_no, c_src%segment_no, [c_src%num_atoms])
    call c2f_int_array(f_dst%residue_no, c_src%residue_no, [c_src%num_atoms])
    call c2f_int_array(f_dst%residue_c_no, c_src%residue_c_no, [c_src%num_atoms])
    call c2f_string_array(f_dst%residue_name, c_src%residue_name, c_src%num_atoms)
    call c2f_string_array(f_dst%atom_name, c_src%atom_name, c_src%num_atoms)
    call c2f_string_array(f_dst%atom_cls_name, c_src%atom_cls_name, c_src%num_atoms)
    call c2f_int_array(f_dst%atom_cls_no, c_src%atom_cls_no, [c_src%num_atoms])
    call c2f_double_array(f_dst%charge, c_src%charge, [c_src%num_atoms])
    call c2f_double_array(f_dst%mass, c_src%mass, [c_src%num_atoms])
    call c2f_double_array(f_dst%inv_mass, c_src%inv_mass, [c_src%num_atoms])
    call c2f_int_array(f_dst%imove, c_src%imove, [c_src%num_atoms])
    call c2f_double_array(f_dst%stokes_radius, c_src%stokes_radius, [c_src%num_atoms])
    call c2f_double_array(f_dst%inv_stokes_radius, c_src%inv_stokes_radius, [c_src%num_atoms])

    call c2f_string_array(f_dst%chain_id, c_src%chain_id, c_src%num_atoms)
    call c2f_double_array(f_dst%atom_coord, c_src%atom_coord, [3, c_src%num_atoms])
    call c2f_double_array(f_dst%atom_occupancy, c_src%atom_occupancy, [c_src%num_atoms])
    call c2f_double_array(f_dst%atom_temp_factor, c_src%atom_temp_factor, [c_src%num_atoms])
    call c2f_double_array(f_dst%atom_velocity, c_src%atom_velocity, [3, c_src%num_atoms])
    call c2f_bool_array(f_dst%light_atom_name, c_src%light_atom_name, [c_src%num_atoms])
    call c2f_bool_array(f_dst%light_atom_mass, c_src%light_atom_mass, [c_src%num_atoms])

    call c2f_int_array(f_dst%molecule_no, c_src%molecule_no, [c_src%num_atoms])

    call c2f_int_array(f_dst%bond_list, c_src%bond_list, [2, c_src%num_bonds])
    call c2f_int_array(f_dst%enm_list, c_src%enm_list, [2, c_src%num_enm_bonds])
    call c2f_int_array(f_dst%angl_list, c_src%angl_list, [3, c_src%num_angles])
    call c2f_int_array(f_dst%dihe_list, c_src%dihe_list, [4, c_src%num_dihedrals])
    call c2f_int_array(f_dst%impr_list, c_src%impr_list, [4, c_src%num_impropers])
    call c2f_int_array(f_dst%cmap_list, c_src%cmap_list, [8, c_src%num_cmaps])

    call c2f_int_array(f_dst%molecule_atom_no, c_src%molecule_atom_no, [c_src%num_molecules])
    call c2f_double_array(f_dst%molecule_mass, c_src%molecule_mass, [c_src%num_molecules])
    call c2f_string_array(f_dst%molecule_name, c_src%molecule_name, c_src%num_molecules)
    call c2f_double_array(f_dst%atom_refcoord, c_src%atom_refcoord, [3, c_src%num_atoms])
    call c2f_double_array(f_dst%atom_fitcoord, c_src%atom_fitcoord, [3, c_src%num_atoms])

    f_dst%num_pc_modes = c_src%num_pc_modes
    call c2f_double_array(f_dst%pc_mode, c_src%pc_mode, [c_src%num_pc_modes])

    f_dst%fep_topology       = c_src%fep_topology
    f_dst%num_hbonds_singleA = c_src%num_hbonds_singleA
    f_dst%num_hbonds_singleB = c_src%num_hbonds_singleB

    call c2f_int_array_static(f_dst%num_atoms_fep, c_src%num_atoms_fep, [5])
    call c2f_int_array_static(f_dst%num_bonds_fep, c_src%num_bonds_fep, [5])
    call c2f_int_array_static(f_dst%num_angles_fep, c_src%num_angles_fep, [5])
    call c2f_int_array_static(f_dst%num_dihedrals_fep, c_src%num_dihedrals_fep, [5])
    call c2f_int_array_static(f_dst%num_impropers_fep, c_src%num_impropers_fep, [5])
    call c2f_int_array_static(f_dst%num_cmaps_fep, c_src%num_cmaps_fep, [5])

    call c2f_int_array(f_dst%bond_list_fep, c_src%bond_list_fep, [5, c_src%nbnd_fep_max, 2])
    call c2f_int_array(f_dst%angl_list_fep, c_src%angl_list_fep, [5, c_src%nangl_fep_max, 3])
    call c2f_int_array(f_dst%dihe_list_fep, c_src%dihe_list_fep, [5, c_src%ndihe_fep_max, 4])
    call c2f_int_array(f_dst%impr_list_fep, c_src%impr_list_fep, [5, c_src%nimpr_fep_max, 4])
    call c2f_int_array(f_dst%cmap_list_fep, c_src%cmap_list_fep, [5, c_src%ncmap_fep_max, 8])

    call c2f_int_array(f_dst%id_singleA, c_src%id_singleA, [c_src%size_id_singleA])
    call c2f_int_array(f_dst%id_singleB, c_src%id_singleB, [c_src%size_id_singleB])
    call c2f_int_array(f_dst%fepgrp, c_src%fepgrp, [c_src%size_fepgrp])
    call c2f_int_array_static(f_dst%fepgrp_bond, c_src%fepgrp_bond, [5, 5])
    call c2f_int_array_static(f_dst%fepgrp_angl, c_src%fepgrp_angl, [5, 5, 5])
    call c2f_int_array_static(f_dst%fepgrp_dihe, c_src%fepgrp_dihe, [5, 5, 5, 5])
    call c2f_int_array_static(f_dst%fepgrp_cmap, c_src%fepgrp_cmap, [5*5*5*5*5*5*5*5])
  end subroutine c2f_s_molecule

  subroutine deallocate_s_molecule_c(self) &
      bind(C, name="deallocate_s_molecule_c")
    use, intrinsic :: iso_c_binding
    implicit none
    type(s_molecule_c), intent(inout) :: self
    character(kind=c_char), pointer :: work_string(:,:)
    integer(c_int), pointer :: work_int(:)
    integer(c_int), pointer :: work_int2(:,:)
    integer(c_int), pointer :: work_int3(:,:,:)
    integer(c_int), pointer :: work_int4(:,:,:,:)
    real(c_double), pointer :: work_double(:)
    real(c_double), pointer :: work_double2(:,:)
    logical(c_bool), pointer :: work_bool(:)

    if (c_associated(self%atom_no)) then
      call C_F_pointer(self%atom_no, work_int, [self%num_atoms])
      deallocate(work_int)
      self%atom_no = c_null_ptr
    end if
    if (c_associated(self%segment_name)) then
      call C_F_pointer(self%segment_name, work_int2, [4, self%num_atoms])
      deallocate(work_int2)
      self%segment_name = c_null_ptr
    end if
    if (c_associated(self%segment_no)) then
      call C_F_pointer(self%segment_no, work_int, [self%num_atoms])
      deallocate(work_int)
      self%segment_no = c_null_ptr
    end if
    if (c_associated(self%residue_no)) then
      call C_F_pointer(self%residue_no, work_int, [self%num_atoms])
      deallocate(work_int)
      self%residue_no = c_null_ptr
    end if
    if (c_associated(self%residue_c_no)) then
      call C_F_pointer(self%residue_c_no, work_int, [self%num_atoms])
      deallocate(work_int)
      self%residue_c_no = c_null_ptr
    end if
    if (c_associated(self%residue_name)) then
      call C_F_pointer(self%residue_name, work_string, [6, self%num_atoms])
      deallocate(work_string)
      self%residue_name = c_null_ptr
    end if
    if (c_associated(self%atom_name)) then
      call C_F_pointer(self%atom_name, work_string, [4, self%num_atoms])
      deallocate(work_string)
      self%atom_name = c_null_ptr
    end if
    if (c_associated(self%atom_cls_name)) then
      call C_F_pointer(self%atom_cls_name, work_string, [6, self%num_atoms])
      deallocate(work_string)
      self%atom_cls_name = c_null_ptr
    end if
    if (c_associated(self%atom_cls_no)) then
      call C_F_pointer(self%atom_cls_no, work_int, [self%num_atoms])
      deallocate(work_int)
      self%atom_cls_no = c_null_ptr
    end if
    if (c_associated(self%charge)) then
      call C_F_pointer(self%charge, work_double, [self%num_atoms])
      deallocate(work_double)
      self%charge = c_null_ptr
    end if
    if (c_associated(self%mass)) then
      call C_F_pointer(self%mass, work_double, [self%num_atoms])
      deallocate(work_double)
      self%mass = c_null_ptr
    end if
    if (c_associated(self%inv_mass)) then
      call C_F_pointer(self%inv_mass, work_double, [self%num_atoms])
      deallocate(work_double)
      self%inv_mass = c_null_ptr
    end if
    if (c_associated(self%imove)) then
      call C_F_pointer(self%imove, work_int, [self%num_atoms])
      deallocate(work_int)
      self%imove = c_null_ptr
    end if
    if (c_associated(self%stokes_radius)) then
      call C_F_pointer(self%stokes_radius, work_double, [self%num_atoms])
      deallocate(work_double)
      self%stokes_radius = c_null_ptr
    end if
    if (c_associated(self%inv_stokes_radius)) then
      call C_F_pointer(self%inv_stokes_radius, work_double, [self%num_atoms])
      deallocate(work_double)
      self%inv_stokes_radius = c_null_ptr
    end if
    if (c_associated(self%chain_id)) then
      call C_F_pointer(self%chain_id, work_string, [1, self%num_atoms])
      deallocate(work_string)
      self%chain_id = c_null_ptr
    end if
    if (c_associated(self%atom_coord)) then
      call C_F_pointer(self%atom_coord, work_double2, [3, self%num_atoms])
      deallocate(work_double2)
      self%atom_coord = c_null_ptr
    end if
    if (c_associated(self%atom_occupancy)) then
      call C_F_pointer(self%atom_occupancy, work_double, [self%num_atoms])
      deallocate(work_double)
      self%atom_occupancy = c_null_ptr
    end if
    if (c_associated(self%atom_temp_factor)) then
      call C_F_pointer(self%atom_temp_factor, work_double, [self%num_atoms])
      deallocate(work_double)
      self%atom_temp_factor = c_null_ptr
    end if
    if (c_associated(self%atom_velocity)) then
      call C_F_pointer(self%atom_velocity, work_double2, [3, self%num_atoms])
      deallocate(work_double2)
      self%atom_velocity = c_null_ptr
    end if
    if (c_associated(self%light_atom_name)) then
      call C_F_pointer(self%light_atom_name, work_bool, [self%num_atoms])
      deallocate(work_bool)
      self%light_atom_name = c_null_ptr
    end if
    if (c_associated(self%light_atom_mass)) then
      call C_F_pointer(self%light_atom_mass, work_bool, [self%num_atoms])
      deallocate(work_bool)
      self%light_atom_mass = c_null_ptr
    end if
    if (c_associated(self%molecule_no)) then
      call C_F_pointer(self%molecule_no, work_int, [self%num_atoms])
      deallocate(work_int)
      self%molecule_no = c_null_ptr
    end if
    if (c_associated(self%bond_list)) then
      call C_F_pointer(self%bond_list, work_int2, [2, self%num_bonds])
      deallocate(work_int2)
      self%bond_list = c_null_ptr
    end if
    if (c_associated(self%enm_list)) then
      call C_F_pointer(self%enm_list, work_int2, [2, self%num_enm_bonds])
      deallocate(work_int2)
      self%enm_list = c_null_ptr
    end if
    if (c_associated(self%angl_list)) then
      call C_F_pointer(self%angl_list, work_int2, [3, self%num_angles])
      deallocate(work_int2)
      self%angl_list = c_null_ptr
    end if
    if (c_associated(self%dihe_list)) then
      call C_F_pointer(self%dihe_list, work_int2, [4, self%num_dihedrals])
      deallocate(work_int2)
      self%dihe_list = c_null_ptr
    end if
    if (c_associated(self%impr_list)) then
      call C_F_pointer(self%impr_list, work_int2, [4, self%num_impropers])
      deallocate(work_int2)
      self%impr_list = c_null_ptr
    end if
    if (c_associated(self%cmap_list)) then
      call C_F_pointer(self%cmap_list, work_int2, [8, self%num_cmaps])
      deallocate(work_int2)
      self%cmap_list = c_null_ptr
    end if
    if (c_associated(self%molecule_atom_no)) then
      call C_F_pointer(self%molecule_atom_no, work_int, [self%num_molecules])
      deallocate(work_int)
      self%molecule_atom_no = c_null_ptr
    end if
    if (c_associated(self%molecule_mass)) then
      call C_F_pointer(self%molecule_mass, work_double, [self%num_molecules])
      deallocate(work_double)
      self%molecule_mass = c_null_ptr
    end if
    if (c_associated(self%molecule_name)) then
      call C_F_pointer(self%molecule_name, work_string, [10, self%num_molecules])
      deallocate(work_string)
      self%molecule_name = c_null_ptr
    end if
    if (c_associated(self%atom_refcoord)) then
      call C_F_pointer(self%atom_refcoord, work_double2, [3, self%num_atoms])
      deallocate(work_double2)
      self%atom_refcoord = c_null_ptr
    end if
    if (c_associated(self%atom_fitcoord)) then
      call C_F_pointer(self%atom_fitcoord, work_double2, [3, self%num_atoms])
      deallocate(work_double2)
      self%atom_fitcoord = c_null_ptr
    end if
    if (c_associated(self%pc_mode)) then
      call C_F_pointer(self%pc_mode, work_double, [self%num_pc_modes])
      deallocate(work_double)
      self%pc_mode = c_null_ptr
    end if
    if (c_associated(self%num_atoms_fep)) then
      call C_F_pointer(self%num_atoms_fep, work_int, [5])
      deallocate(work_int)
      self%num_atoms_fep = c_null_ptr
    end if
    if (c_associated(self%num_bonds_fep)) then
      call C_F_pointer(self%num_bonds_fep, work_int, [5])
      deallocate(work_int)
      self%num_bonds_fep = c_null_ptr
    end if
    if (c_associated(self%num_angles_fep)) then
      call C_F_pointer(self%num_angles_fep, work_int, [5])
      deallocate(work_int)
      self%num_angles_fep = c_null_ptr
    end if
    if (c_associated(self%num_dihedrals_fep)) then
      call C_F_pointer(self%num_dihedrals_fep, work_int, [5])
      deallocate(work_int)
      self%num_dihedrals_fep = c_null_ptr
    end if
    if (c_associated(self%num_impropers_fep)) then
      call C_F_pointer(self%num_impropers_fep, work_int, [5])
      deallocate(work_int)
      self%num_impropers_fep = c_null_ptr
    end if
    if (c_associated(self%num_cmaps_fep)) then
      call C_F_pointer(self%num_cmaps_fep, work_int, [5])
      deallocate(work_int)
      self%num_cmaps_fep = c_null_ptr
    end if
    if (c_associated(self%bond_list_fep)) then
      call C_F_pointer(self%bond_list_fep, work_int3, [2, self%nbnd_fep_max, 5])
      deallocate(work_int3)
      self%bond_list_fep = c_null_ptr
    end if
    if (c_associated(self%angl_list_fep)) then
      call C_F_pointer(self%angl_list_fep, work_int3, [3, self%nangl_fep_max, 5])
      deallocate(work_int3)
      self%angl_list_fep = c_null_ptr
    end if
    if (c_associated(self%dihe_list_fep)) then
      call C_F_pointer(self%dihe_list_fep, work_int3, [4, self%ndihe_fep_max, 5])
      deallocate(work_int3)
      self%dihe_list_fep = c_null_ptr
    end if
    if (c_associated(self%impr_list_fep)) then
      call C_F_pointer(self%impr_list_fep, work_int3, [4, self%nimpr_fep_max, 5])
      deallocate(work_int3)
      self%impr_list_fep = c_null_ptr
    end if
    if (c_associated(self%cmap_list_fep)) then
      call C_F_pointer(self%cmap_list_fep, work_int3, [8, self%ncmap_fep_max, 5])
      deallocate(work_int3)
      self%cmap_list_fep = c_null_ptr
    end if
    if (c_associated(self%id_singleA)) then
      call C_F_pointer(self%id_singleA, work_int, [self%size_id_singleA])
      deallocate(work_int)
      self%id_singleA = c_null_ptr
    end if
    if (c_associated(self%id_singleB)) then
      call C_F_pointer(self%id_singleB, work_int, [self%size_id_singleB])
      deallocate(work_int)
      self%id_singleB = c_null_ptr
    end if
    if (c_associated(self%fepgrp)) then
      call C_F_pointer(self%fepgrp, work_int, [self%size_fepgrp])
      deallocate(work_int)
      self%fepgrp = c_null_ptr
    end if
    if (c_associated(self%fepgrp_bond)) then
      call C_F_pointer(self%fepgrp_bond, work_int2, [5, 5])
      deallocate(work_int2)
      self%fepgrp_bond = c_null_ptr
    end if
    if (c_associated(self%fepgrp_angl)) then
      call C_F_pointer(self%fepgrp_angl, work_int3, [5, 5, 5])
      deallocate(work_int3)
      self%fepgrp_angl = c_null_ptr
    end if
    if (c_associated(self%fepgrp_dihe)) then
      call C_F_pointer(self%fepgrp_dihe, work_int4, [5, 5, 5, 5])
      deallocate(work_int4)
      self%fepgrp_dihe = c_null_ptr
    end if
    if (c_associated(self%fepgrp_cmap)) then
      call C_F_pointer(self%fepgrp_cmap, work_int, [5*5*5*5*5*5*5*5])
      deallocate(work_int)
      self%fepgrp_cmap = c_null_ptr
    end if
  end subroutine deallocate_s_molecule_c

  subroutine allocate_s_molecule_c(self) &
      bind(C, name="allocate_s_molecule_c")
    use conv_f_c_util
    implicit none
    type(s_molecule_c), intent(inout) :: self

    self%atom_no = allocate_c_int_array(self%num_atoms)
    self%segment_name = allocate_c_str_array(self%num_atoms, 4)
    self%segment_no = allocate_c_int_array(self%num_atoms)
    self%residue_no    = allocate_c_int_array(self%num_atoms)
    self%residue_c_no  = allocate_c_int_array(self%num_atoms)
    self%residue_name  = allocate_c_str_array(self%num_atoms, 6)
    self%atom_name     = allocate_c_str_array(self%num_atoms, 4)
    self%atom_cls_name = allocate_c_str_array(self%num_atoms, 6)
    self%atom_cls_no   = allocate_c_int_array(self%num_atoms)
    self%charge            = allocate_c_double_array(self%num_atoms)
    self%mass              = allocate_c_double_array(self%num_atoms)
    self%inv_mass          = allocate_c_double_array(self%num_atoms)
    self%imove             = allocate_c_int_array(self%num_atoms)
    self%stokes_radius     = allocate_c_double_array(self%num_atoms)
    self%inv_stokes_radius = allocate_c_double_array(self%num_atoms)
    self%chain_id = allocate_c_str_array(self%num_atoms, 1)
    self%atom_coord = allocate_c_double_array(self%num_atoms * 3)
    self%atom_occupancy = allocate_c_double_array(self%num_atoms)
    self%atom_temp_factor = allocate_c_double_array(self%num_atoms)
    self%atom_velocity = allocate_c_double_array(self%num_atoms * 3)
    self%light_atom_name = allocate_c_bool_array(self%num_atoms)
    self%light_atom_mass = allocate_c_bool_array(self%num_atoms)
    self%molecule_no = allocate_c_int_array(self%num_atoms)
    self%bond_list = allocate_c_int_array(self%num_bonds * 2)
    self%enm_list = allocate_c_int_array(self%num_enm_bonds * 2)
    self%angl_list = allocate_c_int_array(self%num_angles * 3)
    self%dihe_list = allocate_c_int_array(self%num_dihedrals * 4)
    self%impr_list = allocate_c_int_array(self%num_impropers * 4)
    self%cmap_list = allocate_c_int_array(self%num_cmaps * 8)
    self%molecule_atom_no = allocate_c_int_array(self%num_molecules)
    self%molecule_mass = allocate_c_double_array(self%num_molecules)
    self%molecule_name = allocate_c_str_array(self%num_molecules, 10)
    self%atom_refcoord = allocate_c_double_array(self%num_atoms * 3)
    self%atom_fitcoord = allocate_c_double_array(self%num_atoms * 3)
    self%pc_mode = allocate_c_double_array(self%num_pc_modes)

    self%num_atoms_fep = allocate_c_int_array(5)
    self%num_bonds_fep = allocate_c_int_array(5)
    self%num_angles_fep = allocate_c_int_array(5)
    self%num_dihedrals_fep = allocate_c_int_array(5)
    self%num_impropers_fep = allocate_c_int_array(5)
    self%num_cmaps_fep = allocate_c_int_array(5)
    self%bond_list_fep = allocate_c_int_array(5 * self%nbnd_fep_max * 2)
    self%angl_list_fep = allocate_c_int_array(5 * self%nangl_fep_max * 3)
    self%dihe_list_fep = allocate_c_int_array(5 * self%ndihe_fep_max * 4)
    self%impr_list_fep = allocate_c_int_array(5 * self%nimpr_fep_max * 4)
    self%cmap_list_fep = allocate_c_int_array(5 * self%ncmap_fep_max * 8)
    self%id_singleA = allocate_c_int_array(self%size_id_singleA)
    self%id_singleB = allocate_c_int_array(self%size_id_singleB)
    self%fepgrp = allocate_c_int_array(self%size_fepgrp)
    self%fepgrp_bond = allocate_c_int_array(5 * 5)
    self%fepgrp_angl = allocate_c_int_array(5 * 5 * 5)
    self%fepgrp_dihe = allocate_c_int_array(5 * 5 * 5 * 5)
    self%fepgrp_cmap = allocate_c_int_array(5*5*5*5*5*5*5*5)
  end subroutine allocate_s_molecule_c
end module s_molecule_c_mod
