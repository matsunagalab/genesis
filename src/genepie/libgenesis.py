import ctypes
import os
import threading
from .libloader import load_genesis_lib
from .s_molecule_c import SMoleculeC
from .s_trajectories_c import STrajectoriesC

class LibGenesis:
    """Thread-safe singleton for GENESIS library access.

    Note: While the singleton pattern is thread-safe, the underlying
    Fortran library uses module-level save pointers and is NOT thread-safe.
    All GENESIS function calls should be made from a single thread.
    """

    _instance = None
    _lock = threading.Lock()

    def __new__(cls):
        # Double-check locking pattern for thread safety
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    instance = super().__new__(cls)
                    instance.lib = load_genesis_lib()
                    instance._initialized = True
                    instance._setup_function_signatures()
                    cls._instance = instance
        return cls._instance

    def __init__(self):
        # All initialization done in __new__ to avoid race conditions
        pass

    def _setup_function_signatures(self):

        self.lib.define_molecule_from_file.argtypes = [
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.c_char_p,
                ctypes.POINTER(SMoleculeC),
                ]
        self.lib.define_molecule_from_file.restype = ctypes.c_int

        self.lib.deallocate_s_molecule_c.argtypes = [
                ctypes.POINTER(SMoleculeC)]
        self.lib.deallocate_s_molecule_c.restype = None

        # Atom selection using GENESIS selection syntax
        self.lib.selection_c.argtypes = [
                ctypes.POINTER(SMoleculeC),     # molecule_c
                ctypes.c_char_p,                # selection_str
                ctypes.c_int,                   # selection_len
                ctypes.POINTER(ctypes.c_void_p),  # indices (output)
                ctypes.POINTER(ctypes.c_int),   # n_indices (output)
                ctypes.POINTER(ctypes.c_int),   # status (output)
                ctypes.c_char_p,                # msg (output)
                ctypes.c_int]                   # msglen
        self.lib.selection_c.restype = None

        self.lib.deallocate_selection_c.argtypes = []
        self.lib.deallocate_selection_c.restype = None

        self.lib.allocate_s_molecule_c.argtypes = [
                ctypes.POINTER(SMoleculeC)]
        self.lib.allocate_s_molecule_c.restype = None

        self.lib.deallocate_s_trajectories_c.argtypes = [
                ctypes.POINTER(STrajectoriesC)]
        self.lib.deallocate_s_trajectories_c.restype = None

        self.lib.allocate_s_trajectories_c_array.argtypes = [
                ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
        self.lib.allocate_s_trajectories_c_array.restype = None

        self.lib.deallocate_s_trajectories_c_array.argtypes = [
                ctypes.c_void_p, ctypes.POINTER(ctypes.c_int)]
        self.lib.deallocate_s_trajectories_c_array.restype = None

        # crd_convert info (zerocopy phase 1: get trajectory metadata)
        self.lib.crd_convert_info_c.argtypes = [
                ctypes.POINTER(SMoleculeC),       # molecule_c
                ctypes.c_char_p,                  # trj_filenames (packed)
                ctypes.c_int,                     # n_trj_files
                ctypes.c_int,                     # filename_len
                ctypes.c_int,                     # trj_format
                ctypes.c_int,                     # trj_type
                ctypes.POINTER(ctypes.c_void_p),  # frame_counts_ptr (output)
                ctypes.POINTER(ctypes.c_int),     # n_trajs (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.crd_convert_info_c.restype = None

        # crd_convert zerocopy (phase 2: fill pre-allocated arrays)
        self.lib.crd_convert_zerocopy_c.argtypes = [
                ctypes.POINTER(SMoleculeC),       # molecule_c
                ctypes.c_char_p,                  # trj_filenames (packed)
                ctypes.c_int,                     # n_trj_files
                ctypes.c_int,                     # filename_len
                ctypes.c_int,                     # trj_format
                ctypes.c_int,                     # trj_type
                ctypes.c_void_p,                  # selected_indices
                ctypes.c_int,                     # n_selected
                ctypes.c_int,                     # fitting_method
                ctypes.c_void_p,                  # fitting_indices
                ctypes.c_int,                     # n_fitting
                ctypes.c_int,                     # mass_weighted
                ctypes.c_int,                     # do_centering
                ctypes.c_void_p,                  # centering_indices
                ctypes.c_int,                     # n_centering
                ctypes.c_void_p,                  # center_coord (double[3])
                ctypes.c_int,                     # pbcc_mode
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # frame_counts
                ctypes.c_void_p,                  # coords_ptrs
                ctypes.c_void_p,                  # pbc_box_ptrs
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.crd_convert_zerocopy_c.restype = None

        self.lib.deallocate_frame_counts_c.argtypes = [ctypes.c_void_p]
        self.lib.deallocate_frame_counts_c.restype = None

        # TRJ analysis (zerocopy, pre-allocated result arrays)
        self.lib.trj_analysis_c.argtypes = [
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # dist_list_ptr
                ctypes.c_int,                     # n_dist
                ctypes.c_void_p,                  # angl_list_ptr
                ctypes.c_int,                     # n_angl
                ctypes.c_void_p,                  # tors_list_ptr
                ctypes.c_int,                     # n_tors
                ctypes.c_void_p,                  # dist_ptr (pre-allocated)
                ctypes.c_int,                     # dist_size
                ctypes.c_void_p,                  # angl_ptr (pre-allocated)
                ctypes.c_int,                     # angl_size
                ctypes.c_void_p,                  # tors_ptr (pre-allocated)
                ctypes.c_int,                     # tors_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.trj_analysis_c.restype = None

        # TRJ analysis with lazy DCD loading (memory efficient)
        self.lib.trj_analysis_lazy_c.argtypes = [
                ctypes.c_char_p,                  # dcd_filename
                ctypes.c_int,                     # filename_len
                ctypes.c_int,                     # trj_type
                ctypes.c_int,                     # physical DCD atom count
                ctypes.c_void_p,                  # source selection
                ctypes.c_int,                     # selected atom count
                ctypes.c_void_p,                  # dist_list_ptr
                ctypes.c_int,                     # n_dist
                ctypes.c_void_p,                  # angl_list_ptr
                ctypes.c_int,                     # n_angl
                ctypes.c_void_p,                  # tors_list_ptr
                ctypes.c_int,                     # n_tors
                ctypes.c_int,                     # n_atoms (selected)
                ctypes.c_int,                     # ana_period
                ctypes.c_int,                     # n_frame (result columns)
                ctypes.c_void_p,                  # dist_ptr (pre-allocated)
                ctypes.c_int,                     # dist_size
                ctypes.c_void_p,                  # angl_ptr (pre-allocated)
                ctypes.c_int,                     # angl_size
                ctypes.c_void_p,                  # tors_ptr (pre-allocated)
                ctypes.c_int,                     # tors_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # dcd_nframe_out (output)
                ctypes.POINTER(ctypes.c_int),     # dcd_natom_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.trj_analysis_lazy_c.restype = None

        # TRJ analysis with COM (zerocopy, pre-allocated result arrays)
        self.lib.trj_analysis_com_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                # Atom-based measurements
                ctypes.c_void_p,                  # dist_list_ptr
                ctypes.c_int,                     # n_dist
                ctypes.c_void_p,                  # angl_list_ptr
                ctypes.c_int,                     # n_angl
                ctypes.c_void_p,                  # tors_list_ptr
                ctypes.c_int,                     # n_tors
                # COM distance
                ctypes.c_void_p,                  # cdis_atoms_ptr
                ctypes.c_int,                     # n_cdis_atoms
                ctypes.c_void_p,                  # cdis_offsets_ptr
                ctypes.c_int,                     # n_cdis_offsets
                ctypes.c_void_p,                  # cdis_pairs_ptr
                ctypes.c_int,                     # n_cdis
                # COM angle
                ctypes.c_void_p,                  # cang_atoms_ptr
                ctypes.c_int,                     # n_cang_atoms
                ctypes.c_void_p,                  # cang_offsets_ptr
                ctypes.c_int,                     # n_cang_offsets
                ctypes.c_void_p,                  # cang_triplets_ptr
                ctypes.c_int,                     # n_cang
                # COM torsion
                ctypes.c_void_p,                  # ctor_atoms_ptr
                ctypes.c_int,                     # n_ctor_atoms
                ctypes.c_void_p,                  # ctor_offsets_ptr
                ctypes.c_int,                     # n_ctor_offsets
                ctypes.c_void_p,                  # ctor_quads_ptr
                ctypes.c_int,                     # n_ctor
                # Pre-allocated output arrays
                ctypes.c_void_p,                  # dist_ptr
                ctypes.c_int,                     # dist_size
                ctypes.c_void_p,                  # angl_ptr
                ctypes.c_int,                     # angl_size
                ctypes.c_void_p,                  # tors_ptr
                ctypes.c_int,                     # tors_size
                ctypes.c_void_p,                  # cdis_result_ptr
                ctypes.c_int,                     # cdis_size
                ctypes.c_void_p,                  # cang_result_ptr
                ctypes.c_int,                     # cang_size
                ctypes.c_void_p,                  # ctor_result_ptr
                ctypes.c_int,                     # ctor_size
                # Output
                ctypes.POINTER(ctypes.c_int),     # nstru_out
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.trj_analysis_com_c.restype = None

        # TRJ analysis with COM and lazy DCD loading (memory efficient)
        self.lib.trj_analysis_com_lazy_c.argtypes = [
                # Lazy DCD source
                ctypes.c_char_p,                  # dcd_filename
                ctypes.c_int,                     # filename_len
                ctypes.c_int,                     # trj_type
                ctypes.c_int,                     # physical DCD atom count
                ctypes.c_void_p,                  # source selection
                ctypes.c_int,                     # selected atom count
                # Mass + loop control
                ctypes.c_void_p,                  # mass_ptr (selected space)
                ctypes.c_int,                     # n_atoms (selected)
                ctypes.c_int,                     # ana_period
                ctypes.c_int,                     # n_frame (result columns)
                # Atom-based measurements
                ctypes.c_void_p,                  # dist_list_ptr
                ctypes.c_int,                     # n_dist
                ctypes.c_void_p,                  # angl_list_ptr
                ctypes.c_int,                     # n_angl
                ctypes.c_void_p,                  # tors_list_ptr
                ctypes.c_int,                     # n_tors
                # COM distance
                ctypes.c_void_p,                  # cdis_atoms_ptr
                ctypes.c_int,                     # n_cdis_atoms
                ctypes.c_void_p,                  # cdis_offsets_ptr
                ctypes.c_int,                     # n_cdis_offsets
                ctypes.c_void_p,                  # cdis_pairs_ptr
                ctypes.c_int,                     # n_cdis
                # COM angle
                ctypes.c_void_p,                  # cang_atoms_ptr
                ctypes.c_int,                     # n_cang_atoms
                ctypes.c_void_p,                  # cang_offsets_ptr
                ctypes.c_int,                     # n_cang_offsets
                ctypes.c_void_p,                  # cang_triplets_ptr
                ctypes.c_int,                     # n_cang
                # COM torsion
                ctypes.c_void_p,                  # ctor_atoms_ptr
                ctypes.c_int,                     # n_ctor_atoms
                ctypes.c_void_p,                  # ctor_offsets_ptr
                ctypes.c_int,                     # n_ctor_offsets
                ctypes.c_void_p,                  # ctor_quads_ptr
                ctypes.c_int,                     # n_ctor
                # Pre-allocated output arrays
                ctypes.c_void_p,                  # dist_ptr
                ctypes.c_int,                     # dist_size
                ctypes.c_void_p,                  # angl_ptr
                ctypes.c_int,                     # angl_size
                ctypes.c_void_p,                  # tors_ptr
                ctypes.c_int,                     # tors_size
                ctypes.c_void_p,                  # cdis_result_ptr
                ctypes.c_int,                     # cdis_size
                ctypes.c_void_p,                  # cang_result_ptr
                ctypes.c_int,                     # cang_size
                ctypes.c_void_p,                  # ctor_result_ptr
                ctypes.c_int,                     # ctor_size
                # Output
                ctypes.POINTER(ctypes.c_int),     # nstru_out
                ctypes.POINTER(ctypes.c_int),     # dcd_nframe_out
                ctypes.POINTER(ctypes.c_int),     # dcd_natom_out
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.trj_analysis_com_lazy_c.restype = None

        # RG analysis (pre-allocated result array)
        self.lib.rg_analysis_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr (pointer to NumPy array)
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # analysis_idx (pointer to int array)
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # mass_weighted (0 or 1)
                ctypes.c_void_p,                  # result_ptr (pre-allocated result)
                ctypes.c_int,                     # result_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rg_analysis_c.restype = None

        # RG analysis with lazy DCD loading (memory efficient)
        self.lib.rg_analysis_lazy_c.argtypes = [
                ctypes.c_char_p,                  # dcd_filename
                ctypes.c_int,                     # filename_len
                ctypes.c_int,                     # trj_type
                ctypes.c_int,                     # physical DCD atom count
                ctypes.c_void_p,                  # source selection
                ctypes.c_int,                     # selected atom count
                ctypes.c_void_p,                  # mass_ptr
                ctypes.c_int,                     # n_atoms
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # analysis_idx_ptr
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # mass_weighted
                ctypes.c_void_p,                  # result_ptr (pre-allocated)
                ctypes.c_int,                     # result_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # dcd_nframe_out (output)
                ctypes.POINTER(ctypes.c_int),     # dcd_natom_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rg_analysis_lazy_c.restype = None

        # RMSD analysis (no fitting, pre-allocated result array)
        self.lib.rmsd_analysis_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr
                ctypes.c_void_p,                  # ref_coord_ptr
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # analysis_idx
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # mass_weighted
                ctypes.c_void_p,                  # result_ptr (pre-allocated)
                ctypes.c_int,                     # result_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rmsd_analysis_c.restype = None

        # RMSD analysis with fitting (pre-allocated result array)
        self.lib.rmsd_analysis_fitting_c.argtypes = [
                ctypes.c_void_p,                  # mass_ptr
                ctypes.c_void_p,                  # ref_coord_ptr
                ctypes.c_int,                     # n_atoms
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # fitting_idx_ptr
                ctypes.c_int,                     # n_fitting
                ctypes.c_void_p,                  # analysis_idx_ptr
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # fitting_method
                ctypes.c_int,                     # mass_weighted
                ctypes.c_void_p,                  # result_ptr (pre-allocated)
                ctypes.c_int,                     # result_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rmsd_analysis_fitting_c.restype = None

        # RMSD analysis with lazy DCD loading (memory efficient)
        self.lib.rmsd_analysis_lazy_c.argtypes = [
                ctypes.c_char_p,                  # dcd_filename
                ctypes.c_int,                     # filename_len
                ctypes.c_int,                     # trj_type
                ctypes.c_int,                     # physical DCD atom count
                ctypes.c_void_p,                  # source selection
                ctypes.c_int,                     # selected atom count
                ctypes.c_void_p,                  # mass_ptr
                ctypes.c_void_p,                  # ref_coord_ptr
                ctypes.c_int,                     # n_atoms
                ctypes.c_int,                     # ana_period
                ctypes.c_void_p,                  # fitting_idx_ptr
                ctypes.c_int,                     # n_fitting
                ctypes.c_void_p,                  # analysis_idx_ptr
                ctypes.c_int,                     # n_analysis
                ctypes.c_int,                     # fitting_method
                ctypes.c_int,                     # mass_weighted
                ctypes.c_void_p,                  # result_ptr (pre-allocated)
                ctypes.c_int,                     # result_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # dcd_nframe_out (output)
                ctypes.POINTER(ctypes.c_int),     # dcd_natom_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.rmsd_analysis_lazy_c.restype = None

        # DRMS analysis (pre-allocated result array)
        self.lib.drms_analysis_c.argtypes = [
                ctypes.c_void_p,                  # contact_list_ptr
                ctypes.c_void_p,                  # contact_dist_ptr
                ctypes.c_int,                     # n_contact
                ctypes.POINTER(STrajectoriesC),   # s_trajes_c
                ctypes.c_int,                     # ana_period
                ctypes.c_int,                     # pbc_correct
                ctypes.c_void_p,                  # result_ptr (pre-allocated)
                ctypes.c_int,                     # result_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.drms_analysis_c.restype = None

        # DRMS analysis with lazy DCD loading (memory efficient)
        self.lib.drms_analysis_lazy_c.argtypes = [
                ctypes.c_char_p,                  # dcd_filename
                ctypes.c_int,                     # filename_len
                ctypes.c_int,                     # trj_type
                ctypes.c_int,                     # physical DCD atom count
                ctypes.c_void_p,                  # source selection
                ctypes.c_int,                     # selected atom count
                ctypes.c_void_p,                  # contact_list_ptr
                ctypes.c_void_p,                  # contact_dist_ptr
                ctypes.c_int,                     # n_contact
                ctypes.c_int,                     # n_atoms
                ctypes.c_int,                     # ana_period
                ctypes.c_int,                     # pbc_correct
                ctypes.c_void_p,                  # result_ptr (pre-allocated)
                ctypes.c_int,                     # result_size
                ctypes.POINTER(ctypes.c_int),     # nstru_out (output)
                ctypes.POINTER(ctypes.c_int),     # dcd_nframe_out (output)
                ctypes.POINTER(ctypes.c_int),     # dcd_natom_out (output)
                ctypes.POINTER(ctypes.c_int),     # status (output)
                ctypes.c_char_p,                  # msg (output)
                ctypes.c_int,                     # msglen
                ]
        self.lib.drms_analysis_lazy_c.restype = None

        self.lib.ma_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.ma_analysis_c.restype = None

        # Diffusion analysis (zerocopy, pre-allocated result arrays)
        self.lib.diffusion_analysis_c.argtypes = [
                ctypes.c_void_p,         # msd_ptr
                ctypes.c_int,            # ndata
                ctypes.c_int,            # ncols
                ctypes.c_double,         # time_step
                ctypes.c_double,         # distance_unit
                ctypes.c_int,            # ndofs
                ctypes.c_int,            # start_step
                ctypes.c_int,            # stop_step
                ctypes.c_void_p,         # out_data_ptr (pre-allocated)
                ctypes.c_int,            # out_data_size
                ctypes.c_void_p,         # diff_coeff_ptr (pre-allocated)
                ctypes.c_int,            # n_sets
                ctypes.POINTER(ctypes.c_int),  # status
                ctypes.c_char_p,              # msg
                ctypes.c_int,                 # msglen
        ]
        self.lib.diffusion_analysis_c.restype = None

        self.lib.hb_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.hb_analysis_c.restype = None

        self.lib.aa_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.aa_analysis_c.restype = None

        self.lib.wa_analysis_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.wa_analysis_c.restype = None

        self.lib.mbar_analysis_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.c_int,                     # return_weights
                ctypes.POINTER(ctypes.c_void_p),  # result_fene
                ctypes.POINTER(ctypes.c_int),     # n_replica
                ctypes.POINTER(ctypes.c_int),     # n_blocks
                ctypes.POINTER(ctypes.c_void_p),  # result_weights
                ctypes.POINTER(ctypes.c_int),     # n_weight_replica
                ctypes.POINTER(ctypes.c_int),     # n_weight_step
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.mbar_analysis_c.restype = None

        self.lib.pmf_analysis_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),  # result_pmf
                ctypes.POINTER(ctypes.c_int),     # n_out1 (rows)
                ctypes.POINTER(ctypes.c_int),     # n_out2 (cols)
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.pmf_analysis_c.restype = None

        self.lib.kc_analysis_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.POINTER(STrajectoriesC),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.kc_analysis_c.restype = None

        self.lib.export_pdb_to_string_c.argtypes = [
                ctypes.POINTER(SMoleculeC),
                ctypes.c_void_p,
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
                ctypes.c_int,
                ]
        self.lib.export_pdb_to_string_c.restype = None

        self.lib.allocate_c_int_array.argtypes = [
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.allocate_c_int_array.restype = ctypes.c_void_p

        self.lib.allocate_c_double_array.argtypes = [
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.allocate_c_double_array.restype = ctypes.c_void_p

        self.lib.allocate_c_double_array2.argtypes = [
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.allocate_c_double_array2.restype = ctypes.c_void_p

        self.lib.deallocate_int.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.deallocate_int.restype = None

        self.lib.deallocate_double.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.deallocate_double.restype = None

        self.lib.deallocate_double2.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ]
        self.lib.deallocate_double2.restype = None

        self.lib.deallocate_c_string.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.deallocate_c_string.restype = None

        self.lib.join_s_trajectories_c.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.join_s_trajectories_c.restype = None

        self.lib.deep_copy_s_trajectories_c.argtypes = [
                ctypes.POINTER(ctypes.c_void_p),
                ctypes.POINTER(ctypes.c_void_p),
                ]
        self.lib.deep_copy_s_trajectories_c.restype = None

        # ATDYN MD/Minimization functions
        self.lib.atdyn_md_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),  # result_energies
                ctypes.POINTER(ctypes.c_int),     # result_nframes
                ctypes.POINTER(ctypes.c_int),     # result_nterms
                ctypes.POINTER(ctypes.c_void_p),  # result_final_coords
                ctypes.POINTER(ctypes.c_int),     # result_natom
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.atdyn_md_c.restype = None

        self.lib.atdyn_min_c.argtypes = [
                ctypes.c_char_p,                  # ctrl_text
                ctypes.c_int,                     # ctrl_len
                ctypes.POINTER(ctypes.c_void_p),  # result_energies
                ctypes.POINTER(ctypes.c_int),     # result_nsteps
                ctypes.POINTER(ctypes.c_int),     # result_nterms
                ctypes.POINTER(ctypes.c_void_p),  # result_final_coords
                ctypes.POINTER(ctypes.c_int),     # result_natom
                ctypes.POINTER(ctypes.c_int),     # result_converged
                ctypes.POINTER(ctypes.c_double),  # result_final_gradient
                ctypes.POINTER(ctypes.c_int),     # status
                ctypes.c_char_p,                  # msg
                ctypes.c_int,                     # msglen
                ]
        self.lib.atdyn_min_c.restype = None

        self.lib.deallocate_atdyn_results_c.argtypes = []
        self.lib.deallocate_atdyn_results_c.restype = None

        # Reset atdyn global state for multiple sequential runs
        self.lib.reset_atdyn_state_c.argtypes = []
        self.lib.reset_atdyn_state_c.restype = None
