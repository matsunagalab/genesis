import ctypes

# ライブラリをロード
lib = ctypes.CDLL('./python_interface.la')

# 関数のプロトタイプを定義
lib.define_molecule_from_pdb.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_void_p), ctypes.POINTER(ctypes.c_void_p)]
lib.define_molecule_from_pdb.restype = None

# 関数を呼び出す
pdb_filename = b"molecule.pdb"
num_atoms = ctypes.c_int()
atom_names_ptr = ctypes.c_void_p()
atom_coords_ptr = ctypes.c_void_p()

lib.define_molecule_from_pdb(pdb_filename, ctypes.byref(num_atoms), ctypes.byref(atom_names_ptr), ctypes.byref(atom_coords_ptr))

# 結果を処理する
# ...

# メモリを解放する
# ...

