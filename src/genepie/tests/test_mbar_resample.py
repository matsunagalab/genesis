"""Tests for MBAR weight-based trajectory resampling.

Two groups:

* **Synthetic** (always run, no downloaded data, no MBAR solve): exercise the
  pure resampling step and the DCD round-trip. These are the CI-safe tests.
* **Reference regression** (skipped unless the T-REMD trialanine tutorial data
  has been downloaded): run MBAR on the real energies, compare the weights and
  free energies against the paper reference, and check that resampling by the
  weights reproduces the reweighted phi/psi distribution.
"""
import os
import tempfile
import unittest

import numpy as np

from ..analysis.resample import MbarResampleResult, resample_indices
from ..exceptions import GenesisValidationError
from ..s_molecule import SMolecule
from ..s_trajectories import STrajectories
from . import conftest as C

try:
    import mdtraj  # noqa: F401
    HAS_MDTRAJ = True
except ImportError:
    HAS_MDTRAJ = False


def _tremd_available() -> bool:
    return (C.TREMD_DIR / "README.md").is_file()


# Umbrella-sampling regression data (source repo only): 61 windows spaced 3 deg
# apart along a periodic dihedral, 100 samples each. Used to check that MBAR
# yields the *bias-free* ensemble weights and that resampling by them removes
# the umbrella bias.
_REG_ROOT = C.TEST_DIR.parent.parent.parent / "tests" / "regression_test"
UMB_DIR = _REG_ROOT / "test_analysis" / "trajectories" / "umbrella_1d"
UMB_FENE_REF = (_REG_ROOT / "test_analysis" / "test_mbar_analysis"
                / "umbrella_1d" / "fene.dat.ref")
UMB_NREPLICA = 61
UMB_CONSTANT = (0.06092,) * UMB_NREPLICA
UMB_REFERENCE = tuple(3.0 * i for i in range(UMB_NREPLICA))


def _umbrella_available() -> bool:
    return (UMB_DIR / "1.dat").is_file() and UMB_FENE_REF.is_file()


# ---------------------------------------------------------------------------
# Synthetic tests (CI-safe)
# ---------------------------------------------------------------------------
class TestResampleIndices(unittest.TestCase):
    """The pure weights -> indices step, checked against theory."""

    def test_flatten_order_is_state_then_step(self):
        """A single nonzero weight must always be the drawn index.

        weights[state, step] flattens in C order, so the only nonzero entry at
        (state=1, step=0) of a (2, 2) array is flat index 1*2 + 0 = 2.
        """
        weights = np.array([[0.0, 0.0], [1.0, 0.0]])
        idx = resample_indices(weights, n_samples=100, seed=0)
        self.assertTrue(np.all(idx == 2))

    def test_empirical_frequencies_match_weights(self):
        """Draw counts must converge to the (normalized) weights."""
        weights = np.array([[0.5, 0.1], [0.3, 0.1]])  # already sums to 1
        p = weights.ravel()
        n = 400_000
        idx = resample_indices(weights, n_samples=n, seed=12345)
        freq = np.bincount(idx, minlength=p.size) / n
        # 400k draws: sampling error on each bin is well under 0.01.
        np.testing.assert_allclose(freq, p, atol=0.01)

    def test_unnormalized_weights_are_normalized(self):
        weights = np.array([2.0, 6.0, 2.0])  # sums to 10
        n = 200_000
        idx = resample_indices(weights, n_samples=n, seed=7)
        freq = np.bincount(idx, minlength=3) / n
        np.testing.assert_allclose(freq, [0.2, 0.6, 0.2], atol=0.01)

    def test_default_n_samples_is_total(self):
        weights = np.ones((3, 4))
        idx = resample_indices(weights, seed=1)
        self.assertEqual(idx.shape[0], 12)

    def test_rejects_empty(self):
        with self.assertRaises(GenesisValidationError):
            resample_indices(np.array([]))

    def test_rejects_negative(self):
        with self.assertRaises(GenesisValidationError):
            resample_indices(np.array([0.5, -0.1, 0.6]))

    def test_rejects_zero_sum(self):
        with self.assertRaises(GenesisValidationError):
            resample_indices(np.zeros(5))

    def test_seed_is_reproducible(self):
        weights = np.array([0.25, 0.25, 0.25, 0.25])
        a = resample_indices(weights, n_samples=1000, seed=99)
        b = resample_indices(weights, n_samples=1000, seed=99)
        np.testing.assert_array_equal(a, b)


@unittest.skipUnless(HAS_MDTRAJ, "mdtraj is required for save_dcd")
class TestSaveDcdRoundTrip(unittest.TestCase):
    """STrajectories.save_dcd -> reload must preserve atom/frame counts."""

    def test_save_and_reload(self):
        mol = SMolecule.from_file(pdb=str(C.BPTI_PDB), psf=str(C.BPTI_PSF),
                                  ref=str(C.BPTI_PDB))
        natom = mol.num_atoms
        nframe = 5
        rng = np.random.default_rng(0)
        coords = rng.random((nframe, natom, 3)) * 20.0  # angstrom
        traj = STrajectories.from_numpy(coords)

        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "roundtrip.dcd")
            traj.save_dcd(out, mol)
            self.assertTrue(os.path.isfile(out))
            reloaded = mdtraj.load(out, top=str(C.BPTI_PDB))
            self.assertEqual(reloaded.n_atoms, natom)
            self.assertEqual(reloaded.n_frames, nframe)
            # mdtraj stores nm; save converted angstrom->nm, so /10 must match.
            np.testing.assert_allclose(reloaded.xyz, coords / 10.0, atol=1e-3)

    def test_atom_count_mismatch_raises(self):
        mol = SMolecule.from_file(pdb=str(C.BPTI_PDB), psf=str(C.BPTI_PSF),
                                  ref=str(C.BPTI_PDB))
        coords = np.zeros((2, mol.num_atoms + 1, 3))
        traj = STrajectories.from_numpy(coords)
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(GenesisValidationError):
                traj.save_dcd(os.path.join(d, "x.dcd"), mol)


# ---------------------------------------------------------------------------
# Reference regression tests (require downloaded T-REMD tutorial data)
# ---------------------------------------------------------------------------
@unittest.skipUnless(
    _tremd_available(),
    "T-REMD tutorial data not found. Run: "
    "python -m genepie.tests.download_tremd_data")
class TestMbarWeightsReference(unittest.TestCase):
    """MBAR weights / free energies must match the paper reference."""

    @classmethod
    def setUpClass(cls):
        from .. import genesis_exe
        cls.result = genesis_exe.mbar_analysis_isolated(
            cvfile=C.TREMD_POT_PATTERN,
            nreplica=C.TREMD_NREPLICA,
            input_type="EneSingle",
            dimension=1,
            temperature=C.TREMD_TEMPERATURES,
            target_temperature=C.TREMD_TARGET_TEMPERATURE,
            tolerance=1e-8,
            self_iteration=100,
            newton_iteration=10,
            return_weights=True,
            timeout=600.0,
        )

    def test_shapes(self):
        r = self.result
        self.assertEqual(r.n_replica, C.TREMD_NREPLICA)
        self.assertEqual(r.weights.shape, (C.TREMD_NREPLICA, r.n_step))
        self.assertEqual(r.fene.shape, (C.TREMD_NREPLICA, 1))

    def test_weights_sum_to_one(self):
        self.assertAlmostEqual(float(self.result.weights.sum()), 1.0, places=6)
        self.assertTrue(np.all(self.result.weights >= 0))

    def test_fene_matches_reference(self):
        ref = np.loadtxt(str(C.TREMD_FENE_REF))
        got = self.result.fene[:, 0]
        # The two solves converge to the same optimum; leave headroom for BLAS.
        np.testing.assert_allclose(got, ref, rtol=1e-4, atol=1e-4)

    def test_weights_match_reference(self):
        for k in range(1, C.TREMD_NREPLICA + 1):
            ref = np.loadtxt(C.TREMD_WEIGHT_PATTERN.format(k))[:, 1]
            got = self.result.weights[k - 1]
            self.assertEqual(got.shape, ref.shape)
            # Weights span ~80 orders of magnitude; compare relatively with a
            # small floor so vanishing hot-state weights do not dominate.
            np.testing.assert_allclose(got, ref, rtol=1e-5, atol=1e-12)


@unittest.skipUnless(
    _tremd_available(),
    "T-REMD tutorial data not found. Run: "
    "python -m genepie.tests.download_tremd_data")
class TestResampleReproducesReweighting(unittest.TestCase):
    """Resampling by the weights must reproduce the reweighted phi/psi PMF.

    Uses the reference weights and precomputed dihedrals directly (no MBAR
    solve, no DCD load), so it is a fast, focused check that the resampling
    turns weighted samples into an equivalent uniform-weight ensemble.
    """

    def _load_phi_and_weights(self):
        phi = []
        w = []
        for k in range(1, C.TREMD_NREPLICA + 1):
            tor = np.loadtxt(C.TREMD_TOR_PATTERN.format(k))
            phi.append(tor[:, 1])  # column layout: index phi psi
            w.append(np.loadtxt(C.TREMD_WEIGHT_PATTERN.format(k))[:, 1])
        return np.concatenate(phi), np.stack(w)  # phi (N,), w (nrep, nstep)

    def test_resampled_histogram_matches_weighted(self):
        phi, weights = self._load_phi_and_weights()
        w_flat = weights.ravel()
        w_flat = w_flat / w_flat.sum()

        bins = np.linspace(-180.0, 180.0, 19)  # 18 bins
        weighted_hist, _ = np.histogram(phi, bins=bins, weights=w_flat)

        idx = resample_indices(weights, n_samples=phi.size, seed=2024)
        plain_hist, _ = np.histogram(phi[idx], bins=bins)
        plain_hist = plain_hist / plain_hist.sum()

        # Plain histogram of the resampled ensemble must reproduce the
        # weight-reweighted histogram of the original samples.
        np.testing.assert_allclose(plain_hist, weighted_hist, atol=0.01)


@unittest.skipUnless(
    _tremd_available() and HAS_MDTRAJ,
    "T-REMD tutorial data and mdtraj are required")
class TestResampleTrajectoryEndToEnd(unittest.TestCase):
    """Full pipeline: MBAR weights -> resample DCD frames -> save DCD."""

    def test_pipeline(self):
        from .. import genesis_exe
        mol = SMolecule.from_file(pdb=str(C.TREMD_PDB), psf=str(C.TREMD_PSF),
                                  ref=str(C.TREMD_PDB))
        dcds = [C.TREMD_DCD_PATTERN.format(k)
                for k in range(1, C.TREMD_NREPLICA + 1)]
        n_samples = 3000
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "resampled.dcd")
            res = genesis_exe.mbar_resample_trajectory(
                molecule=mol,
                dcd_files=dcds,
                cvfile=C.TREMD_POT_PATTERN,
                nreplica=C.TREMD_NREPLICA,
                temperature=C.TREMD_TEMPERATURES,
                target_temperature=C.TREMD_TARGET_TEMPERATURE,
                input_type="EneSingle",
                self_iteration=100,
                newton_iteration=10,
                selection="all",
                n_samples=n_samples,
                seed=42,
                output_dcd=out,
                isolated=True,
            )
            self.assertIsInstance(res, MbarResampleResult)
            self.assertEqual(res.indices.shape[0], n_samples)
            self.assertEqual(res.trajectory.nframe, n_samples)
            self.assertEqual(res.trajectory.natom, mol.num_atoms)
            self.assertEqual(res.molecule.num_atoms, mol.num_atoms)

            reloaded = mdtraj.load(out, top=str(C.TREMD_PDB))
            self.assertEqual(reloaded.n_frames, n_samples)
            self.assertEqual(reloaded.n_atoms, mol.num_atoms)

        # Target temperature is the coldest state, so draws must concentrate in
        # the low-index (cold) states and essentially never hit the hottest.
        state_of_draw = res.indices // res.weights.shape[1]
        counts = np.bincount(state_of_draw, minlength=C.TREMD_NREPLICA)
        self.assertGreater(counts[0], counts[-1])
        self.assertEqual(counts[-1], 0)

    def test_mismatched_dcd_count_raises(self):
        mol = SMolecule.from_file(pdb=str(C.TREMD_PDB), psf=str(C.TREMD_PSF),
                                  ref=str(C.TREMD_PDB))
        from .. import genesis_exe
        with self.assertRaises(GenesisValidationError):
            genesis_exe.mbar_resample_trajectory(
                molecule=mol,
                dcd_files=[C.TREMD_DCD_PATTERN.format(1)],  # only 1, need 20
                cvfile=C.TREMD_POT_PATTERN,
                nreplica=C.TREMD_NREPLICA,
                temperature=C.TREMD_TEMPERATURES,
                target_temperature=C.TREMD_TARGET_TEMPERATURE,
                input_type="EneSingle",
            )


@unittest.skipUnless(
    _umbrella_available(),
    "umbrella_1d regression data not found (source repo only).")
class TestUmbrellaBiasFreeWeights(unittest.TestCase):
    """MBAR must return the *bias-free* ensemble weights for umbrella sampling,
    and resampling by them must reproduce the reweighted (unbiased) CV
    distribution.

    This is the umbrella-sampling counterpart of the T-REMD tests: instead of
    reweighting to a target temperature, MBAR reweights every biased window to
    the restraint-free ensemble.
    """

    @classmethod
    def setUpClass(cls):
        from .. import genesis_exe
        cls.result = genesis_exe.mbar_analysis_isolated(
            cvfile=str(UMB_DIR / "{}.dat"),
            nreplica=UMB_NREPLICA,
            input_type="US",
            dimension=1,
            temperature=300.0,
            target_temperature=300.0,
            tolerance=1e-8,
            rest_function=(1,),
            grids=((-1.0, 181.0, 81),),
            constant=(UMB_CONSTANT,),
            reference=(UMB_REFERENCE,),
            is_periodic=(True,),
            box_size=(360.0,),
            return_weights=True,
            timeout=300.0,
        )
        # Per-sample reaction coordinate (dihedral, deg) in (window, step) order.
        cv = [np.loadtxt(str(UMB_DIR / f"{k}.dat"))[:, 1]
              for k in range(1, UMB_NREPLICA + 1)]
        cls.cv = np.concatenate(cv)

    def test_shapes_sum_nonnegative(self):
        w = self.result.weights
        self.assertEqual(w.shape[0], UMB_NREPLICA)
        self.assertAlmostEqual(float(w.sum()), 1.0, places=6)
        self.assertTrue(np.all(w >= 0))
        self.assertTrue(np.all(np.isfinite(w)))

    def test_fene_matches_reference(self):
        ref = np.loadtxt(str(UMB_FENE_REF))
        if ref.ndim == 2:
            ref = ref[:, 0]
        np.testing.assert_allclose(self.result.fene[:, 0], ref,
                                   rtol=1e-5, atol=1e-5)

    # The reaction coordinate is a periodic dihedral spanning [-180, 180] deg.
    CV_BINS = np.linspace(-180.0, 180.0, 37)  # 10-degree bins

    def test_reweighting_removes_the_bias(self):
        """Each biased window is localized to a few degrees, but the reweighted
        (unbiased) ensemble is far broader - the hallmark of bias removal."""
        w = self.result.weights.ravel()
        w = w / w.sum()
        reweighted, _ = np.histogram(self.cv, bins=self.CV_BINS, weights=w)
        reweighted_bins = int(np.count_nonzero(reweighted > 1e-9))

        widest_window = max(
            int(np.count_nonzero(
                np.histogram(np.loadtxt(str(UMB_DIR / f"{k}.dat"))[:, 1],
                             bins=self.CV_BINS)[0]))
            for k in range(1, UMB_NREPLICA + 1))

        # The unbiased ensemble must populate many more bins than any single
        # biased window (windows here cover ~2-3 bins; reweighted covers ~14).
        self.assertGreater(reweighted_bins, 3 * widest_window)

    def test_resample_reproduces_reweighted_distribution(self):
        """Plain histogram of the resampled CV must match the weighted one."""
        weights = self.result.weights
        w_flat = weights.ravel()
        w_flat = w_flat / w_flat.sum()
        weighted_hist, _ = np.histogram(self.cv, bins=self.CV_BINS,
                                        weights=w_flat)

        idx = resample_indices(weights, n_samples=self.cv.size * 20, seed=7)
        plain_hist, _ = np.histogram(self.cv[idx], bins=self.CV_BINS)
        plain_hist = plain_hist / plain_hist.sum()

        np.testing.assert_allclose(plain_hist, weighted_hist, atol=0.01)
