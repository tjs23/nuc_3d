import numpy as np
import unittest
import sys, os

sys.path.insert(1, os.path.join(os.path.dirname(__file__), ".."))

from os.path import dirname, abspath, exists

# Setup import path

thisDir = dirname(abspath(__file__))
sys.path.remove(thisDir)

nucDir = dirname(thisDir)
sys.path.append(nucDir)

from analyses.overlap import overlap


class TestOverlap(unittest.TestCase):

    def test_identity(self):
        """
        The overlap of something with itself should always be 1.
        """
        trackA = { '1' : np.array([[[1,10], [15,20]]], dtype=np.int32)}
        trackB = { '1' : np.array([[[5,8], [10,16], [18, 22]]], dtype=np.int32)}
        trackC = { '1' : np.array([[[1,10]]], dtype=np.int32)}

        self.assertEqual(overlap(trackA, trackA), 1.0)
        self.assertEqual(overlap(trackB, trackB), 1.0)
        self.assertEqual(overlap(trackC, trackC), 1.0)

    def test_A_in_B(self):
        trackA = { '1' : np.array([[[4,10]]], dtype=np.int32)}
        trackB = { '1' : np.array([[[1,15]]], dtype=np.int32)}
        self.assertEqual(overlap(trackA, trackB), 1.0)

    def test_A_overlap_RHS_B(self):
        trackA = { '1' : np.array([[[10,25]]], dtype=np.int32)}
        trackB = { '1' : np.array([[[1,15]]], dtype=np.int32)}
        trackC = { '1' : np.array([[[1,15], [20,25]]], dtype=np.int32)}

        self.assertAlmostEqual(overlap(trackA, trackB), 0.375)
        self.assertAlmostEqual(overlap(trackA, trackC), 0.75)

    def test_A_overlap_LHS_B(self):
        trackA = { '1' : np.array([[[10,25]]], dtype=np.int32)}
        trackB = { '1' : np.array([[[15,35]]], dtype=np.int32)}

        self.assertAlmostEqual(overlap(trackA, trackB), 0.6875)

    def test_B_in_A(self):
        trackA = { '1' : np.array([[[10,25]]], dtype=np.int32)}
        trackB = { '1' : np.array([[[15,20]]], dtype=np.int32)}
        self.assertAlmostEqual(overlap(trackA, trackB), 0.375)

if __name__ == '__main__':
    unittest.main()
