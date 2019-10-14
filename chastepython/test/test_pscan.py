from __future__ import print_function
import unittest
import numpy as np
from chastepython.util.pscan import Scan


class TestScan(unittest.TestCase):
    def test_error_if_param_val_not_iterable(self):
        pass

    def test_error_if_joint_size_wrong(self):
        pass

    def test_docs_example(self):
        ans = [[0.0, 0.1, 10.0],
               [0.0, 0.1, 10.0],
               [2.5, -0.05, 10.0],
               [2.5, -0.05, 10.0],
               [5.0, -0.2, 10.0],
               [5.0, -0.2, 10.0],
               [7.5, -0.35, 10.0],
               [7.5, -0.35, 10.0],
               [10.0, -0.5, 10.0],
               [10.0, -0.5, 10.0],
               [0.0, 0.1, 15.0],
               [2.5, -0.05, 15.0],
               [5.0, -0.2, 15.0],
               [7.5, -0.35, 15.0],
               [10.0, -0.5, 15.0],
               [0.0, 0.1, 20.0],
               [2.5, -0.05, 20.0],
               [5.0, -0.2, 20.0],
               [7.5, -0.35, 20.0],
               [10.0, -0.5, 20.0], ]
        output = []

        f = lambda a, b, c: output.append([a, b, c])
        p = {}
        p['a'] = np.linspace(0, 10, 5)
        p['b'] = np.linspace(0.1, -0.5, 5)
        p['c'] = np.linspace(10, 20, 3)
        default_count = lambda p: 2
        # None specifies to default to previous value
        big_c_count = lambda p: 1 if p['c'] >= 14 else None

        s = Scan.from_dict(p, joint_lists=[['a', 'b']])
        s.add_count(default_count)  # these will be called
        s.add_count(big_c_count)  # in the order they were added
        s.run_scan(f)  # verify scan validity and run

        num_expanded = len(output)
        num_vars = 3
        self.assertEqual(num_expanded, len(ans), "Length not equal")
        for i in range(num_expanded):
            for j in range(num_vars):
                self.assertAlmostEqual(ans[i][j], output[i][j])
