"""Assembly and provenance checks for the published extension."""
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import unittest
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'independent'/'bit-improvement'))
from assembly_options import evaluate, options


class Extension(unittest.TestCase):
    def test_pinned_sources_unmodified(self):
        manifest = json.loads((ROOT/'upstream-manifest.json').read_text())
        self.assertEqual(manifest['commit'], 'f2176bc1124821bf17eb63725bd366d7bdc020a3')
        for name, digest in manifest['unmodified_files_sha256'].items():
            self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(), digest, name)

    def test_both_assembly_profiles(self):
        high, low = options()
        self.assertEqual(Q(high['gap']), Q(23525148,10**17))
        self.assertEqual(Q(low['gap']), Q(3140869,10**15))
        self.assertTrue(high['ok'] and low['ok'])

    def test_negative_controls(self):
        args = ('negative', Q(993,1000), Q(1,500), Q(3668,10**8), Q(1,10**9))
        self.assertFalse(evaluate(*args, c1=Q(14694))['ok'])
        self.assertFalse(evaluate(*args, beta=Q(1,2))['ok'])
        self.assertFalse(evaluate(*args, x=Q(1))['ok'])
        self.assertFalse(evaluate('too large', Q(993,1000), Q(1,500), Q(1,10000), Q(1,10**9))['ok'])
        self.assertFalse(evaluate('arithmetic too costly', Q(993,1000), Q(1,100), Q(3668,10**8), Q(1,10**9))['ok'])

    def test_explicit_prefactor_bounds(self):
        for b in range(1,1000):
            self.assertLess(128*(b+1),768*b**3)
        self.assertGreater(Q(4)-Q(1,4),1)


if __name__=='__main__': unittest.main()
