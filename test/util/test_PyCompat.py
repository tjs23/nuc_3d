# coding=utf-8

import unittest

# Don't import * since PyCompat redefines `str`
from util import PyCompat

class TestPyCompat(unittest.TestCase):
  b_ascii = b'foo'
  s_ascii = u'foo'
  b_uc = b'Oh no, a m\xc3\xb8\xc3\xb8se!'
  s_uc = u'Oh no, a møøse!'

  def test_bytes(self):
    self.assertTrue(isinstance(self.b_ascii, bytes))
    self.assertTrue(isinstance(self.b_uc, bytes))
    self.assertFalse(isinstance(self.s_ascii, bytes))
    self.assertFalse(isinstance(self.s_uc, bytes))

  def test_unicode(self):
    self.assertTrue(isinstance(self.s_ascii, PyCompat.str))
    self.assertTrue(isinstance(self.s_uc, PyCompat.str))
    self.assertEqual(PyCompat.str(u'Thís'), u'Thís')

  def test_tobytes(self):
    self.assertEqual(PyCompat.tobytes(self.s_uc, 'utf-8'), self.b_uc)
    self.assertEqual(PyCompat.tobytes(self.s_ascii, 'utf-8'), self.b_ascii)
