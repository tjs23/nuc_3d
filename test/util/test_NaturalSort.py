import unittest

from util.NaturalSort import *

class TestNaturalKey(unittest.TestCase):
  def test_simple(self):
    l = ['two', '3', '2', '1']
    self.assertEqual(sorted(l, key=naturalKey()), list(reversed(l)))
    self.assertEqual(sorted(l, key=integerKey), list(reversed(l)))
    self.assertEqual(sorted(l, key=decimalKey), list(reversed(l)))
    
    l = ['asd1000asd', 'asd100baf']
    self.assertEqual(sorted(l, key=naturalKey()), list(reversed(l)))

    l = ['asd1000asd', '100baf']
    self.assertEqual(sorted(l, key=naturalKey()), list(reversed(l)))

  def test_decimal(self):
    decimal = ["the1.9th cat", "the1.80th cat"]
    self.assertEqual(sorted(decimal, key=naturalKey(decimals=True)), list(reversed(decimal)))
    self.assertEqual(sorted(decimal, key=decimalKey), list(reversed(decimal)))
    self.assertEqual(sorted(decimal, key=integerKey), decimal)
  
    decimal = ["the1.9th cat", "the1.80.6th cat"]
    self.assertEqual(sorted(decimal, key=decimalKey), list(reversed(decimal)))
