import sys
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))
import numpy as np

from NucApi import Nucleus

def mergeContacts(nuc_file_a, nuc_file_b, nuc_file_c='Merged.nuc', remove_frac=None):
    
  nuc = Nucleus(nuc_file_a)
  nuc.saveAs(nuc_file_c)
  
  group_a = nuc.getDefaultContactGroup()
  group_b = 'copied'

  conts_b = Nucleus(nuc_file_b).getContacts()
  
  nuc.setContacts(group_b, conts_b, transposed=False)
  
  nuc.mergeContacts(group_a, group_b, 'merged', remove_frac)
  
  nuc.save()

  
if __name__ == '__main__':

  nuc_file_a = 'Nextera_02/Q5.nuc'
  nuc_file_b = 'Nextera_02/Q6.nuc'

  mergeContacts(nuc_file_a, nuc_file_b, remove_frac=0.5)
