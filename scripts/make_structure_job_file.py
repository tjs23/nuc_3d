import sys, time
import numpy as np
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from glob import glob
from NucApi import Nucleus, StrucCalcParams
from util.Structure import superimposeCoordPair

IN_FILES = ['P2-E8-PA_test_contacts_100k_s10.nuc', 'P2-E8-PA_test_contacts_100k_s1.nuc', 'P2-E8-PA_test_contacts_100k_s2.nuc', 'P2-E8-PA_test_contacts_100k_s3.nuc',
            'P2-E8-PA_test_contacts_100k_s4.nuc', 'P2-E8-PA_test_contacts_100k_s5.nuc', 'P2-E8-PA_test_contacts_100k_s6.nuc', 'P2-E8-PA_test_contacts_100k_s7.nuc',
            'P2-E8-PA_test_contacts_100k_s8.nuc', 'P2-E8-PA_test_contacts_100k_s9.nuc', 'P2-E8-PA_test_contacts_20k_s10.nuc', 'P2-E8-PA_test_contacts_20k_s1.nuc',
            'P2-E8-PA_test_contacts_20k_s2.nuc', 'P2-E8-PA_test_contacts_20k_s3.nuc', 'P2-E8-PA_test_contacts_20k_s4.nuc', 'P2-E8-PA_test_contacts_20k_s5.nuc',
            'P2-E8-PA_test_contacts_20k_s6.nuc', 'P2-E8-PA_test_contacts_20k_s7.nuc', 'P2-E8-PA_test_contacts_20k_s8.nuc', 'P2-E8-PA_test_contacts_20k_s9.nuc',
            'P2-E8-PA_test_contacts_30k_s10.nuc', 'P2-E8-PA_test_contacts_30k_s1.nuc', 'P2-E8-PA_test_contacts_30k_s2.nuc', 'P2-E8-PA_test_contacts_30k_s3.nuc',
            'P2-E8-PA_test_contacts_30k_s4.nuc', 'P2-E8-PA_test_contacts_30k_s5.nuc', 'P2-E8-PA_test_contacts_30k_s6.nuc', 'P2-E8-PA_test_contacts_30k_s7.nuc',
            'P2-E8-PA_test_contacts_30k_s8.nuc', 'P2-E8-PA_test_contacts_30k_s9.nuc', 'P2-E8-PA_test_contacts_40k_s10.nuc', 'P2-E8-PA_test_contacts_40k_s1.nuc',
            'P2-E8-PA_test_contacts_40k_s2.nuc', 'P2-E8-PA_test_contacts_40k_s3.nuc', 'P2-E8-PA_test_contacts_40k_s4.nuc', 'P2-E8-PA_test_contacts_40k_s5.nuc',
            'P2-E8-PA_test_contacts_40k_s6.nuc', 'P2-E8-PA_test_contacts_40k_s7.nuc', 'P2-E8-PA_test_contacts_40k_s8.nuc', 'P2-E8-PA_test_contacts_40k_s9.nuc',
            'P2-E8-PA_test_contacts_50k_s10.nuc', 'P2-E8-PA_test_contacts_50k_s1.nuc', 'P2-E8-PA_test_contacts_50k_s2.nuc', 'P2-E8-PA_test_contacts_50k_s3.nuc',
            'P2-E8-PA_test_contacts_50k_s4.nuc', 'P2-E8-PA_test_contacts_50k_s5.nuc', 'P2-E8-PA_test_contacts_50k_s6.nuc', 'P2-E8-PA_test_contacts_50k_s7.nuc',
            'P2-E8-PA_test_contacts_50k_s8.nuc', 'P2-E8-PA_test_contacts_50k_s9.nuc', 'P2-E8-PA_test_contacts_60k_s10.nuc', 'P2-E8-PA_test_contacts_60k_s1.nuc',
            'P2-E8-PA_test_contacts_60k_s2.nuc', 'P2-E8-PA_test_contacts_60k_s3.nuc', 'P2-E8-PA_test_contacts_60k_s4.nuc', 'P2-E8-PA_test_contacts_60k_s5.nuc',
            'P2-E8-PA_test_contacts_60k_s6.nuc', 'P2-E8-PA_test_contacts_60k_s7.nuc', 'P2-E8-PA_test_contacts_60k_s8.nuc', 'P2-E8-PA_test_contacts_60k_s9.nuc',
            'P2-E8-PA_test_contacts_70k_s10.nuc', 'P2-E8-PA_test_contacts_70k_s1.nuc', 'P2-E8-PA_test_contacts_70k_s2.nuc', 'P2-E8-PA_test_contacts_70k_s3.nuc',
            'P2-E8-PA_test_contacts_70k_s4.nuc', 'P2-E8-PA_test_contacts_70k_s5.nuc', 'P2-E8-PA_test_contacts_70k_s6.nuc', 'P2-E8-PA_test_contacts_70k_s7.nuc',
            'P2-E8-PA_test_contacts_70k_s8.nuc', 'P2-E8-PA_test_contacts_70k_s9.nuc', 'P2-E8-PA_test_contacts_80k_s10.nuc', 'P2-E8-PA_test_contacts_80k_s1.nuc',
            'P2-E8-PA_test_contacts_80k_s2.nuc', 'P2-E8-PA_test_contacts_80k_s3.nuc', 'P2-E8-PA_test_contacts_80k_s4.nuc', 'P2-E8-PA_test_contacts_80k_s5.nuc',
            'P2-E8-PA_test_contacts_80k_s6.nuc', 'P2-E8-PA_test_contacts_80k_s7.nuc', 'P2-E8-PA_test_contacts_80k_s8.nuc', 'P2-E8-PA_test_contacts_80k_s9.nuc',
            'P2-E8-PA_test_contacts_90k_s10.nuc', 'P2-E8-PA_test_contacts_90k_s1.nuc', 'P2-E8-PA_test_contacts_90k_s2.nuc', 'P2-E8-PA_test_contacts_90k_s3.nuc',
            'P2-E8-PA_test_contacts_90k_s4.nuc', 'P2-E8-PA_test_contacts_90k_s5.nuc', 'P2-E8-PA_test_contacts_90k_s6.nuc', 'P2-E8-PA_test_contacts_90k_s7.nuc',
            'P2-E8-PA_test_contacts_90k_s8.nuc', 'P2-E8-PA_test_contacts_90k_s9.nuc', 
            'P2-E8-PA_test_trans_3k_s10.nuc', 'P2-E8-PA_test_trans_3k_s1.nuc', 'P2-E8-PA_test_trans_3k_s2.nuc', 'P2-E8-PA_test_trans_3k_s3.nuc',
            'P2-E8-PA_test_trans_3k_s4.nuc', 'P2-E8-PA_test_trans_3k_s5.nuc', 'P2-E8-PA_test_trans_3k_s6.nuc', 'P2-E8-PA_test_trans_3k_s7.nuc',
            'P2-E8-PA_test_trans_3k_s8.nuc', 'P2-E8-PA_test_trans_3k_s9.nuc', 'P2-E8-PA_test_trans_4k_s10.nuc', 'P2-E8-PA_test_trans_4k_s1.nuc',
            'P2-E8-PA_test_trans_4k_s2.nuc', 'P2-E8-PA_test_trans_4k_s3.nuc', 'P2-E8-PA_test_trans_4k_s4.nuc', 'P2-E8-PA_test_trans_4k_s5.nuc',
            'P2-E8-PA_test_trans_4k_s6.nuc', 'P2-E8-PA_test_trans_4k_s7.nuc', 'P2-E8-PA_test_trans_4k_s8.nuc', 'P2-E8-PA_test_trans_4k_s9.nuc',
            'P2-E8-PA_test_trans_5k_s10.nuc', 'P2-E8-PA_test_trans_5k_s1.nuc', 'P2-E8-PA_test_trans_5k_s2.nuc', 'P2-E8-PA_test_trans_5k_s3.nuc',
            'P2-E8-PA_test_trans_5k_s4.nuc', 'P2-E8-PA_test_trans_5k_s5.nuc', 'P2-E8-PA_test_trans_5k_s6.nuc', 'P2-E8-PA_test_trans_5k_s7.nuc',
            'P2-E8-PA_test_trans_5k_s8.nuc', 'P2-E8-PA_test_trans_5k_s9.nuc', 'P2-E8-PA_test_trans_6k_s10.nuc', 'P2-E8-PA_test_trans_6k_s1.nuc',
            'P2-E8-PA_test_trans_6k_s2.nuc', 'P2-E8-PA_test_trans_6k_s3.nuc', 'P2-E8-PA_test_trans_6k_s4.nuc', 'P2-E8-PA_test_trans_6k_s5.nuc',
            'P2-E8-PA_test_trans_6k_s6.nuc', 'P2-E8-PA_test_trans_6k_s7.nuc', 'P2-E8-PA_test_trans_6k_s8.nuc', 'P2-E8-PA_test_trans_6k_s9.nuc',
            'P2-E8-PA_test_trans_7k_s10.nuc', 'P2-E8-PA_test_trans_7k_s1.nuc', 'P2-E8-PA_test_trans_7k_s2.nuc', 'P2-E8-PA_test_trans_7k_s3.nuc',
            'P2-E8-PA_test_trans_7k_s4.nuc', 'P2-E8-PA_test_trans_7k_s5.nuc', 'P2-E8-PA_test_trans_7k_s6.nuc', 'P2-E8-PA_test_trans_7k_s7.nuc',
            'P2-E8-PA_test_trans_7k_s8.nuc', 'P2-E8-PA_test_trans_7k_s9.nuc', 'P2-E8-PA_test_trans_8k_s10.nuc', 'P2-E8-PA_test_trans_8k_s1.nuc',
            'P2-E8-PA_test_trans_8k_s2.nuc', 'P2-E8-PA_test_trans_8k_s3.nuc', 'P2-E8-PA_test_trans_8k_s4.nuc', 'P2-E8-PA_test_trans_8k_s5.nuc',
            'P2-E8-PA_test_trans_8k_s6.nuc', 'P2-E8-PA_test_trans_8k_s7.nuc', 'P2-E8-PA_test_trans_8k_s8.nuc', 'P2-E8-PA_test_trans_8k_s9.nuc',
            'P2-E8-PA_test_trans_9k_s10.nuc', 'P2-E8-PA_test_trans_9k_s1.nuc', 'P2-E8-PA_test_trans_9k_s2.nuc', 'P2-E8-PA_test_trans_9k_s3.nuc',
            'P2-E8-PA_test_trans_9k_s4.nuc', 'P2-E8-PA_test_trans_9k_s5.nuc', 'P2-E8-PA_test_trans_9k_s6.nuc', 'P2-E8-PA_test_trans_9k_s7.nuc',
            'P2-E8-PA_test_trans_9k_s8.nuc', 'P2-E8-PA_test_trans_9k_s9.nuc',
            'P2-E8-PA_test_noise_0k_s10.nuc', 'P2-E8-PA_test_noise_0k_s1.nuc', 'P2-E8-PA_test_noise_0k_s2.nuc', 'P2-E8-PA_test_noise_0k_s3.nuc',
            'P2-E8-PA_test_noise_0k_s4.nuc', 'P2-E8-PA_test_noise_0k_s5.nuc', 'P2-E8-PA_test_noise_0k_s6.nuc', 'P2-E8-PA_test_noise_0k_s7.nuc',
            'P2-E8-PA_test_noise_0k_s8.nuc', 'P2-E8-PA_test_noise_0k_s9.nuc', 'P2-E8-PA_test_noise_10k_s10.nuc', 'P2-E8-PA_test_noise_10k_s1.nuc',
            'P2-E8-PA_test_noise_10k_s2.nuc', 'P2-E8-PA_test_noise_10k_s3.nuc', 'P2-E8-PA_test_noise_10k_s4.nuc', 'P2-E8-PA_test_noise_10k_s5.nuc',
            'P2-E8-PA_test_noise_10k_s6.nuc', 'P2-E8-PA_test_noise_10k_s7.nuc', 'P2-E8-PA_test_noise_10k_s8.nuc', 'P2-E8-PA_test_noise_10k_s9.nuc',
            'P2-E8-PA_test_noise_12k_s10.nuc', 'P2-E8-PA_test_noise_12k_s1.nuc', 'P2-E8-PA_test_noise_12k_s2.nuc', 'P2-E8-PA_test_noise_12k_s3.nuc',
            'P2-E8-PA_test_noise_12k_s4.nuc', 'P2-E8-PA_test_noise_12k_s5.nuc', 'P2-E8-PA_test_noise_12k_s6.nuc', 'P2-E8-PA_test_noise_12k_s7.nuc',
            'P2-E8-PA_test_noise_12k_s8.nuc', 'P2-E8-PA_test_noise_12k_s9.nuc', 'P2-E8-PA_test_noise_14k_s10.nuc', 'P2-E8-PA_test_noise_14k_s1.nuc',
            'P2-E8-PA_test_noise_14k_s2.nuc', 'P2-E8-PA_test_noise_14k_s3.nuc', 'P2-E8-PA_test_noise_14k_s4.nuc', 'P2-E8-PA_test_noise_14k_s5.nuc',
            'P2-E8-PA_test_noise_14k_s6.nuc', 'P2-E8-PA_test_noise_14k_s7.nuc', 'P2-E8-PA_test_noise_14k_s8.nuc', 'P2-E8-PA_test_noise_14k_s9.nuc',
            'P2-E8-PA_test_noise_2k_s10.nuc', 'P2-E8-PA_test_noise_2k_s1.nuc', 'P2-E8-PA_test_noise_2k_s2.nuc', 'P2-E8-PA_test_noise_2k_s3.nuc',
            'P2-E8-PA_test_noise_2k_s4.nuc', 'P2-E8-PA_test_noise_2k_s5.nuc', 'P2-E8-PA_test_noise_2k_s6.nuc', 'P2-E8-PA_test_noise_2k_s7.nuc',
            'P2-E8-PA_test_noise_2k_s8.nuc', 'P2-E8-PA_test_noise_2k_s9.nuc', 'P2-E8-PA_test_noise_4k_s10.nuc', 'P2-E8-PA_test_noise_4k_s1.nuc',
            'P2-E8-PA_test_noise_4k_s2.nuc', 'P2-E8-PA_test_noise_4k_s3.nuc', 'P2-E8-PA_test_noise_4k_s4.nuc', 'P2-E8-PA_test_noise_4k_s5.nuc',
            'P2-E8-PA_test_noise_4k_s6.nuc', 'P2-E8-PA_test_noise_4k_s7.nuc', 'P2-E8-PA_test_noise_4k_s8.nuc', 'P2-E8-PA_test_noise_4k_s9.nuc',
            'P2-E8-PA_test_noise_6k_s10.nuc', 'P2-E8-PA_test_noise_6k_s1.nuc', 'P2-E8-PA_test_noise_6k_s2.nuc', 'P2-E8-PA_test_noise_6k_s3.nuc',
            'P2-E8-PA_test_noise_6k_s4.nuc', 'P2-E8-PA_test_noise_6k_s5.nuc', 'P2-E8-PA_test_noise_6k_s6.nuc', 'P2-E8-PA_test_noise_6k_s7.nuc',
            'P2-E8-PA_test_noise_6k_s8.nuc', 'P2-E8-PA_test_noise_6k_s9.nuc', 'P2-E8-PA_test_noise_8k_s10.nuc', 'P2-E8-PA_test_noise_8k_s1.nuc',
            'P2-E8-PA_test_noise_8k_s2.nuc', 'P2-E8-PA_test_noise_8k_s3.nuc', 'P2-E8-PA_test_noise_8k_s4.nuc', 'P2-E8-PA_test_noise_8k_s5.nuc',
            'P2-E8-PA_test_noise_8k_s6.nuc', 'P2-E8-PA_test_noise_8k_s7.nuc', 'P2-E8-PA_test_noise_8k_s8.nuc', 'P2-E8-PA_test_noise_8k_s9.nuc',
            
            ]
            
IN_FILES = ['P2-E8-PA_test_noise_0k_s10.nuc', 'P2-E8-PA_test_noise_0k_s1.nuc', 'P2-E8-PA_test_noise_0k_s2.nuc', 'P2-E8-PA_test_noise_0k_s3.nuc',
            'P2-E8-PA_test_noise_0k_s4.nuc', 'P2-E8-PA_test_noise_0k_s5.nuc', 'P2-E8-PA_test_noise_0k_s6.nuc', 'P2-E8-PA_test_noise_0k_s7.nuc',
            'P2-E8-PA_test_noise_0k_s8.nuc', 'P2-E8-PA_test_noise_0k_s9.nuc', 'P2-E8-PA_test_noise_10k_s10.nuc', 'P2-E8-PA_test_noise_10k_s1.nuc',
            'P2-E8-PA_test_noise_10k_s2.nuc', 'P2-E8-PA_test_noise_10k_s3.nuc', 'P2-E8-PA_test_noise_10k_s4.nuc', 'P2-E8-PA_test_noise_10k_s5.nuc',
            'P2-E8-PA_test_noise_10k_s6.nuc', 'P2-E8-PA_test_noise_10k_s7.nuc', 'P2-E8-PA_test_noise_10k_s8.nuc', 'P2-E8-PA_test_noise_10k_s9.nuc',
            'P2-E8-PA_test_noise_12k_s10.nuc', 'P2-E8-PA_test_noise_12k_s1.nuc', 'P2-E8-PA_test_noise_12k_s2.nuc', 'P2-E8-PA_test_noise_12k_s3.nuc',
            'P2-E8-PA_test_noise_12k_s4.nuc', 'P2-E8-PA_test_noise_12k_s5.nuc', 'P2-E8-PA_test_noise_12k_s6.nuc', 'P2-E8-PA_test_noise_12k_s7.nuc',
            'P2-E8-PA_test_noise_12k_s8.nuc', 'P2-E8-PA_test_noise_12k_s9.nuc', 'P2-E8-PA_test_noise_14k_s10.nuc', 'P2-E8-PA_test_noise_14k_s1.nuc',
            'P2-E8-PA_test_noise_14k_s2.nuc', 'P2-E8-PA_test_noise_14k_s3.nuc', 'P2-E8-PA_test_noise_14k_s4.nuc', 'P2-E8-PA_test_noise_14k_s5.nuc',
            'P2-E8-PA_test_noise_14k_s6.nuc', 'P2-E8-PA_test_noise_14k_s7.nuc', 'P2-E8-PA_test_noise_14k_s8.nuc', 'P2-E8-PA_test_noise_14k_s9.nuc',
            'P2-E8-PA_test_noise_2k_s10.nuc', 'P2-E8-PA_test_noise_2k_s1.nuc', 'P2-E8-PA_test_noise_2k_s2.nuc', 'P2-E8-PA_test_noise_2k_s3.nuc',
            'P2-E8-PA_test_noise_2k_s4.nuc', 'P2-E8-PA_test_noise_2k_s5.nuc', 'P2-E8-PA_test_noise_2k_s6.nuc', 'P2-E8-PA_test_noise_2k_s7.nuc',
            'P2-E8-PA_test_noise_2k_s8.nuc', 'P2-E8-PA_test_noise_2k_s9.nuc', 'P2-E8-PA_test_noise_4k_s10.nuc', 'P2-E8-PA_test_noise_4k_s1.nuc',
            'P2-E8-PA_test_noise_4k_s2.nuc', 'P2-E8-PA_test_noise_4k_s3.nuc', 'P2-E8-PA_test_noise_4k_s4.nuc', 'P2-E8-PA_test_noise_4k_s5.nuc',
            'P2-E8-PA_test_noise_4k_s6.nuc', 'P2-E8-PA_test_noise_4k_s7.nuc', 'P2-E8-PA_test_noise_4k_s8.nuc', 'P2-E8-PA_test_noise_4k_s9.nuc',
            'P2-E8-PA_test_noise_6k_s10.nuc', 'P2-E8-PA_test_noise_6k_s1.nuc', 'P2-E8-PA_test_noise_6k_s2.nuc', 'P2-E8-PA_test_noise_6k_s3.nuc',
            'P2-E8-PA_test_noise_6k_s4.nuc', 'P2-E8-PA_test_noise_6k_s5.nuc', 'P2-E8-PA_test_noise_6k_s6.nuc', 'P2-E8-PA_test_noise_6k_s7.nuc',
            'P2-E8-PA_test_noise_6k_s8.nuc', 'P2-E8-PA_test_noise_6k_s9.nuc', 'P2-E8-PA_test_noise_8k_s10.nuc', 'P2-E8-PA_test_noise_8k_s1.nuc',
            'P2-E8-PA_test_noise_8k_s2.nuc', 'P2-E8-PA_test_noise_8k_s3.nuc', 'P2-E8-PA_test_noise_8k_s4.nuc', 'P2-E8-PA_test_noise_8k_s5.nuc',
            'P2-E8-PA_test_noise_8k_s6.nuc', 'P2-E8-PA_test_noise_8k_s7.nuc', 'P2-E8-PA_test_noise_8k_s8.nuc', 'P2-E8-PA_test_noise_8k_s9.nuc',
            
            ]


#IN_FILES = ['P2-E8-PA_test_contacts_100k_s1.nuc', 'P2-E8-PA_test_contacts_100k_s2.nuc']

def make_job_file_400(nuc_file_path, num_models=10):

  chromosomes = [str(x) for x in range(1,20)] + ['X']
  
  seed = int(time.time())
  
  print('Setup calc params')

  calcParams = StrucCalcParams(distPowerLaw=-0.33, distLower=0.8, distUpper=1.2,
                               bboneLower=0.1, bboneUpper=1.1, seqScale=10000,
                               randSeed=seed, randStart=False,  randWalk=False, randRad=1000.0)

  resolution = int(8e6)
  
  calcParams.addAnnealStage(domainDict={}, bboneSep=[resolution, resolution],
                            tempStart=5000.0, tempEnd=10.0,
                            tempSteps=500, dynSteps=20, timeStep=0.001, useTrans=True)

  resolution = int(4e6)
  
  calcParams.addAnnealStage(domainDict={}, bboneSep=[resolution, resolution],
                            tempStart=5000.0, tempEnd=10.0,
                            tempSteps=500, dynSteps=50, timeStep=0.001, useTrans=True)

  resolution = int(2e6)
  
  calcParams.addAnnealStage(domainDict={}, bboneSep=[resolution, resolution],
                            tempStart=5000.0, tempEnd=10.0,
                            tempSteps=500, dynSteps=100, timeStep=0.001, useTrans=True)
  
  resolution = int(4e5)
  
  calcParams.addAnnealStage(domainDict={}, bboneSep=[resolution, resolution],
                            tempStart=5000.0, tempEnd=10.0,
                            tempSteps=1000, dynSteps=100, timeStep=0.001, useTrans=True)
  
  nuc = Nucleus(nuc_file_path)
  
  group_name = nuc.getContactGroups()[0][0]
  
  nuc.removeIsolatedContacts(group_name, threshold=int(2e6))
  
  file_root, file_ext = splitext(nuc_file_path)
  
  new_nuc_file_path = '%s_%dx_%dkb.nuc' % (file_root, num_models, int(resolution/1000))
  
  nuc.saveAs(new_nuc_file_path)
  
  save_file_path = new_nuc_file_path[:-4] + '_job.npz'
  
  nuc.setupStrucCalcJobFile(save_file_path, group_name, chromosomes=chromosomes,
                            numModels=num_models, calcParams=calcParams, structure='0')
  

def make_job_file_100(nuc_file_path, num_models=10, viol_threshold=6.0, resolution=int(1e5)):
  
  dir_name, file_name = split(nuc_file_path)
  chromosomes = [str(x) for x in range(1,20)] + ['X']
  
  seed = int(time.time())
  
  print('Setup calc params')

  calcParams = StrucCalcParams(distPowerLaw=-0.33, distLower=0.8, distUpper=1.2,
                               bboneLower=0.1, bboneUpper=1.1, seqScale=10000,
                               randSeed=seed, randStart=False,  randWalk=False, randRad=1000.0)
  
  calcParams.addAnnealStage(domainDict={}, bboneSep=[resolution, resolution],
                            tempStart=5000.0, tempEnd=10.0,
                            tempSteps=500, dynSteps=100, timeStep=0.001, useTrans=True)
  
  nuc = Nucleus(nuc_file_path)
  in_num_models = nuc.getNumModels()
  
  if num_models < in_num_models:
  
    names = file_name.split('_')
    names[-2] = '%dx' % num_models
    new_nuc_file_path = join(dir_name, '_'.join(names)) 
    nuc.saveAs(new_nuc_file_path)
    
    coords, chr_ranges = nuc.getAllCoords()
    rmsd_weight_scale = 10.0
    
    model_rmsds = [[] for m in range(in_num_models)]
      
    for i in range(in_num_models-1):
      for j in range(i+1, in_num_models):
        t_coords, rotn, rmsd, atom_rmsds = superimposeCoordPair(coords[i], coords[j], rmsd_weight_scale)
        model_rmsds[i].append(rmsd)
        model_rmsds[j].append(rmsd)
     
    model_rmsds = np.array([np.median(x) for x in model_rmsds])
    
    idx = np.argsort(model_rmsds)
    print ['%.2f' % x for x in model_rmsds[idx]]
    
    idx = idx[:num_models]
   
    
    nuc.setAllCoords(coords[idx], sorted(chromosomes), structure='0')      
    nuc.modelAlign(chromosomes=chromosomes, structure='0')
    nuc.save()
    
      
  group_name = nuc.getContactGroups()[0][0]
  
  if viol_threshold is not None:
    nuc.removeViolatedContacts(group_name, threshold=viol_threshold)
   

  
  names = file_name.split('_')[:-2]
  names += ['%dx' % num_models, '%dkb.nuc' % int(resolution/1000)]
 
  new_nuc_file_path = join(dir_name, '_'.join(names)) 
  
  if new_nuc_file_path == nuc_file_path:
    nuc.save()
  else:  
    nuc.saveAs(new_nuc_file_path)
  
  save_file_path = new_nuc_file_path[:-4] + '_job.npz'
  
  nuc.setupStrucCalcJobFile(save_file_path, group_name, chromosomes=chromosomes,
                            numModels=num_models, calcParams=calcParams, structure='0')

"""
python scripts/make_structure_job_file.py 0 &
python scripts/make_structure_job_file.py 1 &
python scripts/make_structure_job_file.py 2 &
python scripts/make_structure_job_file.py 3 &
python scripts/make_structure_job_file.py 4 &
python scripts/make_structure_job_file.py 5 &
python scripts/make_structure_job_file.py 6 &
python scripts/make_structure_job_file.py 7 &
python scripts/make_structure_job_file.py 8 &
python scripts/make_structure_job_file.py 9 &
python scripts/make_structure_job_file.py 10 &
python scripts/make_structure_job_file.py 11 &
"""

IN_DIR = 'P2-E8-PA_test/'

IN_FILES2 = ['P2E8.nuc', 'P2F8.nuc', 'P2I5.nuc', 'P2J8.nuc', 'P30E18.nuc',
             'P30E4.nuc', 'P30E8.nuc', 'P30F18.nuc', 'Q5.nuc', 'Q6.nuc']

IN_DIR2 = 'paper_structs/'

if __name__ == '__main__':
  
  import sys
  
  if len(sys.argv) > 1:
    job_id = int(sys.argv[1])
  else:
    job_id = None
  
  n_cpu = 12
  num_models = 10
  #num_models = 20
  
  #for i, file_name in enumerate(IN_FILES):
  #  if (job_id is None) or (i % n_cpu == job_id):
  #    print i, file_name
  #    make_job_file_400(IN_DIR + file_name, num_models)
  #print "Total files:", len(IN_FILES)
  
  #in_files = sorted(glob('/home/tjs23/nucleus/paper_structs/*_20x_400kb.nuc'))
  #in_files = sorted(glob('/home/tjs23/nucleus/paper_structs/*_10x_200kb.nuc'))
  in_files = sorted(glob('/home/tjs23/nucleus/P2-E8-PA_test/P2-E8-PA_test_noise_*_10x_200kb.nuc'))
  #in_files = sorted(glob('/home/tjs23/nucleus/P2-E8-PA_test/P2-E8-PA_test_*_10x_200kb.nuc'))
  #in_files = sorted(glob('/home/tjs23/nucleus/paper_structs/*_ambig_10x_100kb.nuc'))
  for i, file_name in enumerate(in_files):
    if (job_id is None) or (i % n_cpu == job_id):
      print i, file_name
      #make_job_file_100(file_name, num_models, 6.0, 200000)
      make_job_file_100(file_name, num_models, 5.0, 100000)
  #    #make_job_file_100(file_name, num_models, None, 100000)
  
  # Server structure calculation process
  # - Run make_job_file_400 on input nuc files makes 400kb.nuc file
  #   + Transfer .npz files to comp cluster
  #   + Repeat models share same input
  #   + Run njobs-nmodels script on comp cluster - last digit selects model to calculate
  #   + When complete copy back _job_coords_*.npy files
  # - Run combine_structure_job_coords.py  400kb.nuc _job_coords_*.npy
  # - Run make_job_file_100 on 400kb.nuc file makes 100kb.nuc file
  #   + Transfer .npz files to comp cluster
  #   + Run njobs-nmodels script on comp cluster
  #   + When complete compy back _job_coords_*.npy files
  # - Run combine_structure_job_coords.py 100kb.nuc _job_coords_*.npy
  # - Run make_job_file_100 on 100kb.nuc with viol_threshold=5.0 file makes new 100kb.nuc file
  #   + Transfer .npz files to comp cluster
  #   + Run njobs-nmodels script on comp cluster
  #   + When complete compy back _job_coords_*.npy files
  # - Run combine_structure_job_coords.py 100kb.nuc _job_coords_*.npy
  
  
  
  
  
  
  
  
  
  
