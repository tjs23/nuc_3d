from cUtil.samread import readPairedSam

fileName = 'data/SiCUP/sample_1438/sample_1438_L001_AGGCAGAA_ACTGCATA_SiCUPPED.bam'

print(readPairedSam(fileName, 'rb'))
