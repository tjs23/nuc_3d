from NucApi import Nucleus, EXTERNAL, DERIVED
import numpy as np
from cUtil.overlapHelper import getPercOL 

def overlap(trackA, trackB):
    """
    Gets the proportion of A that overlaps with B.
    """
    percOL = {} 
    for k in trackA:
        if k in trackB:
            trackA = np.sort(trackA[k][0].astype(np.int32))
            trackB = np.sort(trackB[k][0].astype(np.int32))
            percOL[k] = getPercOL(trackA, trackB)
        else:
            percOL[k] = 0.0 

    totalOL = np.average(np.array([percOL[x] for x in percOL]))

    return(totalOL)


if __name__ == "__main__":
    import itertools as it
    nuc = Nucleus("/Users/Liam/ProjCam/projects/nucleus/data/examples/NXT-33_50k-chipSeq.nuc")

    tracks = ["chipSeq_broad_K36_3_HAP",
              "chipSeq_broad_Nanog_HAP",
              "chipSeq_broad_K27_3_HAP",
              "chipSeq_broad_K4_3_HAP",
              "chipSeq_narrow_K36_3_HAP",
              "chipSeq_narrow_Nanog_HAP",
              "chipSeq_narrow_K27_3_HAP",
              "chipSeq_narrow_K4_3_HAP",
              "LaminB",
              "H3K4me3"]
    for a, b in it.combinations(tracks, 2):
        print("{} vs {}".format(a, b))
        trackA = nuc.getDataTrack(source=EXTERNAL, code=a)
        trackB = nuc.getDataTrack(source=EXTERNAL, code=b)

        print("{}%".format(100 * overlap(trackA, trackB)))

            


        
