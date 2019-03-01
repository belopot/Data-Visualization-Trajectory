
import _pickle as pickle

from base import *
from distance import *
from greedy import *
from base import *
import datetime

import pandas as pd

def readTrajsFromTxtFile(fName):
    df = pd.read_csv(fName,
                     sep=',',
                     header=None,
                     names=['trajID', 't', 'lat', 'lon'],
                     comment='#')
    df['t'] = pd.to_datetime(df.t)
    df = df.sort_values(by='t')
    df = df.groupby(by='trajID')
    return df

if __name__ == "__main__":

    print("Loading trajectories ...")
    trajs = readTrajsFromTxtFile("data/small.txt")
    # print(("trajectories are",trajs))
    rmin, rmax = 0.5, 1

    print("Computing Frechet distances ...")
    distPairs1 = process(trajs, rmin, rmax)
    # print(("Output of Frechet distance", distPairs1))
    distPairs2 = {}
    for k, v in distPairs1.items():
        pth, trID, dist, straj = k[0], k[1].trajID, v, k[1]
        if (pth, trID) in distPairs2:
            distPairs2[(pth, trID)].append((straj, dist))
        else:
            distPairs2[(pth, trID)] = [(straj, dist)]
    print("Computing prerequisite data structures for greedy algorithm ...")
    (strajCov, ptStraj, strajPth, trajCov) = preprocess(trajs, distPairs2)
    #print("Pre-requiste data structure are:")
    # print(("strajCov",strajCov))
    # print(("ptStraj",ptStraj))
    # print(("strajPth",strajPth))
    # print(("trajCov",trajCov))

    c1, c2, c3 = 1, 1, 1

    print("Running greedy algorithm ...")
    retVal = Greedy(trajs, distPairs2, strajCov, ptStraj,
                    strajPth, trajCov, c1, c2, c3)
    #print(("Bundle assigned and unassigned points are",retVal))
