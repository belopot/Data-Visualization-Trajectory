import sys
import _pickle as pickle

from base import *
from heapdict import heapdict


def preprocess(trajs, distPairs):
    """ function responsible for calculating data structures required in greedy algorithm
        Input: trajs ({int : traj})
               distPairs ({(pathlet, int) : [(subTraj, float)])
        Output:
             strajCov:  dict storing point of all subtrajs in distpairs.
            ptStraj :  for each point it stores, the set of subtrajs,distpair containing it
            strajPth  : It stores the list of pathlets for each pathlet in distPairs.
            trajCov : it stores each trajectory points.

"""

    strajCov, ptStraj, strajPth = {}, {}, {}
    for key, value in list(distPairs.items()):
        pth, trID = key[0], key[1]
        for i in range(len(value)):
            straj, dist = value[i][0], value[i][1]
            #print(("straj", straj))
            strajCov[straj] = straj.bounds[1] - straj.bounds[0] + 1
            if straj in strajPth:
                strajPth[straj].append(pth)
            else:
                strajPth[straj] = [pth]
            for j in range(int(straj.bounds[0]), int(straj.bounds[1] + 1)):
                p = tuple(trajs.get_group(trID).iloc[j])
                #print(("point", p))
                if p in ptStraj:
                    ptStraj[p].add(straj)
                else:
                    # print("straj:",set([straj])
                    ptStraj[p] = set([straj])

    trajCov = {}
    for trID, tra in trajs:
        trajCov[trID] = tra.shape[0]
        for i in range(tra.shape[0]):
            p = tuple(tra.iloc[i])
            if p not in ptStraj:  # this means p can't be assigned to any pathlet
                ptStraj[p] = set()

    return (strajCov, ptStraj, strajPth, trajCov)


def processPoint(p, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue):
    """ Responsible for processing a point in each iteration of greedy algorithm.
        Input:
            p : Point that is to be processed.
            ptStraj : for each point it stores, the set of subtrajs,distpair containing it
            strajCov: dict storing point of all subtrajs in distpairs.
            strajPth: It stores the list of pathlets for each pathlet in distPairs.
            trajs : dict mapping IDs to trajectories.
            trajCov  : it stores each trajectory points.
            distPairs : It contaoins mapping of
                    pathlet-trajID pair to a list of subtraj-float pairs where
                    subtraj belongs to its respective traj, and Frechet distance is used to compute float values.
            numUnprocessedPts (int) :Number of unprocessed points
            queue : priority queue

        Output:
            Set having pathlet-trajID pairs to which point
            can be assigned.
    """

    retVal = set()
    for straj in ptStraj[p]:
        trID = straj.trajID
        strajCov[straj] = strajCov[straj] - 1
        for i in range(len(strajPth[straj])):
            pth = strajPth[straj][i]
            retVal.add((pth, trID))
    ptStraj[p] = None  # this is also marking that p is processed.
    numUnprocessedPts[0] = numUnprocessedPts[0] - 1

    # Check if we also have singleton sets.
    if queue is not None:
        id = p[0]
        trajCov[id] -= 1
        # Change priority of trajID to zero if its coverage is 0.
        if trajCov[id] == 0:
            queue[id] = 0
    return retVal


def processSubtraj(straj, strajCov, trajs, trajCov, ptStraj, strajPth, distPairs, numUnprocessedPts, queue):
    """ responsible for processing the points of a subtrajjectory picked in an interation of greedy algorithm.
    """
    trID = straj.trajID
    tr = trajs.get_group(trID)
    retVal = set()
    for i in range(straj.bounds[0], straj.bounds[1] + 1):
        p = tuple(tr.iloc[i])
        # Check if point p is unprocessed.
        if ptStraj[p] is not None:
            retVal = retVal.union(processPoint(
                p, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue))
    if strajCov[straj] != 0:
        print(("Error!! Coverage should have been 0, instead of %d" %
               strajCov[straj]))

    return retVal


def processTraj(trID, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue, unassignedPts):
    """ Responsible for processing the unprocessed points of a trajectory

        This function is called when an interation of the greedy algorithm decides to leave the unprocessed
        points of the trajectory unassigned.

        Args:
            trID (int): ID of the trajectory whose points are input to the function.
            ptStraj
            strajCov
            strajPth
            trajs
            trajCov
            distPairs
            numUnprocessedPts (int) :
            queue : priority queue

        Returns:
            Set having pathlet-trajID pairs to which point can be processed.
    """
    points = trajs[trID].pts
    retVal = set()
    for i in range(len(points)):
        p = points[i]
        # Check if p is unprocessed.
        if ptStraj[p] is not None:
            unassignedPts.append(p)
            retVal = retVal.union(
                processPoint(p, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue))

    if trajCov[trID] != 0:
        print("Error!! Coverage should have been zero.")

    return retVal


def findCovCostRatio(trajStrajDist, c1, c3, strajCov):
    """ Calculates coverage cost ratio for pathlet.

        Args:
            trajStrajDist: dict having subtrajectory for
                    a trajectory ID with distance from the pathlet.
            c1, c3 (float): parameters of the greedy algorithm.
            strajCov : dict storing points of subtrajs.

        Returns:
           coverage-cost ratio of the pathlet.
    """
    curCov, curCost = 0, c1
    for k, v in list(trajStrajDist.items()):
        straj, dist = v[0], v[1]
        if straj is None:
            continue
        curCov += strajCov[straj]
        curCost += 1.0 * c3 * dist
    if curCov == 0:
        return 0
    else:
        return (1.0 * curCov) / (1.0 * curCost)


def optStrajAdvHelp(strajDists, strajCov, r, c3):
    """ Computes subtrajectory of a trajectory with maximum coverage cost ratio
        Args:
            strajDists: list having subtraj-distance pairs of a trajectory.
            strajCov: dict storing points of subtrajs.
            r (float): guess coverage-cost ratio.
            c3 : parameter of greedy algorithm.

        Returns:
            (subtraj, float) or (None,None).
    """
    temp, stra, dista = 0, None, None
    for i in range(len(strajDists)):
        straj, dist = strajDists[i][0], strajDists[i][1]
        cov = strajCov[straj]
        if temp < cov - c3 * r * dist:
            temp = cov - c3 * r * dist
            stra, dista = straj, dist

    return (stra, dista)


def computeStrajsAdvanced(pth, distPairs, pthOptStrajs, strajCov, c1, c3, m, affectedTrajs):
    """ Find subtrajectory for a pathlet with optimal coverage-cost ratio.

        Args:
            pth (pathlet) : pathlet of concern.
            distPairs
            pthOptStrajs: dict that stores optimal subtrajectory of each trajectory with distance for a pathlet.
            strajCov : dict storing point of subtrajs.
            c1, c3 (float) : parameters of the greedy algorithm.
            m (int) : total no. of points.
            affectedTrajs [int] : list of trajIDs with newly covered points.

        Returns:
            Nothing, but updates pthOptStrajs.
    """
    optStrajs = pthOptStrajs[pth]
    ret, temp = {}, {}
    if c3 == 0:
        for trID in affectedTrajs:
            strajDists = distPairs[(pth, trID)]
            (straj, dist) = optStrajAdvHelp(strajDists, strajCov, 1, c3)
            ret[trID] = (straj, dist)
    else:
        # r is a guess on the coverage-cost ratio.
        # It is used to find the exit condition of iteration on the basis of value r.
        summation = 0
        rmin = 1.0 / (c1 + c3)
        r, rmax = rmin, m * rmin
        while r <= rmax:
            for trID in affectedTrajs:
                strajDists = distPairs[(pth, trID)]
                (straj, dist) = optStrajAdvHelp(strajDists, strajCov, r, c3)
                temp[trID] = (straj, dist)
                if straj is not None:
                    summation += (strajCov[straj] - c3 * r * dist)

            if summation < c1 * r:
                break

            else:
                ret = temp
                temp = {}
                summation = 0
            r *= 2
    pthOptStrajs[pth] = ret


def Greedy(trajs, distPairs, strajCov, ptStraj, strajPth, trajCov, c1, c2, c3):
    """ Run the greedy algorithm for pathlet cover.
    The algorithm do one of the two things either leaves a point unassigned or pics pathlet and subtrajectories set assigned to that pathlet on the basis of maximum coverage-cost ratio.
The points that are picked up in each step are said to
    be "processed".

    Args:
        trajs : dict having mapping of ID to traj objects.
        distPairs : dict containing mapping of
                    pathlet-trajID pair to a list of subtraj-float pairs
        strajCov  : It stores points in all subtrajs
                    in distPairs.
        ptStraj  : It stores set of subtrajs in distpairs for each point.
        strajPth :it stores list of pathlets associated with each subtrajs in distpairs
        trajCov  : dict storing the points in each trajectory.
        c1,c2,c3 (float): parameters of the greedy algorithm.

    Returns:
        Pathlet assignments and unassigned points as determined by the greedy algorithm,
    """
    pthOptCovCost, pthOptStrajs = {}, {}

    for key, value in list(distPairs.items()):
        pth, trID, strajDists = key[0], key[1], value
        if pth not in pthOptStrajs:
            pthOptStrajs[pth] = {}
        pthOptStrajs[pth][trID] = (None, None)

    # Compute pthOptStrajs.
    for key, value in list(distPairs.items()):
        pth = key[0]
        affectedTrajs = list(pthOptStrajs[pth].keys())
        computeStrajsAdvanced(pth, distPairs, pthOptStrajs,
                              strajCov, c1, c3, len(ptStraj), affectedTrajs)
        print((pthOptStrajs[pth]))

    # Compute pthOptCovCost using pthOptStrajs.
    for key, value in list(pthOptStrajs.items()):
        pth = key
        pthOptCovCost[pth] = findCovCostRatio(value, c1, c3, strajCov)
        if pthOptCovCost[pth] == 0:
            print("Error")

    # Initialize a max priority queue of pathlets order by coverage cost ratio.
    queue1 = heapdict()
    for pth, ccratio in list(pthOptCovCost.items()):
        # Need to negate, since heapdict is a min-heap.
        queue1[pth] = -1.0 * ccratio

    # Initialize a priority queue of trajs, with coverage to cost ratio of a singleton set, i.e., |T|/c2.
    queue2 = heapdict()
    for trID, cov in list(trajCov.items()):
        queue2[trID] = -(1.0 * cov) / (1.0 * c2)
    pthAssignments = {}
    pthStats = []

    # Unassigned points
    unassignedPts = []

    count = 0
    numUnprocessedPts = [len(ptStraj)]
    print(("number of Unprocessed Points", numUnprocessedPts))
    while numUnprocessedPts[0] > 0:
        print(("number of poits is %d" % numUnprocessedPts[0]))
        x1, x2 = queue1.peekitem(), queue2.peekitem()
        affectedPths = set()

        if x1[1] <= x2[1]:
            x = queue1.popitem()
            pth = x[0]
            fracThickness = 0
            for trID, pair in list(pthOptStrajs[pth].items()):
                straj = pair[0]
                if straj is None:
                    continue
                if pth in pthAssignments:
                    pthAssignments[pth].append(straj)
                else:
                    pthAssignments[pth] = [straj]

                fracThickness += 1.0 * \
                    strajCov[straj] / trajs.get_group(straj.trajID).shape[0]

                affectedPths = affectedPths.union(processSubtraj(
                    straj, strajCov, trajs, trajCov, ptStraj, strajPth, distPairs, numUnprocessedPts, queue2))

            pthStats.append((pth, pthOptCovCost[pth], count, fracThickness))

        else:
            x = queue2.popitem()
            trID = x[0]
            # Process the unassigned point, and also return pathlets whose optimal coverage-cost ratio changes

            affectedPths = affectedPths.union(
                processTraj(trID, ptStraj, strajCov, strajPth, trajs, trajCov, distPairs, numUnprocessedPts, queue2,
                            unassignedPts))

        # Update coverage-cost ratio of affected pathlets.
        affectedPathlets = {path for (path, traID) in affectedPths}
        affectedTrajs = {traID for (path, traID) in affectedPths}

        for path in affectedPathlets:
            computeStrajsAdvanced(
                path, distPairs, pthOptStrajs, strajCov, c1, c3, len(ptStraj), affectedTrajs)
            pthOptCovCost[path] = findCovCostRatio(
                pthOptStrajs[path], c1, c3, strajCov)
            queue1[path] = -1.0 * pthOptCovCost[path]

        count += 1
    print(("Bundles and subtrajectories in the bundle", pthAssignments))
    #print(("pthStats", pthStats))
    print(("Unassigned Points", unassignedPts))
    return (pthAssignments, pthStats, unassignedPts)
