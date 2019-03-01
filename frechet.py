import numpy as np
import pickle
import sys

from base import *


def sqDist(pt1, pt2):

    return (pt1.lat - pt2.lat)**2 + (pt1.lon - pt2.lon)**2


def distPtSegment(p, seg):
    (q1, q2) = seg
    x, y, x1, y1, x2, y2 = p.lon, p.lat, q1.lon, q1.lat, q2.lon, q2.lat

    if x1 == x2 and y1 == y2:  # Degenerate segment.
        return float(np.sqrt(sqDist(p, q1)))

    if x1 == x2:  # Vertical segment.
        if y1 <= y <= y2 or y2 <= y <= y1:
            return abs(x-x1)
        else:
            return min(float(np.sqrt(sqDist(p, q1))), float(np.sqrt(sqDist(p, q2))))

    elif y1 == y2:  # Horizontal segment.
        if x1 <= x <= x2 or x2 <= x <= x1:
            return abs(y-y1)
        else:
            return min(float(np.sqrt(sqDist(p, q1))), float(np.sqrt(sqDist(p, q2))))

    else:
        # Translate so that (x1,y1) is at the origin.
        x, y, x2, y2 = x-x1, y-y1, x2-x1, y2-y1
        m = y2/x2
        c = y + x/m
        # Projection of (x,y) on line passing through origin and (x2,y2).
        x3, y3 = c/(m+1/m), m*c/(m+1/m)
        if x2*x3 + y2*y3 >= 0:  # (x3,y3) between origin and (x2,y2).
            return float(np.sqrt((x3-x)**2 + (y3-y)**2))
        else:
            return min(float(np.sqrt(sqDist(p, q1))), float(np.sqrt(sqDist(p, q2))))


def frechetDec(trajA, trajB, delta):
    """ Decide if the discrete Frechet distance b/w trajectories is at most delta.

    Uses the classic dynamic programming algorithm to find the optimal correspondence.

    Args:
        trajA, trajB ([pt]): input trajectories, represented as list of pt objects.
        delta (float): guess on the discrete Frechet distance.

    Returns:
        True if the discrete Frechet distance b/w trajA and trajB is at most delta,
        False otherwise.
    """

    ptQueue = [(0, 0)]
    visited = set([])
    while len(ptQueue) > 0:
        current = ptQueue.pop(0)
        if current in visited:
            continue
        if current == (len(trajA) - 1, len(trajB) - 1):
            return True
        visited.add(current)
        i = current[0]
        j = current[1]
        # bounds check, add points which are within delta (using squared distance)
        if i + 1 < len(trajA):
            if sqDist(trajA[i + 1], trajB[j]) <= delta**2:
                ptQueue.append((i+1, j))
        if j + 1 < len(trajB):
            if sqDist(trajA[i], trajB[j + 1]) <= delta**2:
                ptQueue.append((i, j+1))
        if i + 1 < len(trajA) and j + 1 < len(trajB):
            if sqDist(trajA[i + 1], trajB[j + 1]) <= delta**2:
                ptQueue.append((i+1, j+1))
    return False
