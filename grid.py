import numpy as np
import sys

from base import *
from distanceUtils import *


class trajGrid:
    """ A class to represent a grid to index trajectory points.

    Atrributes:
        data : dictionary containing mapping of a grid cell
                to a dictionary containing traj points.
        startLat (float): smallest y-value of indexed points.
        startLon (float): smallest x-value of indexed points.
        delta (float): size of grid cell.
        numX (int): number of cells in the x-direction.
        numY (int): number of cells in the y-direction.
    """

    def __init__(self, data, startLat, startLon, delta, numX, numY):
        """ Simple initialization of class. """
        self.data = data
        self.startLat = startLat
        self.startLon = startLon
        self.delta = delta
        self.xCells = numX
        self.yCells = numY


def gridData(trajs, simpTrajs, delta):
    """ Create an index of trajectory points.

    Args:
        trajs: dictionary of traj objects.
        simpTajs: list of trajectories whose points
                need to be indexed.
        delta (float): size of grid cell.

    Returns:
        trajGrid object containing the points of trajectories.
    """

    simpTraj = simpTrajs[0]
    trID = simpTraj.trajID
    p = trajs.get_group(trID).iloc[simpTraj.indices[0]]
    xMin, yMin = p.lon, p.lat
    xMax, yMax = xMin, yMin

    for simpTraj in simpTrajs:
        trID = simpTraj.trajID
        for index in simpTraj.indices:
            p = trajs.get_group(trID).iloc[index]
            xMin, xMax = min(xMin, p.lon), max(xMax, p.lon)
            yMin, yMax = min(yMin, p.lat), max(yMax, p.lat)
    xCells = int(np.ceil((xMax-xMin)/delta))
    yCells = int(np.ceil((yMax-yMin)/delta))
    grid = {}
    for simpTraj in simpTrajs:
        trID = simpTraj.trajID
        for index in simpTraj.indices:
            p = trajs.get_group(trID).iloc[index]
            xCell = int(np.floor((p.lon - xMin)/delta))
            yCell = int(np.floor((p.lat - yMin)/delta))
            if xCell not in grid and yCell not in grid:
                grid[(xCell, yCell)] = {}
            # Add trajID to cell if not there.
            if trID not in grid[(xCell, yCell)]:
                grid[(xCell, yCell)][trID] = []
            grid[(xCell, yCell)][trID].append(index)

    return trajGrid(grid, yMin, xMin, delta, xCells, yCells)
