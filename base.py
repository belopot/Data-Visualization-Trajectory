class pathlet(object):
    """ A class representing a pathlet (essentially a subtrajectory).

        Similar to traj, however avoids storing the points explicitly by referring to the
        trajectory to which the points belong, along with the starting and ending indices
        of the points in the list of points of the trajectory, sorted by timestamp.

        Attributes:
            trajID (int): ID of the trajectory to which points of the pathlet belong.
            bounds (int,int): start and end indices of points in the timestamp-sorted list of
                              trajectory points.
    """

    def __init__(self, trajID, bounds):
        """ Initialize pathlet with points from a trajectory.

            Args:
                trajID (int): ID of trajectory.
                bounds (int,int): start and end indices.
        """
        self.trajID = trajID
        self.bounds = bounds

    def __str__(self):
        """ Return string to be output while printing a pathlet object. """
        return "pathlet TrajID %d ; bounds (%d, %d)" % (self.trajID, self.bounds[0], self.bounds[1])

    def __eq__(self, other):
        """ Define == operator for pathlet objects. """
        return (self.trajID == other.trajID and self.bounds[0] == other.bounds[0] and self.bounds[1] == other.bounds[1])

    def __hash__(self):
        """ Define a hash function so that pathlets can be used as keys in a dict. """
        return hash((self.trajID, self.bounds[0], self.bounds[1]))


class subTraj(object):
    """ A class representing a subtrajectory.

    Exactly identical to a pathlet class. Defined as a class of its own for conceptual reasons.
    """

    def __init__(self, trajID, bounds):
        self.trajID = trajID
        self.bounds = bounds

    def __str__(self):
        return "Subtraj TrajID %d ; bounds (%d, %d)" % (self.trajID, self.bounds[0], self.bounds[1])

    def __eq__(self, other):
        return (self.trajID == other.trajID and self.bounds[0] == other.bounds[0] and self.bounds[1] == other.bounds[1])

    def __hash__(self):
        return hash((self.trajID, self.bounds[0], self.bounds[1]))
