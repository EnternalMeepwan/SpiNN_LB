# 
#   This is a basic implementation of a lattice Boltzmann method on SpiNNaker.
#   This project is currently being released under the GPL 3.0 license. Use it Free as in
#   Freedom. 
#  
#   The up-to-date information of SpiNNaker Project can be found here:
#   https://spinnakermanchester.github.io/
#  

from pacman.model.graphs.machine import MachineEdge


class LatticeEdge(MachineEdge):
    """
    Used for conjunction with a lattice
    """

    def __init__(self, pre_vertex, post_vertex, compass, label=None):
        MachineEdge.__init__(
            self, pre_vertex, post_vertex, label=label)
        self._compass = compass

    @property
    def compass(self):
        return self._compass

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "LatticeEdge: {}:{}".format(self._compass, self._label)