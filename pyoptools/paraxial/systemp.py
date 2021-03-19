from collections import namedtuple
from collections.abc import MutableSequence
from abc import abstractmethod

from numpy import array, dot
PRay = namedtuple("PRay",["y","u"])


class PSurface:
    """
    Superclas to define paraxial surfaces
    """
    @classmethod
    @abstractmethod
    def fromdata(cls, data):
        """Función auxiliar que permite crear una instancia de la clase a
        partir delos datos contenidos en data. Se usa para facilitar la
        creación de las instancias al llenar un PSystem.
        Se debe sobrecargar, para sacar de data la información relevante.
        """
        pass
    
    def __init__(self, d, n=1):
        self.d = d
        self.n = n

    def p_mat(self):
        return array(((1., self.d), (0, 1)))

    def matrix(self):
        return dot(self.s_mat(), self.p_mat())

    @abstractmethod
    def s_mat(self):
        pass

    def propagate(self, ray):

        nray = dot(self.matrix(), array(ray))
        return PRay(nray)
    

class PPlane(PSurface):
    """
    Class to define a paraxial plane surface
    """
    @classmethod
    def fromdata(cls, data):
        n = data["n"]
        d = data["d"]
        return cls(d,n)

    def s_mat(self, n):
        """Matriz de la superficie correspondiente a la refracción

        """
        return array(((1,0),(0, n/self.n)))


class PSpherical(PSurface):
    @classmethod
    def fromdata(cls, data):
        r = data["r"]
        d = data["d"]
        n = data["n"]
        return cls(r, d, n)

    def __init__(self, r, d, n):
        super().__init__(d, n)
        self.r = r

    def s_mat(self, n):
        return array(((1, 0), ((n-self.n)/(self.r*self.n), n/self.n)))


class PObj(PPlane):

    def s_mat(self, n):
        """Los "objetos" no refractan
        """
        return array(((1, 0), (0, 1)))


surfdict = {'SPH': PSpherical.fromdata,
            'OBJ': PObj.fromdata,
            'PLA': PPlane.fromdata,
            }


class PSystem(MutableSequence):
    def __init__(self, l):

        self.surfaces=[]
        
        for s in l:
            tp, data = l
            if isinstance(tp, str):
                tp = surfdict[tp](data)
                self.surfaces.append(tp)
            else:
                raise ValueError
4
