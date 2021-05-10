import numpy as np
import math

class SimpleMaterial:
    def __init__(self, criterionType):
        self.stress = np.zeros(6)
        self.D      = np.zeros([6,6])
        self.yieldCriterion = yieldCriterion.factory(criterionType)

    def computeYC(self):
        self.yieldCriterion.equivalentStress(self.stress)

class yieldCriterion:
    def factory(type):
        if type == "VonMisses":     return VonMisses()
        assert 0, "Bad shape creation: " + type
    factory = staticmethod(factory)

    def __init__(self):
        None

class VonMisses(yieldCriterion):
    def __init__(self):
        super(VonMisses, self).__init__()
    
    def equivalentStress(self):
        S = self.stress
        eqStress = math.sqrt(0.5*((S[0]-S[1])**2+(S[1]-S[2])**2+(S[2]-S[0])*22)+3*(S[3]+S[4]+S[5]))
        return eqStress

class Laminate:
    def __init__(self, matrx, fibre):
        self.flagsSP = np.zeros(6)
        self.matrx = matrx
        self.fibre = fibre
        self.matrxPart = 0.4
        self.fibrePart = 0.6
        self.stress = np.zeros(6)

    def stressesDistribution(self):
        for flag in self.flagsSP:
            if flag == 1:
                strain = self.computeStrain()
                setattr(self.matrx,"stress", np.dot(getattr(self.matrx,"D"),np.transpose(strain) ) )
                setattr(self.fibre,"stress", np.dot(getattr(self.fibre,"D"),np.transpose(strain) ) )
            else:
                setattr(self.matrx,"stress",self.stress[flagsSP.index(flag)])
                setattr(self.fibre,"stress",self.stress[flagsSP.index(flag)])

    def computeStrain(self):
        compositeD = self.matrxPart*getattr(self.matrx,"D")+self.fibrePart*getattr(self.fibre,"D")
        return np.dot(np.linalg.inv(compositeD),self.stress)

    def checkYield(self):
        self.matrx.computeYC()
        self.fibre.computeYC()

matrx = SimpleMaterial("VonMisses")
setattr(matrx,"D",np.matrix([[1000.0,  50.0,  50.0,   0.0,   0.0,   0.0], \
                             [  50.0,1000.0,  50.0,   0.0,   0.0,   0.0], \
                             [  50.0,  50.0,1000.0,   0.0,   0.0,   0.0], \
                             [   0.0,   0.0,   0.0, 200.0,   0.0,   0.0], \
                             [   0.0,   0.0,   0.0,   0.0, 200.0,   0.0], \
                             [   0.0,   0.0,   0.0,   0.0,   0.0, 200.0]]))
fibre = SimpleMaterial("VonMisses")
setattr(fibre,"D",np.matrix([[5000.0, 200.0, 200.0,   0.0,   0.0,   0.0], \
                             [ 200.0,5000.0, 200.0,   0.0,   0.0,   0.0], \
                             [ 200.0, 200.0,5000.0,   0.0,   0.0,   0.0], \
                             [   0.0,   0.0,   0.0,1000.0,   0.0,   0.0], \
                             [   0.0,   0.0,   0.0,   0.0,1000.0,   0.0], \
                             [   0.0,   0.0,   0.0,   0.0,   0.0,1000.0]]))

laminat = Laminate(matrx,fibre)
setattr(laminat,"flagsSP",[1,0,0,0,0,0])
setattr(laminat,"stress",[50,0,0,0,0,0])
laminat.stressesDistribution()
laminat.checkYield()





