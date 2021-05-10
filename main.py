import numpy as np
import math
import sp

class SimpleMaterial:
    def __init__(self, criterionType):
        self.stress = np.zeros(6)
        self.D      = None
        self.YC = yieldCriterion.factory(criterionType)

    def computeYC(self):
        self.YC.equivalentStress(self.stress)

    def setStressIndividually(self,i, value):
        self.stress[i] = value

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
    
    def equivalentStress(self,stress):
        S = stress
        eqStress = math.sqrt(0.5*((S[0]-S[1])**2+(S[1]-S[2])**2+(S[2]-S[0])**2)+3*(S[3]+S[4]+S[5]))
        return eqStress

class Laminate:
    def __init__(self, matrx, fibre, flagsSP):
        self.flagsSP = flagsSP
        self.matrx = matrx
        self.fibre = fibre
        self.matrxPart = 0.4
        self.fibrePart = 0.6
        self.stress = np.zeros(6)
        self.PP = np.zeros([len(self.flagsSP),len(self.flagsSP)])
        self.PS = np.zeros([len(self.flagsSP),len(self.flagsSP)])

        i = 0
        for i in range(0,len(self.flagsSP)):
            if self.flagsSP[i] == 1:
                self.PP[i,i] = 1
            else:
                self.PS[i,i] = 1
        
        self.CT = sp.SP_ComputeCompositeCT(self)

    def stressesDistribution(self):
        counter = 0
        for flag in self.flagsSP:
            if flag == 1:
                strain = np.dot(np.linalg.inv(self.CT),self.stress)
                self.matrx.setStressIndividually(counter, np.dot(getattr(self.matrx,"D"),np.transpose(strain))[counter] )
                self.fibre.setStressIndividually(counter, np.dot(getattr(self.fibre,"D"),np.transpose(strain))[counter] )
                counter += 1
            else:
                self.matrx.setStressIndividually(counter,self.stress[counter])
                self.fibre.setStressIndividually(counter,self.stress[counter])
                counter += 1

    def checkYield(self):
        self.matrx.computeYC()
        self.fibre.computeYC()

    def getMatrxCT(self):
        CT = getattr(self.matrx,"D")
        return CT
    def getFibreCT(self):
        CT = getattr(self.fibre,"D")
        return CT

matrx = SimpleMaterial("VonMisses")
setattr(matrx,"D",np.array([[1000.0,  50.0,  50.0,   0.0,   0.0,   0.0], \
                             [  50.0,1000.0,  50.0,   0.0,   0.0,   0.0], \
                             [  50.0,  50.0,1000.0,   0.0,   0.0,   0.0], \
                             [   0.0,   0.0,   0.0, 200.0,   0.0,   0.0], \
                             [   0.0,   0.0,   0.0,   0.0, 200.0,   0.0], \
                             [   0.0,   0.0,   0.0,   0.0,   0.0, 200.0]]))
fibre = SimpleMaterial("VonMisses")
setattr(fibre,"D",np.array([[5000.0, 200.0, 200.0,   0.0,   0.0,   0.0], \
                             [ 200.0,5000.0, 200.0,   0.0,   0.0,   0.0], \
                             [ 200.0, 200.0,5000.0,   0.0,   0.0,   0.0], \
                             [   0.0,   0.0,   0.0,1000.0,   0.0,   0.0], \
                             [   0.0,   0.0,   0.0,   0.0,1000.0,   0.0], \
                             [   0.0,   0.0,   0.0,   0.0,   0.0,1000.0]]))

laminat = Laminate(matrx,fibre,[1,0,0,0,0,0])
setattr(laminat,"stress",[50,0,0,0,0,0])
laminat.stressesDistribution()
laminat.checkYield()