#%%
import numpy as np
import math
import sp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import brentq

class SimpleMaterial:
    def __init__(self, criterionType, thresholdValue):
        self.stress = np.zeros(6)
        self.strain = np.zeros(6)
        self.E = None
        self.v = None
        self.S = None
        self.D = None
        self.YC = yieldCriterion.factory(criterionType)
        self.threshold = thresholdValue

    def computeYC(self):
        return self.YC.equivalentStress(self.stress)
        
    def computeYCM(self):
        return self.YC.equivalentStressM(self.stress)

    def computeEquivalentStress(self):
        eq = self.computeYC()
        return eq
        
    def computeEquivalentStressM(self):
        eq = self.computeYCM()
        return eq

    def computeAlpha(self):
        alpha = (self.threshold / self.computeYC())
        return alpha

    def computeAlphaM(self):
        alpha = (self.threshold / self.computeYCM())
        return alpha

    def computeG(self):
        self.G = self.E/(2*(1+self.v))

    def computeS(self):
        E = self.E
        v = self.v
        G = self.G
        self.S = np.array([[  1/E, -v/E, -v/E, 0.0, 0.0, 0.0], \
                           [ -v/E,  1/E, -v/E, 0.0, 0.0, 0.0], \
                           [ -v/E, -v/E,  1/E, 0.0, 0.0, 0.0], \
                           [  0.0,  0.0,  0.0, 1/G, 0.0, 0.0], \
                           [  0.0,  0.0,  0.0, 0.0, 1/G, 0.0], \
                           [  0.0,  0.0,  0.0, 0.0, 0.0, 1/G]])

    def computeD(self):
        self.D = np.linalg.inv(self.S)


class yieldCriterion:
    def factory(type):
        if type == "VonMisses":     return VonMisses()
        if type == "MaxStrain":            return MaxStrain()
        if type == "Nonee":            return Nonee()
        assert 0, "Bad shape creation: " + type
    factory = staticmethod(factory)

    def __init__(self):
        None

class Nonee(yieldCriterion):
    def __init__(self):
        super(Nonee, self).__init__()
    
    def equivalentStress(self,stress):
        return 9e99

class VonMisses(yieldCriterion):
    def __init__(self):
        super(VonMisses, self).__init__()
    
    def equivalentStress(self,stress):
        S = stress
        eqStress = math.sqrt( ( (S[0]-S[1])**2 + \
                                (S[1]-S[2])**2 + \
                                (S[2]-S[0])**2 + \
                                6*(S[3]**2+S[4]**2+S[5]**2) )/2 )
        return eqStress
    
    def equivalentStressM(self,stress):
        S = stress
        S[0] = S[0]*(0.2791220556745182)
        eqStress = math.sqrt( ( (S[0]-S[1])**2 + \
                                (S[1]-S[2])**2 + \
                                (S[2]-S[0])**2 + \
                                6*(S[3]**2+S[4]**2+S[5]**2) )/2 )
        return eqStress

def computeFailedP(stress):
    S = stress
    #dintreLarrel = 3*( - (S[1]**2) + 2*S[1]*S[2] - (S[2]**2) - 4*(S[3]+S[4]+S[5]) )
    #eqStress = 0.5*math.sqrt( dintreLarrel  )
    def f(x):
        return ((2000**2-S[1]**2)/(x-S[1]))-x
    
    eqStress, info = brentq(f, 3000, 100, full_output=True)
    
    S[0] = eqStress

    return S

class MaxStrain(yieldCriterion):
    def __init__(self):
        super(MaxStrain, self).__init__()
    
    def equivalentStress(self,stress):
        eqStress = stress[0]-0.22*(stress[1]+stress[2])
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
        self.strain = np.dot(np.linalg.inv(self.CT),self.stress)
        verificationStressC = []

        [MatrxStrain_t, FibreStrain_t] = sp.computeStrainMatrxFibreSP(self)

        #MatrxStrain_t[0] = MatrxStrain_t[0]*(351/1253)

        setattr(self.matrx, "stress", np.dot(getattr(self.matrx,"D"),np.transpose(MatrxStrain_t)) )
        setattr(self.fibre, "stress", np.dot(getattr(self.fibre,"D"),np.transpose(FibreStrain_t)) )

        setattr(self.matrx, "strain", MatrxStrain_t )
        setattr(self.fibre, "strain", FibreStrain_t )

        verificationStressC.append(getattr(self.matrx,"stress")*self.matrxPart + getattr(self.fibre,"stress")*self.fibrePart)
        # verification

        #error = math.sqrt( np.dot(np.transpose(verificationStressC - stress), verificationStressC - stress) )
        #if error < 1.0e-5:
        #    print("Les tensions a matriu i fibra no son correctes")

    def checkYield(self):
        self.matrx.computeYC()
        self.fibre.computeYC()

    def getMatrxCT(self):
        CT = getattr(self.matrx,"D")
        return CT
    def getFibreCT(self):
        CT = getattr(self.fibre,"D")
        return CT

    def computePoint(self):
        dirX = 0
        dirY = 1

        FibraX = [0.0]
        FibraY = [0.0]
        MatrxX = [0.0]
        MatrxY = [0.0]
        CompoX = [0.0]
        CompoY = [0.0]

        #stress = np.zeros(6)

        #stress[dirX] = 1.0
        #stress[dirY] = 0.0
        #self.stress = stress
        #laminat.stressesDistribution()
        #alphaMatrx = self.matrx.computeAlpha()
        #alphaFibre = self.fibre.computeAlpha()

        #coef = alphaMatrx/alphaFibre

        print("angle   alphaMatrx alphaFibre TU1 TU2 TU3 TU4 TU5 TU6")
        for angle in np.arange(0,91.0,1.0):
            #angle = 90-angle

            stress = np.zeros(6)
            stress[dirX] = math.cos(math.radians(angle))
            stress[dirY] = math.sin(math.radians(angle))
            
            self.stress = stress
            laminat.stressesDistribution()
            
            alphaMatrx = self.matrx.computeAlphaM()
            #AAA = computeFailedP(getattr(self.matrx,"stress")*alphaMatrx)
            #setattr(self.fibre, "stress", AAA)
            alphaFibre = self.fibre.computeAlpha()

            TensionUltimaP = (getattr(self.fibre,"stress")*self.fibrePart+getattr(self.matrx,"stress")*self.matrxPart)
            #TensionUltimaP = (stress*alphaFibre)
            TensionUltimaS = (stress*alphaMatrx)
            TensionUltimaMatrx = (getattr(self.matrx,"stress")*alphaMatrx)
            TensionUltimaFibra = (getattr(self.fibre,"stress")*alphaFibre)
            TensionUltimaMatrx[0] = TensionUltimaMatrx[0]*self.matrxPart
            TensionUltimaFibra[0] = TensionUltimaFibra[0]*self.fibrePart

            TU = np.array([ TensionUltimaP[0], \
                            TensionUltimaS[1], \
                            TensionUltimaS[2], \
                            TensionUltimaS[3], \
                            TensionUltimaS[4], \
                            TensionUltimaS[5], ])

            compo = stress*alphaMatrx

            TUm = (stress*alphaMatrx)
            TUf = (stress*alphaFibre)

            FibraX.append( TensionUltimaFibra[dirX] )
            FibraY.append( TensionUltimaFibra[dirY] )
            MatrxX.append( TensionUltimaMatrx[dirX] )
            MatrxY.append( TensionUltimaMatrx[dirY] )
            CompoX.append( compo[dirX] )
            CompoY.append( compo[dirY] )

            print(str(angle) + "   " + \
                  str(alphaMatrx) + "   " + \
                  str(alphaFibre) + "   " + \
                  #str(TU[0]) + " " + \
                  #str(TU[1]) + " " + \
                  str(stress[dirX]) + " " + \
                  str(stress[dirY]) + " " + \
                  str(getattr(self.matrx,"stress")[0]) + " " + \
                  str(getattr(self.matrx,"stress")[1]) + " " + \
                  str(getattr(self.matrx,"stress")[2]) + " " + \
                  str(getattr(self.fibre,"stress")[0]) + " " + \
                  str(getattr(self.fibre,"stress")[1]) + " " + \
                  str(getattr(self.fibre,"stress")[2]) )

            #X.append( self.stress[0]*alphaMatrx )
            #Y.append( self.stress[3]*alphaMatrx )
            #X.append( self.stress[0]*alphaFibre )
            #Y.append( self.stress[3]*alphaFibre )

            #if alphaMatrx < alphaFibre: # Alpha mas pequeña es que rompe antes
            #    self.stress = stress*alphaMatrx
            #    laminat.stressesDistribution()
            #    alphaMatrx = self.matrx.computeAlpha()
            #    alphaFibre = self.fibre.computeAlpha()
            #    X.append( self.stress[0]*alphaMatrx )
            #    Y.append( self.stress[1]*alphaMatrx )
            #else:
            #    self.stress = stress*alphaFibre
            #    laminat.stressesDistribution()
            #    alphaMatrx = self.matrx.computeAlpha()
            #    alphaFibre = self.fibre.computeAlpha()
            #    X.append( self.stress[0]*alphaFibre )
            #    Y.append( self.stress[1]*alphaFibre )

        plt.figure(figsize=(10,10))
        #plt.xlim([750.0,1000.0])
        #plt.ylim([20.0,40.0])
        plt.plot(FibraX,FibraY,'o',color='green') #,'o'
        plt.plot(MatrxX,MatrxY,'o',color='red'  ) #,'o'
        plt.plot(CompoX,CompoY,'o',color='black') #,'o'
        plt.grid()
        plt.show()        

matrx = SimpleMaterial("VonMisses", 30.0)
setattr(matrx, "E", 4670)
setattr(matrx, "v", 0.38)
matrx.computeG()
matrx.computeS()
matrx.computeD()

fibre = SimpleMaterial("VonMisses", 2000.0)
setattr(fibre, "E", 86900)
setattr(fibre, "v", 0.22)
fibre.computeG()
fibre.computeS()
fibre.computeD()

laminat = Laminate(matrx,fibre,[1,0,0,0,0,0])
laminat.computePoint()
