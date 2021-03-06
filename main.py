#%%
import numpy as np
import math

from numpy import linalg
import sp

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

from scipy.optimize import fsolve
from scipy.optimize import brentq

class SimpleMaterial:
    def __init__(self, criterionType):
        self.stress = np.zeros(6)
        self.strain = np.zeros(6)
        self.E = None
        self.v = None
        self.S = None
        self.D = None
        self.YC = yieldCriterion.factory(criterionType)

    def setYC(self,data):
        self.YC.setup(data)

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
        alpha = self.YC.computeAlpha(self.stress)
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

    def mapStress(self,mappingTensor):
        self.stress[0] = self.stress[0]*mappingTensor[0]
        self.stress[1] = self.stress[1]*mappingTensor[1]
        self.stress[2] = self.stress[2]*mappingTensor[2]
        self.stress[3] = self.stress[3]*mappingTensor[3]
        self.stress[4] = self.stress[4]*mappingTensor[4]
        self.stress[5] = self.stress[5]*mappingTensor[5]


class yieldCriterion:
    def factory(type):
        if type == "VonMisses":     return VonMisses()
        if type == "MaxStrain":     return MaxStrain()
        if type == "Ortotrop":      return Ortotrop()
        if type == "Drucker":       return Drucker()
        assert 0, "Bad shape creation: " + type
    factory = staticmethod(factory)

    def __init__(self):
        None

class VonMisses(yieldCriterion):
    def __init__(self):
        super(VonMisses, self).__init__()
        self.thresshold = None

    def setup(self,data):
        self.thresshold = data[0][0]

    def computeAlpha(self,stress):
        t = self.thresshold
        def fun(x):
            S = stress*x
            alpha = + (t**-2)*(S[0]**2+S[1]**2+S[2]**2) \
                    - (t**-2)*(S[1]*S[2]+S[2]*S[0]+S[0]*S[1]) \
                    + (3*t**-2)*(S[3]**2+S[4]**2+S[5]**2) \
                    - 1
            return alpha
        
        alphaConv, info = brentq(fun, 0.0, 10000, full_output=True)
    
        return alphaConv

class Drucker(yieldCriterion):
    def __init__(self):
        super(Drucker, self).__init__()
        self.thresshold = None

    def setup(self,data):
        self.tt = data[0]
        self.tc = data[1]
        self.ts = data[2]

    def computeAlpha(self,stress):
        tt = self.tt
        tc = self.tc
        ts = self.ts
        def fun(x):
            S = stress*x

            sum1 = (tc[0]+tt[0])/(2*tc[0]*tt[0])
            sum2 = (tc[1]+tt[1])/(2*tc[1]*tt[1])
            sum3 = (tc[2]+tt[2])/(2*tc[2]*tt[2])
            F    = 0.5*(-sum1**2+sum2**2+sum3**2)
            G    = 0.5*(+sum1**2-sum2**2+sum3**2)
            H    = 0.5*(+sum1**2+sum2**2-sum3**2)
            L    = 1/(2*(ts[0])**2)
            M    = 1/(2*(ts[1])**2)
            N    = 1/(2*(ts[2])**2)
            I    = (tc[0]-tt[0])/(2*tc[0]*tt[0])
            J    = (tc[1]-tt[1])/(2*tc[1]*tt[1])
            K    = (tc[2]-tt[2])/(2*tc[2]*tt[2])

            alpha = math.sqrt(F*(S[1]-S[2])**2 + G*(S[2]-S[0])**2 + H*(S[0]-S[1])**2 + \
                    (2*L*S[3]**2) + \
                    (2*M*S[4]**2) + \
                    (2*N*S[5]**2)) + \
                    I*S[0] + J*S[1] + K*S[2] - 1

            return alpha
        
        alphaConv, info = brentq(fun, 0.0, 10000, full_output=True)
    
        return alphaConv

class Ortotrop(yieldCriterion):
    def __init__(self):
        super(Ortotrop, self).__init__()
        self.f = None

    def setup(self,data):
        self.f = data[0]
    
    def computeAlpha(self,stress):
        f = self.f
        def fun(x):
            S = stress*x
            theta = [ 0.5*((1)-(1)), 0.5*(((f[1])/(f[1]))-((f[1])/(f[1]))), 0.5*((1)-(1)) ]

            delta = np.zeros(3)
            delta[0] = (+f[2]**2*f[0]**2-f[2]**2*f[1]**2+f[1]**2*f[0]**2)/(f[1]*f[2]*f[0]**2)
            delta[1] = (-f[2]**2*f[0]**2+f[2]**2*f[1]**2+f[1]**2*f[0]**2)/(f[0]*f[2]*f[1]**2)
            delta[2] = (+f[2]**2*f[0]**2+f[2]**2*f[1]**2-f[1]**2*f[0]**2)/(f[1]*f[0]*f[2]**2)
    
            alpha = + ( (S[0]**2)*(f[0]**-2) + (S[1]**2)*(f[1]**-2) + (S[2]**2)*(f[2]**-2) ) \
                    - delta[0]*((S[1]*S[2])/(f[1]*f[2])) \
                    - delta[1]*((S[2]*S[0])/(f[2]*f[0])) \
                    - delta[2]*((S[0]*S[1])/(f[0]*f[1])) \
                    + ( ((S[3]**2)/(f[3]**2)) + ((S[4]**2)/(f[4]**2)) + ((S[5]**2)/(f[5]**2)) ) \
                    + 2*( (theta[0]*(S[0]/f[0])) + (theta[1]*(S[1]/f[1])) + (theta[2]*(S[2]/f[2])) ) \
                    - 1
            return alpha
        
        alphaConv, info = brentq(fun, 0, 1000, full_output=True)
    
        return alphaConv

class MaxStrain(yieldCriterion):
    def __init__(self):
        super(MaxStrain, self).__init__()

    def setup(self,data):
        None
    
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

        [MatrxStrain_t, FibreStrain_t] = sp.computeStrainMatrxFibreSP(self)

        setattr(self.matrx, "stress", np.dot(getattr(self.matrx,"D"),MatrxStrain_t) )
        setattr(self.fibre, "stress", np.dot(getattr(self.fibre,"D"),FibreStrain_t) )

        setattr(self.matrx, "strain", MatrxStrain_t )
        setattr(self.fibre, "strain", FibreStrain_t )

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
        dirZ = 2

        FibraX = []
        FibraY = []
        FibraZ = []
        MatrxX = []
        MatrxY = []
        MatrxZ = []
        CompoX = []
        CompoY = []
        CompoZ = []

        stress = np.zeros(6)

        plt.figure(figsize=(10,10))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.grid()
        plt.xlabel("X axis label")
        plt.ylabel("Y axis label")

        for angle1 in np.arange(0,90.05,2.0):
             
            stress = np.zeros(6)
            module  = 10
            stress[dirZ] = module*math.sin(math.radians(angle1))
            moduleP = module*math.cos(math.radians(angle1))

            for angle2 in np.arange(0,90.05,2.0):
                #angle = 90-angle
    
                stress[dirX] = moduleP*math.cos(math.radians(angle2))
                stress[dirY] = moduleP*math.sin(math.radians(angle2))
                
                self.stress = stress
                laminat.stressesDistribution()
                
                alphaMatrx = self.matrx.computeAlpha()
    
                # Portem el laminat al punt de rotura de la matriu i ho verifiquem mirant lalpha daquesta
                self.stress = stress*alphaMatrx
                laminat.stressesDistribution()
                alphaMatrx = self.matrx.computeAlpha()
    
                ## IDEA FERMIN
                #[MatrxStrain_t, FibreStrain_t] = sp.computeStrainMatrxFibreSPonlyFibre(laminat,np.dot(self.PS, alphaMatrx*self.matrx.stress))
                ##setattr(self.matrx, "stress", np.dot(getattr(self.matrx,"D"),MatrxStrain_t) )
                #setattr(self.fibre, "stress", np.dot(getattr(self.fibre,"D"),FibreStrain_t) )
                ##setattr(self.matrx, "strain", MatrxStrain_t )
                #setattr(self.fibre, "strain", FibreStrain_t )
                #alphaFibre = self.fibre.computeAlpha()
                #fibreStress = getattr(self.fibre, "stress")
                #setattr(self.fibre, "stress", fibreStress*alphaFibre )
    
                ####setattr(self.matrx, "E", 0.0001)
                ####setattr(self.matrx, "v", 0.0)
                ####self.matrx.computeG()
                ####self.matrx.computeS()
                ####self.matrx.computeD()
                ####laminat.stressesDistribution()
                ####setattr(self.matrx, "E", 4670)
                ####setattr(self.matrx, "v", 0.38)
                ####self.matrx.computeG()
                ####self.matrx.computeS()
                ####self.matrx.computeD()
                ## ACABA AQUI EL FERMIN
    
                alphaFibre = self.fibre.computeAlpha()
    
                TensionUltimaP = (getattr(self.fibre,"stress")*self.fibrePart + \
                                  getattr(self.matrx,"stress")*self.matrxPart)
                #TensionUltimaP = (stress*alphaFibre)
                TensionUltimaS = (self.stress*alphaMatrx)
                #TensionUltimaMatrx = (getattr(self.matrx,"stress")*alphaMatrx)
                #TensionUltimaMatrx[0] = TensionUltimaMatrx[0]*self.matrxPart
                #TensionUltimaFibra = (getattr(self.fibre,"stress")*alphaFibre)
                #TensionUltimaFibra[0] = TensionUltimaFibra[0]*self.fibrePart
                TensionUltimaMatrx = self.matrx.stress
                TensionUltimaFibra = self.fibre.stress
    
                TU = np.array([ TensionUltimaP[0], \
                                TensionUltimaS[1], \
                                TensionUltimaS[2], \
                                TensionUltimaS[3], \
                                TensionUltimaS[4], \
                                TensionUltimaS[5], ])
    
                FibraX.append( TensionUltimaFibra[dirX] )
                FibraY.append( TensionUltimaFibra[dirY] )
                FibraZ.append( TensionUltimaFibra[dirZ] )
                MatrxX.append( TensionUltimaMatrx[dirX] )
                MatrxY.append( TensionUltimaMatrx[dirY] )
                MatrxZ.append( TensionUltimaMatrx[dirZ] )
                CompoX.append( TU[dirX] )
                CompoY.append( TU[dirY] )
                CompoZ.append( TU[dirZ] )
    
                print(str(angle1) + "   " + \
                      str(angle2) + "   " + \
                      str(alphaMatrx) + "   " + \
                      str(alphaFibre) + "   " + \
                      str(stress[dirX]) + " " + \
                      str(stress[dirY]) )

            ax.scatter(CompoX,CompoY,CompoZ,linewidth=0.2, antialiased=True)
            plt.show() 

        #ax.scatter(CompoX,CompoY,CompoZ,linewidth=0.2, antialiased=True)
        #plt.grid()
        #plt.xlabel("X axis label")
        #plt.ylabel("Y axis label")
        #plt.show() 

        #ax.scatter(CompoX,CompoY,CompoZ,linewidth=0.2, antialiased=True)
        ###plt.xlim([750.0,1000.0])
        ###plt.ylim([20.0,40.0])
        ##plt.plot_surface(FibraX,FibraY,FibraZ,color='green') #,'o'
        ##plt.plot_surface(MatrxX,MatrxY,MatrxZ,color='red'  ) #,'o'
        ##plt.plot_surface(CompoX,CompoY,CompoZ,'-.',color='black') #,'o'
        #plt.grid()
        #plt.xlabel("X axis label")
        #plt.ylabel("Y axis label")
        #plt.show()   

epsilonYC = 800/86900
sigmaM = epsilonYC*4670     

#matrx = SimpleMaterial("Ortotrop")
#matrx.setYC([[sigmaM,30.0,30.0,30.0,30.0,30.0]])
#matrx = SimpleMaterial("VonMisses")
#matrx.setYC([[30.0]])
matrx = SimpleMaterial("Drucker")
matrx.setYC([[30.0,30.0,30.0],[-30.0,-30.0,-30.0],[30.0,30.0,30.0]])
setattr(matrx, "E", 4670)
setattr(matrx, "v", 0.38)
matrx.computeG()
matrx.computeS()
matrx.computeD()

fibre = SimpleMaterial("VonMisses")
fibre.setYC([[800.0]])
setattr(fibre, "E", 86900)
setattr(fibre, "v", 0.22)
fibre.computeG()
fibre.computeS()
fibre.computeD()

laminat = Laminate(matrx,fibre,[1,0,0,0,0,0])
laminat.computePoint()
