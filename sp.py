import numpy as np

def SP_ComputeCompositeCT(laminat):
      # Calculate matrix serial/parallel matrices
      PP = getattr(laminat, "PP")
      PS = getattr(laminat, "PS")
      
      CTMatrx = laminat.getMatrxCT()
      CTFibre = laminat.getFibreCT()

      CTMatrxPP = np.dot(              PP,  np.dot( CTMatrx, np.transpose(PP) ))
      CTMatrxPS = np.dot(              PP,  np.dot( CTMatrx,              PS  ))
      CTMatrxSP = np.dot( np.transpose(PS), np.dot( CTMatrx, np.transpose(PP) ))
      CTMatrxSS = np.dot( np.transpose(PS), np.dot( CTMatrx,              PS  ))

      # Calculate fibre serial/parallel matrices
      CTFibrePP = np.dot(              PP,  np.dot( CTFibre, np.transpose(PP) ))
      CTFibrePS = np.dot(              PP,  np.dot( CTFibre,              PS  ))
      CTFibreSP = np.dot( np.transpose(PS), np.dot( CTFibre, np.transpose(PP) ))
      CTFibreSS = np.dot( np.transpose(PS), np.dot( CTFibre,              PS  ))

      kMatrx = getattr(laminat,"matrxPart")
      kFibre = getattr(laminat,"fibrePart")
            
      # Calculate A array
      arrayA = computeArrayA( laminat )

      # Compute composite PP PS SP SS

      CTCompoPP = kMatrx*CTMatrxPP + kFibre*CTFibrePP  +  \
                  kMatrx*kFibre*np.dot( (CTFibrePS - CTMatrxPS), np.dot( arrayA, (CTMatrxSP - CTFibreSP) ))
      
      CTCompoPS = kFibre * np.dot( CTFibrePS, np.dot( arrayA, CTMatrxSS ) )  +  \
                  kMatrx * np.dot( CTMatrxPS, np.dot( arrayA, CTFibreSS ) )
      
      CTCompoSP = kMatrx * np.dot( CTFibreSS, np.dot( arrayA, CTMatrxSP ) )  +  \
                  kFibre * np.dot( CTMatrxSS, np.dot( arrayA, CTFibreSP ) )
      
      CTCompoSS = 0.5 *  ( np.dot( CTMatrxSS, np.dot( arrayA, CTFibreSS ) )  +  \
                           np.dot( CTFibreSS, np.dot( arrayA, CTMatrxSS ) )  )
    
      #Obtain the final Composite Stiffness Matrix
      C = np.dot( np.transpose(PP), np.dot( CTCompoPP,              PP  ) )  +  \
          np.dot( np.transpose(PP), np.dot( CTCompoPS, np.transpose(PS) ) )  +  \
          np.dot(              PS , np.dot( CTCompoSP,              PP  ) )  +  \
          np.dot(              PS , np.dot( CTCompoSS, np.transpose(PS) ) ) 

      return C
      
def computeArrayA(laminat):
    #allocate( AI(a%ns,a%ns) )
    PP = getattr(laminat, "PP")
    PS = getattr(laminat, "PS")

    CTMatrx = laminat.getMatrxCT()
    CTFibre = laminat.getFibreCT()
    
    # Calculate matrix serial/parallel matrices
    CTMatrxSS = np.dot( np.transpose(PS), np.dot( CTMatrx, PS  ))

    # Calculate fibres serial/parallel matrices
    CTFibreSS = np.dot( np.transpose(PS), np.dot( CTFibre, PS  ))

    kMatrx = getattr(laminat,"matrxPart")
    kFibre = getattr(laminat,"fibrePart")

    # Calculate A matrix
    AI = kFibre*CTMatrxSS + kMatrx*CTFibreSS

    for i in range(0,len(getattr(laminat,"flagsSP"))):
          if PS[i][i] == 0.0:
              AI = np.delete(AI, i, 0)
              AI = np.delete(AI, i, 1)

    arrayAp = np.linalg.inv(AI)

    arrayA = np.zeros([len(getattr(laminat,"flagsSP")),len(getattr(laminat,"flagsSP"))])

    conectivitats = []
    counter = 0
    for flag in getattr(laminat,"flagsSP"):
        if flag != 1:
            conectivitats.append(counter)
        counter += 1

    for i in range(0,np.shape(arrayAp)[0]):
        for j in range(0,np.shape(arrayAp)[0]):
            arrayA[conectivitats[i],conectivitats[j]] = arrayAp[i,j]

    return arrayA

def computeStrainMatrxFibreSP(laminat):
    PP = getattr(laminat, "PP")
    PS = getattr(laminat, "PS")

    CTMatrx = laminat.getMatrxCT()
    CTFibre = laminat.getFibreCT()
                  
    kMatrx = getattr(laminat,"matrxPart")
    kFibre = getattr(laminat,"fibrePart")
      
    # Obtenemos las deformaciones del compuesto
    LayerStrain_t = getattr(laminat,"strain")
    
    ## Obtain Previous converged strains and constitutive tensors
    MatrxStrain_0 = np.zeros(6) #a%GP(igaus)%MatrxStrain_conv ! Matrix strains
    FibreStrain_0 = np.zeros(6) #a%GP(igaus)%FibreStrain_conv ! Fibre  strains
    
    # Eval parallel fiber strains before loop
    LayerStrainP_t = np.dot(             PP , LayerStrain_t)
    LayerStrainS_t = np.dot(np.transpose(PS), LayerStrain_t) 
    #MatrxStrainP_t = LayerStrainP_t 
    FibreStrainP_t = LayerStrainP_t 
    
    ## Obtain matrix and fiber parallel and serial strains
    MatrxStrainP_0 = np.dot(             PP, MatrxStrain_0)
    MatrxStrainS_0 = np.dot(np.transpose(PS),MatrxStrain_0)
    
    FibreStrainP_0 = np.dot(             PP, FibreStrain_0)
    FibreStrainS_0 = np.dot(np.transpose(PS),FibreStrain_0)
    
    ## Obtain previous strains in the composite
    CompoStrainP_0 = MatrxStrainP_0
    CompoStrainS_0 = kMatrx*MatrxStrainS_0 + kFibre*FibreStrainS_0
    CompoStrain_0  = np.dot(np.transpose(PP),CompoStrainP_0) + \
                     np.dot(             PS ,CompoStrainS_0)
    
    # Obtain the strain's increment in current step    
    CompoStrain_n  = LayerStrain_t - 0.0 # LayerStrain_conv
    CompoStrainP_n = np.dot(             PP, CompoStrain_n)
    CompoStrainS_n = np.dot(np.transpose(PS),CompoStrain_n)
    
    # Obtain some parallel and serial matrices, 
    # required to calculate matrix strain increment
    CTMatrxSS = np.dot( np.transpose(PS), np.dot( CTMatrx, PP  ))
    CTFibreSS = np.dot( np.transpose(PS), np.dot( CTFibre, PS  ))
    CTMatrxSP = np.dot( np.transpose(PS), np.dot( CTMatrx, np.transpose(PP) ))
    CTFibreSP = np.dot( np.transpose(PS), np.dot( CTFibre, np.transpose(PP) ))
    
    ArrayA = computeArrayA(laminat)
    # Calculate matrix parallel and serial strain increment
    MatrxStrainP_n = CompoStrainP_n
    MatrxStrainS_n = np.dot( ArrayA, \
                     ( np.dot(CTFibreSS,CompoStrainS_n) + \
                     ( kFibre * np.dot( CTFibreSP - CTMatrxSP,CompoStrainP_n)) ))
                    
    # Recompose strains to obtain prediction of matrix strain
    MatrxStrain_n = np.dot( np.transpose(PP), MatrxStrainP_n ) + \
                    np.dot(              PS , MatrxStrainS_n )
     
    MatrxStrain_t = MatrxStrain_0 + MatrxStrain_n
       
    # Starts loop to find correct matrix and fibre serie deformation
        
    MatrxStrainP_t = np.dot(             PP ,MatrxStrain_t) # Matrix parallel strains
    MatrxStrainS_t = np.dot(np.transpose(PS),MatrxStrain_t) # Matrix serial strains
    
    FibreStrainS_t = ( (1/kFibre) * LayerStrainS_t ) - ( (kMatrx/kFibre) * MatrxStrainS_t ) # Fibre serial strains
     
    FibreStrain_t  = np.dot( np.transpose(PP), FibreStrainP_t ) + \
                     np.dot(              PS , FibreStrainS_t )  
     
    # Compute MATRIX stresses (in local direction)
    MatrxStrain_n = MatrxStrain_t - MatrxStrain_0
    
    #call VoigtOpe%GetGradDispFromGpStrain( MatrxStrain_t, gradDispMatrx )         
    #call a%EMDs(1)%p%ComputeHistoryAndConstitutiveTensor( igaus, gradDispMatrx )
    #call a%EMDs(1)%p%ComputeStress( igaus, gradDispMatrx )
    #call a%EMDs(1)%p%GetStressTensorPointer( igaus, MatrxStress_t )
     
    # Compute FIBRE stresses (in local direction)
    FibreStrain_n = FibreStrain_t - FibreStrain_0
    
    #call VoigtOpe%GetGradDispFromGpStrain( FibreStrain_t, gradDispFibre )         
    #call a%EMDs(2)%p%ComputeHistoryAndConstitutiveTensor( igaus, gradDispFibre )
    #call a%EMDs(2)%p%ComputeStress( igaus, gradDispFibre )
    #call a%EMDs(2)%p%GetStressTensorPointer( igaus, FibreStress_t )

    return [MatrxStrain_t, FibreStrain_t]

def computeStrainMatrxFibreSPonlyFibre(laminat,tensionsSerie):
    PP = getattr(laminat, "PP")
    PS = getattr(laminat, "PS")

    CTMatrx = laminat.getMatrxCT()
    CTFibre = laminat.getFibreCT()
                  
    kMatrx = getattr(laminat,"matrxPart")
    kFibre = getattr(laminat,"fibrePart")
      
    LayerStrain_t = getattr(laminat,"strain")
    
    LayerStrainP_t = np.dot(PP, LayerStrain_t)
    LayerStrainS_t = np.dot(PS, LayerStrain_t) 
    
    FibreStrainP_t = LayerStrainP_t 

    # Reduim la matriu S i ompliu la diagonal de 1 per invertirla
    FibreS_P = np.dot(np.dot(PP, laminat.fibre.S),PP)
    for i in range(0,6):
        if FibreS_P[i][i] == 0.0:
            FibreS_P[i][i] = 1.0

    FibreStressP_t = np.dot(np.linalg.inv(FibreS_P),LayerStrainP_t-np.dot(laminat.fibre.S,tensionsSerie))

    FibreStress_t = tensionsSerie + np.dot(PP,FibreStressP_t)
    FibreStrain_t = np.dot(laminat.fibre.S,(FibreStress_t))

    FibreStrainP_t = np.dot(PP, FibreStrain_t)
    FibreStrainS_t = np.dot(PS, FibreStrain_t)

    MatrxStrainP_t = FibreStrainP_t
    MatrxStrainS_t = ( LayerStrainS_t - kFibre*FibreStrainS_t )/ kMatrx

    MatrxStrain_t = MatrxStrainP_t + MatrxStrainS_t

    return [MatrxStrain_t, FibreStrain_t]