import numpy as np

def SP_ComputeCompositeCT(laminat):
      # Subroutine auxiliar variables
      
      #Real(8) :: CTCompoPP(a%np,a%np),  & ! Composite stiffness - parallel-parallel
      #           CTCompoPS(a%np,a%ns),  & ! Composite stiffness - parallel-serial
      #           CTCompoSP(a%ns,a%np),  & ! Composite stiffness - serial-parallel 
      #           CTCompoSS(a%ns,a%ns)     ! Composite stiffness - serial-serial
      #
      #Real(8) :: CTFibrePP(a%np,a%np),  & ! Fibres stiffness - parallel-parallel
      #           CTFibrePS(a%np,a%ns),  & ! Fibres stiffness - parallel-serial
      #           CTFibreSP(a%ns,a%np),  & ! Fibres stiffness - serial-parallel 
      #           CTFibreSS(a%ns,a%ns)     ! Fibres stiffness - serial-serial
      #
      #Real(8) :: CTMatrxPP(a%np,a%np),  & ! Matrix stiffness - parallel-parallel
      #           CTMatrxPS(a%np,a%ns),  & ! Matrix stiffness - parallel-serial
      #           CTMatrxSP(a%ns,a%np),  & ! Matrix stiffness - serial-parallel 
      #           CTMatrxSS(a%ns,a%ns)     ! Matrix stiffness - serial-serial
      #
      #Real(8), allocatable :: arrayA(:,:), AI(:,:) ! Auxiliar Matrices

      #allocate( arrayA(a%ns,a%ns) )

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
    #real(rp), pointer :: CTMatrx(:,:)
    #real(rp), pointer :: CTFibre(:,:)                           
                 
    #! Subroutine auxiliar variables

    #Real(8) :: kMatrx, & ! Matrix volumetric participation
    #           kFibre, & ! Fibres volumetric particiation
    #           deter     ! Dummy argument for matrix inversion subroutine

    #Real(8) :: CTMatrxSS(a%ns,a%ns), & ! Matrix stiffness - serial-serial
    #           CTFibreSS(a%ns,a%ns)    ! Fibres stiffness - serial-serial
    #
    #Real(8), allocatable :: arrayA(:,:), AI(:,:)

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