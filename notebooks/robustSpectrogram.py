# Returns: xEst, freq, time, numIter
# data: 1D time series vector
#
def robustSpectrogram(data, fs, window, alpha, tol=0.005, maxIter=30):
    ## Initialize parameters
    numSamples = len(data)
    W = window # window size of 1 second
    K = W    # indice over freq. bands
    N = numSamples//W # number of windows along signal time

    F = np.zeros((K,K))
    Q = np.eye(K) * 0.001
    signal = data
    
    # Fourier along 1000 time points and 1000 frequency points
    for l in range(0, W):
        for k in range(0, K//2):
            F[l,k] = np.cos(2*np.pi*l*(k-1)/K)
            F[l,k+K//2] = np.sin(2*np.pi*l * (k-1)/K)
    # run iterations
    numIter = 1
    while(numIter <= maxIter):
        ## Initialize and Kalman Filter
        xPrior = np.zeros((K,N+1))     # x(n+1|n)
        xPosterior = np.zeros((K,N+1)) # x(n|n)
        pPrior = np.zeros((K,K,N+1))   # p(n+1|n)
        pPosterior = np.zeros((K,K,N+1)) # p(n|n)
        pPrior[:,:,0] = np.eye(K) # set initial priors to have variance of 1, cov = 0

        ## Step 1: Kalman Filter
        for n in range(0, N):
            y = signal[n*W:(n+1)*W]
            
            ## 01: Update priors
            xPrior[:,n+1] = xPosterior[:,n]
            pPrior[:,:,n+1] = pPosterior[:,:,n] + Q
            
            ## 02: Compute Kalman gain
            kGain = np.dot(pPrior[:,:,n+1], F.T).dot(np.linalg.inv(np.dot(F,pPrior[:,:,n+1]).dot(F.T) + np.eye(K)))
            
            ## 03: Compute Posteriors
            xPosterior[:,n+1] = xPrior[:,n+1] + np.dot(kGain, y-np.dot(F,xPrior[:,n+1]))
            pPosterior[:,:,n+1] = pPrior[:,:,n+1] - np.dot(kGain, F).dot(pPrior[:,:,n+1])
        
        
        ## Remove initial conditions
        xPrior = xPrior[:,1:N+1]
        xPosterior = xPostrior[:, 1:N+1]
        pPrior = pPrior[:,:,1:N+1]
        pPosterior = pPosterior[:,:,1:N+1]
        
        ## Step 2: Kalman Smoother
        xSmooth = xPrior
        pSmooth = pPrior
        
        for n in range(N-1,-1,-1):
            # compute Kalman gain smoother
            B = np.dot(pPosterior[:,:,n], np.linalg.inv(pPrior[:,:,n+1]))
            xSmooth[:,n] = xPrior[:,n] + np.dot(B, (xSmooth[:,n+1] - xPrior[:,n+1]))
            pSmooth[:,:,n] = pPosterior[:,:,n] + np.dot(B, (pSmooth[:,:,n+1] - pPrior[:,:,n+1]).dot(B.T))
            
        ## Step 3: Check Breakpoint
        if numIter > 1 and np.linalg.norm(xSmooth-xPrev, 'fro')/np.linalg.norm(xPrev,'fro') < tol:
            break
            
        ## Step 4: Update Q
        Q = np.zeros((K,K))
        for k in range(0, K):
            qTemp = 0
            # perform summation of adjacent state's squared error
            for n in range(1, N):
                qTemp += (xSmooth[k,n] - xSmooth[k,n-1])**2
            Q[k,k] = ((qTemp + np.finfo(float).eps**2)**(1/2)) / alpha
        
        xPrev = xSmooth
        numIter += 1
        print "iteration: ", numIter
        
    xEst = xSmooth[0:K//2,:] - xSmooth[K//2:W,:]*1j
    freq = np.arange(0, K//2, 1.)*fs/K + 1 # frequency to half of sampling freq
    time = np.arange(0, N)*W/fs
    
    return xEst, freq, time, numIter