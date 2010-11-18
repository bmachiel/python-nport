import munkres
import numpy as np

def eigenshuffle(Asequence):
    """
    Consistent sorting for an eigenvalue/vector sequence
    
    Based on eigenshuffle.m for MATLAB by John D'Errico
    http://www.mathworks.com/matlabcentral/fileexchange/22885
    
    Python adaptation by Brecht Machiels
        <brecht.machiels@esat.kuleuven.be>
    
    Requires NumPy (http://numpy.scipy.org)
    and munkres.py by Brian M. Clapper
    (http://www.clapper.org/software/python/munkres/)
    
    Parameters
    ----------
    Asequence : sequence of square arrays
        An array of eigenvalue problems. If Asequence is a 3-d numeric array, 
        then each plane of Asequence must contain a square matrix that will be 
        used to call numpy.linalg.eig.
        
        numpy.linalg.eig will be called on each of these matrices to produce a 
        series of eigenvalues/vectors, one such set for each eigenvalue problem.
        
    Returns
    -------
    Dseq : 
        A pxn array of eigen values, sorted in order to be consistent with each 
        other and with the eigenvectors in Vseq.
    Vseq :
        A 3-d array (pxpxn) of eigenvectors. Each plane of the array will be 
        sorted into a consistent order with the other eigenvalue problems. The 
        ordering chosen will be one that maximizes the energy of the consecutive
        eigensystems relative to each other.
        
    See Also
    --------
    numpy.linalg.eig
    
    Example
    -------
        import numpy as np
        from eigenshuffle import eigenshuffle
    
        def Efun(t):
            return np.array([
                [1,     2*t+1 , t**2 ,   t**3],
                [2*t+1, 2-t   , t**2 , 1-t**3],
                [t**2 , t**2  , 3-2*t,   t**2],
                [t**3 , 1-t**3, t**2 ,  4-3*t]])
        
        Aseq = np.zeros( (21, 4, 4) )
        for i in range(21):
            Aseq[i] = Efun((i + 1 - 11)/10.0)

        [Dseq,Vseq] = eigenshuffle(Aseq)
        
    To see that eigenshuffle has done its work correctly, look at the
    eigenvalues in sequence, after the shuffle.
    
        t = np.arange(-1,1,0.1)
        print np.hstack([np.asarray([t]).T, Dseq])
        
        ... TODO ...
    """
    
    # Is Asequence a 3-d array?
    Ashape = np.shape(Asequence)
    if Ashape[-1] != Ashape[-2]:
        raise Exception, "Asequence must be a (nxpxp) array of eigen-problems, each of size pxp"
    p = Ashape[-1]
    if len(Ashape) < 3:
        n = 1
        Asequence = np.asarray([Asequence], dtype=complex)
    else:
        n = Ashape[0]

    # the initial eigenvalues/vectors in nominal order
    Vseq = np.zeros( (n, p, p), dtype=complex )
    Dseq = np.zeros( (n, p), dtype=complex )

    for i in range(n):
        D, V = np.linalg.eig( Asequence[i] )
        # initial ordering is purely in decreasing order.
        # If any are complex, the sort is in terms of the
        # real part.
        tags = np.argsort(D.real, axis=0)[::-1]
        
        Dseq[i] = D[:, tags]
        Vseq[i] = V[:, tags]
    
    # now, treat each eigenproblem in sequence (after the first one.)
    m = munkres.Munkres()
    for i in range(1, n):
        # compute distance between systems
        D1 = Dseq[i - 1]
        D2 = Dseq[i]
        V1 = Vseq[i - 1]
        V2 = Vseq[i]
        dist = (1 - np.abs( np.dot( np.transpose(V1), V2 ) ) ) \
            * np.sqrt( 
                distancematrix(D1.real, D2.real)**2 +
                distancematrix(D1.imag, D2.imag)**2
            )

#        print i, np.transpose(dist), '\n'

        # Is there a best permutation? use munkres.
        reorder = m.compute(np.transpose(dist))
        reorder = [coord[1] for coord in reorder]
#        print i+1
#        print reorder
#        print "Vseq"
#        print Vseq[i]

        Vs = Vseq[i]
#        Vseq[i] = Vseq[i, :, reorder]
        Vseq[i] = Vseq[i][:, reorder]
#        print "Vseq reorder"
#        print Vseq[i]
        Dseq[i] = Dseq[i, reorder]
#        print "Dseq reorder"
#        print Dseq[i], '\n'

        # also ensure the signs of each eigenvector pair
        # were consistent if possible
        S = np.squeeze( np.sum( Vseq[i - 1] * Vseq[i], 0 ).real ) < 0
       
        Vseq[i] = Vseq[i] * (-S * 2 - 1)
#        print "sign"
#        print Vseq[i]

    return Dseq, Vseq


def distancematrix(vec1, vec2):
    """
    simple interpoint distance matrix
    """
    v1, v2 = np.meshgrid(vec1, vec2)
    return np.abs(v1 - v2)
