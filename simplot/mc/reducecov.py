import numpy as np

################################################################################

class ReducedGaussian:
    def __init__(self, parameter_names, mu, cov, marginalise=[], conditional={}):
        mindices = self._get_indices(marginalise, parameter_names)
        cindices = self._get_indices(conditional.keys(), parameter_names)
        #input indices must be sorted
        mindices.sort()
        cindices.sort()
        observed = [conditional[parameter_names[ci]] for ci in cindices]
        self._parameter_names = [p for p in parameter_names if p not in marginalise and p not in conditional]
        self._mu = ReduceMeanVector()(mu, cov, marginalise=mindices, conditional=cindices, observed=observed)
        self._cov = ReduceCovarianceMatrix()(cov, marginalise=mindices, conditional=cindices)
        return

    def _get_indices(self, l, keylist):
        return [keylist.index(i) for i in l]

    @property
    def mu(self):
        return self._mu

    @property
    def cov(self):
        return self._cov

    @property
    def parameter_names(self):
        return self._parameter_names

################################################################################

class ReduceCovarianceMatrix:
    def __call__(self, cov, marginalise=[], conditional=[]):
        return self._reduce_covariance_matrix(cov, marginalise, conditional)

    def _reduce_covariance_matrix(self, cov, marginalise=[], conditional=[]):
        marginalise = list(sorted(marginalise))
        conditional = list(sorted(conditional))
        cov = np.copy(cov)
        if marginalise:
            cov = _strip_multiple_row_and_col(marginalise, cov)
            conditional = _update_rows_and_columns(marginalise, conditional)
        if conditional:
            invcov = _strip_multiple_row_and_col(conditional, np.linalg.inv(cov))
            cov = np.linalg.inv(invcov)
        return cov

################################################################################

class ReduceMeanVector:
    def __call__(self, mu, cov, marginalise=[], conditional=[], observed=[]):
        return self._reduce_mean_vector(mu, cov, marginalise, conditional, observed)

    def _reduce_mean_vector(self, mu, cov, marginalise, conditional, observed):
        if marginalise:
            mu = [x for ii, x in enumerate(mu) if ii not in marginalise]
            conditional = _update_rows_and_columns(marginalise, conditional)
        if conditional:
            indices1 = [ii for ii in xrange(len(mu)) if ii not in conditional]
            indices2 = conditional
            mu1 = np.array([mu[ii] for ii in indices1])
            mu2 = np.array([mu[ii] for ii in indices2])
            x = np.array(observed)
            M22 = _strip_multiple_row_and_col(indices1, cov)
            M12 = np.array([[cov[i1, i2] for i1 in indices1] for i2 in indices2])
            mu = mu1 +  np.dot(M12, np.dot(np.linalg.inv(M22),(x - mu2)))
        return mu

################################################################################

def _strip_multiple_row_and_col(nset, m):
    nset = list(nset)
    nset.sort()
    for ii in xrange(len(nset)):
        if nset[ii] >= 0:
            m = _strip_row_and_col(nset[ii], m)
        #shift every subsequent index down.
        for jj in xrange(len(nset)):
            if nset[jj] > nset[ii]:
                nset[jj] -= 1
            if nset[jj] == nset[ii]:
                nset[jj] = -1
    return m  

def _strip_row_and_col(n, m):
    N = m.shape[0]
    if n < 0 or n >= N:
        raise ValueError("bad index given, required: 0 < n < N, given: n=%s, N=%s" % (n, N))
    out = np.zeros((N-1, N-1))
    out[0:n,0:n] = m[0:n, 0:n]
    out[n:N-1, n:N-1] = m[n+1:N, n+1:N]
    out[0:n, n:N-1] = m[0:n, n+1:N]
    out[n:N-1, 0:n] = m[n+1:N, 0:n]
    return out

def _update_rows_and_columns(removed, indices):
    removed = list(removed)
    indices = list(indices)
    while len(removed) > 0:
        r = removed.pop()
        #update indices
        for ii in xrange(len(indices)):
            if indices[ii] > r:
                #shift index
                indices[ii] -= 1 
            elif indices[ii] == r:
                raise Exception("invalid index, this index was already removed")
        #update removed
        for ii in xrange(len(removed)):
            if removed[ii] > r:
                #shift index
                removed[ii] -= 1 
            elif removed[ii] == r:
                raise Exception("invalid index, this index was already removed")
    return indices
