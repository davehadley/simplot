import hashlib
import os
import numpy
import cPickle as pickle
import ROOT

_DEFAULT_TMPDIR = "/tmp/simplot_cache"

###############################################################################

def cache(uniquestr, callable_, filelist=None, tmpdir=_DEFAULT_TMPDIR):
    '''A simple interface for CacheNumpy object.
    If a cache already exists it reads it,
    otherwise is evaluates the callable_ (with no arguments)
    and expects it to return data is a form that can be written to disk. 
    If a filelist is given, compares the file modification times to the cache file.
    If he files have been modified since the cache was written, then the cache is overriden. 
    '''
    cn = CachePickle(uniquestr, tmpdir=tmpdir)
    if cn.exists() and (filelist is None or cn.newerthan(*filelist)):
        data = cn.read()
    else:
        data = callable_()
        cn.write(data)
    return data    

def cache_numpy(uniquestr, callable_, filelist=None, tmpdir=_DEFAULT_TMPDIR):
    '''As cache but for numpy arrays. 
    '''
    cn = CacheNumpy(uniquestr)
    if cn.exists() and (filelist is None or cn.newerthan(*filelist)):
        data = cn.read()
    else:
        data = callable_()
        cn.write(data)
    return data

def cache_root(uniquestr, callable_, filelist=None, tmpdir=_DEFAULT_TMPDIR):
    cn = CacheRoot(uniquestr)
    if cn.exists() and (filelist is None or cn.newerthan(*filelist)):
        data = cn.read()
    else:
        data = callable_()
        cn.write(data)
    return data

###############################################################################

class Cache(object):
    def __init__(self, uniquestr, prefix="tmp", postfix=".bin", tmpdir=_DEFAULT_TMPDIR):
        self._uniquestr = uniquestr
        self._tmpdir = tmpdir
        self._prefix = prefix
        self._postfix = postfix
        
    def tmpfilename(self):
        uniquestr = self._uniquestr
        hashstr = hashlib.md5(uniquestr).hexdigest()
        tmpdir = self._tmpdir
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
        return "".join([tmpdir,
                        os.sep,
                        "_".join((self._prefix, hashstr)),
                        self._postfix,
                        ])
    
    def exists(self):
        return os.path.exists(self.tmpfilename()) 
    
    def newerthan(self, *filelist):
        result = False
        if self.exists():
            result = True
            cachetime = os.path.getmtime(self.tmpfilename())
            for fname in filelist:
                if os.path.getmtime(fname) > cachetime:
                    #recently modified 
                    result = False
        return result

###############################################################################

class CacheNumpy(Cache):
    def __init__(self, uniquestr, prefix="tmp", tmpdir=_DEFAULT_TMPDIR):
        super(CacheNumpy, self).__init__(uniquestr=uniquestr, prefix=prefix, tmpdir=tmpdir)
    
    def read(self):
        fname = self.tmpfilename()
        infile = file(fname, "rb")
        return numpy.load(infile)
    
    def write(self, data):
        fname = self.tmpfilename()
        outfile = file(fname, "wb")
        numpy.save(outfile, data)
        return

###############################################################################

class CacheRoot(Cache):
    def __init__(self, uniquestr, prefix="tmp", tmpdir=_DEFAULT_TMPDIR, usepickle=True):
        super(CacheRoot, self).__init__(uniquestr=uniquestr, prefix=prefix, tmpdir=tmpdir, postfix=".root")
        self._usepickle = usepickle
    
    def read(self):
        fname = self.tmpfilename()
        f = ROOT.TFile(fname, "read")
        data = {}
        for key in f.GetListOfKeys():
            print "DEBUG", key.GetName()
            obj = key.ReadObj()
            data[key.GetName()] = obj
        #unpickle
        result = []
        for key in sorted(data.keys()):
            d = data[key]
            if "picklestring_" in key:
                d = pickle.loads(str(d))
            result.append(d)
        return result
    
    def write(self, data):
        try:
            iterable = iter(data)
        except TypeError:
            iterable = [data]
        fname = self.tmpfilename()
        outfile = ROOT.TFile(fname, "recreate")
        for num, obj in enumerate(iterable):
            name = "object_" + str(num)
            try:
                outfile.WriteObject(obj, name)
            except TypeError:
                if not self._usepickle:
                    raise
                name = "picklestring_" + str(num)
                #must be a python object, pickle the 
                string = pickle.dumps(obj, protocol=0)
                cstring = ROOT.TObjString(string)
                outfile.WriteObject(cstring, name)
        outfile.Close()
        return

###############################################################################

class CachePickle(Cache):
    def __init__(self, uniquestr, prefix="tmp", tmpdir=_DEFAULT_TMPDIR):
        super(CachePickle, self).__init__(uniquestr=uniquestr, prefix=prefix, tmpdir=tmpdir)
    
    def read(self):
        fname = self.tmpfilename()
        with open(fname, "rb") as infile:
            data = pickle.load(infile)
        return data
    
    def write(self, data):
        fname = self.tmpfilename()
        with open(fname, "wb") as outfile:
            pickle.dump(data, outfile, protocol=pickle.HIGHEST_PROTOCOL)
        return

###############################################################################

def _test_numpy_cache():
    test_dict = { "A":1, "B":2, "C":3 }
    test_numpy = numpy.ones(shape=(3, 2))
    def func_dict():
        return test_dict
    def func_numpy():
        return test_numpy
    assert cache("_test_simplot_cache_1", func_dict) == test_dict
    assert numpy.array_equal( cache_numpy("_test_simplot_cache_2", func_numpy), test_numpy)
    print "_test_numpy_cache success!"
    return

###############################################################################

def main():
    _test_numpy_cache()
    _test_numpy_cache()
    return

if __name__ == "__main__":
    main()
