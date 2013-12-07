import unittest
import numpy

import simplot.binning as binning

class BinningTests(unittest.TestCase):
    def testComparisonOperator(self):
        #create objects
        b1 = binning.Binning(10, 0.0, 10.0)
        b1clone = binning.Binning(10, 0.0, 10.0)
        b1different = binning.Binning(10, 0.0, 11.0)
        #test cases
        self.assertEqual(b1, b1clone)
        self.assertNotEqual(b1, b1different)
        
    def testConstructor(self):
        #test construction with different methods gives the same result
        b1 = binning.Binning(10, 0.0, 10.0)
        b2 = binning.Binning(numpy.arange(0.0, 11.0, 1.0))
        b3 = binning.Binning([float(i) for i in xrange(0, 11)])
        self.assertEquals(b1, b2)
        self.assertEquals(b1, b3)
        self.assertEquals(b2, b3)
        
    def testVerifyNonSequential(self):
        self.assertRaises(ValueError, binning.Binning, ([0.0, 2.0, 1.0, 3.0]))
    
    def testIteration(self):
        b = binning.Binning(4, 0.0, 10.0)
        r = [0, 1, 2, 3]
        self.assertEquals(list(b), r)
        
    def testRanges(self):
        b = binning.Binning(4, 0.0, 10.0)
        r = [0, 1, 2, 3]
        #length
        self.assertEquals(len(b), len(r))
        #min and max
        self.assertEquals(min(b), min(r))
        self.assertEquals(max(b), max(r))
        
    def testFindBinNumber(self):
        b = binning.Binning(4, 0.0, 4.0)
        for x, n in [(-1.0, -1), # check underflow
                     (0.0, 0), # check low edge of range
                     (0.5, 0), # check middle of bin
                     (1.0, 1), # check up edge of a bin
                     (2.0001, 2), # check some random numbers
                     (3.9999, 3),
                     (4.0, 4), # check up edge of range
                     (5.0, 4), # check the overflow
                     ]:
            self.assertEquals(b.binnumber(x), n)
        return
    
    def testBinEdges(self):
        b = binning.Binning(4, 0.0, 4.0)
        for binnum, lowedge, upedge, bincentre in [#bins in range
                                                   (0, 0.0, 1.0, 0.5),
                                                   (1, 1.0, 2.0, 1.5),
                                                   (2, 2.0, 3.0, 2.5),
                                                   (3, 3.0, 4.0, 3.5),
                                                   #overflow bins
                                                   (-1, None, 0.0, None),
                                                   (4, 4.0, None, None),
                                                   ]:
            if lowedge is not None:
                self.assertEquals(b.binlow(binnum), lowedge)
            else:
                self.assertRaises(IndexError, b.binlow, [binnum])
            if upedge is not None:
                self.assertEquals(b.binhigh(binnum), upedge)
            else:
                self.assertRaises(IndexError, b.binhigh, [binnum])
            if bincentre is not None:
                r = b.bincentre(binnum)
                self.assertAlmostEqual(r, bincentre, delta=1.0e-6,
                                       msg="binnum={0}, expected bincentre={1} but received bincentre={2}".format(binnum, bincentre, r),
                                       )
            else:
                self.assertRaises(IndexError, b.bincentre, [binnum])
        return
    
    
        
        
    

def main():
    unittest.main()

if __name__ == "__main__":
    main()



