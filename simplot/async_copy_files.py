from time import sleep
from multiprocessing import Queue, Process
import shutil

class ListAsyncOperation:
    """Asynchronously applied "operation" to the elements of a list "l".
       Iterating over the class to get the result. If the next "operation" is
       not complete it will block until the next result is available.

       This can be used, for example, to perform IO and then process the result
       while IO is run for the next process. See `async_copy_files` for an
       example.
    """
    def __init__(self, l, operation, timeout=None):
        self._q = Queue()
        self._timeout = timeout
        self._l = list(l)
        self._thread = Process(target=self._apply, args=(operation, self._l, self._q,))
        self._thread.start()

    def __iter__(self):
        n = 0
        while n < len(self._l):
            r = self._q.get(timeout=self._timeout)
            yield r
            n += 1

    def _apply(self, f, l, q):
        for x in l:
            q.put(f(x))

def async_copy_files(flist, timeout=None):
    """Asynchronously copies files from one location to another.
    Takes as a list of pairs of files [(source, destination)...].
    Returns an iterable that iterates over (source, destiation). The
    iterable will block when next() is called if the corresponding files have
    not finished copying.
    """
    def operation(f):
        src, dst = f
        if (not os.path.exists(dst)) or os.path.getmtime(src) > os.path.getmtime(dst):
            shutil.copy(src, dst)
        return src, dst
    return ListAsyncOperation(flist, operation, timeout=timeout)

def _test_list_async_operation():
    def operation(x):
        print "AsyncListOperation", str(x)
        sleep(2)
        print "AsyncListOperation", str(x), "done."
        return x
    def process(f):
        print "Process", f
        sleep(3)
        print "Process", f, "complete"
    for x in ListAsyncOperation(range(3), operation):
        process(x)
    return

def _test():
    for x in async_copy_files(range(3)):
        process(x)

def _main():
    _test_list_async_operation()
    return

if __name__ == "__main__":
    _main()
