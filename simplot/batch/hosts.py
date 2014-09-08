import socket

class KnownHosts:
    #list of known hosts
    WARWICK_CLUSTER_LOGIN = "WARWICK_CLUSTER_LOGIN"
    WARWICK_CLUSTER_WORKER = "WARWICK_CLUSTER_WORKER"
    CSC_DESKTOP = "CSC_DESKTOP"
    OTHER = "OTHER"

def ishost(host):
    return host == _thishostid

def _gethostid(hostname):
    result = None
    if "epp-ui" in hostname:
        result = KnownHosts.WARWICK_CLUSTER_LOGIN
    elif "comp" in hostname:
        result = KnownHosts.WARWICK_CLUSTER_WORKER
    elif ".warwick.ac.uk" in hostname:
        result = KnownHosts.CSC_DESKTOP
    else:
        result = KnownHosts.OTHER
    return result

_thishostname = socket.getfqdn()
_thishostid = _gethostid(_thishostname)
