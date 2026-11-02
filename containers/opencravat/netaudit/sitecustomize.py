# Network audit for the OPENCRAVAT task of noncodingdrivers (containers/opencravat).
#
# Python imports this file at start-up when its folder is on PYTHONPATH (the module sets that
# for the `oc run` call only). With OC_NETAUDIT_LOG set, every IP socket connect() and every
# getaddrinfo() lookup made by OpenCRAVAT or one of its child processes is appended to that
# file as one tab-separated line: pid, call, host, port. Unix-domain sockets (the local
# multiprocessing listener, thousands of calls per run) are not network traffic and are not
# recorded. Nothing else is changed: the calls go ahead as usual, so the audit can only record,
# never alter, a result.
import os
import socket

_log = os.environ.get("OC_NETAUDIT_LOG")

if _log:
    _connect = socket.socket.connect
    _connect_ex = socket.socket.connect_ex
    _getaddrinfo = socket.getaddrinfo

    def _record(call, address):
        if isinstance(address, tuple) and len(address) >= 2:
            host, port = address[0], address[1]
        else:
            host, port = address, ""
        try:
            with open(_log, "a") as fh:
                fh.write(f"{os.getpid()}\t{call}\t{host}\t{port}\n")
        except OSError:
            pass

    _unix = getattr(socket, "AF_UNIX", None)

    def connect(self, address):
        if self.family != _unix:
            _record("connect", address)
        return _connect(self, address)

    def connect_ex(self, address):
        if self.family != _unix:
            _record("connect_ex", address)
        return _connect_ex(self, address)

    def getaddrinfo(host, port, *args, **kwargs):
        _record("getaddrinfo", (host, port))
        return _getaddrinfo(host, port, *args, **kwargs)

    socket.socket.connect = connect
    socket.socket.connect_ex = connect_ex
    socket.getaddrinfo = getaddrinfo
