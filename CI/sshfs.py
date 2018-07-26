from fs.sshfs import SSHFS as SSHFSBase
from fs.base import FS
from fs import errors
import paramiko
import socket

class SSHFS(SSHFSBase):
    def __init__(self, host, user, passwd, port=22, compress=False, timeout=10, keepalive=10):
        FS.__init__(self)

        self._user = user
        self._host = host
        self._port = port
        self._client = client = paramiko.SSHClient()
        self._locale = None

        try:
            client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

            client.connect(
                socket.gethostbyname(host), port, user, passwd,
                compress=compress, timeout=timeout,
                allow_agent=False,look_for_keys=False
            )
            if keepalive > 0:
                client.get_transport().set_keepalive(keepalive)
            self._sftp = client.open_sftp()
            self._platform = None

        except (paramiko.ssh_exception.SSHException,            # protocol errors
                paramiko.ssh_exception.NoValidConnectionsError, # connexion errors
                socket.gaierror, socket.timeout) as e:          # TCP errors

            message = "Unable to create filesystem: {}".format(e)
            raise errors.CreateFailed(message)
