# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 17:11:03 2012

@author: a92549
"""

import paramiko


class Connection:
    """Used to manage connection and work with sipsik
    (or basically, with some limitation, any other remote server, that
    supports ssh).
    """
    def __init__(self, username='mishin1991', password='Max7120',
                 hostname='sipsik.chem.ut.ee'):
        """Set variables."""
        self.username = username
        self.password = password
        self.hostname = hostname
        self.ssh = None

    def connect(self):
        """Connection.
        """
        self.ssh = paramiko.SSHClient()
        self.ssh.load_system_host_keys()
        self.ssh.connect(self.hostname, username=self.username,
                         password=self.password)

    def disconnect(self):
        self.ssh.close()

    def execute_command(self, command, verbose=True):
        """Executes command on connected server and prints Stderr. In verbose
        mode prints Stdout as well.
        Commands should be given as string.
        """
        _, stdout, stderr = self.ssh.exec_command(command)
        stdout = stdout.read()
        stderr = stderr.read()
        if stderr:
            print stderr
        if verbose:
            if stdout:
                print stdout

    def change_path(self, local_path):
        """Changes given (absolulte) path into path, relative to mikro
        home directory.

        Is suitable to locate files on sipsik, as all servers have one
        storage per user."""
        local_home = local_path.split('/', 2)[2]
        return '/home2/' + local_home
