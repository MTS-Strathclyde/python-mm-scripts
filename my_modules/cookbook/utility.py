# -*- coding: utf-8 -*-
"""
Created on Fri May 16 09:49:13 2014

@author: max
"""

def is_float(number):
    try:
        float(number)
        return True
    except ValueError:
        return False
        
import os

def try_to_make_dir(dirname):
    """Tries to make directory, passes in case it exists."""
    try:
        os.mkdir(dirname)
    except OSError, e:
        if e.errno == 17:
            pass # Directory already exists, all is well
        else:
            raise e


import subprocess
import threading
import signal

class RunCmd(threading.Thread):
    """ Will only work on Linux. And sometimes it will not work even on Linux. """
    def __init__(self, cmd, timeout, cwd='.'):
        """ Run subprocess for a fixed ammount of time and then kill it.
        
        Parameters
        ----------
        cmd : list
            Command to execute.
        
        timeout : float
            Time in minutes after which process will be killed
            
        cwd : string, default .
            Directory in which process will be run
        """
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout*60
        self.cwd = '.'

    def run(self):
        self.p = subprocess.Popen(self.cmd, preexec_fn=os.setsid,
                           stdout=subprocess.PIPE,
                           cwd=self.cwd)
        self.p.wait()
    
    def run_and_timeout(self):
        """Run subprocess
        
        Returns
        -------
        out : tuple
            Tuple has stout and stderr outputs        
        """
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            os.killpg(self.p.pid, signal.SIGTERM)
            self.join()
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print self.cmd
            print 'Has run out of time'
            print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        return self.p.communicate()


#import paramiko
#
#class Connection:
#    """Used to manage connection and work with sipsik
#    (or basically, with some limitation, any other remote server, that
#    supports ssh).
#    """
#    def __init__(self, username='xpb13212', password='',
#                 hostname='node004'):
#        """Set variables."""
#        self.username = username
#        self.password = password
#        self.hostname = hostname
#        self.ssh = None
#
#    def connect(self):
#        """Connection.
#        """
#        self.ssh = paramiko.SSHClient()
#        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
#        #self.ssh.load_system_host_keys()
#        self.ssh.connect(self.hostname, username=self.username,
#                         password=self.password)
#
#    def disconnect(self):
#        self.ssh.close()
#
#    def execute_command(self, command, verbose=True):
#        """Executes command on connected server and prints Stderr. In verbose
#        mode prints Stdout as well.
#        Commands should be given as string.
#        """
#        _, stdout, stderr = self.ssh.exec_command(command)
#        stdout = stdout.read()
#        stderr = stderr.read()
#        if stderr:
#            print stderr
#        if verbose:
#            if stdout:
#                print stdout


def list_ext(ext, folder='.'):
    """Returns absolute pathes of files, which have given extension.
    By default looks in working directory.
    Returns them as list.
    """
    if folder[-1] == '/':
        folder = folder[:-1]
    folder_path = os.path.join(os.getcwd(), folder)
    all_files = os.listdir(folder_path)
    ext_list = []
    for filename in all_files:
        if os.path.splitext(filename)[1] == ext:
            ext_list.append(os.path.join(folder_path, filename))
    return ext_list            
            
def change_ext(path, new_ext):
    """Accepts filenames as well. Extension should be supplied with dot.
    Example:
        >>> change_extension('f.ext', '.o')
        f.o
    """
    return os.path.splitext(path)[0] + new_ext
    
def get_filename(path):
    """Retuns filename, remove path and extension."""
    with_ext = os.path.split(path)[1]
    return os.path.splitext(with_ext)[0]


def remove_path(path):
    """Returns filename without path"""
    return os.path.split(path)[1]


def remove_extenison(path):
    """Returns path without extension"""
    return os.path.splitext(path)[0]


def add_to_name_but_keep_ext(path, str_to_add):
    """Returns path with appended string between end of old file name and 
    extension."""
    name, ext = os.path.splitext(path)
    return name + str_to_add + ext
    
def file_exist(path):
    try:
        with open(path) as f:
            return True
    except IOError:
       return False


