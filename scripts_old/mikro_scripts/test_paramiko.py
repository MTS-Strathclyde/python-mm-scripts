#!/usr/bin/python

import sys
import paramiko as pm
import os

HOST = 'sipsik.chem.ut.ee'
USER = 'mishin1991'
PASSWORD = 'Max7120'

client = pm.SSHClient()
client.load_system_host_keys()
client.connect(HOST, username=USER, password=PASSWORD)

channel = client.invoke_shell()
stdin = channel.makefile('wb')
stdout = channel.makefile('rb')

stdin.write('''
cd log
ls
exit
''')
print stdout.read()

stdout.close()
stdin.close()
client.close()
