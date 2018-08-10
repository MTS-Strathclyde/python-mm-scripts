#!/usr/bin/python2



import argparse
import shutil
import os
import sys
import subprocess



runtxt = """for i in {pdbs}
do
    name=${{i%.pdb}}
    rism3d_pc.py $name.pdb -p $name.prmtop -x {xvv} {args}
done
"""


def process_command_line(argv):
    """Processes arguments

    Parameters
    ----------
    argv : list
        Command line arguments.

    Returns
    -------
    out : argparse.Namespace
        Namespace of command line arguments.
    """
    parser = argparse.ArgumentParser(description=""" Create rism3d_pc run directory.""")
    #Positional args
    parser.add_argument('files', metavar='molec.p*',
                        help="""Input files. Should be equal number of pdbs and prmtops""",
                        nargs='+')
    #Optional args
    parser.add_argument('-x', '--xvv',
                        help="""rism3d_pc.py xvv file relative to pwd""")
    parser.add_argument('-a', '--args',
                        help="""rism3d_pc.py extra args surrounded by quotes""",
                        default='')
    parser.add_argument('-n', '--name',
                        help="""Name of the directory [calc]""", default='calc')
    parser.add_argument('--nohup',
                        help="""Put job under nohup.""", action='store_true')
    return parser.parse_args(argv)


def main(argv):
    args = process_command_line(argv)
    try:
        os.mkdir(args.name)
    except OSError:
        pass
    pdbs = []
    for f in args.files:
        if f.endswith('.pdb'):
            pdbs.append(os.path.split(f)[1])
        shutil.copy(f, args.name)
    pdb_names = ' '.join(pdbs)
    orig_dir = os.getcwd()
    os.chdir(args.name)
    xvv_path = os.path.relpath(args.xvv, args.name)
    with open('run.sh', 'w') as f:
        f.write(runtxt.format(pdbs=pdb_names, xvv=xvv_path, args=args.args))
    if args.nohup:
        subprocess.Popen(['nohup', 'bash', 'run.sh', '&'])
    else:
        subprocess.call(['bash', 'run.sh'])
    os.chdir(orig_dir)

if __name__ == '__main__':
    main(sys.argv[1:])

