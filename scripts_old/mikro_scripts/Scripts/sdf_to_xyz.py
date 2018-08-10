import pybel
import sys


def main():
    f = sys.argv[1]

    mols = pybel.readfile("sdf", f)

    i = 0
    for mol in mols:
        mol.write("xyz", "mol" + str(i) + ".xyz")
        i += 1

if __name__ == "__main__":
    main()
