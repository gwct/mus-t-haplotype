import sys, os, argparse

print()
print("###### Build site pages ######");
print("PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])))
print("# Script call: " + " ".join(sys.argv) + "\n----------");

parser = argparse.ArgumentParser(description="Gets stats from a bunch of abyss assemblies.");
parser.add_argument("--all", dest="all", help="Build all pages", action="store_true", default=False);
parser.add_argument("--index", dest="index", help="Without --all: build index.html. With --all: exlude index.html", action="store_true", default=False);
parser.add_argument("--cactus", dest="cactus", help="Without --all: build cactus.html. With --all: exlude cactus.html", action="store_true", default=False);
parser.add_argument("--windows", dest="windows", help="Without --all: build windows.html. With --all: exlude windows.html", action="store_true", default=False);
parser.add_argument("--trees", dest="trees", help="Without --all: build trees.html. With --all: exlude trees.html", action="store_true", default=False);
parser.add_argument("--inversions", dest="inversions", help="Without --all: build inversions.html. With --all: exlude inversions.html", action="store_true", default=False);
args = parser.parse_args();
# Input options.

#cwd = os.getcwd();
os.chdir("generators");

pages = {
    'index' : args.index,
    'cactus' : args.cactus,
    'windows' : args.windows,
    'trees' : args.trees,
    'inversions' : args.inversions
}

if args.all:
    pages = { page : False if pages[page] == True else True for page in pages };

if pages['index']:
    os.system("python index_generator.py");

if pages['cactus']:
    os.system("Rscript cactus_generator.r");

if pages['windows']:
    os.system("Rscript windows_generator.r");

if pages['trees']:
    os.system("Rscript trees_generator.r");

if pages['inversions']:
    os.system("Rscript inversions_generator.r");

print("----------\nDone!");


