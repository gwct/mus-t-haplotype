import sys, os, argparse

print()
print("###### Build site pages ######");
print("PYTHON VERSION: " + ".".join(map(str, sys.version_info[:3])))
print("# Script call: " + " ".join(sys.argv) + "\n----------");

parser = argparse.ArgumentParser(description="Gets stats from a bunch of abyss assemblies.");
parser.add_argument("--all", dest="all", help="Build all pages", action="store_true", default=False);
parser.add_argument("--index", dest="index", help="Without --all: build index.html. With --all: exlude index.html", action="store_true", default=False);
# parser.add_argument("--assemblies", dest="assemblies", help="Without --all: build assemblies.html. With --all: exlude assemblies.html", action="store_true", default=False);
# parser.add_argument("--trees", dest="trees", help="Without --all: build trees.html. With --all: exlude trees.html", action="store_true", default=False);
# parser.add_argument("--annotations", dest="annotations", help="Without --all: build annotations.html. With --all: exlude annotations.html", action="store_true", default=False);
parser.add_argument("--cactus", dest="cactus", help="Without --all: build cactus.html. With --all: exlude cactus.html", action="store_true", default=False);
# parser.add_argument("--cactustests", dest="cactustests", help="Without --all: build cactus_tests.html. With --all: exlude cactus_tests.html", action="store_true", default=False);
parser.add_argument("--windows", dest="windows", help="Without --all: build windows.html. With --all: exlude windows.html", action="store_true", default=False);
# parser.add_argument("--analyses", dest="analyses", help="Without --all: build analyses.html. With --all: exlude analyses.html", action="store_true", default=False);
# parser.add_argument("--people", dest="people", help="Without --all: build people.html. With --all: exlude people.html", action="store_true", default=False);
# parser.add_argument("--links", dest="links", help="Without --all: build links.html. With --all: exlude links.html", action="store_true", default=False);
args = parser.parse_args();
# Input options.

#cwd = os.getcwd();
os.chdir("generators");

pages = {
    'index' : args.index,
    # 'assemblies' : args.assemblies,
    # 'trees' : args.trees,
    # 'annotations' : args.annotations,
    'cactus' : args.cactus,
    # 'cactustests' : args.cactustests,
    'windows' : args.windows,
    # 'analyses' : args.analyses,
    # 'people' : args.people,
    # 'links' : args.links
}

if args.all:
    pages = { page : False if pages[page] == True else True for page in pages };

if pages['index']:
    os.system("python index_generator.py");

# if pages['annotations']:
#     os.system("Rscript annotations_generator.r");

# if pages['assemblies']:
#     os.system("Rscript assemblies_generator.r");

# if pages['trees']:
#     os.system("Rscript trees_generator.r");

if pages['cactus']:
    os.system("Rscript cactus_generator.r");

# if pages['cactustests']:
#     os.system("Rscript cactus_tests_generator.r");

if pages['windows']:
    os.system("Rscript windows_generator.r");

# if pages['analyses']:
#     os.system("Rscript analyses_generator.r");

# if pages['people']:
#     os.system("python people_generator.py");

# if pages['links']:
#     os.system("python links_generator.py");
    
print("----------\nDone!");


