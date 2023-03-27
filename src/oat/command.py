#!/usr/bin/env python3
import sys

def main():
    if sys.argv[1] == 'gui':
        import oat.gui.gui as gui
        gui.main()
    # elif sys.argv[1] == 'cli':
    else:
        import oat.cli.cli as cli
        cli.main(sys.argv[1:])
        
if __name__ == "__main__":
    for dist in __import__('pkg_resources').working_set:
        print(dist.project_name.replace('Python', ''))
    main()
