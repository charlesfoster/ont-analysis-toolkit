#!/usr/bin/env python
import sys

def main():
    if sys.argv[1] == 'gui':
        import oat.gui.gui as gui
        gui.main()
    else:
        import oat.cli.cli as cli
        cli.main(sys.argv[1:])
        
if __name__ == "__main__":
    main()
