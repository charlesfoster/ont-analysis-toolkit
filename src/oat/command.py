#!/usr/bin/env python
import sys
from oat import __version__

def main():
    if len(sys.argv) == 1:
        print(
            """\n\033[95m
                                            ,d
                                            88
                   ,adPPYba,   ,adPPYYba, MM88MMM
                  a8"     "8a ""     `Y8    88
                  8b       d8 ,adPPPPP88    88
                  "8a,   ,a8" 88,    ,88    88,
                   `"YbbdP"'  `"8bbdP"Y8   "Y888

            OAT: ONT Analysis Toolkit (version {})\033[0m
        """.format(
                __version__
            )
        )

        print("\n'oat' can be deployed using either a graphical user interface (GUI) or using a command line interface (CLI)")
        print("\n\t* To run the GUI: oat gui")
        print("\n\t* To run the CLI: oat [options] input_spreadsheet.csv")
        print("\n\t* To see the CLI options: oat --help\n")
        sys.exit(1)
    if sys.argv[1] == 'gui':
        import oat.gui.gui as gui
        gui.main()
    else:
        import oat.cli.cli as cli
        cli.main(sys.argv[1:])

if __name__ == "__main__":
    main()
