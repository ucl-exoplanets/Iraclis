"""Usage: pipeline.py -p PARFILE
          pipeline.py -d DATADIR
          pipeline.py -D DATADIR [PROCEDURE] [PARSTRING]
          pipeline.py -T

"""


def console():
    from iraclis.pipeline import process_visit, run_test
    import docopt

    arguments = docopt.docopt(__doc__)

    if arguments['-p'] or arguments['-d'] or arguments['-D']:
        process_visit(arguments['PARFILE'], arguments['DATADIR'], arguments['PROCEDURE'], arguments['PARSTRING'])

    # for developers use only

    elif arguments['-T']:
        run_test()


if __name__ == '__main__':
    console()
