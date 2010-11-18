import touchstone, citi

def main():
    import sys
    from optparse import OptionParser

    usage = "usage: %prog [options] <infile> <outfile (without extension)>"
    parser = OptionParser(usage=usage)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="be verbose")
    parser.add_option("-f", "--format",
                      action="store", type="string", dest="format",
                      metavar="FORMAT", default="tstone",
                      help="set output format (tstone or citi)")
    parser.add_option("-r", "--recombine", 
                      action="store", type="string", dest="recombine",
                      metavar="PORTLIST", help="recombine ports")

    (options, args) = parser.parse_args()

    if len(args) != 2:
        print("ERROR: Incorrect number of arguments supplied")
        parser.print_help()
        sys.exit(1)

    for format in (touchstone, citi):
        try:
            nport = format.read(args[0], options.verbose)
            break
        except:
            nport = None
    
    if nport is None:
        print("ERROR: Unsupported input file format")
        sys.exit(1)

    if options.recombine:
        nport = nport.recombine(eval("[" + options.recombine + "]"))

    outformats = {'tstone': touchstone, 'citi': citi}
    outformat = outformats[options.format.lower()]
    outformat.write(nport, args[1])

