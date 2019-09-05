#!/Users/anselmo/anaconda3/bin/python3 python3

from optparse import OptionParser
import os
import sys
import inspect
import time
import urllib
import getpass
import csv
import datetime
import re

from modules.getEntrezHgncByEnsg import getEntrezHgncByEnsg



def reannotateGtfWithEntrezHgnc(gtf,outfile):

    #print(gtf)
    #print(outfile)

    with open(gtf, 'r+') as fin, open(outfile, 'w') as fout:
        geneid_current = " "
        geneid_returned_current = " "
        allrows = csv.reader(fin,delimiter="\t")
        for row in allrows:
            if row[0].startswith('#'):                #print out header info
                fout.write("/t".join(row)+"\n")
                continue
            else:                                     #start replacing EnsemblID with EnsemblID:EntrezID:Hgnc:HgncID:HgncSymbol
                infolist = row[8].split(";")
                #print(infolist)
                geneid_index = 0
                for i, j in enumerate(infolist):
                    if j.startswith('gene_id') == True:
                        geneid_index = i
                        #print(geneid_index) 
                geneid = infolist[geneid_index].replace('gene_id','').replace('"','').replace(' ','')  #geneid refers to EnsemblID
                #print(geneid)

                geneid_returned = " "                 #geneid_returned will be: EnsemblID:EntrezID:Hgnc:HgncID:HgncSymbol
                if geneid == geneid_current:
                    geneid_returned = geneid_returned_current
                else:
                    geneid_current = geneid
                    geneid_returned = getEntrezHgncByEnsg(geneid,skipHgncApi=False)
                    time.sleep(1.0/15.0)                            #requests only at 15 times per second
                    geneid_returned_current = geneid_returned
                    #print(geneid_returned)

                p = re.compile(geneid)
                #print(row)
                row_reconstituted = "\t".join(row)
                row_replaced = p.sub(geneid_returned,row_reconstituted)
                #print(row_replaced)
                fout.write(row_replaced+"\n")
 
        fin.close()
        fout.close()
     
    return None

def main():

    #############################################################################
    ########################### Command Line Flags ##############################
    parser = OptionParser()
    parser.add_option("-i", "--gtf", dest="gtf", default="none",
                      help="gtf file to reannotate with HGNC symbol, id and Entrez gene id.")
    parser.add_option("-o", "--outfile", dest="outfile", default="./gtf.reannotated.out",
                      help="reannotated gtf file with full file path and name. ex. /home/GRCh38.p12.reannot.gtf")
    parser.add_option("-v", "--verbose", action="store_true", dest="verb", default=False,
                      help="Option to report all output. Useful for debugging.  TRUE when set. DEFAULT=FALSE when not set.")

    (options, args) = parser.parse_args()

    usage = "usage: %prog [options]"

    ###############################################################################
    gtf = options.gtf
    outfile = options.outfile
    verb = options.verb

    reannotateGtfWithEntrezHgnc(gtf,outfile)

if __name__ == '__main__':
    main()

