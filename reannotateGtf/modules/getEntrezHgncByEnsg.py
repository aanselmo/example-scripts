#!/usr/bin/env python3

from optparse import OptionParser
import os
import sys
import inspect
import time
import urllib
import getpass
import csv
import datetime
import requests

from modules.queryHgncApi import queryHgncApi

### Query Ensembl API with EnsemblGeneID 
###  and return Concatenated String containing EnsemblGeneID:EntrezID:HGNC_ID:HGNC_Symbol

def getEntrezHgncByEnsg(ensgID,skipHgncApi):

    # Defaults if information not available in Ensembl API
    entrez_id = "."
    hgnc_id = "HGNC:."
    hgnc_symbol = "."
    fillInWithHgncApi = True

    # Query arguments    
    server = "http://rest.ensembl.org"
    ext_entrez = "/xrefs/id/"+ensgID+"?external_db=EntrezGene"
    ext_hgnc = "/xrefs/id/"+ensgID+"?external_db=HGNC"

    ## Query the Ensembl API
    # Entrez ID
    r_entrez = requests.get(server+ext_entrez, headers={ "Content-Type" : "application/json"})
    if not r_entrez.ok:
        r_entrez.raise_for_status()
        sys.exit()

    # HGNC ID and Symbol
    r_hgnc = requests.get(server+ext_hgnc, headers={ "Content-Type" : "application/json"})
    if not r_hgnc.ok:
        r_hgnc.raise_for_status()
        sys.exit()

    # Parse Returned Entrez entry 
    entrez_all = r_entrez.json()
    #print(entrez_all)
    if len(entrez_all) == 0:
        entrez_id = "."
    else:    
        entrez_id = entrez_all[0]['primary_id']

    # Parse Returned HGNC entry
    hgnc_all = r_hgnc.json()
    #print(hgnc_all)
    if len(hgnc_all) == 0:
        hgnc_id = "HGNC:."
        hgnc_symbol = "."
    else:
        hgnc_id = hgnc_all[0]['primary_id']
        hgnc_symbol = hgnc_all[0]['display_id']
   
    if skipHgncApi == True:
        fillInWithHgncApi = False
    if fillInWithHgncApi == True:
        #### Query failures against HGNC API
        ## if no entrez or hgnc data available from Ensembl API
        ## use data from HGNC API
        entrez_id_sub = "."
        hgnc_id_sub = "HGNC:."
        hgnc_symbol_sub = "."

        if entrez_id == "." and hgnc_id == "HGNC:.":
            entrez_id_sub, hgnc_id_sub, hgnc_symbol_sub = queryHgncApi("ensembl_gene_id",ensgID,verbose=False)
            time.sleep(1.0/15.0)
        elif entrez_id == "." and hgnc_id != "HGNC:.":
            entrez_id_sub, hgnc_id_sub, hgnc_symbol_sub = queryHgncApi("hgnc_id",hgnc_id,verbose=False)
            time.sleep(1.0/15.0)
        elif hgnc_id != "." and hgnc_id == "HGNC:.":
            entrez_id_sub, hgnc_id_sub, hgnc_symbol_sub = queryHgncApi("entrez_id",entrez_id,verbose=False)
            time.sleep(1.0/15.0)

        # Use HGNC API values if none available from EnsemblAPI
        if entrez_id == ".":
            entrez_id = entrez_id_sub
        if hgnc_id == "HGNC:.":
            hgnc_id = hgnc_id_sub
        if hgnc_symbol == ".":
            hgnc_symbol = hgnc_symbol_sub 
        
    ens_entrezid_hgnc = ensgID+":"+entrez_id+":"+hgnc_id+":"+hgnc_symbol

    #print(ens_entrezid_hgnc)

    return ens_entrezid_hgnc


def main():
    #############################################################################
    ########################### Command Line Flags ##############################
    parser = OptionParser()
    parser.add_option("-i", "--ensemblGeneID", dest="ensgID", default="none",
                      help="Ensembl gene id to search.")
    parser.add_option("-f", "--skipHgncApi", action="store_true", dest="skipHgncApi", default=False,
                      help="Option to skip using HGNC API to fill in what EnsemblAPI misses. TRUE when set. DEFAULT=False when not set.")
    parser.add_option("-v", "--verbose", action="store_true", dest="verb", default=False,
                      help="Option to report all output. Useful for debugging.  TRUE when set. DEFAULT=FALSE when not set.")

    (options, args) = parser.parse_args()

    usage = "usage: %prog [options]"

    ###############################################################################
    ensgID = options.ensgID                # ensembl gene id  example: ENSG0000000000001
    verb = options.verb
    skipHgncApi = options.skipHgncApi

    getEntrezHgncByEnsg(ensgID,skipHgncApi)

if __name__ == '__main__':
    main()
