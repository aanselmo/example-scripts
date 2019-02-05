#!/Users/aanselmo/anaconda3/bin/python

import httplib2 as http
import json
import sys
import os

try:
    from urlparse import urlparse
except ImportError:
    from urllib.parse import urlparse

headers = {
 'Accept': 'application/json',
}

uri = 'http://rest.genenames.org'

def queryHgncApi(id_type,id,verbose):
    path = " "
    if id_type == "entrez_id":
        path = '/fetch/entrez_id/'+id
    elif id_type == "hgnc_id":
        path = '/fetch/hgnc_id/'+id
    elif id_type == "ensembl_gene_id":
        path = '/fetch/ensembl_gene_id/'+id
    else:
        print("Query must be EnsGeneID or HGNC_ID or Entrez_ID ..")
        print('First argment is either "ensembl_gene_id" or "entrez_id" or "hgnc_id" ..')
        sys.exit()

    target = urlparse(uri+path)
    method = 'GET'
    body = ''

    h = http.Http()

    response, content = h.request(
        target.geturl(),
        method,
        body,
        headers)

    #print('Response status: ',response['status'])

    entrez_id = "."
    hgnc_id = "HGNC:."
    hgnc_symbol = "."

    if response['status'] == '200':
     # assume that content is a json reply
     # parse content with the json module 
        data = json.loads(content)
        #print(data)

        try:
            entrez_id = data['response']['docs'][0]['entrez_id']
        except Exception as e:
            if verbose == True:
                print('Caught entrez_id exception: ',e)
        try:
            hgnc_id = data['response']['docs'][0]['hgnc_id']
        except Exception as e:
            if verbose == True:
                print('Caught hgnc_id exception: ',e)
        try:
            hgnc_symbol = data['response']['docs'][0]['symbol']
        except Exception as e:
            if verbose == True:
                print('Caught hgnc_symbol exception : ',e)
    else:
        print('Error detected: ' + response['status'])

    return entrez_id,hgnc_id,hgnc_symbol


def main():
    ## sys.argv[1] should be either hgnc_id or entrez_id (number only)
    entrez_id, hgnc_id, hgnc_symbol = queryHgncApi(sys.argv[1],sys.argv[2],sys.argv[3])

    print(entrez_id+":"+hgnc_id+":"+hgnc_symbol)
    #print(hgnc_symbol)

if __name__ == '__main__':
    main()    
