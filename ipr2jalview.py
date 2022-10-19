#!/bin/env python

import collections
import requests
import json
import sys
import xml.etree.cElementTree as ET
from argparse import ArgumentParser
from requests.adapters import HTTPAdapter, Retry

ipr_ws='https://www.ebi.ac.uk/interpro/api/entry/all/protein/Uniprot/'
up_ws='https://rest.uniprot.org/uniprotkb/'

# The colour palette is colour-blind friendly, selected from the microshades cvd palette
# see https://github.com/KarstensLab/microshades
palette=[ "4E7705", "9D654C", "098BD9", "148F77", "7D3560", "6D9F06", 
		  "C17754", "F09163", "56B4E9", "009E73", "A1527F", "97CE2F", 
		  "FCB076", "7DCCFF", "43BA8F", "CC79A7", "BDEC6F", "FFD5AF", 
		  "BCE1FF", "48C9B0", "EFB6D6" ]

def get_ipr(accession):

	"""
	Retrieves interpro features for a specified uniprot accession

	Required parameters:
		accession(str): uniprot accession

	Returns:
		ipr_json(dict): dict of returned data
	"""

	s=requests.Session()
	retries=Retry(total=5, backoff_factor=1,status_forcelist=[500,502,503])
	s.mount('https://',HTTPAdapter(max_retries=retries))

	url="{}/{}".format(ipr_ws,accession)
	r=requests.get(url)
	if r.status_code==204:
		print('Accession {} not found in InterPro'.format(accession))
		sys.exit(1)

	ipr_json=json.loads(r.text)
	
	return(ipr_json)

def get_uniprot_id(accession):

	"""
	Finds uniprot name for for accession through uniprot API

	required parameters:
		accession(str): uniprot accession

	returns:
		name(str): uniprot name
	"""

	s=requests.Session()
	retries=Retry(total=5, backoff_factor=1,status_forcelist=[500,502,503])
	s.mount('https://',HTTPAdapter(max_retries=retries))

	url="{}/{}.xml".format(up_ws, accession)
	r=requests.get(url)
	if r.status_code==204:
		print('Accession {} not found in UniProt'.format(accession))
		sys.exit(1)

	tree=ET.ElementTree(ET.fromstring(r.text))
	root=tree.getroot()
	name=root.find('./{http://uniprot.org/uniprot}entry/{http://uniprot.org/uniprot}name').text

	return(name)

def format_features(ipr,include_list,unip_id):

	"""
	Generates (mostly) formatted lines appropriate for producing jalview features file
	These will be lacking the colour mapping since this is only generated one we know
	what all the feature types are.

	Required parameters:
		ipr(dict): parsed data for individual interpro record
		include_list(list): list of interpro entries to include
		unip_id(str): Uniprot entry name, which will be substituted for the accession if provided

	Returns:
		ipr_features(list): list of strings representing each feature
		ipr_type(str): type of feature from metadata->type in parsed ipr data
	"""

	ipr_acc=ipr['metadata']['accession']
	ipr_name=ipr['metadata']['name']
	ipr_type=ipr['metadata']['type']

	ipr_features=list()

	if (include_list and ipr_acc in include_list) or include_list==None:
		for protein in ipr['proteins']:
			
			if unip_id:
				prot_acc=unip_id
			else:
				prot_acc=protein['accession'].upper()

			for loc in protein['entry_protein_locations']:
				for fragment in loc['fragments']:
					start=fragment['start']
					end=fragment['end']
					ipr_features.append("{}\t{}\t-1\t{}\t{}\t{}".format( ipr_name, prot_acc, start, end, ipr_type))

	return(ipr_features, ipr_type)

def generate_output(formatted, outfile):

	"""
	Takes (mostly) formatted lines feature lines, picks colours for each feature type
	and generates appropriate output file

	Required params:
		formatted(list): list of strings representing each feature 
		outfile(str): Path of output file to create

	Returns:
		None
	"""

	feat_colours=dict()
	outlines=list()
	
	# allocate a colour to each feature type - these should be first in the output
	for ipr_type in formatted.keys():
		feat_col=palette.pop(0)
		feat_colours[ipr_type]=feat_col
			
		outlines.append("{}\t{}\n".format(ipr_type,feat_col))

	outlines.append('startgroup	interpro\n')
	# iterate through formatted lines appending appropriate colour...
	for ipr_type in formatted.keys():
		features=formatted[ipr_type]
		for group in features:
			for feature in group:
				feat_type=feature.split("\t")[-1]
				line="{}\t{}\n".format(feature,feat_colours[feat_type])
				outlines.append(line)
	outlines.append('endgroup	interpro\n')

	with open(outfile,'w') as out_handle:
		out_handle.writelines(outlines)

if __name__=='__main__':

	parser=ArgumentParser(description='Obtains InterPro annotations for defined accession and formats as Jalview features')
	parser.add_argument('-a', '--accession', dest='accession', help="Uniprot accession to request", type=str, required=True)
	parser.add_argument('-i', '--interpro',  dest='iprs', action="extend",nargs="+", type=str, help='InterPro annotation to include i.e. IPR00001. May be specified multiple times')
	parser.add_argument('-o', '--output', dest='output', type=str, help='Output file to create', required=True)
	parser.add_argument('-u', '--uniprot_id', dest='unip_id', default=False, action="store_true", help='Label features with uniprot name rather than accession')
	args=parser.parse_args()

	unip_id=None

	ipr=get_ipr(args.accession)

	# optionally output features mapped to uniprot name instead of accession
	if args.unip_id:
		unip_id=get_uniprot_id(args.accession)

	formatted=collections.defaultdict(list)

	for i in ipr['results']:
		f, ipr_type=format_features(i, args.iprs,unip_id)
		formatted[ipr_type].append(f)

	if len(formatted) > len(palette):
		print('Your going to need a bigger colour palette...')	
		sys.exit(1)

	generate_output(formatted, args.output)
