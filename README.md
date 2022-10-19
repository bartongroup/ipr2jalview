# ipr2jalview

Reformats InterPro annotations into a Jalview format feature file 

A simple script which retrieves InterPro annotations for a protein sequence and formats them appropriately for loading into Jalview.

At present, the individual member databases hits are not reported, just the overall InterPro entry.

## Installation

Just download the script somewhere which is useful to you...

## Usage

Running `ipr2jalview.py -h` will provide basic usage information.

The following arguments are accepted:

* -a: UniProt accession of sequence of interest [required]
* -o: path to output file to create [required]
* -u: Map features to UniProt entry name instead of accession. This is probably required in most cases to successfully map features to the sequence in Jalview, which uses the uniprot record name.
* -i: InterPro accession to include in outputs. This can be specified multiple times. By default, if this option is not selected, results will be produced for all InterPro accessions

## Colour Selection

The colour palette is a colour-blind friendly palette selected from [microshades](https://karstenslab.github.io/microshades/) CVD palette. There are 21 colours available - if you end up with more feature types than this it will complain and fail, in which case you will need to update the `palette` list at the top of the script.

The values in this list can also be tweaked to suit your preferences


