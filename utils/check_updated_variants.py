#!/usr/bin/env python3

import os
import sys
import argparse

#Read mapping results
def get_mapping(fn):
    var = {}
    f = open(fn)
    header = f.readline().strip().split("\t")
    for line in f:
        tmp = line.strip().split("\t")
        var[tmp[0]] = tmp
    f.close()
    return header, var

#Create update analysis
def create_update_check(fn1, fn2, fn3, header1, header2, header3, var1, var2, var3):
    f1 = open(fn1, "w")
    f2 = open(fn2, "w")
    f3 = open(fn3, "w")
    f1.write(header1 + "\n")
    f2.write(header2 + "\n")
    f3.write(header3 + "\n")
    
    for k, v in var1.items():
        if k not in var2 and k not in var3:
            f1.write("\t".join(v) + "\n")
        elif k in var2:
            if v != var2[k]:
                f2.write("\t".join(v + var2[k]) + "\n")
        elif k in var3:
            if v != var3[k]:
                f3.write("\t".join(v + var3[k]) + "\n")
        else:
            pass
    f1.close()
    f2.close()
    f3.close()

parser = argparse.ArgumentParser(description='Check result files for update analysis')
parser.add_argument('-cd', '--curr_dna', type=str, required=True, help='current dna mapping file')
parser.add_argument('-cp', '--curr_prot', type=str, required=True, help='current prot mapping file')
parser.add_argument('-pd', '--prev_dna', type=str, required=True, help='previous dna mapping file')
parser.add_argument('-pp', '--prev_prot', type=str, required=True, help='previous prot mapping file')
parser.add_argument('-o', '--out_dir', type=str, required=True, help='output directory for update analysis files')

args = vars(parser.parse_args())
prev_dna_fn = args['prev_dna']
prev_prot_fn = args['prev_prot']
curr_dna_fn = args['curr_dna']
curr_prot_fn = args['curr_prot']
out_dir = args['out_dir']

if not os.path.exists(out_dir):
	os.makedirs(out_dir)

dna_header, prev_dna = get_mapping(prev_dna_fn)
_, curr_dna = get_mapping(curr_dna_fn)
prot_header, prev_prot = get_mapping(prev_prot_fn)
_, curr_prot = get_mapping(curr_prot_fn)

#Updated
fn1 = out_dir + "/curr_civic_dna_added.tsv"
fn2 = out_dir + "/curr_civic_dna_updateof_prev_civic_dna.tsv"
fn3 = out_dir + "/curr_civic_dna_updateof_prev_civic_prot.tsv"
header1 = "\t".join(dna_header)
header2 = "\t".join(["curr_" + t for t in dna_header] + ["prev_" + t for t in dna_header])
header3 = "\t".join(["curr_" + t for t in dna_header] + ["prev_" + t for t in prot_header])
create_update_check(fn1, fn2, fn3, header1, header2, header3, curr_dna, prev_dna, prev_prot)

fn1 = out_dir + "/curr_civic_prot_added.tsv"
fn2 = out_dir + "/curr_civic_prot_updateof_prev_civic_prot.tsv"
fn3 = out_dir + "/curr_civic_prot_updateof_prev_civic_dna.tsv"
header1 = "\t".join(prot_header)
header2 = "\t".join(["curr_" + t for t in prot_header] + ["prev_" + t for t in prot_header])
header3 = "\t".join(["curr_" + t for t in prot_header] + ["prev_" + t for t in dna_header])
create_update_check(fn1, fn2, fn3, header1, header2, header3, curr_prot, prev_prot, prev_dna)

#Changed
fn1 = out_dir + "/prev_civic_dna_removed.tsv"
fn2 = out_dir + "/prev_civic_dna_updatedby_curr_civic_dna.tsv"
fn3 = out_dir + "/prev_civic_dna_updatedby_curr_civic_prot.tsv"
header1 = "\t".join(dna_header)
header2 = "\t".join(["prev_" + t for t in dna_header] + ["curr_" + t for t in dna_header])
header3 = "\t".join(["prev_" + t for t in dna_header] + ["curr_" + t for t in prot_header])
create_update_check(fn1, fn2, fn3, header1, header2, header3, prev_dna, curr_dna, curr_prot)

fn1 = out_dir + "/prev_civic_prot_removed.tsv"
fn2 = out_dir + "/prev_civic_prot_updatedby_curr_civic_prot.tsv"
fn3 = out_dir + "/prev_civic_prot_updatedby_curr_civic_dna.tsv"
header1 = "\t".join(prot_header)
header2 = "\t".join(["prev_" + t for t in prot_header] + ["curr_" + t for t in prot_header])
header3 = "\t".join(["prev_" + t for t in prot_header] + ["curr_" + t for t in dna_header])
create_update_check(fn1, fn2, fn3, header1, header2, header3, prev_prot, curr_prot, curr_dna)
