"""
@author: Alexey Markin
"""

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from dendropy import Tree
from typing import List
import re
import collections
import subprocess
import os
import sys
from datetime import datetime
import multiprocessing
import statistics

segments = {'PB2': ('1', 'PB2'),
            'PB1': {'2', 'PB1'},
            'PA': {'3', 'PA'},
            'HA': {'4', 'HA', 'H1', 'H3'},
            'NP': {'5', 'NP'},
            'NA': {'6', 'NA', 'N1', 'N2'},
            'MP': {'7', 'MP', 'M'},
            'NS': {'8', 'NS'}}

start_codon = 'ATG'
stop_codons = ['TAA', 'TGA', 'TAG']


def extract_dates(path: str) -> str:
    records = SeqIO.parse(path, 'fasta')
    file_name = path + '.dates.csv'
    dates_file = open(file_name, 'w+')
    dates_file.write('name, date\n')
    for record in records:
        name = record.name
        date_str = name.split('|')[-1]
        date = datetime.strptime(date_str, '%Y-%m-%d')
        dec_date = date.year + ((date.month-1)*30 + date.day)/365.0
        dates_file.write('%s, %.2f\n' % (name.replace('/', '_').replace('|', '_'),
                                         dec_date))
    dates_file.close()
    return file_name


def parse_iav_record(iav: SeqRecord, sep='|'):
    tokens = iav.id.split(sep)
    for token in tokens:
        # First see if token specifies the segment.
        for seg, aliases in segments.items():
            if token.upper() in aliases:
                iav.segment = seg
        
        if re.match(r'A/([a-zA-Z_\-]+/)?[^/]+/[^/]+/[^/]+', token):
            # Found strain name.
            iav.strain_name = token
        elif re.match(r'[\d\-/]{4,}', token) and not re.match(r'\d{5,}', token):
            # Date: consists of digits and '/' or '-' symbols.
            try:
                iav.date_str = token
                if token.count('/') == 2:
                    iav.date = datetime.strptime(token, '%m/%d/%Y')
                elif token.count('/') == 1:
                    iav.date = datetime.strptime(token, '%m/%Y')
                elif token.count('-') == 2:
                    iav.date = datetime.strptime(token, '%Y-%m-%d')
                elif token.count('-') == 1:
                    iav.date = datetime.strptime(token, '%Y-%m')
                else:
                    iav.date = datetime.strptime(token, '%Y')
            except Exception as e:
                print(e)
                raise Exception('Cannot parse date %s in header %s' % (token, iav.id))
        elif re.match('swine|human', token.lower()):
            # Host.
            iav.host = token.lower()
        elif re.match(r'H\dN\d', token.upper()):
            # Subtype.
            iav.subtype = token.upper()
        elif token.startswith('EPI_ISL'):
            iav.epi_isl_id = token


def unify_and_classify(fasta_path):
    dates = {}  # Make sure to unify dates.
    for segment in segments.keys():
        seg_path = f'{segment}-{fasta_path}'
        final_path = segment + '.final.aln'
        cds_path = seg_path + '.cds.aln'
        records = list(AlignIO.read(cds_path, 'fasta'))
        for r in records:
            parse_iav_record(r)
            date = r.date.strftime('%Y-%m-%d')
            if dates.get(r.strain_name):
                date = dates.get(r.strain_name)
            else:
                dates[r.strain_name] = date
            r.id = '|'.join([r.strain_name, date])
            r.description = ''
        SeqIO.write(records, final_path, 'fasta')


def concatenate_wgs(fasta_path: str):
    segment_by_name = {}
    dates = {}
    wgs_path = f'{fasta_path}.wgs.aln'
    for segment in segments.keys():
        seg_path = f'{segment}-{fasta_path}'
        cds_path = seg_path + '.cds.aln'
        records = list(AlignIO.read(cds_path, 'fasta'))
        for r in records:
            parse_iav_record(r)
            date = r.date.strftime('%Y-%m-%d')
            if not dates.get(r.strain_name):
                dates[r.strain_name] = date
            if not segment_by_name.get(r.strain_name, None):
                segment_by_name[r.strain_name] = [r]
            else:
                segment_by_name[r.strain_name].append(r)

    wgs_records = []
    for name, seg_list in segment_by_name.items():
        wgs_record = seg_list[0]
        for record in seg_list[1:]:
            wgs_record.seq += record.seq
        wgs_record.id = '|'.join([name, dates[name]])
        wgs_record.description = ''
        wgs_records.append(wgs_record)
    SeqIO.write(wgs_records, wgs_path, 'fasta')


def find_start_codon_pos(record: SeqRecord):
    return record.seq.upper().find(start_codon)


def find_stop_codon_pos(record: SeqRecord):
    return max((record.seq.upper().rfind(stop_codon) for stop_codon in stop_codons))


def prune_to_cds(records: List[SeqRecord]):
    starts = [find_start_codon_pos(r) for r in records]
    ends = [find_stop_codon_pos(r) + 3 for r in records]
    start = int(statistics.median(starts))
    end = int(statistics.median(ends))
    print(start, end)
    
    for r in records:
        r.seq = r.seq[start:end]    
        

def align_and_check(fasta_path: str):
    all_records = []
    cds_all_path = fasta_path + '.cds.fasta'
    for segment in segments.keys():
        seg_path = f'{segment}-{fasta_path}'
        aln_path = seg_path + '.aln'
        cds_path = seg_path + '.cds.aln'
        with open(aln_path, 'w') as aln:
            subprocess.call(['mafft', '--thread', str(multiprocessing.cpu_count()),
                             seg_path], stdout=aln)
        records = list(AlignIO.read(aln_path, 'fasta'))
        prune_to_cds(records)
        SeqIO.write(records, cds_path, 'fasta')
        all_records.extend(records)
    for r in all_records:
        parse_iav_record(r)
        r.id = r.strain_name + '_' + r.segment
        r.description = ''
        r.seq = r.seq.replace('-', '')
    SeqIO.write(all_records, cds_all_path, 'fasta')


def process_wgs(fasta_path: str):
    records = list(SeqIO.parse(fasta_path, 'fasta'))
    segs_by_name = {}
    for record in records:
        name_match = re.search(r'A/([^/]+/)?[^/]+/[^/]+/[^/|]+', record.id)
        if not name_match:
            print('No valid name for', record.id)
            continue
        parse_iav_record(record)
        record.name = record.strain_name
        if not record.segment:
            print('Cannot determine the segment for', record.id)
            continue
        segs = segs_by_name.get(record.strain_name, set())
        segs.add(record.segment)
        segs_by_name[record.strain_name] = segs
    
    wgs_names = [name for name in segs_by_name.keys() if
                 len(segs_by_name.get(name, set())) == 8]
    
    wgs_iavs = {seg: [] for seg in segments.keys()}
    added_sequences = set()
    for record in records:
        if record.name in wgs_names and record.segment and\
           (record.name.upper(), record.segment) not in added_sequences:
            wgs_iavs[record.segment].append(record)
            added_sequences.add((record.name.upper(), record.segment))
    for segment in segments.keys():
        SeqIO.write(wgs_iavs[segment], f'{segment}-{fasta_path}', 'fasta')
    align_and_check(fasta_path)
    unify_and_classify(fasta_path)

    concatenate_wgs(fasta_path)
    

if __name__ == '__main__':
    args = sys.argv[1:]
    process_wgs(args[0])
