#!/usr/bin/env python3

import os
import csv
import json
import subprocess
from collections import namedtuple, defaultdict
from Bio import Entrez
from Bio.Blast import NCBIXML
from dig_up_primers import DigProject

_hit_columns = ['microsatellite', 'seq_id', 'start', 'end', 'strand', 'align_length', 'identities', 'ssrs']
_Hit = namedtuple('_Hit', _hit_columns)


def _fetch_microsatellites(species, retmax=18000):
    def _entrez(handle):
        data = Entrez.read(handle)
        handle.close()
        return data

    db = 'nucleotide'
    term = f'("{species}"[Organism] AND microsatellite[Title]) AND tandem[All Fields] AND biomol_genomic[PROP]'
    print(f'Search: {db} "{term}"')
    search = _entrez(Entrez.esearch(db=db, term=term, usehistory='y', retmax=1, rettype='json'))
    data = []
    if (count := int(search['Count'])):
        n_pages = ((count - 1) // retmax) + 1
        args = dict(db=db, query_key=search['QueryKey'], WebEnv=search['WebEnv'], retmax=retmax, rettype='json')
        for page in range(n_pages):
            print(f'    Fetch: page {page + 1} / {n_pages}, records {count}')
            data += _entrez(Entrez.efetch(**args))
            args['retstart'] = (page + 1) * retmax
    return data


def _get_microsatellites(out_dir, species):
    m_filename = os.path.join(out_dir, 'microsatellites.json')
    if os.path.isfile(m_filename):
        with open(m_filename, 'r') as _in:
            return json.load(_in)
    #
    m_data = _fetch_microsatellites(species)
    with open(m_filename, 'w', encoding='utf-8') as _out:
        json.dump(m_data, _out, indent=2)
    return m_data


def _blast_microsatellites(out_dir, first_assembly, params):
    b_xml = 'blast_result.xml'
    b_csv = 'blast_result.csv'
    blast_filename = os.path.join(out_dir, b_xml)
    if not os.path.isfile(blast_filename):
        cmd = ['blastn', '-db', first_assembly.get_blast_db(), '-query', 'microsatellites.fa', '-outfmt', '5', '-out', b_xml,
               '-perc_identity', str(params.perc_identity), '-qcov_hsp_perc', str(params.qcov_hsp_perc),
               '-num_threads', str(params.num_threads)]
        print(f"Executing: {' '.join(cmd)}  (cwd: {out_dir})")
        subprocess.run(cmd, cwd=out_dir)

    # Collect results
    blast_csv = os.path.join(out_dir, b_csv)
    if not os.path.isfile(blast_csv):
        b_data = []
        with open(blast_filename, 'r') as result:
            for record in NCBIXML.parse(result):
                if record.alignments:  # Primer side (left, rigth)
                    microsatellite = record.query
                    for align in record.alignments:
                        seq_id = align.hit_id
                        # If primer_idx starts with NC_, than Blast hit_id can be like: ref|NC_041794.1| or gb|MPGU01000001.1|
                        if seq_id.startswith('ref|') and seq_id.endswith('|'):
                            seq_id = seq_id[4:-1]
                        elif seq_id.startswith('gb|') and seq_id.endswith('|'):
                            seq_id = seq_id[3:-1]
                        assert not seq_id.endswith('|'), seq_id

                        for hsp in align.hsps:
                            assert hsp.strand[0] == 'Plus'  # Just to be sure!
                            b_data.append(_Hit(
                                microsatellite, seq_id, hsp.sbjct_start, hsp.sbjct_end, int(hsp.strand[1] == 'True'), hsp.align_length, hsp.identities, []))
        with open(blast_csv, 'w') as _csv:
            _csv.write(';'.join(_hit_columns[:-1]) + '\n')
            _csv.write('\n'.join(';'.join(map(str, row[:-1])) for row in b_data))
    else:
        with open(blast_csv, 'r') as _csv:
            next(_csv)
            b_data = [_Hit(*fs[:2], *map(int, fs[2:]), []) for line in _csv if (fs := line.strip().split(';'))]

    mic_2_hits = defaultdict(list)        
    for h in b_data:
        mic_2_hits[h.microsatellite].append(h)
    return mic_2_hits


def check_microsatellites(params):
    # Project settings
    settings_file = os.path.join(params.project_working_dir, 'ssrs_settings.json')
    with open(settings_file, 'r') as _in:
        project = DigProject(json.load(_in))
    first_assembly = project.assembly_objs[0]

    # Ensure output directory
    out_dir = params.output_subdirectory or f"microsatellites_{params.species.replace(' ', '_')}"
    out_dir = os.path.join(params.project_working_dir, out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Fetch (or load) microsatellie data
    m_data = _get_microsatellites(out_dir, params.species)

    # Create query file
    query_filename = os.path.join(out_dir, 'microsatellites.fa')
    with open(query_filename, 'w', encoding='utf-8') as _out:
        _out.write(''.join(f">{r['GBSeq_locus']}\n{r['GBSeq_sequence'].upper()}\n" for r in m_data))

    # Blast microsatellites on the first assembly
    mic_2_hits = _blast_microsatellites(out_dir, first_assembly, params)

    # Read project data
    seq_id_2_ssrs = defaultdict(list)
    ssrs_misa = first_assembly.get_misa_ssrs()
    for ssr in ssrs_misa:  # ssr_idx='0_1', seq_id='NC_041789.1', start=25911, end=25940
        seq_id_2_ssrs[ssr.seq_id].append(ssr)
    ssrs_misa = set(s.ssr_idx for s in ssrs_misa)
    ssrs_p3 = set(r.ssr_idx for r in first_assembly.get_primer_ssrs())
    ssrs_rm = set(r.ssr_idx for r in project.read_primers(first_assembly.rm_ssrs_csv))
    ssrs_amps = [set(s[0].rsplit('_', 1)[0] for s in a.get_aplified_primers()) for a in project.assembly_objs]
    ssrs_final = set(r.ssr_idx for r in project.read_results())

    lines = []
    not_found = []
    no_hit_with_ssrs = []
    num_in_final = 0
    prod_min, prod_max = map(int, project.params['product_size_range'].split('-'))
    _ran = lambda _l: '+' if prod_min <= _l <= prod_max else '-'
    _pm = lambda b: '+' if b else '-'
    for m in sorted(m_data, key=lambda m: m['GBSeq_locus']):
        microsatellite = m['GBSeq_locus']
        m_len = int(m['GBSeq_length'])
        if m_hits := mic_2_hits.get(microsatellite):
            in_final = False
            mm = microsatellite
            for hit in m_hits:
                for ssr in seq_id_2_ssrs.get(hit.seq_id, []):
                    if hit.start < ssr.start and ssr.end < hit.end:
                        si = ssr.ssr_idx
                        if si in ssrs_final:
                            in_final = True
                        _amp = '  '.join(_pm(si in s) for s in ssrs_amps)
                        lines.append(f'    {mm:<10} : {m_len:>4}  {_ran(m_len)}    {_pm(si in ssrs_misa)}  {_pm(si in ssrs_p3)}' +
                                     f'  {_pm(si in ssrs_rm)}  {_amp}  {_pm(si in ssrs_final)}')
                        mm = ''
            if mm:
                no_hit_with_ssrs.append(microsatellite)
            if in_final:
                num_in_final += 1
        else:
            not_found.append(microsatellite)
            # lines.append(f'    {microsatellite:<10} : not found')
    #
    l_amps = '\n'.join(f'      Amplification {i} : {len(ssrs):>6}' for i, ssrs in enumerate(ssrs_amps))
    mics = '\n'.join(lines)
    c_size = 10
    nf = '\n'.join(('    ' + ', '.join(not_found[i:i + c_size])) for i in range(0, len(not_found), c_size))
    nhs = '\n'.join(('    ' + ', '.join(no_hit_with_ssrs[i:i + c_size])) for i in range(0, len(no_hit_with_ssrs), c_size))
    

    summary = f"""Summary

Microsatellites:
    Number        : {len(m_data)}
    Length range  : {min(len(m['GBSeq_sequence']) for m in m_data)} - {max(len(m['GBSeq_sequence']) for m in m_data)}

Project data:
    MISA min repeats  : {project.params['misa_repeats']}
    Product range     : {prod_min} - {prod_max}
    MISA SSRs         : {len(ssrs_misa):>6}
    Contigs with SSRs : {len(seq_id_2_ssrs):>6}
    Primer 3 SSRS     : {len(ssrs_p3):>6}
    RepeatMasker SSRS : {len(ssrs_rm):>6}
    Amplifications
{l_amps}
    Final SSRS        : {len(ssrs_final):>6}

Blast data:
    Mapped microsatellites : {len(mic_2_hits):>6}
    Overall hits           : {sum(len(v) for v in mic_2_hits.values()):>6}

Not blasted ({len(not_found)}):
{nf}

Blast hits without our SSRs ({len(no_hit_with_ssrs)}):
{nhs}

In final: {num_in_final}

Process:          Len In   MI P3 RM {' '.join(f'A{i}' for i in range(len(ssrs_amps)))} FI
{mics}
"""

    print(summary)
    with open(os.path.join(out_dir, 'summary.txt'), 'w') as _out:
        _out.write(summary)



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Checks existing markers against dig_up_primer's results")
    parser.add_argument('species', help="Species name to collecting microsatellite sequences")
    parser.add_argument('-p', '--project-working-dir', default='.', help="Directory of project with results.")
    parser.add_argument('-o', '--output-subdirectory', help="Subdirectory of project's directory to store results. Default is species name")

    parser.add_argument('-P', '--perc-identity', default=70, type=int, help="Blast perc_identity argument")
    parser.add_argument('-Q', '--qcov-hsp-perc', default=70, type=int, help="Blast qcov_hsp_perc argument")
    parser.add_argument('-T', '--num-threads', default=1, type=int, help="Number of threads to use")

    check_microsatellites(parser.parse_args())

