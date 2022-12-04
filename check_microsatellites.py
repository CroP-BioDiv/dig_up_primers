#!/usr/bin/env python3

import os
import csv
import json
import subprocess
import pandas
from collections import defaultdict
from Bio import Entrez
from Bio.Blast import NCBIXML
from dig_up_primers import DigProject, _misa_ini


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
            m_data = json.load(_in)
    else:
        m_data = _fetch_microsatellites(species)
        with open(m_filename, 'w', encoding='utf-8') as _out:
            json.dump(m_data, _out, indent=2)

    def _repeat_region(m):
        for d in m.get('GBSeq_feature-table', []):
            if d.get('GBFeature_key') == 'repeat_region' and (loc := d.get('GBFeature_location')):
                start, end = loc.split('..')
                start, end = (int(start) - 1), int(end)
                # Sanity check
                seq = m['GBSeq_sequence'][start:end] if end - start < len(m['GBSeq_sequence']) / 2 else '???'
                return dict(start=start, end=end, sequence=seq)

    return dict((m['GBSeq_locus'],
                  dict(seq_id=m['GBSeq_locus'], sequence=m['GBSeq_sequence'], length=len(m['GBSeq_sequence']),
                      repeats=[], hits=[], repeat_region=_repeat_region(m)))
              for m in m_data)


def _blast_microsatellites(out_dir, m_data, first_assembly, params):
    # Create query file
    query_filename = os.path.join(out_dir, 'microsatellites.fa')
    with open(query_filename, 'w', encoding='utf-8') as _out:
        _out.write(''.join(f">{m['seq_id']}\n{m['sequence'].upper()}\n" for m in m_data.values()))

    b_xml = 'blast_result.xml'
    blast_filename = os.path.join(out_dir, b_xml)
    if not os.path.isfile(blast_filename):
        cmd = ['blastn', '-db', first_assembly.get_blast_db(), '-query', 'microsatellites.fa', '-outfmt', '5', '-out', b_xml,
               '-perc_identity', str(params.perc_identity), '-qcov_hsp_perc', str(params.qcov_hsp_perc),
               '-num_threads', str(params.num_threads)]
        print(f"Executing: {' '.join(cmd)}  (cwd: {out_dir})")
        subprocess.run(cmd, cwd=out_dir)

    # Collect results
    with open(blast_filename, 'r') as result:
        for record in NCBIXML.parse(result):
            if record.alignments:
                m_hits = m_data[record.query]['hits']  # This is a list!
                for align in record.alignments:
                    seq_id = align.hit_id
                    # If primer_idx starts with known prefixes (NC_, ..), than Blast hit_id can be like: ref|<id>|, gb|<id|, emb|<id>|
                    if seq_id.endswith('|'):
                        seq_id = seq_id[seq_id.index('|') + 1:-1]
                    m_hits.extend(dict(seq_id=seq_id, hsp=hsp, repeats=[], ssrs=[]) for hsp in align.hsps)


def _misa_on_all(project, out_dir, m_data):
    misa_f = 'mics_hits.fa'
    misa_result_filename = os.path.join(out_dir, f'{misa_f}.misa')
    if not os.path.isfile(misa_result_filename):
        # Create query file
        query_filename = os.path.join(out_dir, misa_f)
        with open(query_filename, 'w', encoding='utf-8') as _out:
            for m in m_data.values():
                _out.write(''.join(f">{m['seq_id']}\n{m['sequence'].upper()}\n"))
                for idx, hit in enumerate(m['hits']):
                    _out.write(''.join(f">{m['seq_id']}_{idx}\n{hit['hsp'].sbjct.upper()}\n"))

        with open(os.path.join(out_dir, 'misa.ini'), 'w') as out:
            out.write(_misa_ini.format(repeats=project.params['misa_repeats']))
        cmd = ['misa.pl', misa_f]
        print(f"Executing: {' '.join(cmd)}  (cwd: {out_dir})")
        subprocess.run(cmd, cwd=out_dir)

    # Read MISA output
    with open(misa_result_filename, 'r') as _misa:
        header = next(_misa)
        for line in _misa:
            mic_id, _, ssr_type, ssr, size, start, end = line.split()
            ssr = dict(ssr_type=ssr_type, ssr=ssr, size=int(size), start=int(start), end=int(end))
            hit_idx = None
            if mic := m_data.get(mic_id):
                mic['repeats'].append(ssr)
            else:
                mic_id, hit_idx = mic_id.rsplit('_', 1)
                m_data[mic_id]['hits'][int(hit_idx)]['repeats'].append(ssr)


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

    # Blast microsatellites on the first assembly
    _blast_microsatellites(out_dir, m_data, first_assembly, params)

    # Find repeats on microsatellites and blasted hits
    _misa_on_all(project, out_dir, m_data)

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

    for m in m_data.values():
        for m_hit in m['hits']:
            hsp = m_hit['hsp']
            for ssr in seq_id_2_ssrs.get(m_hit['seq_id'], []):
                if hsp.sbjct_start < ssr.start and ssr.end < hsp.sbjct_end:
                    si = ssr.ssr_idx
                    m_hit['ssrs'].append(dict(ssr=ssr,
                                              in_misa=(si in ssrs_misa),
                                              in_p3=(si in ssrs_p3),
                                              in_rm=(si in ssrs_rm),
                                              in_amplifies=[(si in s) for s in ssrs_amps],
                                              in_final=(si in ssrs_final)))

    # Results in excel
    _pm = lambda b: '+' if b else '-'
    _R = lambda repeats: '; '.join(rp['ssr'] for rp in repeats) if repeats else None
    _MC = lambda m: [m['seq_id'], None, m['length'], _R(m['repeats']),
                     rp['sequence'] if (rp := m.get('repeat_region')) else None]
    _HI = lambda h: [None, h['seq_id'], len(h['hsp'].sbjct), _R(h['repeats'])]  # h['hsp'].identities,
    _S = lambda s: [None, None, None, f"({s['ssr'].motif}){s['ssr'].repeats}", None, _pm(s['in_p3']), _pm(s['in_rm'])] + \
        [_pm(a) for a in s['in_amplifies']] + [_pm(s['in_final'])]

    sheets = []
    columns = ['Microsatellite', 'Sequence', 'Length', 'Repeats', 'Declared repeats', 'Primer3', 'RepeatMasker'] + \
        [f'Amplification_{i + 1}' for i in range(len(project.assembly_objs))] + ['Final']
    s_m_data = sorted(m_data.values(), key=lambda m: m['seq_id'])

    if rows := [_MC(m) for m in s_m_data if not m['hits']]:
        sheets.append(('Not found', columns[:5], rows))

    rows = []  # Todo: dodati header-e
    for m in s_m_data:
        if m['hits'] and all(not h['ssrs'] for h in m['hits']):
            rows.append(_MC(m))
            rows.extend(_HI(h) for h in m['hits'])
    if rows:
        sheets.append(('Without SSRs', columns[:5], rows))


    rows = []
    for m in s_m_data:
        if m['hits'] and any(h['ssrs'] for h in m['hits']):
            rows.append(_MC(m))
            for h in m['hits']:
                rows.append(_HI(h))
                rows.extend(_S(s) for s in h['ssrs'])
    if rows:
        sheets.append(('With SSRs', columns, rows))

    writer = pandas.ExcelWriter(os.path.join(out_dir, 'results.xlsx'))
    for name, columns, rows in sheets:
        pandas.DataFrame(rows, columns=columns).to_excel(writer, sheet_name=name, index=None, header=True)
    writer.save()

    # Summary
    l_amps = '\n'.join(f'      Amplification {i} : {len(ssrs):>6}' for i, ssrs in enumerate(ssrs_amps))
    c_size = 10
    not_found = [i for i, m in sorted(m_data.items()) if not m['hits']]
    nf = '\n'.join(('    ' + ', '.join(not_found[i:i + c_size])) for i in range(0, len(not_found), c_size))
    no_hit_with_ssrs = [i for i, m in sorted(m_data.items()) if m['hits'] and all(not h['ssrs'] for h in m['hits'])]
    nhs = '\n'.join(('    ' + ', '.join(no_hit_with_ssrs[i:i + c_size])) for i in range(0, len(no_hit_with_ssrs), c_size))

    summary = f"""Summary

Microsatellites:
    Number        : {len(m_data)}
    Length range  : {min(m['length'] for m in m_data.values())} - {max(m['length'] for m in m_data.values())}

Project data:
    MISA min repeats  : {project.params['misa_repeats']}
    Product range     : {project.params['product_size_range']}
    MISA SSRs         : {len(ssrs_misa):>6}
    Contigs with SSRs : {len(seq_id_2_ssrs):>6}
    Primer 3 SSRS     : {len(ssrs_p3):>6}
    RepeatMasker SSRS : {len(ssrs_rm):>6}
    Amplifications
{l_amps}
    Final SSRS        : {len(ssrs_final):>6}

Blast data:
    Mapped microsatellites : {sum(int(bool(m['hits'])) for m in m_data.values()):>6}
    Overall hits           : {sum(len(m['hits']) for m in m_data.values()):>6}

Not blasted ({sum(int(not m['hits']) for m in m_data.values())}):
{nf}

Blast hits without our SSRs ({len(no_hit_with_ssrs)}):
{nhs}

In final: {sum(1 for m in m_data.values() if any(any(s['in_final'] for s in h['ssrs']) for h in m['hits']))}
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

