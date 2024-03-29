#!/usr/bin/env python3

import os
import shutil
import re
import json
import time
import subprocess
from collections import namedtuple, defaultdict
from itertools import chain, product
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from Bio.Blast import NCBIXML


class TimeSteps:
    def __init__(self, filename):
        self.filename = filename
        if os.path.isfile(filename):
            with open(filename, 'r') as _in:
                self.steps = json.load(_in)
        else:
            self.steps = []  # dicts
        self.steps_in_progress = []

    def store(self):
        with open(self.filename, 'w', encoding='utf-8') as _out:
            json.dump(self.steps, _out, indent=4)

    def start(self, name):
        self.steps_in_progress.append(dict(name=name, start=time.time()))

    def _readable(self, lasted):
        fs = [0, 0, 0]
        if lasted >= 3600:
            fs[0] = f'{int(lasted // 3600)}h'
            lasted = lasted % 3600
        if lasted >= 60:
            fs[1] = f'{int(lasted // 60)}m'
            lasted = lasted % 60
        fs[2] = f'{int(round(lasted))}s'
        f_i = 0 if fs[0] else (1 if fs[1] else 2)
        return ' '.join(fs[f_i:])

    def end(self):
        d = self.steps_in_progress.pop()
        d['end'] = time.time()
        d['lasted'] = d['end'] - d['start']
        d['readable'] = self._readable(d['lasted'])
        for old in self.steps:
            if old['name'] == d['name']:
                if old['lasted'] < d ['lasted']:  # Take longer
                    old.update((a, d[a]) for a in ('start', 'end', 'lasted', 'readable'))
                    self.store()
                return
        self.steps.append(d)
        self.store()

    def start_method(self, name, method, *args, **kwargs):
        self.start(name)
        r = method(*args, **kwargs)
        self.end()


class _Primer(namedtuple('_Primer', 'p_idx, left, right, length, gc, tm')):
    def __new__(cls, p_idx, left, right, length, gc, tm):
        return super(_Primer, cls).__new__(cls, int(p_idx), left, right, int(length), float(gc), float(tm))

    @classmethod
    def columns(cls):
        return cls._fields


class _SSR(namedtuple('_SSR', 'ssr_idx, seq_id, start, end, motif, repeats, seq_part, primers')):
    def __new__(cls, ssr_idx, seq_id, start, end, motif, repeats, seq_part, primers=None):
        if primers is None:
            primers = []
        assert isinstance(primers, list), primers
        return super(_SSR, cls).__new__(
            cls, ssr_idx, seq_id, int(start), int(end), motif, int(repeats), seq_part, primers)

    def list_no_primers(self):
        return self[:-1]

    @classmethod
    def columns(cls):
        return cls._fields[:-1]


class _AssemblyProtocol:
    def __init__(self, project, assembly_idx, assembly_file, assembly_dir):
        self.project = project
        self.assembly_idx = assembly_idx
        self.assembly_file = assembly_file
        self.assembly_dir = assembly_dir
        self.misa_dir = os.path.join('1_MISA', assembly_dir)
        self.misa_ssrs_csv = os.path.join(self.misa_dir, 'ssrs.csv')
        self.rm_ssrs_csv = os.path.join(self.rm_dir, 'ssrs.csv')
        self.primer3_primers_csv = os.path.join(self.primer3_dir, 'primers.csv')
        self.amplify_dir = os.path.join('5_amplify', assembly_dir)
        self.amplify_primers_csv = os.path.join(self.amplify_dir, 'primers.csv')
        self.blast_db_dir = os.path.join('BlastDB', assembly_dir)
        #
        self._misa_ssrs = None
        self._rm_ssrs = None
        self._primer3_ssrs = None
        self._amplified_primers = None

    def read_assembly_dict(self):
        return dict((rec.id, str(rec.seq)) for rec in SeqIO.parse(self.assembly_file, 'fasta'))

    #
    def store_misa_ssrs(self, output_ssrs):
        bps = self.project.params['space_around']
        store_ssrs = []
        _parse = SeqIO.parse(self.assembly_file, 'fasta')
        _current_seq_id = None
        for seq_id, start, end, motif, repeats in output_ssrs:
            # Find sequence record
            if _current_seq_id != seq_id:
                for record in _parse:
                    if record.id == seq_id:
                        _current_seq_id = seq_id
                        break
                else:
                    assert False, seq_id
            #
            if start > bps and end + bps < len(record.seq):  # Is SSR on sequence end
                seq_part = record.seq[start - 1 - bps:end + bps]  # Stores sequence part
                store_ssrs.append([f'{self.assembly_idx}_{len(store_ssrs)}', seq_id, start, end, motif, repeats, seq_part])
        #
        self.project.write_csv(self.misa_ssrs_csv, store_ssrs, _SSR.columns())
        self._misa_ssrs = [_SSR(*args) for args in store_ssrs]
        return self._misa_ssrs

    def get_misa_ssrs(self):
        if self._misa_ssrs is None:
            self._misa_ssrs = [_SSR(*args) for args in self.project.read_csv(self.misa_ssrs_csv)]
        return self._misa_ssrs

    #
    def store_rm_ssrs(self, rm_ssrs):
        self._rm_ssrs = rm_ssrs
        self.project.write_csv(self.rm_ssrs_csv, [sn.list_no_primers() for sn in rm_ssrs], _SSR.columns())

    def get_rm_ssrs(self):
        if self._rm_ssrs is None:
            self._rm_ssrs = [_SSR(*args) for args in self.project.read_csv(self.rm_ssrs_csv)]
        return self._rm_ssrs

    #
    def store_primer_ssrs(self, p_ssrs):
        self._primer3_ssrs = p_ssrs
        self.project.write_primers(self.primer3_primers_csv, p_ssrs)

    def get_primer_ssrs(self):
        if self._primer3_ssrs is None:
            self._primer3_ssrs = self.project.read_primers(self.primer3_primers_csv)
        return self._primer3_ssrs

    #
    @staticmethod
    def merge_ssrs(assembly_objs, filename):
        # Just concatenate content of Primer3 files
        with open(filename, 'wb') as wfd:
            for ao in assembly_objs:
                with open(ao.get_ssrs_filename_for_merging(), 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

    #
    def store_aplified_primers(self, ampl_res):
        self.project.write_csv(self.amplify_primers_csv, ampl_res, ['primer_idx', 'num_repeats', 'ampl_length'])
        self._amplified_primers = ampl_res

    def get_aplified_primers(self):
        if self._amplified_primers is None:
            self._amplified_primers = list(self.project.read_csv(self.amplify_primers_csv))
        return self._amplified_primers

    #
    def get_blast_db(self):
        if not os.path.isdir(self.blast_db_dir):
            self.project.ensure_dir(self.blast_db_dir)
            cmd = ['makeblastdb', '-dbtype', 'nucl', '-parse_seqids',
                   '-in', self.assembly_file, '-out', 'db', '-title', 'db']
            self.project._exe_prog(cmd, self.blast_db_dir)
        return os.path.abspath(os.path.join(self.blast_db_dir, 'db'))


class _Assembly_RM_Primer3(_AssemblyProtocol):
    def __init__(self, project, assembly_idx, assembly_file, assembly_dir):
        self.rm_dir = os.path.join('2_RepeatMasker', assembly_dir)
        self.primer3_dir = os.path.join('3_Primer3', assembly_dir)
        super().__init__(project, assembly_idx, assembly_file, assembly_dir)

    def get_subdirs(self):
        return ['1_MISA', '2_RepeatMasker', '3_Primer3', '5_amplify']

    def get_ssrs_for_repeat_masker(self):
        return self.get_misa_ssrs()

    def get_ssrs_for_primer3(self):
        return self.get_rm_ssrs()

    def get_ssrs_filename_for_merging(self):
        return self.primer3_primers_csv


class _Assembly_Primer3_RM(_AssemblyProtocol):
    def __init__(self, project, assembly_idx, assembly_file, assembly_dir):
        self.primer3_dir = os.path.join('2_Primer3', assembly_dir)
        self.rm_dir = os.path.join('3_RepeatMasker', assembly_dir)
        super().__init__(project, assembly_idx, assembly_file, assembly_dir)

    def get_subdirs(self):
        return ['1_MISA', '2_Primer3', '3_RepeatMasker', '5_amplify']

    def get_ssrs_for_primer3(self):
        return self.get_misa_ssrs()

    def get_ssrs_for_repeat_masker(self):
        return self.get_primer_ssrs()

    def store_rm_ssrs(self, rm_ssrs):
        self._rm_ssrs = rm_ssrs
        self.project.write_primers(self.rm_ssrs_csv, rm_ssrs)

    def get_ssrs_filename_for_merging(self):
        return self.rm_ssrs_csv


#
class DigProject:
    _primer_line_start = '  '
    _merged_primers_csv = '4_merged_primers.csv'
    _final_primers_csv = '6_final_primers.csv'

    def __init__(self, params):
        self.params = params
        self.assembly_cls = _Assembly_RM_Primer3 if params.get('workflow', 'rp') == 'rp' else _Assembly_Primer3_RM
        self._merged_primers = None
        self.num_threads = params.get('num_threads', 1)
        # MISA
        misa_repeats = params['misa_repeats'].split(',')
        self.used_repeats = [int(r.split('-')[0]) for r in misa_repeats]
        self.min_repeats = dict((int(r.split('-')[0]), int(r.split('-')[1])) for r in misa_repeats)
        self.misa_repeats = ' '.join(misa_repeats)
        #
        zero_pad = '{:0>' + str(len(str(len(self.params['assemblies'])))) + '}'
        self.assembly_objs = [self.assembly_cls(self, idx, a_f, f'assembly_{zero_pad.format(idx)}')
                              for idx, a_f in enumerate(params['assemblies'])]

    def get_time_steps(self):
        return TimeSteps('duration_of_steps.json')

    def merge_ssrs(self):
        self.assembly_cls.merge_ssrs(self.assembly_objs, self._merged_primers_csv)

    def merge_ssrs_first(self):
        os.symlink(self.assembly_objs[0].get_ssrs_filename_for_merging(), self._merged_primers_csv)

    def get_merged_primers(self):
        if self._merged_primers is None:
            self._merged_primers = self.read_primers(self._merged_primers_csv)
        return self._merged_primers

    def write_primers(self, csv_filename, p_ssrs):
        pls = self._primer_line_start
        with open(csv_filename, 'w') as _csv:
            _csv.write(';'.join(_SSR.columns()) + '\n')
            _csv.write(pls + ';'.join(_Primer.columns()) + '\n')
            for sn in p_ssrs:
                _csv.write(';'.join(map(str, sn.list_no_primers())) + '\n')
                assert sn.primers
                for p in sn.primers:
                    _csv.write(pls + ';'.join(map(str, p)) + '\n')

    def read_primers(self, csv_filename):
        pls = self._primer_line_start
        ssrs = []
        with open(csv_filename, 'r') as _csv:
            next(_csv)  # Skip header
            next(_csv)
            for line in _csv:
                fields = line.strip().split(';')
                if line.startswith(pls):
                    ssrs[-1].primers.append(_Primer(*fields))
                else:
                    ssrs.append(_SSR(*fields))
        return ssrs

    #
    def write_results(self, final_primers):
        final_primers.sort()  # Sort results by the score
        amp_cs = tuple(chain(*((f'a{idx}_repeats', f'a{idx}_length') for idx in range(1, len(self.assembly_objs) + 1))))
        self.write_csv(self._final_primers_csv,
                       [ssr[:-2] + primer + tuple(chain(*amp)) + score for score, ssr, primer, amp in final_primers],
                        _SSR.columns()[:-1] + _Primer.columns() + amp_cs + ('score', 'max_score'))

    def read_results(self):
        amp_cs = tuple(chain(*((f'a{idx}_repeats', f'a{idx}_length') for idx in range(1, len(self.assembly_objs) + 1))))
        columns = _SSR.columns()[:-1] + _Primer.columns() + amp_cs + ('score', 'max_score')
        res_cls = namedtuple('_Result', columns)
        # Set right data types
        _n = lambda x: x
        f_ms = [_n, _n, int, int, _n, int, int, _n, _n, int, float, float] + [int] * (2 * len(self.assembly_objs) + 2)
        with open(self._final_primers_csv, 'r') as _csv:
            next(_csv)  # Skip header
            return [res_cls(*(f(v) for f, v in zip(f_ms, line.strip().split(';')))) for line in _csv]

    #
    def delete_from(self, delete_from):
        # No more than 9 steps :-)
        for f in [self._merged_primers_csv, self._final_primers_csv]:
            if f[0].isdigit() and int(f[0]) >= delete_from and os.path.isfile(f):
                print(f'Removing file: {f}')
                os.remove(f)

        for _dir in self.assembly_objs[0].get_subdirs():
            f = os.path.basename(_dir)
            if f[0].isdigit() and int(f[0]) >= delete_from and os.path.isdir(_dir):
                print(f'Removing directory: {_dir}')
                shutil.rmtree(_dir)

    #
    def _exe_prog(self, cmd, cwd):
        print(f"Executing: {' '.join(cmd)}  (cwd: {cwd})")
        subprocess.run(cmd, cwd=cwd)

    def ensure_dir(self, _dir):
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    def write_csv(self, csv_filename, data, header):
        with open(csv_filename, 'w') as _csv:
            _csv.write(';'.join(header) + '\n')
            _csv.write('\n'.join(';'.join(map(str, row)) for row in data))

    def read_csv(self, csv_filename):
        with open(csv_filename, 'r') as _in:
            next(_in)  # Skip header
            yield from (line.strip().split(';') for line in _in)

    def write_text(self, filename, text):
        with open(filename, 'w') as _out:
            _out.write(text)


# -------------------------------------------------------------------
# MISA
# -------------------------------------------------------------------
_misa_ini = """
definition(unit_size,min_repeats):             {repeats}
interruptions(max_difference_between_2_SSRs):  1
GFF:                                           false
"""
_misa_ssr = re.compile(r'^\(([ATCG]+)\)([0-9]+)$')


def _misa(a_data, proj):
    if os.path.isfile(a_data.misa_ssrs_csv):
        return

    # Run MISA (if needed)
    misa_result_filename = os.path.join(a_data.misa_dir, 'input.fa.misa')
    if not os.path.exists(misa_result_filename):
        proj.ensure_dir(a_data.misa_dir)

        # Set input data and parameters
        i_f = os.path.join(a_data.misa_dir, 'input.fa')
        if not os.path.exists(i_f):
            os.symlink(a_data.assembly_file, i_f)

        with open(os.path.join(a_data.misa_dir, 'misa.ini'), 'w') as out:
            out.write(_misa_ini.format(repeats=proj.misa_repeats))

        proj._exe_prog(['misa.pl', 'input.fa'], a_data.misa_dir)

    # Filter only simple SSRs
    id_2_ssrs = defaultdict(list)  # id -> list of SSRs
    num_misa_ssrs = 0
    num_composite = 0
    num_motif_1 = 0  # Note: newer version take care about shorter
    with open(misa_result_filename, 'r') as _misa:
        header = next(_misa)
        for line in _misa:
            num_misa_ssrs += 1
            _id, _, ssr_type, ssr, size, start, end = line.split()
            # Remove composite repeats
            if ssr_type.startswith('c'):
                num_composite += 1
                continue
            assert ssr_type.startswith('p'), ssr_type
            _s, _rep = _misa_ssr.findall(ssr)[0]

            # Remove repeats with motif of length 1 if that length is not specified
            if 1 not in proj.used_repeats and all(x == _s[0] for x in _s):
                num_motif_1 += 1
                continue
            #
            id_2_ssrs[_id].append((_id, int(start), int(end), _s, int(_rep)))

    # Remove SSRs that are close
    output_ssrs = []
    num_close = 0
    space_around = proj.params['space_around']
    # Note: it is important to iterate in input sequences order. It is guaranteed in python 3.7+!
    for _id, ssrs in id_2_ssrs.items():
        ssrs.sort()  # Inline sort. Note: newer MISA version take a care about sorting!
        to_remove = set(chain.from_iterable((idx, idx + 1) for idx, ssr2 in enumerate(ssrs[1:])
                                            if ssr2[1] - ssrs[idx][2] <= space_around))
        num_close += len(to_remove)
        output_ssrs.extend([ssr for idx, ssr in enumerate(ssrs) if idx not in to_remove])

    # Store SSRs in CSV file
    long_ssrs = ''
    if sml := proj.params['ssr_max_length']:
        _num_before = len(output_ssrs)
        output_ssrs = [ssr for ssr in output_ssrs if len(ssr[3]) * ssr[4] <= sml]
        long_ssrs = f'\nLong SSRs             : {_num_before - len(output_ssrs)}'
    if mms := proj.params.get('misa_max_ssrs'):  # For testing!
        output_ssrs = output_ssrs[:mms]
    num_before_store = len(output_ssrs)
    output_ssrs = a_data.store_misa_ssrs(output_ssrs)

    motif_1 = f'SSRs of motif with 1  : {num_motif_1}\n' if 1 not in proj.used_repeats and num_motif_1 else ''
    proj.write_text(os.path.join(a_data.misa_dir, 'report.txt'), f"""REPORT

SSRs found by MISA    : {num_misa_ssrs}
Composite SSRs        : {num_composite}{long_ssrs}
Removed close to SSRs : {num_close}
Removed close to ends : {num_before_store - len(output_ssrs)}
{motif_1}
Good SSRs             : {len(output_ssrs)}
""")


# -------------------------------------------------------------------
# RepeatMasker
# -------------------------------------------------------------------
def _repeat_masker(a_data, proj):
    if os.path.isfile(a_data.rm_ssrs_csv):
        return

    # Input file for RepeatMasker
    ssrs = dict((sn.ssr_idx, sn) for sn in a_data.get_ssrs_for_repeat_masker())
    proj.ensure_dir(a_data.rm_dir)
    query_fa = os.path.join(a_data.rm_dir, 'query.fa')
    with open(query_fa, 'w') as _query:
        for sn in ssrs.values():
            _query.write(f'>{sn.ssr_idx}\n')
            _query.write(f'{sn.seq_part}\n')

    # Call RepeatMasker
    if not os.path.isfile(query_fa + '.masked'):
        # cmd = ['RepeatMasker', 'query.fa']
        cmd = ['RepeatMasker', '-e', 'hmmer', 'query.fa']
        if proj.num_threads > 1:
            # ToDo: add an argument. For now suppose nhmmer - 2 cores
            if (_nt := proj.num_threads // 2) > 1:
                cmd.extend(['-pa', str(_nt)])
        proj._exe_prog(cmd, a_data.rm_dir)

    # Process output
    close_to_coding = set()
    with open(query_fa + '.out', 'r') as _in:
        for _ in range(3):  # Skip header
            next(_in)
        for line in _in:
            fields = line.split()
            if fields[10] not in ('Low_complexity', 'Simple_repeat'):
                close_to_coding.add(fields[4])

    max_repeat_diff = proj.params['max_repeat_diff']
    rm_ssrs = []
    for sec_record in SeqIO.parse(query_fa + '.masked', 'fasta'):
        sn = ssrs[sec_record.id]
        num_diff = sum(1 for a, b in zip(sn.seq_part.upper(), str(sec_record.seq).upper()) if a != b)
        if (num_diff - len(sn.motif) * sn.repeats) <= max_repeat_diff:
            rm_ssrs.append(sn)

    #
    a_data.store_rm_ssrs(rm_ssrs)
    proj.write_text(os.path.join(a_data.rm_dir, 'report.txt'), f"""REPORT

Input SSRs       : {len(ssrs)}
Good SSRs        : {len(rm_ssrs)}

Close to coding  : {len(close_to_coding)}
Low complexity   : {len(ssrs) - len(rm_ssrs) - len(close_to_coding)}
""")


# -------------------------------------------------------------------
# Primer3
# -------------------------------------------------------------------
_sequence_params = """SEQUENCE_ID={seq_id}
SEQUENCE_TEMPLATE={seq}
SEQUENCE_TARGET={target}
=
"""

_primer_params = """Primer3 File - http://primer3.org
P3_FILE_TYPE=settings

P3_FILE_ID=SSR primers
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_NUM_RETURN={primer3_num_return}
PRIMER_MIN_SIZE={min_size}
PRIMER_MAX_SIZE={max_size}
PRIMER_MIN_GC={min_gc}
PRIMER_MAX_GC={max_gc}
PRIMER_MIN_TM={min_tm}
PRIMER_MAX_TM={max_tm}
PRIMER_MAX_POLY_X={max_poly_x}
PRIMER_PRODUCT_SIZE_RANGE={product_size_range}
PRIMER_EXPLAIN_FLAG=1
=
"""


def _exec_primer3(a_data, proj, ssrs, input_f, output_f, error_f):
    num_around = 10
    na_2 = 2 * num_around
    target = proj.params['space_around'] + 1 - num_around
    with open(os.path.join(a_data.primer3_dir, input_f), 'w') as _out:
        for sn in ssrs:
            _out.write(_sequence_params.format(
                seq_id=sn.ssr_idx, seq=sn.seq_part, target=f"{target},{len(sn.motif) * sn.repeats + na_2}"))

    # Execute command
    # ToDo: max returned, default is 5
    proj._exe_prog(['primer3_core', '--p3_settings_file=primer3.set', '--strict_tags',
                    f'--output={output_f}', f'--error={error_f}', input_f], a_data.primer3_dir)


def _primer3(a_data, proj):
    if os.path.isfile(a_data.primer3_primers_csv):
        return

    ssrs = dict((sn.ssr_idx, sn) for sn in a_data.get_ssrs_for_primer3())

    # Write settings file
    proj.ensure_dir(a_data.primer3_dir)
    with open(os.path.join(a_data.primer3_dir, 'primer3.set'), 'w') as _out:
        _out.write(_primer_params.format(**proj.params))

    # Write input file
    output_files = []
    if proj.num_threads == 1:
        output_files.append('primer3.output')
        if not os.path.isfile(os.path.join(a_data.primer3_dir, output_files[-1])):
            _exec_primer3(a_data, proj, ssrs.values(), 'primer3.input', 'primer3.output', 'primer3.error')
    else:
        with ThreadPoolExecutor(max_workers=proj.num_threads) as executor:
            l_ssrs = list(ssrs.values())
            k, m = divmod(len(l_ssrs), proj.num_threads)
            for i in range(proj.num_threads):
                output_files.append(f'primer3_{i}.output')
                if not os.path.isfile(os.path.join(a_data.primer3_dir, output_files[-1])):
                    executor.submit(_exec_primer3, a_data, proj, l_ssrs[i * k + min(i, m):(i + 1) * k + min(i + 1, m)],
                                    f'primer3_{i}.input', output_files[-1], f'primer3_{i}.error')

    # Process output
    num_primers = 0
    for o_file in output_files:
        with open(os.path.join(a_data.primer3_dir, o_file), 'r') as _in:
            seq_data = dict()
            for line in _in:
                line = line.strip()
                if line == '=':
                    if num := int(seq_data['PRIMER_PAIR_NUM_RETURNED']):
                        ssr_idx = seq_data['SEQUENCE_ID']
                        for p_idx in range(num):
                            num_primers += 1
                            ssrs[ssr_idx].primers.append(_Primer(
                                p_idx,
                                seq_data[f'PRIMER_LEFT_{p_idx}_SEQUENCE'],
                                seq_data[f'PRIMER_RIGHT_{p_idx}_SEQUENCE'],
                                seq_data[f'PRIMER_PAIR_{p_idx}_PRODUCT_SIZE'],
                                seq_data[f'PRIMER_INTERNAL_{p_idx}_GC_PERCENT'],
                                seq_data[f'PRIMER_INTERNAL_{p_idx}_TM']))
                    #
                    seq_data = dict()
                else:
                    key, value = line.split('=', 1)
                    seq_data[key] = value

    # Store data
    num_input_ssr = len(ssrs)
    ssrs = [sn for _, sn in sorted(ssrs.items()) if sn.primers]
    a_data.store_primer_ssrs(ssrs)
    proj.write_text(os.path.join(a_data.primer3_dir, 'report.txt'), f"""REPORT

Input SSRs              : {num_input_ssr}
Number of primers       : {num_primers}
SSRs with primers       : {len(ssrs)}
Regions without primers : {num_input_ssr - len(ssrs)}
""")


# -------------------------------------------------------------------
# Amplify primers
# -------------------------------------------------------------------
def _amplify(a_data, proj):
    if os.path.isfile(a_data.amplify_primers_csv):
        return

    ssrs = proj.get_merged_primers()
    ssr_idx_2_ssr = dict((ssr.ssr_idx, ssr) for ssr in ssrs)

    # Run blastn, if needed
    res_xml = os.path.join(a_data.amplify_dir, 'results.xml')
    if not os.path.isfile(res_xml):
        proj.ensure_dir(a_data.amplify_dir)
        with open(os.path.join(a_data.amplify_dir, 'query.fa'), 'w') as _query:
            for sn in ssrs:
                for p in sn.primers:
                    _query.write(f'>{sn.ssr_idx}_{p.p_idx}_L\n{p.left}\n')
                    _query.write(f'>{sn.ssr_idx}_{p.p_idx}_R\n{p.right}\n')

        cmd = ['blastn', '-db', a_data.get_blast_db(), '-query', 'query.fa', '-outfmt', '5', '-out', 'results.xml',
               '-task', 'blastn-short', '-perc_identity', '100', '-evalue', '1e-1']
        if proj.num_threads > 1:
            cmd.extend(['-num_threads', str(proj.num_threads)])
        # evalue=1e-1 filters sequences shorter than 15char. I hope :-)
        # ToDo: calculate evalue?
        # E = m x n  / 2^(bit-score) ; m - query sequence length ; n - total database length (sum of all sequences)
        proj._exe_prog(cmd, a_data.amplify_dir)

    # Collect data from result file
    blast_res = defaultdict(lambda: defaultdict(lambda: ([], [])))  # primer_idx -> (seq_id -> pair of left/right lists)
    with open(res_xml, 'r') as result:
        for record in NCBIXML.parse(result):
            if record.alignments:  # Primer side (left, rigth)
                primer_idx_side = record.query
                primer_idx = primer_idx_side[:-2]
                primer_side = int(primer_idx_side[-1] == 'R')
                query_length = record.query_length

                for align in record.alignments:
                    seq_id = align.hit_id
                    # If primer_idx starts with known prefixes (NC_, ..), than Blast hit_id can be like: ref|<id>|, gb|<id|, emb|<id>|
                    if seq_id.endswith('|'):
                        seq_id = seq_id[seq_id.index('|') + 1:-1]

                    for hsp in align.hsps:
                        # Match has to be exact! Controlled, not 100%, with arguments perc_identity=100, evalue=1e-1
                        if hsp.align_length != query_length or hsp.identities != query_length:
                            continue
                        # Just to be sure!
                        assert hsp.strand[0] == 'Plus'
                        assert hsp.query_start == 1, hsp.query_start
                        assert hsp.query_end == query_length, hsp.query_end
                        blast_res[primer_idx][seq_id][primer_side].append((hsp.sbjct_start, hsp.sbjct_end, hsp.strand[1]))

    # Find amplifications
    more_amplifications = 0
    zero_amplifications = 0
    one_wrong_amplifications = 0
    num_blast_results = 0
    max_ampl = proj.params['max_amplification_length']
    min_p_length, max_p_length = map(int, proj.params['product_size_range'].split('-'))
    _seqs = a_data.read_assembly_dict()

    ampl_res = []
    for primer_idx, seq_data in blast_res.items():
        num_blast_results += 1
        amplifications = 0
        ampl_length = None
        ampl_region = None
        for seq_id, (left, right) in seq_data.items():
            if left and right:
                # Note: lists are short, not needed anything clever, like sorting and bisect()!
                for (l_start, l_end, l_strand), (r_start, r_end, r_strand) in product(left, right):
                    if l_strand != r_strand:  # Can't be in the same direction!
                        dist = (r_start - l_start) if l_strand == 'Plus' else (l_start - r_start)
                        if 0 <= dist <= max_ampl:
                            amplifications += 1
                            if amplifications > 1:
                                more_amplifications += 1
                                break
                            ampl_length = dist
                            _all_ls = (l_start, l_end, r_start, r_end)
                            ampl_region = [seq_id, min(_all_ls), max(_all_ls)]
                if amplifications > 1:
                    break
        #
        if amplifications == 1:
            # Check is amplified region in good range
            ssr_idx, p_idx = primer_idx.rsplit('_', 1)  # ToDo: split() is better
            if min_p_length <= ampl_length <= max_p_length and \
               (desc := _ampl_num_motifs(proj, ssr_idx_2_ssr[ssr_idx], int(p_idx), _seqs[ampl_region[0]], *ampl_region[1:])):
                ampl_res.append([primer_idx, *desc])
            else:
                one_wrong_amplifications += 1
        elif not amplifications:
            zero_amplifications += 1
    #
    a_data.store_aplified_primers(ampl_res)
    proj.write_text(os.path.join(a_data.amplify_dir, 'report.txt'), f"""REPORT

Blast results            : {num_blast_results}
One amplification        : {len(ampl_res)}

More amplifications      : {more_amplifications}
One wrong amplifications : {one_wrong_amplifications}
Zero amplifications      : {zero_amplifications}
""")


def _ampl_num_motifs(proj, ssr, p_idx, contig, start, end):
    # return [1, 1]
    primer = ssr.primers[p_idx]
    inner = contig[max(0, start + len(primer.left) - 1):(end - len(primer.right))]
    min_r = max(3, proj.min_repeats[len(ssr.motif)] - 3)
    if (m := re.search(rf'({ssr.motif}){{{min_r},}}', inner, re.IGNORECASE)):
        s, e = m.span()
        return [(e - s) // len(ssr.motif), end - start + 1]


# -------------------------------------------------------------------
# Finilize
# -------------------------------------------------------------------
def _finalize(proj):
    if os.path.isfile(proj._final_primers_csv):
        return

    # Collect set of primer idxs
    amplifications = [dict((p_idx, [n, _l]) for p_idx, n, _l in a.get_aplified_primers()) for a in proj.assembly_objs]
    primers = set.intersection(*(set(d.keys()) for d in amplifications))

    # Score methods
    _center = lambda x, c, r: abs(x - c) / r
    def _asm_c(ssr, primer, amp):
        min_r = proj.min_repeats[len(ssr.motif)]
        score = sum(int(n < min_r) for n, _ in asm) * 100

    #
    _mehtods = []
    if pp_len := proj.params['prefer_primer_length']:
        _mehtods.append(lambda _, primer, _a: _center(len(primer.left), *pp_len) + _center(len(primer.right), *pp_len))
    if pp_idx := proj.params['prefer_primer_idx']:
        _mehtods.append(lambda _, primer, _a: _center(primer.p_idx, *pp_idx))
    if p_gc := proj.params['prefer_gc']:
        _mehtods.append(lambda _, primer, _a: _center(primer.gc, *p_gc))
    if p_tm := proj.params['prefer_tm']:
        _mehtods.append(lambda _, primer, _a: _center(primer.tm, *p_tm))

    if _mehtods:
        def _calc_score(ssr, primer, amp):  # Less is better
            # int(s*100) removes digits, precision is futile :-)
            scores = [int(100 * m(ssr, primer, amp)) for m in _mehtods]
            return sum(scores), max(scores)
    else:
        def _calc_score(ssr, primer, amp):
            return 0, 0

    # Collect data
    final_primers = []
    current_assembly = None
    last_ssr_idx = None
    for p_idx in sorted(primers):  # Assemblies are sorted, and primers should be sorted if less than 10 are used :-)
        ssr_idx, primer_index = p_idx.rsplit('_', 1)
        if last_ssr_idx == ssr_idx:  # SSR already processed
            continue

        amp = [d[p_idx] for d in amplifications]
        if len(amp) > 1 and all(a == amp[0] for a in amp[1:]):  # All apmplifications are the same!
            continue

        last_ssr_idx = ssr_idx
        primer_index = int(primer_index)
        assembly_idx = ssr_idx.split('_', 1)[0]
        #
        if assembly_idx != current_assembly:
            current_assembly = assembly_idx
            a_data = proj.assembly_objs[int(assembly_idx)]
            a_ssrs = dict((ssr.ssr_idx, ssr) for ssr in a_data.get_primer_ssrs())
        #
        ssr = a_ssrs[ssr_idx]
        primer = ssr.primers[primer_index]
        final_primers.append((_calc_score(ssr, primer, amp), ssr, primer, amp))
    #
    proj.write_results(final_primers)
    proj.write_text('report.txt', f"""REPORT

Final results            : {len(final_primers)}
""")


# -------------------------------------------------------------------
# Workflow
# -------------------------------------------------------------------
def find_params(params):
    settings_file = os.path.join(params.working_dir, 'ssrs_settings.json')
    if os.path.isfile(settings_file):
        with open(settings_file, 'r') as _in:
            settings = json.load(_in)

    else:
        # ToDo: check data?
        assert params.assembly
        settings = dict(assemblies=[os.path.abspath(a) for a in params.assembly],
                        workflow=params.workflow,
                        misa_max_ssrs=params.misa_max_ssrs,
                        space_around=params.space_around,
                        misa_repeats=params.misa_repeats,
                        ssr_max_length=params.ssr_max_length,
                        max_repeat_diff=params.max_repeat_diff,
                        primer3_num_return=params.primer3_num_return,
                        product_size_range=params.product_size_range,
                        min_size=params.min_size, max_size=params.max_size,
                        min_gc=params.min_gc, max_gc=params.max_gc,
                        min_tm=params.min_tm, max_tm=params.max_tm,
                        max_poly_x=params.max_poly_x,
                        max_amplification_length=params.max_amplification_length,
                        prefer_primer_length=_prefered(params.prefer_primer_length, [23, 3]),
                        prefer_primer_idx=_prefered(params.prefer_primer_idx, [0, 3]),
                        prefer_gc=_prefered(params.prefer_gc, [50, 10]),
                        prefer_tm=_prefered(params.prefer_gc, [60, 2]))
        #
        if not os.path.exists(params.working_dir):
            os.makedirs(params.working_dir)
        with open(settings_file, 'w', encoding='utf-8') as _out:
            json.dump(settings, _out, indent=4)
        #
        os.chdir(params.working_dir)
    #
    settings['num_threads'] = params.num_threads
    return settings


def _prefered(desc, default):
    if not desc:
        return default
    fields = desc.split(':')
    if len(fields) == 2:
        return [float(fields[0]), float(fields[1])]


def process_project(params, delete_from=None):
    proj = DigProject(params)
    time_steps = proj.get_time_steps()

    if delete_from:
        proj.delete_from(delete_from)

    if not os.path.isfile(proj._merged_primers_csv):
        a_data = proj.assembly_objs[0]
        # for a_data in proj.assembly_objs:  # To locate SSRs in all assemblies
        time_steps.start_method('1_MISA', _misa, a_data, proj)
        if proj.params['workflow'] == 'rp':
            time_steps.start_method('2_RepeatMasker', _repeat_masker, a_data, proj)
            time_steps.start_method('3_Primer3', _primer3, a_data, proj)
        else:
            time_steps.start_method('2_Primer3', _primer3, a_data, proj)
            time_steps.start_method('3_RepeatMasker', _repeat_masker, a_data, proj)

        time_steps.start_method('4_Merge', proj.merge_ssrs_first)
        # proj.merge_ssrs()                  # In case of locating SSRs in all assemblies

    for idx, a_data in enumerate(proj.assembly_objs):
        time_steps.start_method(f'5_Amplify_{idx}', _amplify, a_data, proj)

    time_steps.start_method('6_Finalize', _finalize, proj)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""Finds SSRs primers. Works on one or more assemblies of releate species.""")

    #
    parser.add_argument(
        '-D', '--delete-from', type=int, help="Delete project's data from given step.")
    parser.add_argument(
        '-T', '--num-threads', default=1, type=int, help="Number of threads to use")

    #
    parser.add_argument(
        '-w', '--working-dir', default='.', help="Project's directory. All data will be stored in this directory.")
    parser.add_argument(
        '-W', '--workflow', default='pr', choices=('pr', 'rp'),
        help="Workflow to use. Primer3 -> RepeatMasker or RepeatMasker -> Primer3")

    # Assembly data
    parser.add_argument('-a', '--assembly', action='append', help='Assembly')

    # General stuff
    parser.add_argument(
        '-s', '--space-around', default=200, type=int,
        help='Number of bp to take around SSR to check for other SSRs, repeats, and find primers.')

    # Misa params
    parser.add_argument(
        '-m', '--misa-repeats', default='2-10,3-7',
        help='MISA min repeat definition. Simiar format as in MISA, except space is replaced with comma.')
    parser.add_argument('--ssr-max-length', default=80, type=int, help='Limit SSR length')
    parser.add_argument('--misa-max-ssrs', type=int, help='Limit number of SSRs MISA outputs. For testing.')

    # RepeatMasker params
    parser.add_argument(
        '--max-repeat-diff', type=int, default=5,
        help='How many bps can RepeatMasker repeat be longer than MISA repeat')

    # Primer3 params
    parser.add_argument('--primer3-num-return', type=int, default=1, help="Primer3: PRIMER_NUM_RETURN")
    parser.add_argument('--product-size-range', default='100-250', help="Primer3: PRIMER_PRODUCT_SIZE_RANGE")
    parser.add_argument('--min-size', type=int, default=19, help="Primer3: PRIMER_MIN_SIZE")
    parser.add_argument('--max-size', type=int, default=23, help="Primer3: PRIMER_MAX_SIZE")
    parser.add_argument('--min-gc', type=int, default=40, help="Primer3: PRIMER_MIN_GC")
    parser.add_argument('--max-gc', type=int, default=60, help="Primer3: PRIMER_MAX_GC")
    parser.add_argument('--min-tm', type=int, default=58, help="Primer3: PRIMER_MIN_TM")
    parser.add_argument('--max-tm', type=int, default=62, help="Primer3: PRIMER_MAX_TM")
    parser.add_argument('--max-poly-x', type=int, default=4, help="Primer3: PRIMER_MAX_POLY_X")

    # Amplification params
    parser.add_argument(
        '--max-amplification-length', type=int, default=2000, help="Max length to check for primer amplification")

    # Primer preferences. Preferences are set with center and unit radius. Format is str: center,unit_radius.
    parser.add_argument(
        '--prefer-primer-length', help="Prefer primers with length closer to given values. Example 23:3")
    parser.add_argument(
        '--prefer-primer-idx', help="Prefer primers with length closer to given index. Example 0:1")
    parser.add_argument(
        '--prefer-gc', help="Prefer primers with GC closer to given values. Example 50:10")
    parser.add_argument(
        '--prefer-tm', help="Prefer primers with length closer to given values. Example 60:2")

    #
    cmd_params = parser.parse_args()
    params = find_params(cmd_params)
    process_project(params, delete_from=cmd_params.delete_from)
