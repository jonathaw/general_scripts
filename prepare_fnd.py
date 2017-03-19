#!/usr/bin/env python3.5
"""
a script to preapre fold & dock runs
"""
import argparse
import os
import sys
import re
import shutil
import random
import string
import time
import subprocess
from collections import OrderedDict
from AASeq import AASeq, read_seq, read_seqs
from Logger import Logger
import seq_funcs
from MyPDB import parse_PDB


rosetta_path = '/home/labs/fleishman/jonathaw/Rosetta/'
rosetta_bin_path = '/home/labs/fleishman/jonathaw/Rosetta/main/source/bin/'
rosetta_scripts_default = "rosetta_scripts.default.linuxgccrelease"
rosetta_scripts_static = "rosetta_scripts.static.linuxgccrelease"
fnd_protocol = '/home/labs/fleishman/jonathaw/elazaridis/protocols/fnd_7Dec.xml'
design_protocol = '/home/labs/fleishman/jonathaw/elazaridis/protocols/design_13Nov.xml'
# rosetta_path = '/home/labs/fleishman/sarel/rosetta/'
# rosetta_bin_path = '/home/labs/fleishman/sarel/rosetta/main/source/bin/'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', type=str, help='prtoein name')
    parser.add_argument('-mode', type=str, default='fnd')
    parser.add_argument('-fasta', type=str, default=None,
                        help='path to fasta file')
    parser.add_argument('-seq', type=str, default=None,
                        help='sequence to model')
    parser.add_argument('-path', type=str, default=os.getcwd()+'/',
                        help='path to create files in')
    parser.add_argument('-topo_string', type=str, default=None,
                        help='topology string')
    parser.add_argument('-topo_string_file', type=str,
                        default=None, help='topo string file')
    parser.add_argument('-some_pdb', default=None,
                        help='pdb for fragment picker')
    parser.add_argument('-native_pdb', type=str,
                        default='/home/labs/fleishman/jonathaw/elazaridis/gpa_modeling/gpa_wt_28Jun/data/1afo_AB.pdb', help='native PDB file')
    parser.add_argument('-span_file', type=str)
    parser.add_argument('-symm_type', type=str, default='cn',
                        help='symmetry type (cn)')
    parser.add_argument('-symm_num', type=int, default=2,
                        help='number of units in symmetry')
    parser.add_argument('-is_symm', type=bool, default=True,
                        help='is this symmetric?')
    parser.add_argument('-fnd_protocol', type=str,
                        default=fnd_protocol, help='rosetta script for fold & dock')
    parser.add_argument('-design_protocol', type=str, default=design_protocol)
    parser.add_argument('-nstruct', type=int, default=100,
                        help='nstruct argument to pass to rosetta')
    parser.add_argument('-num_jobs', type=int, default=1000,
                        help='number of jobs to create')
    parser.add_argument('-queue', type=str, default='fleishman',
                        help='queue to run on')
    parser.add_argument('-avoid_psipred', type=bool, default=True,
                        help='whether to use psipred')
    parser.add_argument('-database',
                        default='/home/labs/fleishman/jonathaw/Rosetta/main/database')
    parser.add_argument('-energy_func', default=['beta', 'mpframework'],
                        type=str, nargs='+')
    parser.add_argument('-ss2_file', default=None,
                        help='whther an external ss2 file is provided')
    parser.add_argument('-1st_span_orient', default=None,
                        help='predetermine the 1st spans orientation')
    parser.add_argument('-2nd_span_orient', default=None)
    parser.add_argument('-chain', default='A',
                        help='if only a pdb is given, which chains to use')
    parser.add_argument('-hetero', type=bool, default=False, help='is hetero?')
    parser.add_argument('-rosetta', default='default')
    args = vars(parser.parse_args())
    if args['some_pdb'] is None:
        args['some_pdb'] = args['native_pdb']
    if args['rosetta'] == 'default':
        args['rosetta'] = rosetta_scripts_default
    elif args['rosetta'] == 'static':
        args['rosetta'] = rosetta_scripts_static

    args['time'] = time.strftime("%d.%0-m")
    args['logger'] = Logger('log_file_%s.log' % args['time'])

    if args['mode'] == 'fnd':
        if args['hetero']:
            fnd_hetero(args)
        else:
            fnd(args)


    elif args['mode'] == 'from_design':
        from_design(args)

    else:
        print('no mode')


def fnd_hetero(args: dict):
    args['is_symm'] = False
    args['symm_num'] = 0
    args['seq'] = read_seqs(args['fasta'])
    test_input_new(args)

    make_data_new(args)

    make_bin(args)
    make_fnd_jobs(args)
    args['logger'].close()


def test_input_new(args):
    assert args['fasta'] is not None, "no fasta file provided"
    args['seqs'] = read_seqs(args['fasta'])
    assert args['topo_string_file'] is not None, "no topo string file provided"
    args['topo_strings'] = OrderedDict()
    for p in open(args['topo_string_file'], 'r').read().split('>')[1:]:
        s = p.split('\n')
        args['topo_strings'].update({s[0]: s[1]})
    print(args['seqs'])
    print(args['topo_strings'])
    for k, v in args['seqs'].items():
        assert len(args['topo_strings'][k]) == len(v.get_seq()), "for %s not same length" % k


def from_design(args):
    pass


def fnd(args):
    test_input(args)

    make_data(args)
    make_bin(args)
    make_fnd_jobs(args)

    args['logger'].close()


def test_input(args):
    if args['seq'] is None:
        pdb = parse_PDB(args['native_pdb'])
        seq_dict = pdb.get_seq()
        if args['hetero']:
            args['seq'] = seq_dict
        else:
            args['seq'] = seq_dict[args['chain']]
    args['logger'].log('NO SEQ GIVEN, using %s, %i AAs' % (args['seq'],
                                                           len(args['seq'])))
    if args['topo_string'] is None:
        if args['hetero']:
            args['topo_string'] = {}
            for name, seq in args['seq'].items():
                args['topo_string'][name] = 'H'*len(seq)
        else:
            args['topo_string'] = 'H' * len(args['seq'])
        args['logger'].log('NO TOPO-STRING GIVEN, using helix, creating file')
        with open('%s.ts' % args['name'], 'w+') as fout:
            if args['hetero']:
                for name, ts in args['topo_string'].items():
                    fout.write('>%s\n%s' % (name, ts))
            else:
                fout.write('>%s\n' % args['name'])
                fout.write('%s' % args['topo_string'])
        args['logger'].log('NO SEQ GIVEN, using %s, %i AAs' %
                           (args['seq'], len(args['seq'])))

    if args['hetero']:
        ok = all(len(args['seq'][name]) == len(args['topo_string'][name])
                 for name in args['seq'].keys())
    else:
        ok = len(args['seq']) == len(args['topo_string'])
    if args['mode'] == 'fnd':
        args['xml'] = fnd_protocol
    elif args['mode'] == 'design':
        args['xml'] = design_protocol
    else:
        sys.exit('mode was neither design nor fnd...')
    if not ok:
        sys.exit('sequence and topo string not same length')


def make_fnd_jobs(args):
    for energy_func in args['energy_func']:
        args['fnd_path_%s' % energy_func] = args['path'] + 'fnd_%s_%s/' % (args['time'], energy_func)
        os.mkdir(args['fnd_path_%s' % energy_func])
        os.chdir(args['fnd_path_%s' % energy_func])
        args['logger'].log('creating jobs for energy function: %s' % energy_func)
        job_args = {'@%s%s' % (args['path_bin'], args['flags_file_%s' % energy_func]): None}
        for i in range(args['num_jobs']):
            job_args['-out:prefix'] = '%i_' % i
            job_args['-out:file:silent'] = 'ab_%i.out' % i

            make_jobs(args, job_args)
        os.chdir('../')


def make_jobs(args, job_args=None):
    if job_args is None:
        job_args = {}
    pwd = os.getcwd() + '/'
    num = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(8))
    jobname = '%sjob.%s' % (pwd, num)
    outname = '%sout.%s' % (pwd, num)
    errname = '%serr.%s' % (pwd, num)
    cmdname = '%scommand' % pwd
    with open(jobname, 'w+') as job:
        job.write('#!/bin/bash\n')
        job.write('. /usr/share/lsf/conf/profile.lsf\n')
        job.write('cd ' + pwd + '\n')
        if args['queue'] == 'new-all.q':
            job.write(
                '/apps/RH6U4/blcr/0.8.5/bin/cr_run %s%s ' % (rosetta_bin_path, args['rosetta']))
        else:
            job.write('%s%s ' % (rosetta_bin_path, args['rosetta']))
        for k, v in job_args.items():
            if k == 'pdb_dump':
                continue
            if v is None:
                job.write('%s ' % k)
            else:
                job.write('%s %s ' % (k, v))

        job.write('\n')
    with open(cmdname, 'a+') as cmd:
        if args['queue'] == 'new-all.q':
            cmd.write(str(
                'bsub -C 1024 -u /dev/null -N -u /dev/null -R rusage[mem=1024] -G fleishman-wx-grp-lsf -q ' +
                args['queue'] + ' -o ' + outname + ' -e ' + errname + ' /apps/RH6U4/blcr/0.8.5/bin/cr_run ' +
                jobname + '\n'))
        else:
            cmd.write(str('bsub -N -u /dev/null -G fleishman-wx-grp-lsf -q ' + args['queue'] + ' -o ' +
                          outname + ' -e ' + errname + ' ' + jobname + '\n'))
    subprocess.call(['chmod', '+x', jobname])


def make_bin(args):
    args['path_bin'] = args['path'] + 'bin/'
    os.mkdir(args['path_bin'])
    os.chdir(args['path_bin'])

    create_fnd_flags(args, True)
    create_fnd_flags(args, False)


def create_fnd_flags(args, fnd) -> None:
    if args['hetero']:
        spans = get_hetero_spans(args)
    else:
        spans = get_spans_from_topo_string(args)

    number_of_units = 1
    if args['is_symm']:
        number_of_units = args['symm_num']
    for energy_func in args['energy_func']:
        args['flags_file_%s' % energy_func] = 'fnd_%s_%s_%s.flags' % (args['name'], args['time'], energy_func)
        args['logger'].log('creating flag file at %s' % args['flags_file_%s' % energy_func])
        if fnd:
            file_name = args['flags_file_%s' % energy_func]
        else:
            file_name = 'just_score_%s.flags' % energy_func
        with open(file_name, 'w+') as fout:
            fout.write('# general data\n')
            if fnd:
                if args['hetero']:
                    fout.write('-parser:protocol %s\n' % '/home/labs/fleishman/jonathaw/elazaridis/protocols/fnd_hetero_5Mar.xml')
                else:
                    fout.write('-parser:protocol %s\n' % args['xml'])
            else:
                fout.write('-parser:protocol %s\n' % "/home/labs/fleishman/jonathaw/elazaridis/protocols/just_score.xml")
            fout.write('-database %s\n' % args['database'])
            if fnd:
                fout.write('-in:file:fasta %s%s\n' % (args['path_data'], args['fasta']))
                fout.write('-in:file:native %s%s\n\n' % (args['path_data'], args['native_pdb'].split('/')[-1]))
            # fout.write('# membrane adjustments\n')
            # fout.write('-mp:scoring:hbond\n')
            # if energy_func not in  ['ResSolv', 'beta']:
                # fout.write('-parser:script_vars steepness=10\n')
                # fout.write('-parser:script_vars membrane_core=15\n')
                # fout.write('-in:file:spanfile %s%s.span\n' % (args['path_data'], args['name']))
            # else:
                # fout.write('-parser:script_vars steepness=4\n')
                # fout.write('-parser:script_vars membrane_core=10\n')
            # fout.write('\n')

            # fout.write('-mp:setup:spanfiles %s%s.span\n' % (args['path_data'], args['name']))
            # fout.write('-in:file:spanfile %s%s.span\n' % (args['path_data'], args['name']))
            if fnd:
                fout.write('-nstruct %i\n' % args['nstruct'])
                fout.write('-mute all\n')
            fout.write('-overwrite\n\n')
            # fout.write('-parser:script_vars energy_function=%s\n' % energy_func)
            # fout.write('-parser:script_vars score_func=mpframework_docking_cen_ELazaridis\n')
            if fnd:
                fout.write('# fragment stuff\n')
                if  args['hetero']:
                    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n trying to name the frags correctly, LOOK HERE!!!!! \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                    fout.write('-parser:script_vars frags9mersA=%sfrags_9%s.200.9mers\n' % (args['path_data'], args['seqs'].keys()[0]))
                    fout.write('-parser:script_vars frags3mersA=%sfrags_3%s.200.3mers\n' % (args['path_data'], args['seqs'].keys()[0]))
                    fout.write('-parser:script_vars frags9mersB=%sfrags_9%s.200.9mers\n' % (args['path_data'], args['seqs'].keys()[1]))
                    fout.write('-parser:script_vars frags3mersB=%sfrags_3%s.200.3mers\n' % (args['path_data'], args['seqs'].keys()[1]))
                else:
                    fout.write('-parser:script_vars frags9mers=%sfrags.200.9mers\n' % args['path_data'])
                    fout.write('-parser:script_vars frags3mers=%sfrags.200.3mers\n' % args['path_data'])

                if args['is_symm']:
                    fout.write('-parser:script_vars symm_file=%s\n\n' % (args['path_data'] + 'C%i.symm' % args['symm_num']))
            # fout.write('-parser:script_vars start_res=%i\n' % 1)
            # fout.write('-parser:script_vars end_res=%i\n' % len(args['AASeq'])*number_of_units)
            fout.write('# membrane spans:\n')
            for i, span in enumerate(spans):
                fout.write('-parser:script_vars span_start_%i=%i\n' % (i+1, span['start']))
                fout.write('-parser:script_vars span_end_%i=%i\n' % (i+1, span['end']))
                fout.write('-parser:script_vars span_orientation_%s=%s\n' % (i+1, span['orientation']))
            if fnd:
                if energy_func == 'mpframework':
                    for i in [0, 1, 2, 3, 5]:
                        fout.write('-parser:script_vars score_func_%s=score%s_membrane\n' % (str(i), str(i) if i == 0 else ''))
                else:
                    for i in [0, 1, 2, 3, 5]:
                        fout.write('-parser:script_vars score_func_%i=score%i_%s\n' %
                                (i, i, 'elazaridis' if energy_func in ['ResSolv', 'beta'] else ''))
            fout.write('\n')
            # if energy_func != 'ResSolv':
                # fout.write('-parser:script_vars steepness=10\n')
                # fout.write('-parser:script_vars membrane_core=15\n')
                # fout.write('-in:file:spanfile %s%s.span\n' % (args['path_data'], args['name']))
            # else:
                # fout.write('-parser:script_vars steepness=10\n')
                # fout.write('-parser:script_vars membrane_core=15\n')
            fout.write('# energy function stuff:\n')
            if args['mode'] == 'fnd':
                fout.write('-parser:script_vars energy_function=%s\n\n' % energy_func)
            fout.write('# membrane adjustments\n')
            fout.write('-mp:scoring:hbond\n')
            if energy_func not in  ['ResSolv', 'beta']:
                fout.write('-parser:script_vars steepness=10\n')
                fout.write('-parser:script_vars membrane_core=15\n')
                fout.write('-in:file:spanfile %s%s.span\n' % (args['path_data'], args['name']))
            else:
                fout.write('-parser:script_vars steepness=4\n')
                fout.write('-parser:script_vars membrane_core=10\n')
                fout.write('\n')
            if energy_func in ['ResSolv', 'beta']:
                fout.write('-score::elec_memb_sig_die\n')
                fout.write('-corrections::beta_nov16\n')
                fout.write('-score:memb_fa_sol\n')
                fout.write('\n')


def get_hetero_spans(args):
    spans_list = []
    spans_list.append({'start': 1,
                       'end': len(list(args['seqs'].items())[0][1].get_seq()),
                       'orientation': args['1st_span_orient']})
    spans_list.append({'start': len(list(args['seqs'].items())[0][1].get_seq())+1,
                       'end': len(list(args['seqs'].items())[0][1].get_seq()) + len(list(args['seqs'].items())[1][1].get_seq()),
                       'orientation': args['2nd_span_orient']})
    return spans_list


def get_spans_from_topo_string(args) -> list:
    topology_list = []
    topo_string = args['topo_string']
    hhh = re.compile('[hH]*')
    h_list = [(a.start(), a.end()) for a in hhh.finditer(topo_string) if a.end() - a.start() > 1
              and a.end() - a.start() >= 1]

    if args['1st_span_orient'] is None:
        direction = 'fwd' if topo_string[h_list[0][0] - 1] == '1' else 'rev'
    else:
        direction = args['1st_span_orient']
    for h in h_list:
        topology_list.append((h[0]+1, h[1], direction))
        direction = 'fwd' if direction == 'rev' else 'rev'
    spans = []
    for span in topology_list:
        ori = ''
        if span[2] == 'fwd':
            ori = 'in2out'
        elif span[2] =='rev':
            ori = 'out2in'
        else:
            sys.exit('illegal orientation found for span', span)
        spans.append({'start': span[0], 'end': span[1], 'orientation': ori})
    new_spans = spans.copy()
    for i in range(args['symm_num']-1):
        for span in spans:
            new_span = span.copy()
            new_span['start'] += len(args['seq'])
            new_span['end'] += len(args['seq'])
            new_spans.append(new_span)
    if 'logger' in args.keys():
        args['logger'].log('found the following spans for symmetry %i' % args['symm_num'])
        for span in new_spans:
            args['logger'].log('start: %i end: %i orientation: %s' % (span['start'], span['end'], span['orientation']))
    return new_spans


def make_data_new(args):
    args['path_data'] = args['path'] + 'data/'
    os.mkdir(args['path_data'])
    os.chdir(args['path_data'])

    for name, aaseq in args['seqs'].items():
        with open("%s.fasta" % name, 'w+') as fout:
            fout.write(aaseq.write())
        args['logger'].log('creating helix only ss2 for %s' % name)
        with open(args['path_data']+name+'.ss2', 'w+') as fout:
            fout.write('# PSIPRED VFORMAT (PSIPRED V3.5)\n\n')
            for i, aa in aaseq.enumerate():
                fout.write('%4i %s %s   0.000  1.000  0.000\n' % (i, aa, 'H'))

    run_fragment_picker_new(args)
    create_span(args)

    shutil.copy(args['native_pdb'], args['path_data'])
    os.chdir(args['path'])


def make_data(args):
    args['path_data'] = args['path'] + 'data/'
    os.mkdir(args['path_data'])
    os.chdir(args['path_data'])

    make_or_read_fasta(args)
    create_psipred(args)

    run_fragment_picker(args)
    create_span(args)
    if args['is_symm']:
        create_symm(args)
    shutil.copy(args['native_pdb'], args['path_data'])
    os.chdir(args['path'])


def create_psipred(args):
    if not args['avoid_psipred']:
        seq_funcs.run_psipred(args['fasta'])
    elif args['ss2_file'] != None:
        args['logger'].log('using provided ss2 file, provided here: %s' % args['ss2_file'])
        shutil.copy(args['ss2_file'], args['path_data']+args['name']+'.ss2')
    else:
        args['logger'].log('creating helix only ss2')
        with open(args['path_data']+args['name']+'.ss2', 'w+') as fout:
            fout.write('# PSIPRED VFORMAT (PSIPRED V3.5)\n\n')
            for i, aa in args['AASeq'].enumerate():
                fout.write('%4i %s %s   0.000  1.000  0.000\n' % (i, aa, 'H'))


def create_symm(args):
    cmd = 'python2.7 %s%s -symm_type %s -nsub %i > %s' % (rosetta_path,
                                           'main/source/src/apps/public/symmetry/make_symmdef_file_denovo.py',
                                           args['symm_type'], args['symm_num'], 'C%i.symm' % args['symm_num'])
    args['logger'].log('creating symmetry file using command: %s' % cmd)
    os.system(cmd)


def create_span(args):
    """
    wither copy span file, or create one from topo string file / entry
    """
    args['logger'].log('creating span file')
    if args['hetero']:
        with open(args['topo_string_file'], 'r') as fin:
            args['ts_dict'] = {}
            s = fin.read().split('>')
            ts = ''
            for a in s[1:]:
                s1 = a.split('\n')
                args['ts_dict'][s1[0]] = s1[1]
                ts += s1[1]
        ok = all([len(args['ts_dict'][name]) == len(args['seq'][name]) for name in args['seq'].keys()])
        assert ok, "topo_strings not same length as seqs"
        args['topo_string'] = ts
    else:
        if args['span_file'] is not None:
            args['logger'].log('copying span file from %s' % args['span_file'])
            shutil.copy(args['span_file'], args['path_data']+args['name']+'.span')
            return
        if args['topo_string'] is None and args['topo_string_file'] is not None:
            with open(args['topo_string_file'], 'r') as fin:
                args['logger'].log('reading topo string from %s' % args['topo_string'])
                args['topo_string'] = fin.read().rstrip()
        assert len(args['AASeq']) == len(args['topo_string']), "topo string and seq not same length"
    args['logger'].log('found topo string: %s' % args['topo_string'])
    topology_list = []
    hhh = re.compile('[hH]*')
    h_list = [(a.start(), a.end()) for a in hhh.finditer(args['topo_string']) if a.end() - a.start() > 1
              and a.end() - a.start() >= 1]
    for h in h_list:
        direction = 'fwd' if args['topo_string'][h[0] - 1] == '1' else 'rev'
        topology_list.append((h[0], h[1] - 1, direction))
    args['logger'].log('topology list %r' % topology_list)

    args['logger'].log('creating span file at %s' % args['path_data']+args['name']+'.span')
    with open(args['path_data']+args['name']+'.span', 'w+') as fout:
        fout.write('span file generated by prepare_fnd\n')
        if args['hetero']:
            fout.write('%i %i\n' % (len(topology_list), sum([len(a) for a in args['seqs'].values()])))
        else:
            fout.write('%i %i\n' % (len(topology_list), len(args['AASeq'].get_seq())))
        fout.write('antiparallel\n')
        fout.write('%s\n' % ('n2c' if topology_list[0][2] == 'fwd' else 'c2n'))
        for h in topology_list:
            fout.write('\t%i\t%i\n' % (h[0]+1, h[1]+1))


def make_or_read_fasta(args):
    if args['hetero']:
        args['seqs'] = read_seqs(args['fasta'])
        with open('%s.fasta' % args['name'], 'w+') as fout:
            for name, seq in args['seqs'].items():
                fout.write('>%s\n%s\n' % (name, seq.get_seq()))
            args['fasta'] = '%s.fasta' % args['name']
            return
    if args['fasta'] is None and args['seq'] is not None:
        args['AASeq'] = AASeq(args['seq'], args['name'])
        args['AASeq'].write_file(args['name']+'.fasta')
        args['fasta'] = args['name']+'.fasta'
        args['logger'].log('wrote fasta file to %s, the seq is %s' % (args['fasta'], args['AASeq'].get_seq()))
    elif args['fasta'] is not None and args['seq'] is None:
        args['AASeq'] = read_seq(args['fasta'])
        args['AASeq'].write_file(args['name'] + '.fasta')
        args['logger'].log('read seq from fasta at %s, seq is %s' % (args['fasta'], args['AASeq'].get_seq()))
    else:
        assert False, "neither fasta file nor seq provided"


def run_fragment_picker_new(args):
    args['logger'].log('creating fragment_picker_simple.wghts')
    with open('fragment_picker_simple.wghts', 'w+') as fout:
        fout.write('# score name          priority  wght   max_allowed  extras\n')
        fout.write('RamaScore               400    2.0     -       predA\n')
        fout.write('SecondarySimilarity     350     1.0     -       predA\n')
        fout.write('FragmentCrmsd             0     0.0     -\n')

    for name, aaseq in args['seqs'].items():
        for n in [3, 9]:
            args['logger'].log('creating fragment_picker flags for %s %i' % (name, n))
            with open('fragment_picker_%s_%i.flags' % (name, n), 'w+') as fout:
                fout.write('#input databases\n')
                fout.write('-database   %s\n' % args['database']) #'(rosetta_path+'main/database/'))
                fout.write('-in::file::vall  %stools/fragment_tools/vall.jul19.2011.gz\n' % rosetta_path)

                fout.write('# Query-related input files\n')
                fout.write('-in::file::fasta    %s\n' % (args['path_data']+name+'.fasta'))
                fout.write('-frags::ss_pred %s predA\n' % (args['path_data']+name+'.ss2'))
                fout.write('-s  %s\n' % args['some_pdb'])

                fout.write('# Weights file\n')
                fout.write('-frags::scoring::config %sfragment_picker_simple.wghts\n' % args['path_data'])
                fout.write('# What should we do?\n')
                fout.write('-frags::bounded_protocol\n')

                fout.write('# three-mers only, please\n')
                fout.write('-frags::frag_sizes %i\n' % n)
                fout.write('-frags::n_candidates 200\n')
                fout.write('-frags::n_frags 200\n')

                fout.write('# Output\n')
                fout.write('-out::file::frag_prefix frags_%i%s\n' % (n, name))
                fout.write('-frags::describe_fragments frags.fsc\n')

                for n in [3, 9]:
                    args['logger'].log('running fragment_picker for %s %i' % (name, n))
                    cmd = '%sfragment_picker.default.linuxgccrelease @fragment_picker_%s_%i.flags &' % (rosetta_bin_path, name, n)
                    args['logger'].log('running command:\n%s' % cmd)
                    os.system(cmd)


def run_fragment_picker(args):
    args['logger'].log('creating fragment_picker_simple.wghts')
    with open('fragment_picker_simple.wghts', 'w+') as fout:
        fout.write('# score name          priority  wght   max_allowed  extras\n')
        fout.write('RamaScore               400    2.0     -       predA\n')
        fout.write('SecondarySimilarity     350     1.0     -       predA\n')
        fout.write('FragmentCrmsd             0     0.0     -\n')

    for n in [3, 9]:
        args['logger'].log('creating fragment_picker flags for %i' % n)
        with open('fragment_picker_%i.flags' % n, 'w+') as fout:
            fout.write('#input databases\n')
            fout.write('-database   %s\n' % args['database']) #'(rosetta_path+'main/database/'))
            fout.write('-in::file::vall  %stools/fragment_tools/vall.jul19.2011.gz\n' % rosetta_path)

            fout.write('# Query-related input files\n')
            fout.write('-in::file::fasta    %s\n' % (args['path_data']+args['name']+'.fasta'))
            fout.write('-frags::ss_pred %s predA\n' % (args['path_data']+args['name']+'.ss2'))
            fout.write('-s  %s\n' % args['some_pdb'])

            fout.write('# Weights file\n')
            fout.write('-frags::scoring::config		%sfragment_picker_simple.wghts\n' % args['path_data'])

            fout.write('# What should we do?\n')
            fout.write('-frags::bounded_protocol\n')

            fout.write('# three-mers only, please\n')
            fout.write('-frags::frag_sizes 		%i\n' % n)
            fout.write('-frags::n_candidates		200\n')
            fout.write('-frags::n_frags		200\n')

            fout.write('# Output\n')
            fout.write('-out::file::frag_prefix		frags\n')
            fout.write('-frags::describe_fragments 	frags.fsc\n')
    for n in [3, 9]:
        args['logger'].log('running fragment_picker for %i' % n)
        cmd = '%sfragment_picker.default.linuxgccrelease @fragment_picker_%i.flags &' % (rosetta_bin_path, n)
        args['logger'].log('running command:\n%s' % cmd)
        os.system(cmd)

if __name__ == '__main__':
    main()
