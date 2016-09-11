#!/usr/bin/env python
import os
import sys


prefices = ['out.', 'err.', 'job.']

def main():
    if len(sys.argv) < 2:
        pwd = './'
    else:
        pwd = os.getcwd() + '/' + sys.argv[1]
    pwd += '/' if pwd[-1] != '/' else ''
    all_files = os.listdir(pwd)
    results = {a: 0 for a in prefices}
    dirs = []
    for file_ in all_files:
        if file_ == '': continue
        if os.path.isdir(pwd + file_):
            dirs.append(file_)
            continue
        s = file_.split('.')
        # print(file_, s)
        if len(s) == 2:
            pref = s[0] + '.'
            suff = '.' + s[1]
            if pref in prefices:
                results[ pref ] += 1
            else:
                if suff not in results.keys():
                    results[ suff ] = 0
                results[ suff ] += 1
        else:
            if s[0] not in results.keys():
                results[ s[0] ] = 0
            results[ s[0] ] += 1

    i = 1.
    msg = ''
    for k, v in results.items():
        msg += '%s: %i\t' % (k, v)
        if i / 3 == 0:
            msg += '\n'
        i += 1.0
    print(msg)



if __name__ == '__main__':
    main()


# def summarize_dir(dir_):
    # import re
    # from my_colors import my_colors
    # my_colors = my_colors()
    # file_list = os.listdir(dir_)
    # results = {'jobs': 0, 'outs': 0, 'errs': 0, 'pdbs': 0, 'fastas': 0}
    # dirs = []
    # otrs = []
    # for entry in file_list:
        # if re.match('job\..*', entry):
            # results['jobs'] += 1
        # elif re.match('out\..*', entry):
            # results['outs'] += 1
        # elif re.match('err\..*', entry):
            # results['errs'] += 1
        # elif re.match('.*\.pdb', entry):
            # results['pdbs'] += 1
        # elif re.match('.*\.fasta', entry):
            # results['fastas'] += 1
        # elif os.path.isdir(entry):
            # dirs.append(entry)
        # else:
            # otrs.append(entry)
    # for k, v in results.items():
        # if v != 0:
            # print my_colors[my_colors.keys()[results.keys().index(k)]].format('Found %i %s' % (v, k))
    # print 'Found dirs: ' + '\n'.join(x for x in dirs)
    # print 'Found files: ' + '\n'.join(x for x in otrs)



# if __name__ == '__main__':
    # import os
    # summarize_dir(os.getcwd())
