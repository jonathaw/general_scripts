#!/usr/bin/env python
def summarize_dir(dir_):
    import re
    from my_colors import my_colors
    my_colors = my_colors()
    file_list = os.listdir(dir_)
    results = {'jobs': 0, 'outs': 0, 'errs': 0, 'pdbs': 0, 'fastas': 0}
    dirs = []
    otrs = []
    for entry in file_list:
        if re.match('job\..*', entry):
            results['jobs'] += 1
        elif re.match('out\..*', entry):
            results['outs'] += 1
        elif re.match('err\..*', entry):
            results['errs'] += 1
        elif re.match('.*\.pdb', entry):
            results['pdbs'] += 1
        elif re.match('.*\.fasta', entry):
            results['fastas'] += 1
        elif os.path.isdir(entry):
            dirs.append(entry)
        else:
            otrs.append(entry)
    for k, v in results.items():
        if v != 0:
            print my_colors[my_colors.keys()[results.keys().index(k)]].format('Found %i %s' % (v, k))
    print 'Found dirs: ' + '\n'.join(x for x in dirs)
    print 'Found files: ' + '\n'.join(x for x in otrs)



if __name__ == '__main__':
    import os
    summarize_dir(os.getcwd())