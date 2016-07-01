def main(args):
    with open(args['job_maker_log'], 'r') as log_file:
        log_cont = log_file.read()
    log_split = log_cont.split('name.for.files')
    results = []
    for log in log_split:
        coh_name = ''
        chosens = []
        for line in log.split('\n'):
            line_split = line.split(' ')
            if len(line_split) <= 1:
                continue
            if line_split[0] == 'coh':
                if line_split[1] in [a.keys()[0] for a in results]:
                    break
                coh_name = line_split[1]
            if line_split[1] == 'chosen':
                chosens.append((line_split[0], int(line_split[5])))
        if coh_name != '':
            results.append({coh_name: chosens})
    return results



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-job_maker_log')
    args = vars(parser.parse_args())
    results = main(args)
    for ent in results:
        if len(ent.values()[0]) > 1:
            print ent.keys()[0] + ' ' + ' '.join(['%s %i' % (a[0], a[1]) for a in ent.values()[0]])
        else:
            print ent.keys()[0] + ' HAD NO BBs !!!!!7'