import subprocess

users = ['gideonla', 'benezer', 'ravitn', 'cnorn', 'jonatha', 'assafa']
results = {}

proc = subprocess.Popen(['bjobs', '-u', 'all'], stdout=subprocess.PIPE)
bjobs = proc.stdout.read()

bjobs_split = bjobs.split('\n')

others_fleish_run = 0
others_fleish_pnd = 0
others_newall_run = 0
others_newall_pnd = 0

for user in users:
    new_all_pend = 0
    new_all_run = 0
    fleish_pend = 0
    fleish_run = 0

    for line in bjobs_split:
        split = line.split()
        print split
        if user in split:
            if 'fleishman' in split:
                if 'RUN' in split:
                    fleish_run += 1
                elif 'PEND' in split:
                    fleish_pend += 1
            if 'new-all.q' in split:
                if 'RUN' in split:
                    new_all_run += 1
                elif 'PEND' in split:
                    new_all_pend += 1
        else:
            if 'fleishman' in split:
                if 'RUN' in split:
                    others_fleish_run += 1
                elif 'PEND' in split:
                    others_fleish_pnd += 1
            if 'new-all.q' in split:
                if 'RUN' in split:
                    others_newall_run += 1
                elif 'PEND' in split:
                    others_newall_pnd += 1

    results[user] = {'fleish_run': fleish_run, 'fleish_pend': fleish_pend, 'new-all_run': new_all_run,
                     'new-all_pend': new_all_pend}

results['others'] = {'fleish_run': others_fleish_run, 'fleish_pend': others_fleish_pnd,
                     'new-all_run': others_newall_run, 'new-all_pend': others_newall_pnd}

print '%13s%13s%13s%13s%13s' % ('NAME', 'Fleish run', 'Fleish pnd', 'New-all run', 'New-all pnd')
for user in results.keys():
    print '%13s%13s%13s%13s%13s' % (user, results[user]['fleish_run'], results[user]['fleish_pend'],
                                    results[user]['new-all_run'], results[user]['new-all_pend'], )