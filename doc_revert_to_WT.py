#!/usr/bin/env python
"""
a script to take a design and revert it to WT in any position not allowed for design due to Ca binding loop etc.
"""
def main():
    # import sys
    wt_dict = {
        '4uyq': {'seq': 'KVVPKFIYGDVDGNGSVRINDAVLIRDYVLGKINEFPYEYGMLAADVDGNGSIKSIDAVLVRDYVLGKIFLFPVEEKEE',
                 'design': {22: 'V', 25: 'R', 26: 'D', 29: 'L', 58: 'V', 61: 'R', 65: 'L'},
                 'pssm': {18: 'I', 29: 'L', 54: 'S', 65: 'L'}},
        '4uyp': {'seq': 'KVVPKFIYGDVDGNGSVRSIDAVLIRDYVLGKINEFPYEYGMLAADVDGNGSIKINDAVLVRDYVLGKIFLFPVEEKEE',
                 'design': {18: 'S', 22: 'V', 25: 'R', 26: 'D', 29: 'L', 58: 'V', 61: 'R', 65: 'L'},
                 'pssm': {18: 'S', 29: 'L', 54: 'I', 65: 'L'}},
        '5new': {'seq': 'VYGDLDGDGEVDVFDLILMRKAVENGDTERFEAADLNCDGVIDSDDLTYHSEYLHGIRKTLPVEY',
                 'design': {13: 'F', 16: 'I', 19: 'R', 20: 'K', 23: 'E', 44: 'D', 47: 'T', 48: 'Y', 50: 'S',
                            51: 'E'}, 'pssm': {12: 'V', 43: 'S', 54: 'H'}},
        '1ohz': {'seq': 'ESSSVLLGDVNGDGTINSTDLTMLKRSVLRAITLTDDAKARADVDKNGSINSTDVLLLSRYLLRVIDKFPVAENP',
                 'design': {21: 'T', 24: 'K', 25: 'R', 52: 'T', 55: 'L', 58: 'S', 59: 'R'},
                 'pssm': {17: 'S', 28: 'L', 51: 'S', 62: 'L'}},
        '2vn5': {'seq': 'VIVYGDYNNDGNVDALDFAGLKKYIMAADHAYVKNLDVNLDNEVNSTDLAILKKYLLGMVSKLPSN',
                 'design': {15: 'L', 18: 'A', 21: 'K', 22: 'K', 49: 'A', 52: 'K', 53: 'K', 56: 'L'},
                 'pssm': {14: 'A', 25: 'M', 45: 'S'}},
        '3ul4': {'seq': 'VVLNGDLNRNGIVNDEDYILLKNYLLRGNKLVIDLNVADVNKDGKVNSTDCLFLKKYILGLITI',
                 'design': {15: 'E', 18: 'I', 21: 'K', 22: 'N', 51: 'L', 54: 'K', 55: 'K'},
                 'pssm': {14: 'D', 17: 'Y', 25: 'L', 47: 'S', 58: 'L'}},
        '4dh2': {'seq': 'WNKAVIGDVNADGVVNISDYVLMKRYILRIIADFPADDDMWVGDVNGDNVINDIDCNYLKRYLLHMIREFPKNSYNSAPTF',
                 'design': {17: 'S', 20: 'V', 23: 'K', 24: 'R', 56: 'N', 59: 'K', 60: 'R'},
                 'pssm': {16: 'I', 19: 'Y', 27: 'L', 52: 'D', 63: 'L'}},
    }

    for k, v in wt_dict.items():
        # print k, v
        for d1, d2 in v['design'].items():
            # print 'design', d1, d2
            assert d2 == v['seq'][d1], "problem at design %s %s %s %i" % (k, d2, v['seq'][d1], d1)
        for d1, d2 in v['pssm'].items():
            # print 'pssm', d1, d2
            assert d2 == v['seq'][d1], "problem at pssm %s %s %s %i" % (k, d2, v['seq'][d1], d1)

    design_dict = {
        'ac21': {'seq': 'KVVPKFIYGDVNGNGRVTADDAAAVREYYLGKINEFPYEYGMLAADVDGNGSIKINDADLIAEYAAGSITLFPVEEKEE', 'wt': '4uyp'},
        'ac26': {'seq': 'KVVPKFIYGDVNGNGRVTSDDAALILAYVLGWINEFPYEYGMLAADVDGNGSIKLADADLVAKYAQGRLTLFPVEEKEE', 'wt': '4uyq'},
        'ac31': {'seq': 'VYGDLDGDGEINSFDLKLTKEAVIGKDTERFEAADLACEGTINALDAALHLAYLSGQLTSLPYRY', 'wt': '5new'},
        'ac39': {'seq': 'VIVYGDYNNDGNVDALDYAGLKKFILTSDHAYVKNADTNNNNSVNAVDLANLRSYLLGMVSKLPSN', 'wt': '2vn5'},
        'ac34': {'seq': 'VYGDLDGDGEVDVFDLILTKKATEGRDTERFEAADLACDGTINSLDLALHRAYLSGHLTSLPTAY', 'wt': '5new'},
        'ac53': {'seq': 'KVVPKFIYGDVDGNGSVRSIDAVLVKKYVLGKITEFPYEYGMLAADVNGSGKVNSSDAALIEEYVLGKIFLFPVEEKEE', 'wt': '4uyp'},
        'ac44': {'seq': 'WNKAVIGDVNANGRVDDSNLVLVQQYLLKQIADFPADDDMWVGDVNGDNVINSIDLAYVEKYLLHLITEFPKNSYNSAPTF', 'wt': '4dh2'},
        'ac42': {'seq': 'VVLNGDLNRNGIVNDEDLLLLKKYLLQGNNNVIDLNVADVTNNGKIDASNVAALRQYILGLITI', 'wt': '3ul4'},
        'ac36': {'seq': 'ESSSVLLGDVNGDGTVNAADLALLAKHVARKRTLTDDAKARADVNKNGKINSADIAALAAILANLRDKFPVAENP', 'wt': '1ohz'},
        'ac50': {'seq': 'VYGDLDGDGEIDEFDLILTRKAVEGRDTERFEAADLACEGTINALDAALHAAYLSGQLTTLPYRY', 'wt': '5new'},
    }
    for design, details in design_dict.items():
        wt_seq = str(wt_dict[details['wt']]['seq'][:])
        ds_seq = str(details['seq'][:])
        temp_seq = wt_seq
        pssm = parse_pssm(details['wt'], wt_seq)
        print 'ps seq', ''.join([i['type'] for i in pssm.values()])
        print 'wt_seq', wt_seq
        print 'ds_seq', ds_seq, design
        print pssm
        assert len(wt_seq) == len(ds_seq), 'wt %i, ds %i not equal' % (len(wt_seq), len(ds_seq))
        assert len(pssm.keys()) == len(ds_seq), 'at %s pssm has %i, design seq has %i' % (design, len(pssm.keys()), len(ds_seq))
        changes = []
        differences_tot = 0
        for i, aa in enumerate(wt_seq):
            # if the position is allowed for design, replace
            if i in wt_dict[details['wt']]['design'].keys() and ds_seq[i] != aa:
                print 'making a change', i, aa, details['seq'][i]
                temp_seq = replace_by_index(temp_seq, i, ds_seq[i])
                changes.append(i)
            # if the position is allowed for pssm, check:
            if i in wt_dict[details['wt']]['pssm'].keys():
                # if it's the same as in the WT, no change
                if wt_seq[i] == ds_seq[i]:
                    continue
                # if its PSSM is non-negative, allow replacement (design)
                elif pssm[i][ds_seq[i]] >= 0:
                    temp_seq = replace_by_index(temp_seq, i, ds_seq[i])
                    changes.append(i)
                else:
                    continue
                    # elif i in wt_dict[design]['pssm']
            differences_tot += 1 if aa != ds_seq[i] else 0
        print "for design", design
        print 'WT', wt_dict[details['wt']]['seq']
        print 'ne', temp_seq
        print 'dz', details['seq']

        changes_string = ''.join(['-'] * len(temp_seq))
        for i in changes:
            changes_string = replace_by_index(changes_string, i, '^')
        for i in wt_dict[details['wt']]['pssm'].keys():
            changes_string = replace_by_index(changes_string, i, 'p')
        print 'ch', changes_string
        print 'changes located at', changes, len(changes)
        print 'there were %i differences design~WT. %i were accepted' % (differences_tot, len(changes))

def parse_pssm(wt_name, seq):
    with open('/home/labs/fleishman/jonathaw/doc_analysis_16Aug/' + wt_name + '.pssm') as f:
  #  with open('/Volumes/labs/fleishman/jonathaw/doc_analysis_16Aug/' + wt_name + '.pssm') as f:
        cont = f.read().split('\n')
    print 'reading %s pssm' % wt_name
    result = {}
    aas = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    for l in cont:
        s = l.split()
        if s == []:
            continue
        # print s
        try:
            num = int(s[0])-1
        except:
            continue
        result[num] = {}
        result[num]['type'] = s[1]
        assert s[1] == seq[num], "at pos %i, pssm %s, but wt_sew %s" % (num, s[1], seq[num])
        for i, sc in enumerate(s[2:22]):
            result[num][aas[i]] = int(sc)
    return result


def replace_by_index(s, i, r):
    """
    :param s:string
    :param i:index
    :param r:char to replace at i
    :return:s with s[i] as r
    >>> s = 'abcdefg'
    >>> i = 3
    >>> r = 'X'
    >>> replace_by_index(s, i, r)
    'abcXefg'
    """
    return str(s[:i] + r + s[i + 1:])


if __name__ == '__main__':
    main()
    # parse_pssm('4uyp')
