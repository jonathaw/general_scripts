#!/usr/bin/env python3.5
"""
scripts to test proteins for MS experiments
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from seq_funcs import read_multi_fastas


flanks = {'coh': {'cbm': 'MANTPVSGNLKVEFYNSNPSDTTNSINPQFKVTNTGSSAIDLSKLTLRYYYTVDGQKDQTFWCDHAAIIGSNGSYNGITSNVKGTFVKMSSST'
                         'NNADTYLEISFTGGTLEPGAHVQIQGRFAKNDWSNYTQSNDYSFKSASQFVEWDQVTAYLNGVLVWGKEPGGSVVPSTQPVTTPPATTKPPAT'
                         'TIPPSDDPNAGS',
                  '1ohz': {'start': 'D', 'end': 'NAT'}},
          'doc': {'xyn': 'MSHHHHHHKNADSYAKKPHISALNAPQLDQRYKNEFTIGAAVEPYQLQNEKDVQMLKRHFNSIVAENVMKPISIQPEEGKFNFEQADRIVKFA'
                         'KANGMDIRFHTLVWHSQVPQWFFLDKEGKPMVNECDPVKREQNKQLLLKRLETHIKTIVERYKDDIKYWDVVNEVVGDDGKLRNSPWYQIAGI'
                         'DYIKVAFQAARKYGGDNIKLYMNDYNTEVEPKRTALYNLVKQLKEEGVPIDGIGHQSHIQIGWPSEAEIEKTINMFAALGLDNQITELDVSMY'
                         'GWPPRAYPTYDAIPKQKFLDQAARYDRLFKLYEKLSDKISNVTFWGIADNHTWLDSRADVYYDANGNVVVDPNAPYAKVEKGKGKDAPFVFGP'
                         'DYKVKPAYWAIIDHKVVP',
                  '1ohz': {'start': 'ESSSVLL', 'end': 'RVIDKFPVAENP'},
                  '2vn5': {'start': 'V', 'end': 'SKLPSN'},
                  '3ul4': {'start': 'V', 'end': ''},
                  '4dh2': {'start': 'WNK', 'end': 'NSAPTF'},
                  '5new': {'start': '', 'end': 'Y'}}
          }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-coh_files', nargs='+')
    parser.add_argument('-doc_files', nargs='+')
    parser.add_argument('-threshold', type=float, default=30.0)
    args = vars(parser.parse_args())

    cohs = OrderedDict({k: v for a in args['coh_files'] for k, v in read_multi_fastas(a).items()})
    docs = OrderedDict({k: v for a in args['doc_files'] for k, v in read_multi_fastas(a).items()})

    mw_df = pd.DataFrame(columns=['coh_seq', 'coh_weight']+list(docs.keys()))
    for i, coh in enumerate(cohs.values()):
        coh.add_prefix = flanks['coh']['cbm']
        coh_weight = coh.calc_molecular_weight()
        print(coh.get_seq(), coh_weight)
        weights_combined = []
        for doc in docs.values():
            doc.add_prefix(flanks['doc']['xyn'])
            weights_combined.append(coh_weight+doc.calc_molecular_weight())
        mw_df.loc[coh.name] = [coh, coh_weight] + weights_combined
    print(mw_df)

    diffs = []
    for coh1 in cohs.keys():
        for doc1 in docs.keys():
            coh_doc_1 = mw_df[doc1][coh1]
            for coh2 in cohs.keys():
                for doc2 in docs.keys():
                    if coh1 != coh2 and doc1 != doc2:
                        diff = abs(coh_doc_1 - mw_df[doc2][coh2])
                        diffs.append(diff)
                        if diff <= args['threshold']:
                            print('%s %s and %s %s have a weight difference of only %f' % (coh1, doc1, coh2, doc2, diff))

    plt.boxplot(diffs)
    plt.show()


if __name__ == '__main__':
    main()
