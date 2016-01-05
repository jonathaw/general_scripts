"""
binding data from cohesin microarray
"""
def binding_data() -> dict:
    """
    :return: {'coh_name': {'doc_name': True/False}}
    """
    data = {'ca_a5': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ca_a2': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ca_a1': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'cc_ox': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },

            'ct_s1': {'ac_9b': False, 'ac_scaa': True, 'ac_scab': False, 'af_doc': False, 'bc_48a': True, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': True, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ct_2p2': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': True, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ct_ob4': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': True, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ct_ob1': {'ac_9b': False, 'ac_scaa': True, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': True, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },

            'cc_a8': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': True, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'cc_a1': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': True, 'ct_11b': True, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'rf_e': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': True, },
            'rf_c': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },

            'ct_oa': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': True, 'ct_48b': True, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ct_a3': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': True, 'ct_48b': True, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ct_a2': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': True, 'ct_48b': True, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ct_a1': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': True, 'ct_48b': True, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },

            'rf_b6': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': True, 'rf_scab': False, },
            'rf_b1': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': True, 'rf_scab': False, },
            'rf_a3': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'af_76': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': True, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },

            'bc_b3': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': True,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'bc_a11': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': True, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'bc_a5': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': True, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ac_d3': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },

            'af_75': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': True, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },

            'ac_d1': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ac_c3': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': True, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ac_b1': {'ac_9b': False, 'ac_scaa': True, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, },
            'ac_a5': {'ac_9b': False, 'ac_scaa': False, 'ac_scab': False, 'af_doc': False, 'bc_48a': False, 'bc_scaa': False,
                      'ca_9e': False, 'cc_5a': False, 'ct_11b': False, 'ct_48b': False, 'ct_cipa': False, 'ctcipadoc': False,
                      'rf_44b': False, 'rf_scaa': False, 'rf_scab': False, }
    }
    return data

if __name__ == '__main__':
    data = binding_data()
    count_trues = 0
    count_all = 0
    for k, v in data.items():
        for k1, v1 in v.items():
            count_trues +=1 if v1 else 0
            count_all += 1
    print(list(data.keys()))
    print('found %i trues' % count_trues)
    print('found a total of %i pairs' % count_all)
