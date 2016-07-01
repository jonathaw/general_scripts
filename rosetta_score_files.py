#!/usr/bin/env python
def score_passes_thresholds(score_line, fields_dict, thresholds_dict):
    score_line = score_line.split()
    if float(score_line[fields_dict['ddg']]) < thresholds_dict['ddg']\
            and float(score_line[fields_dict['sasa']]) > thresholds_dict['sasa'] \
            and float(score_line[fields_dict['pack']]) > thresholds_dict['pack']\
            and float(score_line[fields_dict['shape']]) > thresholds_dict['shape']\
            and float(score_line[fields_dict['buried']]) <= thresholds_dict['buried']:
        return True
    else:
        return False


def how_many_purples_in_file(file_name):
    thresholds = {'ddg': -20.0, 'sasa': 1300.0, 'pack': 0.6, 'shape': 0.5, 'buried': 1.0}
    filter_fields = {}
    purples_num = 0
    score_file = open(file_name, 'r')
    for file_line in score_file:
        line_split = file_line.split()
        if 'SCORE:' in line_split:
            if line_split[1] == 'score':
                filter_fields['ddg'] = line_split.index('a_ddg')
                filter_fields['sasa'] = line_split.index('a_sasa')
                filter_fields['pack'] = line_split.index('a_packstat')
                filter_fields['shape'] = line_split.index('a_shape')
                filter_fields['buried'] = line_split.index('a_buried_2')
            else:
                if score_passes_thresholds(file_line, filter_fields, thresholds):
                    purples_num += 1
    score_file.close()
    return purples_num


def score2dict(file_name):
    thresholds = {'ddg': -20.0, 'sasa': 1300.0, 'pack': 0.6, 'shape': 0.5, 'buried': 1.0}
    filter_fields = {}
    results = {}
    score_file = open(file_name, 'r')
    for file_line in score_file:
        line_split = file_line.split()
        if 'SCORE:' in line_split:
            if line_split[1] == 'score':
                filter_fields['ddg'] = line_split.index('a_ddg')
                filter_fields['sasa'] = line_split.index('a_sasa')
                filter_fields['pack'] = line_split.index('a_packstat')
                filter_fields['shape'] = line_split.index('a_shape')
                filter_fields['buried'] = line_split.index('a_buried_2')
                filter_fields['description'] = line_split.index('description')
            else:
                temp = {line_split[filter_fields['description']]:
                            {'ddg': float(line_split[filter_fields['ddg']]),
                            'sasa': float(line_split[filter_fields['sasa']]),
                            'pack': float(line_split[filter_fields['pack']]),
                            'shape': float(line_split[filter_fields['shape']]),
                            'buried': float(line_split[filter_fields['buried']]),
                            'description': line_split[filter_fields['description']],
                            'purple': score_passes_thresholds(file_line, filter_fields, thresholds)}}
                results.update(temp)
    score_file.close()
    return results


def scored_dict_passes_thresholds(score_dict_entry, thresholds_dict):
    if score_dict_entry['ddg'] < thresholds_dict['ddg'] and \
        score_dict_entry['sasa'] > thresholds_dict['sasa'] and \
        score_dict_entry['pack'] > thresholds_dict['pack'] and \
        score_dict_entry['shape'] > thresholds_dict['shape'] and \
        score_dict_entry['ddg'] < thresholds_dict['ddg'] and \
        score_dict_entry['buried'] < thresholds_dict['buried']:
        return True
    else:
        return False


if __name__ == '__main__':
    import sys
    print how_many_purples_in_file(sys.argv[1])