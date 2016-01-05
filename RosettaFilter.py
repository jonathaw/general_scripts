class Filter():
    def __init__(self, name: str, typ: str, threshold: float, limits: list, under_over: str,
                 g_name: str=None):
        self.filter_name = name
        self.filter_type = typ
        self.upper_limit = limits[1]
        self.lower_limit = limits[0]
        self.threshold = threshold
        self.under_over = under_over
        self.g_name = g_name if g_name is not None else typ
        self.failed = 0
        self.passed = 0
        self.max_tested = -100000
        self.min_tested = 100000
        self.all_seen = []

    def __str__(self):
        return 'name: %-11s type: %-9s [%i, %i] threshold %f %s %s' % (self.filter_name, self.filter_type, self.lower_limit,
                                                                  self.upper_limit, self.threshold, self.under_over,
                                                                  self.g_name)

    def __repr__(self) -> str:
        return 'name: %-11s type: %-9s [%i, %i] threshold %f %s %s' % (self.filter_name, self.filter_type, self.lower_limit,
                                                                  self.upper_limit, self.threshold, self.under_over,
                                                                  self.g_name)

    def pass_or_fail(self, score: float, verbose: bool=False) -> bool:
        """
        :param score: a score, float or int
        :return: True if passes the threshold. False if not
        >>> filter = Filter('a_sasa', 'sasa', 1300, [0, 10000], 'over')
        >>> filter.pass_or_fail(1000)
        False
        >>> filter.pass_or_fail(1500)
        True
        """
        self.all_seen.append(score)
        self.max_tested = max([self.max_tested, score])
        self.min_tested = min([self.min_tested, score])
        if self.under_over == 'under':
            if self.threshold >= score:
                self.passed += 1
                if verbose:
                    print('passed %s, threhsold %f, with %f' % (self.filter_type, self.threshold, score))
                return True
            else:
                self.failed += 1
                return False
        elif self.under_over == 'over':
            if self.threshold <= score:
                self.passed += 1
                return True
            else:
                self.failed += 1
                return False

    def within_limits(self, score: dict) -> bool:
        """
        :param score: a score
        :return: True if score is within the filters limits
        >>> filter = Filter('a_sasa', 'sasa', 1300, [0, 10000], 'over')
        >>> filter.within_limits(1000)
        True
        >>> filter.within_limits(100000)
        False
        """
        return self.lower_limit <= score <= self.upper_limit

    def set_threshold(self, threshold) -> None:
        self.threshold = threshold


class RunFilters():
    def __init__(self):
        """

        """
        self.filters = {}

    def __repr__(self) -> str:
        return '%s\n'.join([str(a) for a in self.filters])

    def __str__(self) -> str:
        return '%s\n'.join([str(a) for a in self.filters])

    def __len__(self) -> int:
        return len(self.filters)

    def __getitem__(self, item):
        return self.filters[item]

    def items(self) -> (str, Filter):
        for k, v in self.filters.items():
            yield k, v

    def keys(self):
        return self.filters.keys()

    def values(self):
        return self.filters.values()

    def report(self) -> str:
        msg = 'report\n'
        for k, flt in self.items():
            msg += 'Filter %s, threshold %5.f, limits [%5.f, %5.f], passed %i, failed %i, highest %5.f, lowest %5.f\n' % \
                   (flt.filter_name, flt.threshold, flt.lower_limit, flt.upper_limit, flt.passed, flt.failed,
                    flt.max_tested, flt.min_tested)
        return msg

    def append_filter(self, filter: Filter) -> None:
        self.filters[filter.filter_type] = filter

    def set_thresholds(self, thresholds) -> None:
        for flt in self.filters:
            flt.set_threshold(thresholds[flt.filter_type])

    def test_all(self, score: dict, verbose=False) -> (bool, str):
        """
        :param score: a score dictionary {'filter_name': score}
        :return: True/False if all filters pass and lists of passed and failed filters
        """
        tests = []
        msg = False
        for flt in self.values():
            if flt.filter_type == 'rmsd' and 'rmsd' not in score.keys():  # or flt.filter_type == 'hbonds':
                # msg = 'DID NOT CONSIDER %s !!!!! ' % flt.filter_type.upper()
                continue
            tests.append(flt.pass_or_fail(score[flt.filter_type], verbose))
        return all(tests), msg


def score2dict(file_name: str, verbose=False) -> dict:
    have_fields = False
    results = {}
    with open(file_name, 'r') as fin:
        cont = fin.read().split('\n')
    for l in cont:
        s = l.split()
        if len(s) < 2:
            continue
        if (s[1] == 'total_score' or s[1] == 'score') and not have_fields:
            fields = {a if a[:2] != 'a_' else a[2:]: i for i, a in enumerate(s) if a != 'rms'}
            try:
                fields['rmsd'] = s.index('a_rms')
                fields.pop('rms')  # remove the a_rms, experimental
            except:
                if verbose:
                    print('No rmsd')
            have_fields = True
        elif (s[0] == 'SCORE:' or 'SCORE' in s[0]) and 'score' not in s[1]:
            if len(s) != len(list(fields.keys()))+1 if 'rms' in s else 0:  # adding 1 because I remove rms
                continue
            results[s[fields['description']]] = {a: float(s[i]) for a, i in fields.items()
                                                 if a not in ['SCORE:', 'description'] and 'SCORE' not in a}
            if 'score' not in results[s[fields['description']]].keys():
                results[s[fields['description']]]['score'] = results[s[fields['description']]]['total_score']
            elif 'total_score' not in results[s[fields['description']]].keys():
                results[s[fields['description']]]['total_score'] = results[s[fields['description']]]['score']
            results[s[fields['description']]]['description'] = s[fields['description']]
    return results