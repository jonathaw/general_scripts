#!/usr/bin/env python3.5
__author__ = 'jonathan'


class lsf_user():
    def __init__(self, name):
        global QUEUES, MODES
        # fleishman-priority only has the - at the end...
        QUEUES = ['fleishman', 'new-all.q', 'fleishman-']
        MODES = ['RUN', 'PEND', 'SSUSP', 'UNKWN', 'PSUSP', 'USUSP']
        self.name = name
        self.results = {}
        for queue in QUEUES:
            self.results[queue] = {}
            for mode in MODES:
                self.results[queue][mode] = 0

    def add_to_count(self, queue, mode):
        self.results[queue][mode] += 1

    def print_results(self):
        msg = '%-9s' % self.name
        for queue in QUEUES:
            for mode in MODES:
                msg += '%-6s' % self.results[queue][mode]
        return msg


def main():
    import subprocess
    global QUEUES, MODES
    QUEUES = ['fleishman', 'new-all.q', 'fleishman-priority']
    MODES = ['RUN', 'PEND', 'SSUSP', 'UNKWN', 'PSUSP']
    colorred = "\033[01;31m{0}\033[00m"
    coloryellow = "\033[01;20m{0}\033[00m"
    colorgreen = "\033[01;32m{0}\033[00m"
    colorblue = "\033[01;34m{0}\033[00m"
    colorpurple = "\033[01;35m{0}\033[00m"
    colorlightpurple = "\033[01;36m{0}\033[00m"
    proc = subprocess.Popen(['bjobs', '-u', 'all'], stdout=subprocess.PIPE)
    bjobs = proc.stdout.read()
    bjobs_split = bjobs.split('\n')
    all_users = {}
    all_users_names = []
    for line in bjobs_split:
        line_split = line.split()
        if len(line_split) > 3 and line_split[3] in QUEUES:
            if line_split[1] not in all_users_names:
                all_users[line_split[1]] = lsf_user(line_split[1])
                all_users_names.append(line_split[1])
            all_users[line_split[1]].add_to_count(line_split[3], line_split[2])

    # make sure fleishman-priority is named correctly
    print('%-9s' % '' + ''.join("%-30s" % queue if queue != 'fleishman-' else 'fleishman-priority' for queue in QUEUES))
    print('%-9s' % 'NAME' + ''.join(''.join("%-6s"  % mode for mode in MODES) for queue in QUEUES))

    for user in all_users:
        if user == 'jonatha':
            print(colorred.format(all_users[user].print_results()))
        elif user == 'elazara':
            print(colorlightpurple.format(all_users[user].print_results()))
        elif user == 'gideonl':
            print(colorpurple.format(all_users[user].print_results()))
        # elif user == 'orim':
        #     print colorlightpurple.format(all_users[user].print_results())
        elif user == 'longo':
            print(colorpurple.format(all_users[user].print_results()))
        else:
            print(all_users[user].print_results())


if __name__ == '__main__':
    main()