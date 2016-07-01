#!/usr/bin/env python3.5
import os
from collections import OrderedDict

lab_root = '/home/labs/fleishman'


class Friend:
    def __init__(self, name, user_name):
        self.name = name
        self.log = parse_friend_log(user_name, self.name)

class User:
    def __init__(self, name: str, icon: str, friend_names: list, chat_file: str):
        self.name = name
        self.icon = icon
        self.friend_names = friend_names
        self.chat_file = chat_file
        self.friends = []
        self.initialize_friends()

    def __repr__(self):
        return "%s %s\nfriends: %s" % (self.name, self.icon, ', '.join(self.friends))

    def initialize_friends(self):
        for f_name in self.friend_names:
            friend = Friend(name=f_name, user_name=self.name)

    def read_friends(self):
        for friend in self.friends:
            if os.path.isfile('%s/%s/.pytexter_profile' % (lab_root, friend)):


def main():
    user = parse_profile(lab_root + '/' + 'jonathaw/.pytexter_profile')


def parse_friend_log(user_name, friend_name):
    result = OrderedDict()
    for l in open('%s/%s_log' % (friend_name, friend_name), 'r'):
        s = l.split(' ### ')
        if s != []:
            if s[0] == user_name:
                result.update({int(s[1]): s[2]})
    return result

def parse_profile(profile_file):
    result = {}
    for l in open(profile_file, 'r'):
        s = l.split()
        if s != []:
            result[s[0]] = s[1:]
    return User(name=result['name'][0], icon=result['icon'][0], friends=result['friends'])

if __name__ == '__main__':
    print(parse_profile(lab_root + '/' + 'jonathaw/.pytexter_profile'))
    main()