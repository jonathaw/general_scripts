#!/usr/bin/env python3.5
from LSFManager import am_i_running
__author__ = 'jonathan'

if am_i_running():
    print('cron is running')
else:
    print('cron is not running')