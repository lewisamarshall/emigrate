# import emigrate
from Migrate import Migrate
from data_structure import Electrolyte, Electromigration

import sys
import json


class CLI(object):

    electromigration = None

    def __init__(self):
        self.listener()
        return None

    def open(self, file):
        if self.electromigration:
            self.close()
        self.electromigration = Electromigration(filename=file, mode='r')

    def close(self):
        if self.electromigration:
            self.electromigration.hdf5.close()
            self.electromigration = None

    def frame(self, frame):
        serial = self.electromigration[frame].serialize()
        serial['n_electrolytes'] = \
            len(self.electromigration.electrolytes.keys())
        print json.dumps(serial)
        sys.stdout.flush()

    def listener(self):
        for command in iter(sys.stdin.readline, ''):
            command = command.strip().split()
            if command[0] == 'open':
                self.open(command[1])
            elif command[0] == 'frame':
                self.frame(int(command[1]))
            elif command[0] == 'close':
                self.close()
            elif command[0] == 'exit':
                break
            else:
                sys.stderr.write('unknown command')
                sys.stderr.flush()

if __name__ == '__main__':
    CLI()
