"""Command line interface for emigrate."""
# import emigrate
from .Solver import Solver
from .Frame import Frame
from .FrameSeries import FrameSeries

import sys
import json


class CLI(object):
    """Command line interface for emigrate."""

    frames = None

    def __init__(self):
        self.listener()
        return None

    def open(self, open_file):
        if self.frames:
            self.close()
        self.frames = FrameSeries(filename=open_file, mode='r')

    def close(self):
        if self.frames:
            self.frames.hdf5.close()
            self.frames = None

    def frame(self, frame):
        serial = self.frames[frame].serialize()
        serial['n_electrolytes'] = \
            len(self.frames.frames.keys())
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
