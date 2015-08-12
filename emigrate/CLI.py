"""Command line interface for emigrate."""
from .Solver import Solver
from .Frame import Frame
from .FrameSeries import FrameSeries

# import sys
import click


@click.group(name='emigrate', chain=True)
@click.pass_context
def cli(ctx):
    pass


@cli.command()
@click.pass_context
@click.argument('filename', type=click.Path(exists=True))
@click.option('--frame', '-f', prompt=True, type=click.INT)
def open(ctx, filename, frame):
    """Open an emgrate file and return a serialized frame."""
    if ctx.obj:
        close(ctx)
    ctx.obj = FrameSeries(filename=filename, mode='r')
    # click.echo(ctx.obj[frame].serialize())

@cli.command()
@click.pass_context
@click.option('--frame', '-f', prompt=True, type=click.INT)
def plot(ctx, frame):
    print 'plotting'
    print ctx.obj

def close(ctx):
    if ctx.obj:
        ctx.obj.hdf5.close()
        ctx.obj = None

#     def prompt(self):
#         for command in iter(sys.stdin.readline, ''):
#             command = command.strip().split()
#             if command[0] == 'open':
#                 self.open(command[1])
#             elif command[0] == 'frame':
#                 self.frame(int(command[1]))
#             elif command[0] == 'close':
#                 self.close()
#             elif command[0] == 'exit':
#                 break
#             else:
#                 sys.stderr.write('unknown command')
#                 sys.stderr.flush()

#         print json.dumps(serial)
#         sys.stdout.flush()
#
#     def listener(self):
#         for command in iter(sys.stdin.readline, ''):
#             command = command.strip().split()
#             if command[0] == 'open':
#                 self.open(command[1])
#             elif command[0] == 'frame':
#                 self.frame(int(command[1]))
#             elif command[0] == 'close':
#                 self.close()
#             elif command[0] == 'exit':
#                 break
#             else:
#                 sys.stderr.write('unknown command')
#                 sys.stderr.flush()

if __name__ == '__main__':
    cli()
