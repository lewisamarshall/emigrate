"""Command line interface for emigrate."""
from .Solver import Solver
from .Frame import Frame
from .FrameSeries import FrameSeries
from .deserialize import deserialize

import sys
import os
from math import ceil
import click
from matplotlib import pyplot

try:
    import simplejson as json
except:
    import json


@click.group(name='emigrate', chain=True)
@click.pass_context
def cli(ctx):
    pass


@cli.command()
@click.pass_context
@click.argument('filename', type=click.Path(exists=True))
@click.option('--io', is_flag=True)
def load(ctx, filename, io):
    """Open an emgrate file and return a serialized frame."""
    _, file_extension = os.path.splitext(filename)
    if file_extension == '.hdf5':
        ctx.obj['frame_series'] = FrameSeries(filename=filename, mode='r')
    elif file_extension == '.json':
        with open(filename) as f:
            ctx.obj['frame'] = deserialize(f.read())
    else:
        raise RuntimeError("Can't load {} files.".format(file_extension))

    if io:
        for frame in iter(sys.stdin.readline, ''):
            click.echo(ctx.obj['frame_series'][int(frame)].serialize(compact=True))


@cli.command()
@click.pass_context
@click.argument('output', type=click.Path(exists=False))
@click.option('-f', '--frame', type=click.INT, default=None)
def plot(ctx, output, frame):
    if frame:
        ctx.obj['frame'] = ctx.obj['frame_series'][frame]
    if not ctx.obj['frame']:
        n = click.prompt('Frame', default=1, type=click.INT)
        ctx.obj['frame'] = ctx.obj['frame_series'][n]

    frame = ctx.obj['frame']
    for ion, ion_concentration in zip(frame.ions, frame.concentrations):
        pyplot.plot(frame.nodes, ion_concentration, '-', label=ion)
    pyplot.xlabel('x (mm)')
    pyplot.ylabel('concentration (M)')
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmax=frame.nodes[-1])
    pyplot.legend()
    pyplot.savefig(output, bbox_inches='tight')


@cli.command()
@click.pass_context
@click.option('-f', '--frame', type=click.INT, default=None)
def echo(ctx, frame):
    if frame:
        ctx.obj['frame'] = ctx.obj['frame_series'][frame]
    if not ctx.obj['frame']:
        n = click.prompt('Frame', default=1, type=click.INT)
        ctx.obj['frame'] = ctx.obj['frame_series'][n]

    click.echo(ctx.obj['frame'].serialize(compact=True))


@cli.command()
@click.pass_context
@click.option('-i', '--infile', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False), default=None)
def construct(ctx, infile, output):
    infile = click.format_filename(infile)
    with open(infile, 'r') as inputfile:
        constructor = deserialize(inputfile.read())
    ctx.obj['frame'] = Frame(constructor)
    if output:
        with open(output, 'w') as loc:
            loc.write(ctx.obj['frame'].serialize())


@cli.command()
@click.pass_context
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-t', '--time', type=float)
@click.option('-d', '--dt', type=float)
def solve(ctx, output, dt, time):
    solver = Solver(ctx.obj['frame'], filename=output)

    with click.progressbar(solver.iterate(dt, time),
                           length=int(ceil(time/dt)),
                           label='Solving...',
                           ) as bar:
        for frame in bar:
            pass

def close(ctx):
    if ctx.obj.get('frame_series', None):
        ctx.obj['frame_series'].hdf5.close()
        ctx.obj['frame_series'] = None

def main():
    cli(obj={'frame_series': None, 'frame': None})

if __name__ == '__main__':
    main()
