"""Command line interface for emigrate."""
from .Solver import Solver
from .Frame import Frame
from .FrameSeries import FrameSeries

import sys
import click
from matplotlib import pyplot
from math import ceil

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
def open(ctx, filename, io):
    """Open an emgrate file and return a serialized frame."""
    ctx.obj['frame_series'] = FrameSeries(filename=filename, mode='r')

    if io:
        for frame in iter(sys.stdin.readline, ''):
            click.echo(ctx.obj['frame_series'][int(frame)].serialize())


@cli.command()
@click.pass_context
@click.argument('output', type=click.Path(exists=False))
def plot(ctx, output):
    if not ctx.obj['frame']:
        n = click.prompt('Frame', default=1, type=click.INT)
        ctx.obj['frame'] = ctx.obj['frame_series'][n]

    pyplot.figure()
    frame = ctx.obj['frame']
    for ion, ion_concentration in zip(frame.ions, frame.concentrations):
        pyplot.plot(frame.nodes, ion_concentration, '-', label=ion)
    pyplot.xlabel('x (mm)')
    pyplot.ylabel('concentration (M)')
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmax=frame.nodes[-1])
    pyplot.legend()
    pyplot.savefig('output', bbox_inches='tight')


@cli.command()
@click.pass_context
def echo(ctx):
    if not ctx.obj['frame']:
        n = click.prompt('Frame', default=1, type=click.INT)
        ctx.obj['frame'] = ctx.obj['frame_series'][n]

    click.echo(ctx.obj['frame'].serialize())


@cli.command()
@click.pass_context
@click.option('-i', '--input', type=click.Path(exists=True))
@click.option('-o', '--output', type=click.Path(exists=False))
def construct(ctx, input, output):
    constructor = deserialize(input)
    ctx.obj['frame'] = Frame(constructor)
    if output:
        with open(output, 'w') as loc:
            json.dump(ctx.obj['frame'], loc)


@cli.command()
@click.pass_context
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-t', '--time', type=click.Path(exists=True))
@click.option('-d', '--dt', type=click.Path(exists=False))
def solve(ctx, output, dt, time):
    solver = Solver(ctx.obj['frame'], filename=output)

    with click.progressbar(length=time,
                           label='Solving...'
                           ) as bar:
        for frame in solver.iterate(dt, time):
            bar.update(dt)


def close(ctx):
    if ctx.obj.get('frame_series', None):
        ctx.obj['frame_series'].hdf5.close()
        ctx.obj['frame_series'] = None

if __name__ == '__main__':
    cli(obj={'frame_series': None, 'frame': None})
