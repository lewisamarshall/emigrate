"""Command line interface for emigrate."""
from .Solver import Solver
from .Frame import Frame
from .FrameSeries import FrameSeries

import sys
import click
from matplotlib import pyplot


@click.group(name='emigrate', chain=True)
@click.pass_context
def cli(ctx):
    pass


@cli.command()
@click.pass_context
@click.argument('filename', type=click.Path(exists=True))
@click.option('--io', is_flag=True)
@click.option('--frame', '-f', prompt=False, default=None, type=click.INT)
def open(ctx, filename, frame, io):
    """Open an emgrate file and return a serialized frame."""
    ctx.obj['frame_series'] = FrameSeries(filename=filename, mode='r')
    if frame:
        ctx.obj['frame'] = ctx.obj['frame_series'][frame]
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
    for ion_concentration in ctx.obj['frame'].concentrations:
        pyplot.plot(ctx.obj['frame'].nodes, ion_concentration, '-')
    pyplot.xlabel('x (mm)')
    pyplot.ylabel('concentration (M)')
    pyplot.savefig('output')


@cli.command()
@click.pass_context
def echo(ctx):
    if not ctx.obj['frame']:
        n = click.prompt('Frame', default=1, type=click.INT)
        ctx.obj['frame'] = ctx.obj['frame_series'][n]

    click.echo(ctx.obj['frame'].serialize())


def close(ctx):
    if ctx.obj.get('frame_series', None):
        ctx.obj['frame_series'].hdf5.close()
        ctx.obj['frame_series'] = None

if __name__ == '__main__':
    cli(obj={'frame_series': None, 'frame': None})
