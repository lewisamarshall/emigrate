"""Command line interface for emigrate."""
from .Solver import Solver
from .Frame import Frame
from .Sequence import Sequence
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
@click.argument('path', type=click.Path(exists=True))
@click.option('--io', is_flag=True)
def load(ctx, path, io):
    """Open an emgrate file and return a serialized frame."""
    _, file_extension = os.path.splitext(path)
    if file_extension == '.hdf5':
        ctx.obj['sequence'] = Sequence(path=path)
    elif file_extension == '.json':
        with open(path) as f:
            ctx.obj['frame'] = deserialize(f.read())
    else:
        raise RuntimeError("Can't load {} files.".format(file_extension))

    if io:
        for frame_num in iter(sys.stdin.readline, ''):
            frame = ctx.obj['sequence'][int(frame_num)]
            click.echo(frame.serialize(compact=True))


@cli.command()
@click.pass_context
@click.argument('output', type=click.Path(exists=False))
@click.option('-f', '--frame', type=click.INT, default=None)
def plot(ctx, output, frame):

    ensure_frame(ctx, frame)
    frame = ctx.obj['frame']

    for ion, ion_concentration in zip(frame.ions, frame.concentrations):
        pyplot.plot(frame.nodes, ion_concentration, '-', label=ion.name)

    pyplot.xlabel('x (mm)')
    pyplot.ylabel('concentration (M)')
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmax=frame.nodes[-1])
    pyplot.legend()
    pyplot.savefig(output, bbox_inches='tight')
    close(ctx)


@cli.command()
@click.pass_context
@click.option('-f', '--frame', type=click.INT, default=None)
def echo(ctx, frame):
    ensure_frame(ctx, frame)
    click.echo(ctx.obj['frame'].serialize(compact=True))
    close(ctx)


@cli.command()
@click.pass_context
@click.option('-i', '--infile', type=click.Path(exists=True), default=None)
@click.option('-o', '--output', type=click.Path(exists=False), default=None)
@click.option('--io', is_flag=True)
def construct(ctx, infile, output, io):
    if io:
        for constructor in iter(sys.stdin.readline, ''):
            constructor = deserialize(constructor)
            click.echo(Frame(constructor).serialize(compact=True))
    else:
        infile = click.format_filename(infile)
        with open(infile, 'r') as inputfile:
            constructor = deserialize(inputfile.read())
        ctx.obj['frame'] = Frame(constructor)
        if output:
            with open(output, 'w') as loc:
                loc.write(ctx.obj['frame'].serialize())
    close(ctx)


@cli.command()
@click.pass_context
@click.option('-o', '--output', type=click.Path(exists=False))
@click.option('-t', '--time', type=float)
@click.option('-d', '--dt', type=float)
def solve(ctx, output, dt, time):
    solver = Solver(ctx.obj['frame'])

    with click.progressbar(solver.iterate(output, dt, time),
                           length=int(ceil(time/dt)),
                           label='Solving...',
                           ) as bar:
        for frame in bar:
            pass


def close(ctx):
    if ctx.obj['sequence']:
        ctx.obj['sequence'].close()
        ctx.obj['sequence'] = None


def ensure_frame(ctx, frame):
    if frame:
        ctx.obj['frame'] = ctx.obj['sequence'][frame]
    if not ctx.obj['frame']:
        n = click.prompt('Frame', default=1, type=click.INT)
        ctx.obj['frame'] = ctx.obj['sequence'][n]


def main():
    cli(obj={'sequence': None, 'frame': None})

if __name__ == '__main__':
    main()
