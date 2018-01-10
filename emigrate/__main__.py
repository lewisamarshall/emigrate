"""Command line interface for emigrate."""
from .__version__ import __version__
from .Solver import Solver
from .Frame import Frame
from .Sequence import Sequence
from .deserialize import deserialize

import sys
import os
from math import ceil
import click
import numpy as np
from matplotlib import pyplot
import matplotlib.animation as manimation
FFMpegWriter = manimation.writers['ffmpeg']

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
    ctx.obj['path'] = path
    _, file_extension = os.path.splitext(path)
    if file_extension == '.hdf5':
        ctx.obj['sequence'] = Sequence(path=path)
    elif file_extension == '.json':
        with open(path) as f:
            ctx.obj['frame'] = deserialize(f.read())
    else:
        raise RuntimeError("Can't load {} files.".format(file_extension))

    if io:
        sequence = ctx.obj['sequence']
        header = {'length': len(sequence),
                  'version': sequence.version()}
        click.echo(json.dumps(header))
        for frame_index in iter(sys.stdin.readline, ''):
            if 'exit' in frame_index:
                break
            try:
                frame = sequence[int(frame_index)]
                click.echo(frame.serialize(compact=True))
            except (IndexError, ValueError) as e:
                msg = {'error': repr(e)}
                click.echo(json.dumps(msg))


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
@click.option('--field', '-f', is_flag=True)
def movie(ctx, field):
    metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
    writer = FFMpegWriter(fps=15, metadata=metadata)

    sequence = ctx.obj['sequence']

    fig = pyplot.figure()
    lines = dict()
    frame = sequence[0]
    if field:
        line, = pyplot.plot(frame.nodes, frame.field, '-')
        pyplot.xlabel('x (mm)')
        pyplot.ylabel('electric field (V/m)')
        pyplot.xlim(xmax=frame.nodes[-1])
        pyplot.ylim([0, 10000])
        savename = os.path.splitext(ctx.obj['path'])[0]+'_field.mp4'
        with writer.saving(fig, savename, 100):
            for frame in sequence:
                line.set_data(frame.nodes, frame.field)
                writer.grab_frame()
        return

    for ion, ion_concentration in zip(frame.ions, frame.concentrations):
        lines[ion.name], = pyplot.plot([], [], '-', label=ion.name)

    pyplot.xlabel('x (mm)')
    pyplot.ylabel('concentration (M)')
    pyplot.ylim(ymin=0)
    pyplot.xlim(xmax=frame.nodes[-1])
    pyplot.legend()

    savename = os.path.splitext(ctx.obj['path'])[0]+'.mp4'
    with writer.saving(fig, savename, 100):
        for frame in sequence:
            for ion, ion_concentration in zip(frame.ions, frame.concentrations):
                lines[ion.name].set_data(frame.nodes, ion_concentration)
            writer.grab_frame()

@cli.command()
@click.pass_context
@click.option('--red', '-r', type=str, default=None)
@click.option('--green', '-g', type=str, default=None)
@click.option('--blue', '-b', type=str, default=None)
def spacetemp(ctx, red, green, blue):
    n = 1000
    sequence = ctx.obj['sequence']
    frame0 = sequence[0]
    nodes = np.linspace(frame0.nodes[0], frame0.nodes[-1], n)
    extent = [0, nodes[-1], 0, sequence[-1].time]

    slices = dict()
    for ion in frame0.ions:
        slices[ion.name] = np.zeros((len(sequence), n))

    for idx, frame in enumerate(sequence):
        for ion, concentration in zip(frame.ions, frame.concentrations):
            new_data = np.interp(nodes, frame.nodes, concentration)
            slices[ion.name][idx, :] += new_data

    # for name, data in slices.items():
    #     pyplot.imshow(data, origin='lower', extent=extent, aspect='auto')
    #     pyplot.xlabel('distance (m)')
    #     pyplot.ylabel('time (s)')
    #     pyplot.title(name)
    #     pyplot.savefig(ctx.obj['path']+'_{}_.png'.format(name))
    #     pyplot.clf()

    if all([red, green, blue]):
        print(slices.keys())
        red_slice = slices[red]
        green_slice = slices[green]
        blue_slice = slices[blue]
        color = np.zeros(red_slice.shape + (3,))
        color[:, :, 0] = red_slice
        color[:, :, 1] = green_slice
        color[:, :, 2] = blue_slice
        color /= color.max()/5
        pyplot.imshow(color, origin='lower', extent=extent, aspect='auto')
        pyplot.xlabel('distance (m)')
        pyplot.ylabel('time (s)')
        pyplot.savefig(ctx.obj['path']+'_rgb.png')
        pyplot.clf()

@cli.command()
@click.pass_context
# @click.option('--ion', '-i', type=str, default=None)
@click.option('--location', '-l', type=float, default=None)
def gram(ctx, location):
    sequence = ctx.obj['sequence']
    frame0 = sequence[0]
    times = [f.time for f in sequence]
    # nodes = np.linspace(frame0.nodes[0], frame0.nodes[-1], n)
    if location is None: location = frame0.nodes[-1]

    slices = dict()
    for ion in frame0.ions:
        slices[ion.name] = np.zeros((len(sequence), ))

    for idx, frame in enumerate(sequence):
        for ion, concentration in zip(frame.ions, frame.concentrations):
            new_data = np.interp(location, frame.nodes, concentration)
            slices[ion.name][idx] += new_data

    for name, data in slices.items():
        pyplot.plot(times, data)
        pyplot.xlabel('time (s)')
        pyplot.ylabel('concentration (M)')
        pyplot.title(name)
        pyplot.savefig(ctx.obj['path']+'_{}_electropherogram.png'.format(name))
        pyplot.clf()

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
            save = False
            if 'save' in constructor.keys():
                save = constructor.pop('save')
            ctx.obj['frame'] = Frame(constructor)
            click.echo(ctx.obj['frame'].serialize(compact=True))
            if save:
                with open(save, 'w') as loc:
                    loc.write(ctx.obj['frame'].serialize())
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
@click.option('-t', '--time', type=float, default=None)
@click.option('-d', '--dt', type=float, default=1)
@click.option('--io', is_flag=True)
def solve(ctx, output, dt, time, io):
    solver = Solver(ctx.obj['frame'])

    if io:
        for frame in solver.iterate(output, dt, time):
            click.echo(frame.serialize(compact=True))
    else:
        if time is not None:
            with click.progressbar(solver.iterate(output, dt, time),
                                   length=int(ceil(time/dt)),
                                   label='Solving...',
                                   ) as bar:
                for frame in bar:
                    pass
        else:
            with click.progressbar(length=10,
                                   label='Solving (Ctrl-C to stop)...',
                                   show_eta=False, show_percent=False,
                                   ) as bar:
                i = 0
                for frame in solver.iterate(output, dt):
                    i +=1
                    if i == 10:
                        i = 0
                        bar.update(-10)
                    bar.update(1)



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
    cli(obj={'sequence': None, 'frame': None, 'filename': None})

if __name__ == '__main__':
    main()
