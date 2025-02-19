#!/bin/env python3
# -*- coding: utf-8 -*-


import mdtraj as md
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import os


valid_tops = set(['pdb', 'pdb.gz', 'h5', 'lh5', 'prmtop', 'parm7', 'prm7',
                  'psf', 'mol2', 'hoomdxml', 'gro', 'arc', 'hdf5', 'gsd'])
valid_trajs = set(['arc', 'dcd', 'binpos', 'xtc', 'trr', 'hdf5', 'h5', 'ncdf',
                   'netcdf', 'nc', 'pdb.gz', 'pdb', 'lh5', 'crd', 'mdcrd',
                   'inpcrd', 'restrt', 'rst7', 'ncrst', 'lammpstrj', 'dtr',
                   'stk', 'gro', 'xyz.gz', 'xyz', 'tng', 'xml', 'mol2',
                   'hoomdxml', 'gsd'])


def parse_arguments():
    desc = ('\n Simple MD trajectories analysis')
    usage = '%(prog)s -traj trajectory_file [options]'
    parser = argparse.ArgumentParser(prog='md_analysis',
                                     description=desc,
                                     add_help=True,
                                     allow_abbrev=False,
                                     usage=usage)
    # -- Arguments: Trajectory ------------------------------------------------
    traj = parser.add_argument_group(title='Trajectory input options')
    traj.add_argument('-traj', dest='trajectory_file', action='store',
                      help='Path to the trajectory file [required]', type=str,
                      required=True, default=None, metavar='trajectory')
    traj.add_argument('-top', dest='topology_file', action='store',
                      help='Path to the topology file', type=str,
                      required=False, default=None, metavar='topology')

    # -- Arguments: Analysis --------------------------------------------------
    analysis = parser.add_argument_group(title='Analiysis options')
    analysis.add_argument('-sel', dest='selections', action='store',
                          help='Atom selections (MDTraj syntax) separated by'
                          ' ":" (e.g. "all:backbone") [default: %(default)s]',
                          required=False, default='all', metavar='selections')
    analysis.add_argument('-labels', dest='labels', action='store',
                          help='Name of the selections for the legend of the'
                          ' rmsd plot separated by ":" (e.g. "all:backbone")',
                          required=False, type=str, default=None,
                          metavar='labels')
    analysis.add_argument('-i_time', dest='i_time', action='store',
                          help='time of the first frame (ns), note is the'
                          ' first in the trajectory and not the one in the -f'
                          ' option, if none is provided frame number will be'
                          ' used instead [default: %(default)s]',
                          required=False, default=None, metavar='initial_time',
                          type=float)
    analysis.add_argument('-tstep', dest='tstep', action='store',
                          help='Time step between frames in the original'
                          ' trajectory (ns), if none is provided'
                          ' frame number will be use instead'
                          ' [default: %(default)s]', metavar='time_step',
                          default=None, type=float, required=False)

    # -- Arguments: Output ----------------------------------------------------
    out = parser.add_argument_group(title='Output options')
    out.add_argument('-odir', dest='outdir', action='store',
                     help='Output directory to store the analysis'
                     ' [default: %(default)s]',
                     type=str, required=False, default='./', metavar='[path]')

    args = parser.parse_args()
    return args


def is_valid_traj(traj, valid_trajs):
    traj_ext = os.path.basename(traj).split('.')[-1]
    if traj_ext not in valid_trajs:
        raise ValueError('\n\n>>> Trajectory with extension "{}" is not valid.'
                         .format(traj_ext)
                         + '\nOptions are: {}'.format(valid_trajs))

    return True


def need_for_top(traj_file):
    traj_ext = os.path.basename(traj_file).split('.')[-1]
    if traj_ext in ['h5', 'lh5', 'pdb']:
        return False

    return True


def is_valid_top(top_file, valid_tops):
    top_ext = os.path.basename(top_file).split('.')[-1]
    if top_ext not in valid_tops:
        raise ValueError('\n\nERROR: Topology with extension "{}" is not valid.'
                         .format(top_ext)
                         + '\nOptions are: {}'.format(valid_tops))

    return True


def load_traj(traj_file, top_file, valid_trajs, valid_tops):
    if need_for_top(traj_file) and not top_file:
        traj_ext = os.path.basename(traj_file).split('.')[-1]
        raise ValueError('\n\nERROR: Trajectory files with extension {} need'
                         'a topology file'.format(traj_ext)
                         + '\nOptions are {}'.format(valid_tops))

    if is_valid_traj(traj_file, valid_trajs) and need_for_top(traj_file):
        if is_valid_top(top_file, valid_tops):
            return md.load(traj_file, top=top_file)

    if is_valid_traj(traj_file, valid_trajs) and not need_for_top(traj_file):
        return md.load(traj_file)


def atoms_sel(traj, selection):
    try:
        sel_idx = traj.topology.select(selection)
    except Exception:
        raise ValueError('\n\nERROR: The provided selection "{}"'
                         'is not valid in MDTraj'.format(selection))

    if sel_idx.size == 0:
        raise ValueError('\n\nERROR: The provided selection "{}"'
                         'correspond to no atoms'.format(selection))

    return sel_idx

def general_canvas(figsize, dpi):
    """
    Customization of plots

    Returns:
        None
    """
    mpl.rc('figure', figsize=figsize, dpi=dpi)
    mpl.rc('xtick', direction='in', top=False)
    mpl.rc('xtick.major', top=False)
    mpl.rc('xtick.minor', top=False)
    mpl.rc('ytick', direction='in', right=True)
    mpl.rc('ytick.major', right=False)
    mpl.rc('ytick.minor', right=False)
    mpl.rc('axes', labelsize=20)
    plt.rcParams['axes.autolimit_mode'] = 'data'
    mpl.rc('lines', linewidth=3, color='k')
    mpl.rc('font', family='monospace', size=20)
    mpl.rc('grid', alpha=0.5, color='gray', linewidth=1, linestyle='--')

    return


if __name__ == '__main__':
    args = parse_arguments()

    # =========================================================================
    # Checking out dir
    # =========================================================================
    if args.outdir != './':
        try:
            os.makedirs(args.outdir)
        except FileExistsError:
            raise Exception('\n\nERROR: Output dir already exist, please especify'
                            'a different one')

    # =========================================================================
    # Loading and processing the trajectory and selections
    # =========================================================================
    print('\n** Loading trajectory **')
    traj = load_traj(args.trajectory_file, args.topology_file, valid_trajs,
                     valid_tops)
    total_frames = traj.n_frames

    print(f'Trajectory with {total_frames} frames loaded')

    selections = args.selections.split(':')  # splitting different selections
    if not args.labels:
        print('\nINFO: No labels provided for the rmsd legend.\n'
              'Using generic one.')
        labels = ['sel_{}'.format(sel) for sel in
                  range(1, len(selections) + 1)]
    if args.labels:
        if len(labels) != len(selections):
            print('\nWARNING: Number of labels does not match with the number\n'
                  'of selections. Using generic labels')
            labels = ['sel_{}'.format(sel) for sel in
                      range(1, len(selections) + 1)]
        else:
            labels = args.labels.split(':')  # splitting different labels

    # =========================================================================
    #  RMSD
    # =========================================================================
    general_canvas([12, 8], 300)
    fig, ax = plt.subplots()

    print('\n** Calculating RMSD **')

    if (args.tstep == None) or (args.i_time == None):
        print('\nINFO: Using frame number in the x axis')
        x_axis = np.arange(1, total_frames+1, 1)
        x_axis_label = r'Frame'

    elif (args.tstep != None) and (args.i_time != None):
        print('\nINFO: Using time in the x axis')
        x_axis = np.arange(args.i_time,
                           args.i_time + total_frames * args.tstep,
                           args.tstep)
        x_axis_label = r'Time $(ns)$'

    rmsd_per_sel = np.zeros((len(selections), traj.n_frames))
    for idx, selection in enumerate(selections):
        traj_atoms = atoms_sel(traj, selection)
        rmsd = md.rmsd(traj, traj, frame=0, atom_indices=traj_atoms)
        rmsd *= 10
        rmsd_per_sel[idx] = rmsd

        ax.plot(x_axis, rmsd, label=labels[idx])

    ax.set_ylabel(r'RMSD $(\AA)$')
    ax.set_xlabel(x_axis_label)
    ax.legend(loc='lower right')
    fig.savefig('{}.png'.format(os.path.join(args.outdir, 'RMSD')))
    plt.close(fig)

    with open(os.path.join(args.outdir, 'rmsd.dat'), 'w') as rmsd_out:
        rmsd_out.writelines(['{},'.format(sel) for sel in labels])
        rmsd_out.write('Frame\n')
        for frame in range(traj.n_frames):
            rmsd_out.writelines(['{:.2f},'.format(sel[frame])
                                 for sel in rmsd_per_sel])
            rmsd_out.write('{}\n'.format(frame + 1))

    # =========================================================================
    #  RMSF per residue
    # =========================================================================
    print('\n** Calculating RMSF per residue **')

    traj_atoms = atoms_sel(traj, 'name CA')
    rmsf = md.rmsf(traj, traj, frame=0, atom_indices=traj_atoms)

    fig, ax = plt.subplots()
    ax.plot(np.arange(1, len(rmsf)+1, 1), rmsf, color='black')
    ax.set_xlabel(r'Residue')
    ax.set_ylabel(r'RMSF $(\AA)$')
    fig.savefig('{}.png'.format(os.path.join(args.outdir, 'RMSF')))
    plt.close(fig)

    with open(os.path.join(args.outdir, 'rmsf.dat'), 'w') as rmsf_out:
        rmsf_out.write('RMSF,ResId\n')
        for res, value in enumerate(rmsf):
            rmsf_out.write('{:.2f},{}\n'.format(rmsf[res], res + 1))

