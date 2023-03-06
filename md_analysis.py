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
    usage = '%(prog)s -traj trajectory [options]'
    parser = argparse.ArgumentParser(prog='md_analysis',
                                     description=desc,
                                     add_help=True,
                                     allow_abbrev=False,
                                     usage=usage)
    # -- Arguments: Trajectory ------------------------------------------------
    traj = parser.add_argument_group(title='Trajectory input options')
    traj.add_argument('-traj', dest='trajectory_file', action='store',
                      help='Path to the trajectory file', type=str,
                      required=True, default=None, metavar='trajectory')
    traj.add_argument('-top', dest='topology_file', action='store',
                      help='Path to the topology file', type=str,
                      required=False, default=None, metavar='topology')
    traj.add_argument('-f', dest='first', action='store',
                      help='First frame to analyze', type=int,
                      required=False, default=0, metavar='first_frame')
    traj.add_argument('-l', dest='last', action='store',
                      help='Last frame to analyze', type=int,
                      required=False, default=None, metavar='last_frame')
    traj.add_argument('-s', dest='stride', action='store',
                      help='Stride of frames to analyze', type=int,
                      required=False, default=1, metavar='stride')
    traj.add_argument('-sel', dest='selections', action='store',
                      help='Atom selection (MDTraj syntax)',
                      required=False, default='all', metavar='selections')
    # -- Arguments: Output ----------------------------------------------------
    out = parser.add_argument_group(title='Output options')
    out.add_argument('-odir', dest='outdir', action='store',
                     help='Output directory to store the analysis',
                     type=str, required=False, default='./', metavar='.')

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
        raise ValueError('\n\n>>> Topology with extension "{}" is not valid.'
                         .format(top_ext)
                         + '\nOptions are: {}'.format(valid_tops))

    return True


def load_traj(traj_file, top_file, valid_trajs, valid_tops):
    if need_for_top(traj_file) and not top_file:
        traj_ext = os.path.basename(traj_file).split('.')[-1]
        raise ValueError('\n\n>>> Trajectory files with extension {} need'
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
        raise ValueError('\n\n>>> The provided selection "{}"'
                         ' is not valid in MDTraj'.format(selection))

    if sel_idx.size == 0:
        raise ValueError('\n\n>>> The provided selection "{}"'
                         ' correspond to no atoms'.format(selection))

    return sel_idx


def range_traj(traj, first, last, stride):
    n_frames = traj.n_frames
    first_range = range(0, n_frames - 1)
    last_range = range(first + 1, n_frames)

    if last is not None:
        stride_range = range(1, last - first)
    else:
        stride_range = range(1, n_frames - first)

    if first not in first_range:
        raise ValueError('\n\n>>> First frame must be in the interval'
                         ' [{}, {}]'.format(first_range.start,
                                            first_range.stop))
    if last and (last not in last_range):
        raise ValueError('\n\n>>> Last frame must be in the interval'
                         ' [{}, {}]'.format(last_range.start, last_range.stop))
    if stride not in stride_range:
        raise ValueError('\n\n>>> Stride must be in the interval'
                         ' [{}, {}]'.format(stride_range.start,
                                            stride_range.stop))
    sliced = slice(first, last, stride)
    if sliced not in [slice(0, n_frames, 1), slice(0, None, 1)]:
        return traj.slice(sliced)

    return traj


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

    # -- Load arguments -------------------------------------------------------
    args = parse_arguments()

    # -- Load and process trajectory ------------------------------------------
    print('\n** Loading trajectory **')
    traj = load_traj(args.trajectory_file, args.topology_file, valid_trajs,
                     valid_tops)
    if traj.n_frames != 1:
        traj = range_traj(traj, args.first, args.last, args.stride)

    selections = args.selections.split(':')

    # -- Outdir process -------------------------------------------------------
    if args.outdir != './':
        try:
            os.makedirs(args.outdir)
        except FileExistsError:
            raise Exception('\n\n>>> Output dir already exist, please especify'
                            ' a different one')

    # -- RMSD -----------------------------------------------------------------
    general_canvas([12, 8], 300)
    fig, ax = plt.subplots()

    print('\n** Calculating RMSD **')

    rmsd_per_sel = np.zeros((len(selections), traj.n_frames))
    for idx, selection in enumerate(selections):
        traj_atoms = atoms_sel(traj, selection)
        rmsd = md.rmsd(traj, traj, frame=0, atom_indices=traj_atoms)
        rmsd *= 10
        rmsd_per_sel[idx] = rmsd

        ax.plot(np.arange(1, traj.n_frames+1, 1), rmsd, label=selection)

    ax.set_ylabel(r'RMSD $(\AA)$')
    ax.set_xlabel(r'Frame')
    ax.legend(loc='lower right')
    fig.savefig('{}.png'.format(os.path.join(args.outdir, 'RMSD')))
    plt.close(fig)

    rmsd_out = np.vstack((np.arange(1, traj.n_frames+1, 1), rmsd_per_sel))
    np.savetxt('rmsd.dat', rmsd_out.T, delimiter=',')

    # -- RMSF per residue -----------------------------------------------------
    print('\n** Calculating RMSF per residue **')
    residues = np.arange(0, traj.n_residues)
    rmsd_per_res = np.zeros((len(residues), traj.n_frames))

    backbone = 'name CA or name C or name N or name O or name H'  # for ncaa

    for res in residues:
        traj_atoms = atoms_sel(traj, 'resid {} and ({})'.format(str(res),
                                                                backbone))
        rmsd_per_res[res] = md.rmsd(traj, traj, frame=0,
                                    atom_indices=traj_atoms)
    rmsf = [np.mean(rmsd)*10 for rmsd in rmsd_per_res]

    fig, ax = plt.subplots()
    ax.plot(np.arange(1, traj.n_residues + 1, 1), rmsf, color='black')
    ax.set_xlabel(r'Residue')
    ax.set_ylabel(r'RMSF $(\AA)$')
    fig.savefig('{}.png'.format(os.path.join(args.outdir, 'RMSF')))
    plt.close(fig)

    rmsf_out = np.vstack((np.arange(1, traj.n_residues + 1, 1), rmsf))
    np.savetxt('rmsf.dat', rmsf_out.T, delimiter=',')
