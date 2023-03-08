# md_analysis
> A simple script for MD analysis. Several selections can be use for the RMSD
> analysis (following MDTraj syntaxis). In the case of RMSF, md_analysis select
> automatically the backbone of the protein. Please provide an input structure
> without ions or solvent molecules.

## Download and use

Before using md_analysis you need to have MDTraj for RMSD and RMSF
calculations. It is recomended to install MDTraj using conda:
`conda install -c conda-forge mdtraj`

After installing MDTraj you can clone the repository:
```
git clone https://github.com/CAMDgraz/md_analysis
cd md_analysis
```
Then, you should be able to see the help by typing:
`python md_analysis.py -h`

```
$ python md_analysis.py -h

usage: md_analysis -traj trajectory_file [options]

Simple MD trajectories analysis

optional arguments:
  -h, --help            show this help message and exit

Trajectory input options:
  -traj trajectory      Path to the trajectory file [required]
  -top topology         Path to the topology file
  -f first_frame        First frame to analyze starting at 0 [default: 1]
  -l last_frame         Last frame to analyze (counting from 0) [default: last frame]
  -s stride             Stride of frames to analyze [default: 1]
  -sel selections       Atom selections (MDTraj syntax) separated by ":" (e.g.
                        "all:backbone") [default: all]

Analiysis options:
  -labels LABELS [LABELS ...]
                        Name of the selections for the legend of the rmsd plot

Output options:
  -odir [path]          Output directory to store the analysis

```
