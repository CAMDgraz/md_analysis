# md_analysis
> A simple script for MD analysis. Several selections can be use for the RMSD
> analysis (following MDTraj syntaxis). In the case of RMSF, md_analysis select
> automatically the backbone atoms. Please provide an input structure
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
  -f first_frame        First frame to analyze starting at 1 [default: 1]
  -l last_frame         Last frame to analyze [default: last frame]
  -s stride             Stride of frames to analyze [default: 1]

Analiysis options:
  -sel selections       Atom selections (MDTraj syntax) separated by ":" (e.g. "all:backbone") [default: all]
  -labels labels        Name of the selections for the legend of the rmsd plot separated by ":" (e.g. "all:backbone")
  -i_time initial_time  time of the first frame (ns), note is the first in the trajectory and not the one in the -f option, if none is provided frame number will be used instead
                        [default: None]
  -tstep time_step      Time step between frames in the original trajectory (ns), if none is provided frame number will be use instead [default: None]

Output options:
  -odir [path]          Output directory to store the analysis [default: ./]
```
