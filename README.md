# dyn2rel
`dyn2rel` is a set of MATLAB scripts to convert particle and tomogram data between two subtomogram averaging softwares: [Dynamo](https://www.dynamo-em.org/) and [RELION-4.0](https://relion.readthedocs.io/en/release-4.0/).

This version includes:
- `dyn2rel/get_particles_star.m` - to export a `Dynamo`-style particle table as `.tbl` file to the `RELION-4.0`-style `.star` particle file;
- `dyn2rel/get_tomograms_star.m` - to export a `Dynamo`-style tomogram list as `.vll` file to the `RELION-4.0`-style `.star` tomogram file;

Additionally, it includes:
- `dyn2rel/get_defocus_txt.m` - to prepare of an `IMOD`-style defocus `.txt` file from `GCTF` log files for import in `RELION-4.0`. 
