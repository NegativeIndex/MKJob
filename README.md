# MKJob
generate all kinds of job files for Argon cluster.

The job file is a Linux Bash program. It provides two functions:
resource request and job logging.



---

## Comsol Multiphysics job files

### v1.0 (10/22/2019)

I wrote two files to implement the function.

* **dwt-comsol-job-file.job**: This is the template file

* **mkjob-comsol.py**: This python code will modify the template job
    file to create real one.

---
## Meep job files

### v1.0 (10/22/2019)

I used five files to demonstrate how to generate pyhon job files.

* **rad.py**: This is the actual Python simulation file and stays in
    the root of thw working directory.

* **create_3D_jobs.py**: This command is placed in the same folder as
    **rad.py**. It will generate a lot of sub-folder; each has a
    standalong simulation. It modifies the parameters in
    **rad.py**. It also calles other files to generate appropriate job
    files.

* **dwt-mpi-job-file.job**: This is the template job file.

* **mkjob-mpi**: This command will copy the template file to the
   working folder and modify it,

* **modify-meep-jobfile.py**: This command will change the cluster
    resource.

