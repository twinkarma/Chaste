# Parameter Sweeping Tutorial

In this tutorial we will go through how to setup chaste for parameter sweeping on the HPC cluster. 


## Development Environment, pulling a a Singularity image

ShARC supports the [Singularity](https://sylabs.io/docs/) containerisation technology. If you're not developing on Linux, it is also possible to use docker instead. You will want to install either [Singularity](https://sylabs.io/docs/) or [docker](https://docker.com) on your local machine. The use of containers allow us to standardise the development environment and build dependencies on your local machine and your cluster. 

To get the latest Singularity image of Chaste, run the following:

```
singularity pull docker://chaste/chaste-docker:latest
```

You should see a file `chaste-docker_latest.sif`. 

Run the following to get into the image's bash shell:

```
singularity exec chaste-docker_latest.sif /bin/bash
```

You can run `exit` to exit from the image.

You will be using the image you've just pulled to do the building and running of the code.

Note: As you will only be using the container for replying on build dependency, you will not need to modify the image during development which needs admin permission.  


## Initial code

You now need your own version of Chaste code. The Chaste repository can be found at [https://github.com/Chaste/Chaste](https://github.com/Chaste/Chaste) and this can be forked or cloned directly depending on your requirements. We will just clone directly in this example:

```
git clone https://github.com/Chaste/Chaste.git ~/Chaste
```

The above command clones Chaste code into the `Chaste` folder in our home directory. You will now need create a custom Chaste app or project with a main() function to be able to accept parameters.

The the instructions on how to create  Chaste executable apps and user project can be followed below:
  * Make executable apps: https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/BuildingExecutableApps
  * User projects: https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/UserProjects
	

## Modify your main() function to accept parameters

The parameter sweeper follows the boost::program_options format, where each parameter is passed in the format of `--param_name param_value` e.g. if your simulation has parameters a and b, the following format: 

```
my_chaste_app --output_dir "/outtput/path/[simulation id]" --a [param a's value] --b [param b's value] 
```

The `output_dir` parameter is always passed to your app and all outputs of the app should be enforced to only write to this directory.


### Automatic main() generator

A `main.cpp` can be automatically generated for your with by calling XXXX 

```
python run.py param1,param2,param3,param4 main.cpp
```

## Declaring parameters and generating the batch script

A python API has been created, the class `ParamSweeper` has been created to facilitate this parameter sweeping process. 


We create a python script `sweep.py` for defining our parameter and creating a batch job:

```
import numpy as np
from chastepython.util.parametersweep import ParamSweeper
p = {}
p['a'] = np.linspace(0, 10, 5)
p['b'] = np.linspace(0.1, -0.5, 5)
p['c'] = [10, 20, 30]

# Feed to param sweeper
sweeper = ParamSweeper()

exec_cmd = "chastepython/test/test_params.sh"
output_dir = "/tmp/myoutdir_sge"
scheduler = ParamSweeper.SGE

# Test output for SGE
sweeper.generate_batch_output(output_dir=output_dir,
                              exec_cmd=exec_cmd,
                              parameters=p,
                              batch_params=["-M myemail@mydomain.com", "-m bes"])
```

Notice that we can additional parameters to the batch script by adding to `batch_params`. In our example we've set the scheduler to alert us of the job status `-m bes` and provided our e-mail `-M myemail@mydomain.com`.

We then run the `sweep.py` script:

```
python sweep.py
```

In the output directory `XXX` that you specified, you will see the `params.json`, `runsimulation.py` and `batch.sge.sh` file. 

  * `params.json` Contains an expanded list of parameters that will be explored.
  * `batch.sge.sh` Batch script containing a task array for running through all of the parameters to be explored.
  * `runsimulation.py` Simulation runner script for running individual instances of the simulation.   


You can submit this file to the SGE scheduler to start your parameter sweeping task:

```
qsub outputdir/batch.sge.sh
```



## Running the sweep locally 

Sweep can also be run locally on your machine with the `perform_serial_sweep` command:

```
import numpy as np
from chastepython.util.parametersweep import ParamSweeper
p = {}
p['a'] = np.linspace(0, 10, 5)
p['b'] = np.linspace(0.1, -0.5, 5)
p['c'] = [10, 20, 30]

# Feed to param sweeper
sweeper = ParamSweeper()

exec_cmd = "chastepython/test/test_params.sh"
output_dir = "/tmp/myoutdir_sge"

sweeper.perform_serial_sweep(output_dir=output_dir, exec_cmd=exec_cmd, parameters=p)
```



```
# Pulls a singularity image
singularity pull docker://chaste/chaste-docker:latest


# Building a standard target
singularity exec https://github.com/Chaste/Chaste.git chaste
cd chaste/build
cmake ..

# Specify the location of files written out
export CHASTE_TEST_OUTPUT="~/"

```


