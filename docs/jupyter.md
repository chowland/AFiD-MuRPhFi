# Post-processing: Jupyter notebooks on HPCs

Jupyter notebooks are a great way to view and analyse the outputs from the simulations.
One useful feature is that you can set a Jupyter notebook running on a cluster and interact with it from your local machine.
This is extremely helpful when running large simulations so you don't need to download the whole simulation output to analyse it.
Here's a step-by-step guide to setting up access to Jupyter on a cluster.

**Warning:** This guide will show you how to run Jupyter notebooks on the login node of a cluster.
This is usually fine, but should **not** be used to perform computationally intensive analysis.
You will annoy other users and the HPC staff if you hog the resources of a login node.
If you need to run heavy post-processing, best to submit it as a full job on the compute nodes.

## Step 1: Check for a Jupyter module!
Since Jupyter as a system is now so widely used, the HPC may offer Jupyter as a loadable module.
This will be the easiest to use, and will likely provide a lot of helpful Python packages up front.
A quick search of the HPC documentation for Jupyter should show this up.

### Step 1a: Install Jupyter yourself in a virtual environment
If Jupyter is not available, then it is best to create a new virtual environment to keep it in.
Preferably `conda` (via `miniconda`) will be available on the HPC, and you can create a new environment with the command
```
conda create -n jupyter_env -c conda-forge notebook
```
This should intall a new conda environment called `jupyter_env` on the system with the `notebook` package and all its prerequisites.
Activate the environment by calling `conda activate jupyter_env`.

If you want to install more packages, this can be achieved using `conda install ...`

### Step 1b: What if my cluster can't access the internet?
Quite often, for security reasons, an HPC cluster does not make outgoing connections to the internet.
Unfortunately this means it cannot directly download software from conda or github or whatever source you want.

Luckily, we can get around this by using a HTTP proxy!
When you ssh onto the HPC, you need to use a command like
```
ssh -R 34567 username@my_hpc.nl
```
where we have chosen a specific 5-digit port as `34567`.
Then, add the following lines to your `.bashrc` or `.bash_profile` on the HPC system:
```
export http_proxy=socks5://localhost:34567
export https_proxy=socks5://localhost:34567
```
where the number must match the port from the ssh command.

Instead of writing `-R 34567` each time you want to access the HPC with internet access, you can create an entry in `.ssh/config` with the line
```
RemoteForward 34567
```

## Step 2: Run a notebook on the cluster and access it

Now that a Jupyter installation has been found or created, we should be ready to run a notebook on the HPC and link to it from your local machine.

Firstly, you need to set up port forwarding as you access the HPC.
Choose a 4 digit number to use as the port you want to forward (here we choose 8192), then on your *local machine*, run
```
ssh -L 8192:localhost:8192 username@myhpc.nl
```

In the same terminal, which should now be connected to the cluster, set up the Jupyter environment (by loading whatever modules you need or activating your virtual environment).
Once the system is ready, start the notebook with the following options:
```
jupyter notebook --no-browser --port=8192 --ip=127.0.0.1
```

You should now be able to use the link that pops up in the terminal to access the notebook from your web browser.