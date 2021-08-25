# Specify parent image. Please select a fixed tag here.
ARG BASE_IMAGE=registry.git.rwth-aachen.de/jupyter/profiles/rwth-courses:2020-ss.1
FROM ${BASE_IMAGE}

# Create new Environment with Python3.9
RUN conda create -n py39 python=3.9

# Activate py39 Environment
SHELL ["conda", "run", "-n", "py39", "/bin/bash", "-c"]

# Install iPyKernel
# RUN pip install ipykernel && ipython kernel install --user --name=py39 && conda deactivate
RUN conda activate py39 && pip install ipykernel && ipython kernel install --user --name=py39 && conda deactivate 

# Add modified Notebook Config. Sets py39 as default environment
ADD ./jupyter_notebook_config.py /etc/jupyter/jupyter_notebook_config.py

# Install packages via requirements.txt
ADD requirements.txt .
RUN pip install -r requirements.txt
