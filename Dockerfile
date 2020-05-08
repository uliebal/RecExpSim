# Specify parent image. Please select a fixed tag here.
ARG BASE_IMAGE=registry.git.rwth-aachen.de/jupyter/profiles/rwth-courses:2020-ss.1
FROM ${BASE_IMAGE}

# Install packages via requirements.txt
ADD requirements.txt .
RUN pip install -r requirements.txt
