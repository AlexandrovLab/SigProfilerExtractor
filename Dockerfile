# Start with a base Ubuntu image and install Python
FROM ubuntu:22.04

# Avoid prompts from apt
ARG DEBIAN_FRONTEND=noninteractive

# Install Python, git, and other dependencies, and apply updates
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y python3-pip python3-dev git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Argument to specify the commit, defaults to the 'master' branch
ARG COMMIT_SHA=master

# Set the working directory in the container
WORKDIR /usr/src/app

# Install the package directly from the specific commit on GitHub
RUN pip3 install 'git+https://github.com/AlexandrovLab/SigProfilerExtractor.git@'${COMMIT_SHA}

# Create a non-root user
RUN useradd -m -s /bin/bash spm_user

# Change the ownership of the working directory
RUN chown -R spm_user:spm_user /usr/src/app

# Switch to the non-root user for security
USER spm_user