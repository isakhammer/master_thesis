
FROM ubuntu:20.04
LABEL maintainer="Isak Hammer"

SHELL ["/bin/bash", "-c"]

ENV TZ=Europe/Oslo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Create a new user which is not root
ARG ID=1000
ENV USERNAME=kahuna
RUN groupadd -g $ID $USERNAME && \
    useradd -r -u $ID -m -g $USERNAME -G sudo -s /bin/bash $USERNAME

ENV HOME_DIR /home/$USERNAME

# Set passwords:
# User: root, Password: $USERNAME
# User: $USERNAME, Password: $USERNAME
RUN echo "root:${USERNAME}" | chpasswd
RUN echo "${USERNAME}:${USERNAME}" | chpasswd
ARG DEBIAN_FRONTEND=noninteractive

# Essential
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -yq apt-utils dialog
RUN apt-get install -yq build-essential software-properties-common
RUN apt-get update

# potential dependencies
RUN apt-get install -yq git cmake make
RUN apt-get install -yq bc libblas-dev liblapack-dev

# Deep bug somwhere when installing matplotlib on ubuntu20
# https://stackoverflow.com/questions/25674612/ubuntu-14-04-pip-cannot-upgrade-matplotllib
RUN apt-get install -yq libfreetype6-dev libxft-dev

# python3.8
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get -y update
RUN apt-get install -y python3.8 python3.8-tk python3.8-dev
RUN apt-get -y install python3-pip
RUN /usr/bin/python3.8 -m pip install --upgrade pip
RUN pip3 install --upgrade pip

RUN mkdir -p $HOME_DIR/project_thesis
ENV WORK_DIR $HOME_DIR/project_thesis
WORKDIR $WORK_DIR

COPY . .

RUN  . install-ngsuite.sh
CMD ["bash"]
