
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
ENV WORKSPACE_DIR $HOME_DIR/catkin_rdv
ENV PIPELINE_DIR $WORKSPACE_DIR/src/rdv_pipeline

# Set passwords:
# User: root, Password: $USERNAME
# User: $USERNAME, Password: $USERNAME
RUN echo "root:${USERNAME}" | chpasswd
RUN echo "${USERNAME}:${USERNAME}" | chpasswd

CMD ["bash"]
