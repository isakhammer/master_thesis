FROM dealii/dealii:master-focal-root

RUN mkdir -p /root/code
ENV WORK_DIR /root/code
WORKDIR $WORK_DIR
COPY . .
RUN echo "source ${WORK_DIR}/common_scripts.sh" >> /root/.bashrc

SHELL ["/bin/bash", "-c"]
