FROM isakhammer/workbox:latest
COPY . .
RUN cd latex && pdflatex main.tex

WORKDIR $WORK_DIR/ngsolve
RUN python3 HDG_biharmonic.py
RUN python3 DG_biharmonic.py
RUN python3 DG_poission.py

WORKDIR $WORK_DIR/julia
RUN julia run_tests.jl

