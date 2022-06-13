FROM isakhammer/workbox:latest
COPY . .
RUN cd latex && pdflatex main.tex

WORKDIR $WORK_DIR/ngsolve
RUN python3 HDG_biharmonic.py
RUN python3 DG_biharmonic.py
RUN python3 DG_poission.py

WORKDIR $WORK_DIR/julia
RUN julia DG_poission.jl
RUN julia DG_biharmonic.jl
RUN julia DG_brenner.jl
RUN julia biharmonic_julia_test.jl

