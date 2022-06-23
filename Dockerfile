FROM isakhammer/workbox:latest
COPY . .
RUN cd latex && pdflatex main.tex

WORKDIR $WORK_DIR/ngsolve
RUN python3 HDG_biharmonic.py
RUN python3 DG_biharmonic.py
RUN python3 DG_poission.py

WORKDIR $WORK_DIR/julia
# RUN julia run_tests.jl
RUN julia DG_poission.jl
RUN julia DG_brenner.jl
# RUN julia main.jl
RUN DISPLAY=:0 && xvfb-run -s '-screen 0 1024x768x24' julia --color=yes 'main.jl'
RUN julia cahn_hilliard.jl


