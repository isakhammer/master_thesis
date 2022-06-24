FROM isakhammer/workbox:latest
COPY . .
RUN cd latex && pdflatex main.tex



RUN apt-get update
RUN apt-get install -y xorg-dev mesa-utils xvfb libgl1 freeglut3-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libxext-dev
WORKDIR $WORK_DIR/julia
RUN DISPLAY=:0 xvfb-run -s '-screen 0 1024x768x24' julia --color=yes 'main.jl'
RUN julia cahn_hilliard.jl
# RUN julia run_tests.jl
RUN julia DG_poission.jl
RUN julia DG_brenner.jl
# RUN julia main.jl


WORKDIR $WORK_DIR/ngsolve
RUN python3 HDG_biharmonic.py
RUN python3 DG_biharmonic.py
RUN python3 DG_poission.py

