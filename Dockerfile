FROM isakhammer/workbox:latest
COPY . .
RUN cd latex && pdflatex main.tex

RUN cd ngsolve && python3 HDG_biharmonic.py
RUN cd ngsolve && python3 DG_biharmonic.py
RUN cd ngsolve && python3 DG_poission.py

RUN cd julia && julia DG_poission.jl
RUN cd julia && julia DG_biharmonic.jl

