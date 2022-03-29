FROM isakhammer/workbox:latest
COPY . .
RUN cd latex && pdflatex main.tex

