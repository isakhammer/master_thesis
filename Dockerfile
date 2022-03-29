FROM isakhammer/workbox:latest
COPY . .
RUN cd latex && latexmk main.tex

