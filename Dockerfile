FROM isakhammer/workbox:latest
COPY . .
RUN latexmk latex/main.tex

