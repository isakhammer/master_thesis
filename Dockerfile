FROM isakhammer/workbox:latest
COPY . .
RUN cd latex && pdf main.tex

