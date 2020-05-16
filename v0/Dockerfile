FROM rocker/shiny
WORKDIR /home/shiny
RUN mkdir /home/shiny/www
COPY www/* /home/shiny/www/
COPY installPrerequisites.R .
RUN R -f ./installPrerequisites.R .
COPY packages/nglShiny_0.99.13.tar.gz .
RUN R CMD INSTALL nglShiny_0.99.13.tar.gz .
COPY 4our.R .
COPY *.html ./
CMD ["R", "-e", "shiny::runApp('/home/shiny/4our.R', port=3838, host='0.0.0.0')"]

