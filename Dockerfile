FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates libcurl4-openssl-dev libssl-dev libxml2-dev libuv1-dev \
    libfontconfig1-dev libfreetype6-dev libpng-dev libjpeg-dev libtiff-dev \
    libharfbuzz-dev libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('shiny', 'DT', 'jsonlite', 'dplyr', 'purrr', 'stringr', 'htmltools', 'plotly', 'visNetwork'), repos='https://cloud.r-project.org', Ncpus=2)"

WORKDIR /opt/clinreport
COPY app.R run.R ./
COPY R/ ./R/
COPY www/ ./www/
COPY bin/clinreport /usr/local/bin/clinreport
RUN chmod 0755 /usr/local/bin/clinreport

EXPOSE 3838
ENTRYPOINT ["clinreport"]
CMD ["--help"]
