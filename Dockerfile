FROM rocker/r-ver:4.4.1

ENV LANG=C.UTF-8 \
    LANGUAGE=C.UTF-8 \
    LC_ALL=C.UTF-8

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates libcurl4-openssl-dev libssl-dev libxml2-dev libuv1-dev \
    fontconfig fonts-dejavu-core locales \
    libfontconfig1-dev libfreetype6-dev libpng-dev libjpeg-dev libtiff-dev \
    libharfbuzz-dev libfribidi-dev \
    && sed -i 's/^# *en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && locale-gen en_US.UTF-8 \
    && fc-cache -f \
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
