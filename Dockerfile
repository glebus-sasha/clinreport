FROM rocker/r-ver:4.3.2

# --- Системные зависимости ---
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libglpk-dev \
    dos2unix \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# --- CRAN пакеты ---
RUN R -e "install.packages(c('dplyr','stringr','readr','ggplot2','plotly','htmlwidgets','optparse','RobustRankAggreg'), repos='https://cloud.r-project.org')"

# --- Bioconductor пакеты ---
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install(c('EnhancedVolcano','WebGestaltR'), update=FALSE, ask=FALSE)"

# --- Папка для твоих скриптов, отдельно от системных бинарников ---
RUN mkdir -p /usr/local/my-scripts
COPY scripts/ /usr/local/my-scripts/

# --- Делаем их исполняемыми и убираем возможные CRLF (Windows) ---
RUN chmod +x /usr/local/my-scripts/*.R && dos2unix /usr/local/my-scripts/*.R

# --- Алиасы для удобного вызова без .R ---
RUN for f in /usr/local/my-scripts/*.R; do ln -s "$f" "/usr/local/bin/$(basename ${f%.R})"; done

# --- Рабочая директория ---
WORKDIR /usr/local/my-scripts

# --- По умолчанию интерактивный shell ---
ENTRYPOINT ["/bin/bash"]