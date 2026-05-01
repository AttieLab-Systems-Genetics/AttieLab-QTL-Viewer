# Use a base image that allows specifying R version
FROM rocker/shiny-verse:latest

USER root

# Install necessary *additional* system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    gfortran \
    libblas-dev \
    liblapack-dev \
    cmake \
    libsodium-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages using the script
COPY install_packages.R /tmp/install_packages.R
RUN Rscript /tmp/install_packages.R && rm /tmp/install_packages.R

# v3: belt-and-suspenders install for data-layer packages — this catches anything
# install_packages.R might have missed and fails the build loudly if CRAN is down.
RUN Rscript -e "install.packages(c('fst', 'promises', 'future', 'future.apply', 'paws.storage', 'filelock'), repos='https://cloud.r-project.org'); for (p in c('fst','promises','future','future.apply','paws.storage','filelock')) if (!requireNamespace(p, quietly=TRUE)) stop('missing: ', p)"

# Force app logs to stderr
ENV SHINY_LOG_STDERR=1

# v3: cache mount point for S3-fetched files. Override QTLAPP_CACHE_DIR at
# runtime to point at a host-mounted volume if you want cache persistence
# across container restarts.
RUN mkdir -p /cache && chmod -R 777 /cache
ENV QTLAPP_CACHE_DIR=/cache

# Copy the application files
COPY app.R /srv/shiny-server/
COPY R /srv/shiny-server/R/
COPY utils /srv/shiny-server/utils/

# v3: entrypoint that writes container env vars into /srv/shiny-server/.Renviron
# before launching shiny-server. Required because shiny-server's `su --login`
# strips the environment from the R process.
COPY deployment/latest_app/entrypoint.sh /usr/local/bin/qtlapp-entrypoint.sh
RUN chmod +x /usr/local/bin/qtlapp-entrypoint.sh

# Set ownership of the app directory
RUN chown -R shiny:shiny /srv/shiny-server && \
    chmod -R 755 /srv/shiny-server

# Expose the port Shiny Server runs on
EXPOSE 3838

# Stay as root to access mounted directories, but shiny-server will drop to shiny user
CMD ["/usr/local/bin/qtlapp-entrypoint.sh"]
