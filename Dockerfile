FROM centos:8

RUN set -eux; \
  # PostgreSQL
  dnf -y module disable postgresql; \
  dnf -y install https://download.postgresql.org/pub/repos/yum/reporpms/EL-8-x86_64/pgdg-redhat-repo-latest.noarch.rpm; \
  # EPEL
  dnf -y install epel-release; \
  # DNF plugins
  dnf -y install dnf-plugins-core; \
  # PowerTools
  dnf -y config-manager --set-enabled powertools; \
  # FMI smartmet
  dnf -y install https://download.fmi.fi/smartmet-open/rhel/8/x86_64/smartmet-open-release-21.3.26-2.el8.fmi.noarch.rpm; \
  # Fix smartmet bugs
  sed -i \
    -e 's/proxy/#proxy/g' \
    -e 's|8$basearch|8/$basearch|g' \
    /etc/yum.repos.d/smartmet-open.repo

RUN set -eux; \
  dnf -y install \
    himan-bin \
    himan-lib \
    himan-plugins

# Example code
WORKDIR /app

COPY ./example/seaicing/seaicing.json ./
COPY ./example/seaicing/param-file.txt ./
COPY ./example/seaicing/seaicing.grib ./

CMD himan -f seaicing.json --no-database --param-file param-file.txt seaicing.grib
