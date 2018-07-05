# Use CentOS 7 as baseline
FROM centos:7

# Add Himan reposiry
RUN echo -e "[himan]\nname=Himan\nbaseurl=https://download.fmi.fi/himan/rhel/7/x86_64\nenabled=1\ngpgcheck=0\n" > /etc/yum.repos.d/himan.repo

# Add smartmet-open repository (newbase library)
RUN echo -e "[smartmet-open]\nname=Smartmet Open\nbaseurl=https://download.fmi.fi/smartmet-open/rhel/7/x86_64\nenabled=1\ngpgcheck=0\n" > /etc/yum.repos.d/smartmet-open.repo

# Add epel repository
RUN rpm -ivh https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm

# Add postgres 9.5 repo for libpqxx 5
RUN rpm -ivh https://download.postgresql.org/pub/repos/yum/9.5/redhat/rhel-7-x86_64/pgdg-centos95-9.5-3.noarch.rpm

# Install Himan and dependencies
RUN yum -y install \
	himan-bin \
	himan-lib \
	himan-plugins \
	wget

# Run example from Himan documentation
# Store example input files to container root as we mount /tmp from host

WORKDIR /

RUN wget \
	https://raw.githubusercontent.com/fmidev/himan/master/example/seaicing/seaicing.json \
	https://raw.githubusercontent.com/fmidev/himan/master/example/seaicing/param-file.txt \
	https://github.com/fmidev/himan/raw/master/example/seaicing/seaicing.grib

WORKDIR /tmp

CMD himan -f ../seaicing.json --no-database --param-file ../param-file.txt ../seaicing.grib
