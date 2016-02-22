#!/bin/bash
# configuration script to run at first boot via Openstack
# using the user_data and cloud-init.

echo "userdata running on hostname: $(uname -n)"

# install and setup CVMFS
yum-config-manager --add-repo http://cvmrepo.web.cern.ch/cvmrepo/yum/cvmfs/x86_64/
yum install -y --nogpgcheck cvmfs cvmfs-config-default

cat >> /etc/cvmfs/default.local << EOF
CVMFS_HTTP_PROXY='http://ca-proxy-meyrin.cern.ch:3128;http://ca-proxy.cern.ch:3128;http://ca01.cern.ch:3128|http://ca02.cern.ch:3128|http://ca03.cern.ch:3128|http://ca04.cern.ch:3128|http://ca05.cern.ch:3128|http://ca06.cern.ch:3128'
CVMFS_CACHE_BASE='/var/lib/cvmfs'
CVMFS_FORCE_SIGNING='yes'
CVMFS_REPOSITORIES='atlas-condb.cern.ch,atlas.cern.ch'
EOF

/usr/sbin/setenforce 0
cvmfs_config setup
cvmfs_config probe

# install gcc and g++
yum install -y gcc-4.4.7-16.el6.x86_64
yum install -y gcc-c++-4.4.7-16.el6.x86_64

# setup GitLab CI runner
yum-config-manager --add-repo http://linuxsoft.cern.ch/mirror/packages.gitlab.com/runner/gitlab-ci-multi-runner/el/6/x86_64/
yum install -y --nogpgcheck gitlab-ci-multi-runner

gitlab-runner register -n -u https://gitlab.cern.ch/ci -r JcQLxS39oaAEsUU3ERwq --name ATS-Runner1 --executor shell --tag-list AFS,CVMFS,ATS,slc6

# set some environment variables
echo "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase" >> /home/gitlab-runner/.bashrc
echo 'alias setupATLAS="source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"' >> /home/gitlab-runner/.bashrc
