#!/usr/bin/env python
#
# Helper script to setup a VM using the CERN opencloud service
# and initialise a GitLab CI Runner on it
#
# The shell profile needs to be configured according to
# http://clouddocs.web.cern.ch/clouddocs/tutorial/create_your_openstack_profile.html

import os
import tempfile
import shutil
import sys
import argparse
import logging
import subprocess

# Setup basic logging
logger = logging.getLogger('GLRunnerInit')
hdlr = logging.StreamHandler(sys.stdout)
frmt = logging.Formatter("%(name)s.%(funcName)s %(levelname)s %(message)s")
hdlr.setFormatter(frmt)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)
    
def main():
    parser = argparse.ArgumentParser(description='Initialise a GitLab Runner',
                                     epilog='Example',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('hostname', metavar='HOST',help='hostname of the created VM')
    parser.add_argument('-f','--flavor',metavar='FLAVOR',default='m1.small',choices=['m1.tiny','m1.small','m1.medium','m1.large'],help='flavor of VM (check \'openstack flavor list\')')
    parser.add_argument('-i','--image',metavar='IMAGE',default='8b082f2d-9728-4c3f-8440-cf110f03f745',help='boot image for VM (check \'openstack image list\')')
    parser.add_argument('-k','--key',metavar='KEY',required=True,help='name of public-private key pair registered with openstack')
    parser.add_argument('-n','--name',metavar='NAME',default='custom-runner',help='name of GitLab CI runner')
    parser.add_argument('-t','--tags',metavar='TAGS',nargs='*',default=[],help='list of tags attached to the GitLab CI Runner')
    parser.add_argument('-v','--debug',action="store_true",help="switch logging into DEBUG mode")
    parser.add_argument('--token',metavar='TOKEN',required=True,help='token of GitLab project the runner should be associated to')
    parser.add_argument('--url',metavar='URL',default='https://gitlab.cern.ch/ci',help='URL of GitLab instance where the GitLab project is hosted')

    # Parse and handle initial arguments
    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    logger.debug("parsed input values: " + str(args))

    # make temporary copy of user data file
    temp_dir = tempfile.gettempdir()
    temp_file = os.path.join(temp_dir,'tmp_install_software.sh')
    shutil.copy2('install_software.sh',temp_file)
    # append command to create GitLab runner with correct arguments
    tags = ','.join(args.tags)
    gitlab_cmd = "gitlab-runner register -n -u {0} -r {1} --name {2} --executor shell --tag-list {3}".format(args.url,args.token,args.name,tags)
    logger.debug("Using following command to register GitLab runner:\n{0}".format(gitlab_cmd))
    try:
        append_cmd = "echo '{0}' >> {1}".format(gitlab_cmd,temp_file)
        logger.debug("Calling '{0}'".format(append_cmd))
        output = subprocess.check_output(append_cmd,shell=True)
    except subprocess.CalledProcessError:
        logger.warning("Executing '{0}' failed".format(append_cmd))
        sys.exit(1)

    # create openstack VM
    openstack_cmd = "openstack server create --key-name {0} --flavor {1} --image {2} --user-data {3} {4}".format(args.key,args.flavor,args.image,temp_file,args.hostname)        
    try:
        logger.debug("Calling '{0}'".format(openstack_cmd))
        output = subprocess.check_output(openstack_cmd,shell=True)
    except subprocess.CalledProcessError:
        logger.warning("Executing '{0}' failed".format(openstack_cmd))
    
if __name__ == '__main__':
    main()
