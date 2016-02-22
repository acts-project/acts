#!/usr/bin/env python
#
# Helper script to setup a VM using the CERN opencloud service
# and initialise a GitLab CI Runner on it
#
# The shell profile needs to be configured according to
# http://clouddocs.web.cern.ch/clouddocs/tutorial/create_your_openstack_profile.html

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
    parser.add_argument('--user-data',metavar='USER_DATA',dest='user_data',help='user data file supplied to cloud-init')

    # Parse and handle initial arguments
    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    logger.debug("parsed input values: " + str(args))
    openstack_cmd = "openstack server create --key-name {0} --flavor {1} --image {2} {3}".format(args.key,args.flavor,args.image,args.hostname)
    if args.user_data:
        openstack_cmd += " --user-data {0}".format(args.user_data)
        
    try:
        logger.debug("Calling '{0}'".format(openstack_cmd))
        output = subprocess.check_output(openstack_cmd,shell=True)
    except subprocess.CalledProcessError:
        logger.warning("Executing '{0}' failed".format(openstack_cmd))
    
if __name__ == '__main__':
    main()
