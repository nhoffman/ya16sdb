#!/usr/bin/env python3

"""Given the parameter server_name, defines the variables
EIP_ALLOCATION_ID and EIP_PUBLIC_IP.

Example:

  tasks:
    - name: set variables
      set_eip_id:
        server_name: oltgtest.labmed.uw.edu

"""

import json
import subprocess
import sys
import shutil
import shlex
import socket

import boto3

def parse_args(infile):
    args = {}
    with open(infile) as f:
        for token in shlex.split(f.read()):
            if '=' in token:
                key, val = token.strip().split('=', 1)
                args[key] = val

    return args


args = parse_args(sys.argv[1])

try:
    server_name = args['server_name']
except KeyError:
    sys.exit('the parameter "server_name" is required')

# get the ip address corresponding to server_name
try:
    ip = socket.gethostbyname(server_name)
except socket.gaierror:
    sys.exit('could not find IP address for "{}"'.format(server_name))

# describe the elastic IP associated with this address
ec2 = boto3.client('ec2')
response = ec2.describe_addresses(PublicIps=[ip])
eip = response['Addresses'][0]

output = {
    'changed': True,
    'ansible_facts': {
        'EIP_ALLOCATION_ID': eip['AllocationId'],
        'EIP_PUBLIC_IP': ip,
    }
}

print(json.dumps(output))
