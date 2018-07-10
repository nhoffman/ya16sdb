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
import socket
import sys

import boto3
from ansible.module_utils.basic import AnsibleModule


def main():

    module = AnsibleModule(
        argument_spec={
            'server_name': {'required': False, 'type': 'str'},
            'server_ip': {'required': False, 'type': 'str'}
        })

    server_name = module.params['server_name']
    server_ip = module.params['server_ip']

    if server_name:
        # get the ip address corresponding to server_name
        try:
            ip = socket.gethostbyname(server_name)
        except socket.gaierror:
            sys.exit('could not find IP address for "{}"'.format(server_name))
    elif server_ip:
        ip = server_ip
    else:
        sys.exit('either server_name or server_ip is required')

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


if __name__ == '__main__':
    main()
