#!/usr/bin/env python3

"""Given the parameter stack_name, defines the variable INSTANCE_IP.

Example:

  tasks:
    - name: set instance IP
      set_instance_id:
        stack_name: ya16sdb

"""

import json
import socket
import sys

import boto3
from ansible.module_utils.basic import AnsibleModule


def main():

    module = AnsibleModule(
        argument_spec={
            'stack_name': {'required': False, 'type': 'str'},
        })

    stack_name = module.params['stack_name']

    # describe the elastic IP associated with this address
    ec2 = boto3.client('ec2')
    response = ec2.describe_instances(
        Filters=[{'Name': 'tag:aws:cloudformation:stack-name',
                  'Values': [stack_name]}])

    ip = response['Reservations'][0]['Instances'][0]['PublicIpAddress']

    output = {
        'changed': True,
        'ansible_facts': {
            'INSTANCE_IP': ip,
        }
    }

    module.exit_json(**output)


if __name__ == '__main__':
    main()
