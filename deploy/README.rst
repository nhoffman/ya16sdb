===============================
 Application deployment to AWS
===============================

.. contents::

Requirements
============

* python 3.6+
* Credentials for an AWS IAM user with necessary privileges
* An SSH key (created in the AWS console and saved in ``secrets.yml``)
* An EIP (Elastic IP) address

Initial setup
=============

Virtualenv for deployment
-------------------------

::

   python3 -m venv deploy-env
   source deploy-env/bin/activate
   pip install -U pip wheel
   pip install -r requirements-deploy.txt

NOTE: mitogen will be installed to the virtualenv and so the path to the ansible mitogen plugin in the ansible.cfg is within the venv::

  strategy_plugins = deploy-env/lib/python3.6/site-packages/ansible_mitogen/plugins/strategy
  strategy = mitogen_linear

AWS credentials
---------------

See the following resources:

* https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html
* https://docs.aws.amazon.com/cli/latest/userguide/cli-environment.html

Generate and retrieve a key and secret access key in the AWS console
using instructions in the "getting started" page above. The detailed
instructions are in the "To get the access key ID and secret access
key for an IAM user" section of the first page.

Once you have the keys, set up your AWS profile like this::

  % aws configure
  AWS Access Key ID [None]: ...
  AWS Secret Access Key [None]: ...
  Default region name [None]: us-west-2
  Default output format [None]: json

Be sure to enter ``uw-west-2`` and ``json`` for the second two prompts
as above.

The steps above add credentials to a default profile stored in
``~/.aws``; if you are using multiple AWS profiles, perform the setup
using ``aws configure --profile profilename`` and then set an environment
variable ``AWS_PROFILE=profilename`` in the deployment environment.

Now you should be able to test your credentials::

  aws sts get-caller-identity


Deployment
==========

Deploying the CloudFormation stack
----------------------------------

Deploy the entire stack and application::

  deploy/deploy.yml -i deploy/hosts --vault-password-file vault_pass.txt

The first time you run the deployment, add the following parameters to
update your ssh configuration to install the private key and define
the ``ya16sdb`` ssh alias::

  deploy/deploy.yml -i deploy/hosts --vault-password-file vault_pass.txt -e ssh_config=yes -t ssh-config

See the lines added to ``~/.ssh/config``::

  grep -A 6 ya16sdb ~/.ssh/config

After the above, you should be able to ssh into the instance::

  ssh ya16sdb

Register the domain name
------------------------

After the steps above, determine the public DNS name for the ELB,
either by looking in the console (Services --> Cloudformation -->
ya16sdb --> Resources --> link to load balancer --> DNS name ), or
using the CLI, eg::

  % aws elb describe-load-balancers | jq -r '.LoadBalancerDescriptions[].DNSName' | grep ya16sdb
  ya16sdb-LoadBalanc-9OS62P5C7TIF-608540592.us-west-2.elb.amazonaws.com

Add the DNS name of the ELB as a CNAME record to
``ya16sdb.labmed.uw.edu`` using the self-service DNS registration
portal: https://networks.uw.edu/networks/dns/resources




Notes
=====

Editing the vault (don't commit ``vault_pass.txt``)::

  EDITOR="$EMACSCLIENT -nw" ansible-vault edit --vault-password-file vault_pass.txt deploy/secrets.yml


