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

This will also update ``~/.ssh/config`` to define an ssh alias for the
instance.  After the above, you should be able to ssh into the
instance::

  ssh ya16sdb

See the lines added to ``~/.ssh/config``::

  grep -A 6 ya16sdb ~/.ssh/config

Note that ssh access is via the public IP of the *instance* which will
change if the stack is re-created. In this case, update your ssh
config using::

  deploy/deploy.yml -i deploy/hosts --vault-password-file vault_pass.txt -t ssh-config

Register the domain name
------------------------

After the steps above, Add the DNS name of the ELB as a CNAME record to
``ya16sdb.labmed.uw.edu`` using the self-service DNS registration
portal: https://networks.uw.edu/networks/dns/resources

Determine the public DNS name for the ELB, either by looking in the
console (Services --> Cloudformation --> ya16sdb --> Resources -->
<load balancer url> --> DNS name ), or using the CLI, eg::

  % aws elb describe-load-balancers | jq -r '.LoadBalancerDescriptions[].DNSName' | grep ya16sdb
  ya16sdb-LoadBalanc-9OS62P5C7TIF-608540592.us-west-2.elb.amazonaws.com

It may take a while for the registration of ``ya16sdb.labmed.uw.edu`` to propagate.

Installing Dokku
----------------

Until I figure out how to restart the instance as part of the
deployment above, first restart the webserver::

  ssh ya16sdb "sudo restart"

Wait a few minutes, and make sure that you can ssh in. Now you can install Dokku::

  deploy/deploy-instance.yml -i deploy/hosts

After dokku is installed, install the ``ya16sdb`` app::

  ssh dokku@ya16sdb apps:create ya16sdb
  git remote add ya16sdb dokku@ya16sdb:ya16sdb
  git subtree push --prefix dash ya16sdb master
  ssh dokku@ya16sdb domains:add ya16sdb ya16sdb.labmed.uw.edu
  ssh dokku@ya16sdb storage:mount ya16sdb /var/lib/dokku/data/storage/ya16sdb:/storage
  ssh ya16sdb "sudo mkdir -p -m 0777 /var/lib/dokku/data/storage/ya16sdb"

Install a robots.txt (https://github.com/dwyl/dokku-robots.txt). ssh into the instance and run this::

  sudo dokku plugin:install https://notabug.org/candlewaster/dokku-robots.txt.git robots.txt

Then::

  ssh dokku@ya16sdb robots.txt:disallow ya16sdb

Now you should be able to visit the application at
http://ya16sdb.labmed.uw.edu

Updating the application
------------------------

After initial deployment, you can update the app with a git push::

  git subtree push --prefix dash ya16sdb master

Updating the data
-----------------

Upload a new feather file::

  scp filter_details.feather.gz ya16sdb:/var/lib/dokku/data/storage/ya16sdb

Restart the application::

  ssh dokku@ya16sdb ps:restart ya16sdb



Notes
=====

Editing the vault
-----------------
::

  EDITOR="$EMACSCLIENT -nw" ansible-vault edit --vault-password-file vault_pass.txt deploy/secrets.yml

Don't commit ``vault_pass.txt``!

Identifying an AMI
------------------

Identifying the 16.04 AMI appears to be somewhat tricky: the best
option may be to go to the console and find the default 16.04 AMI in
the first step after selecting the option to launch a new EC2
instance.

It's also possible to list AMIs in chronological order by creation date::

  aws ec2 describe-images --filters "Name=name,Values=ubuntu/images/hvm-ssd/ubuntu-xenial-16.04-amd64-server*" | \
  jq -c '.Images | sort_by(.CreationDate) | .[] | select(.State == "available") | {ImageId, Description}' | uniq | tail -n 10

Replacing the webserver instance
--------------------------------

The easiest way of replacing the EC2 instance without replacing the
entire stack is to change ``AMIId`` in ``CLOUDFORMATION_PARAMS`` and
re-deploy.
