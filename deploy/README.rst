========================================
 Deployment of Bokeh application to AWS
========================================

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

Install Mitogen::

  wget https://github.com/dw/mitogen/archive/stable.zip
  unzip stable.zip


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
using ``aws configure --profile oltg`` and then set an environment
variable ``AWS_PROFILE=oltg`` in the deployment environment.

Now you should be able to test your credentials::

  aws sts get-caller-identity


Deployment
==========

Deploy the entire stack and application::

  deploy/deploy.yml -i deploy/hosts --vault-password-file vault_pass.txt



Notes
=====

Editing the vault (don't commit ``vault_pass.txt``)::

  EDITOR="$EMACSCLIENT -nw" ansible-vault edit --vault-password-file vault_pass.txt deploy/secrets.yml


