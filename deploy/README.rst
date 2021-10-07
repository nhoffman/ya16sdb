===============================
 Application deployment to AWS
===============================

.. contents::

Requirements
============

* python 3.6+
* An existing dokku (or Heroku, maybe) instance
* Credentials for an AWS IAM user with necessary privileges to create
  an S3 bucket and IAM user
* An ssh key for a dokku server
* vault password in ``vault_pass.txt``

Running the application locally
===============================

Create the virtualenv::

  python3 -m venv py3-env
  source py3-env/bin/activate
  pip install -U pip wheel
  pip install -r requirements.txt

Run the application using flask's internal server::

  ./app.py

Or using gunicorn::

  gunicorn app:app.server

If you want to run the application locally but load the data from the S3 bucket::

  DATA_FILE=s3://ya16sdb-data/filter_details.feather.gz ./app.py

Installation to a Dokku serer
=============================

set up an environment
---------------------

::

   python3 -m venv deploy-env
   source deploy-env/bin/activate
   pip install -U pip wheel
   pip install -r requirements-deploy.txt

NOTE: mitogen will be installed to the virtualenv and so the path to the ansible mitogen plugin in the ansible.cfg is within the venv::

  strategy_plugins = deploy-env/lib/python3.6/site-packages/ansible_mitogen/plugins/strategy
  strategy = mitogen_linear

Deployment
----------

This deployment depends on the labmed dokku-stack, in order to provide application
routing as well as an IAM instance profile that has the Ya16sdb policy attached to it.
Please see the dokku-stack documentation for deployment information there.

In particular, you must add the configuration for the target server to
your ssh config. Clone, enter, and set up the dokku-stack repo and run
the following::

  AWS_PROFILE=saml ./deploy-stack.yml -e ACCOUNT=prod -e TARGET=dokku-stack-apps -t ssh-config

Return to this repo. Create the S3 bucket and IAM Policy, and install
and configure the application (only needs to be done once, or to
update dokku configuration)::

  deploy/deploy-dash.yml -e ENV=prod -t deploy-stack

  deploy/deploy-dash.yml -e ENV=prod -t deploy-app -i deploy/hosts

Add a git remote (once for each newly-cloned repo)::

  git remote add dokku-stack-apps-public-prod dokku@dokku-stack-apps-public-prod:ya16sdb

Push the application to the Dokku instance::

  git subtree push --prefix dash dokku-stack-apps-public-prod master

Updating the data file
----------------------

Copy the feather file to the S3 bucket::

  aws s3 cp filter_details.feather.gz s3://ya16sdb-data/filter_details.feather.gz

Note that this above assumes that your default AWS profile has access
to this bucket. If not, add the following to ``~/.aws/credentials``::

  [ya16sdb]
  aws_access_key_id = xxxx
  aws_secret_access_key = xxxx

And indicate the profile in the call to ``aws``::

  aws --profile ya16sdb s3 cp filter_details.feather.gz s3://ya16sdb-data/filter_details.feather.gz

Restart the server to recognize the new file::

  ssh dokku@dokku-stack-prod ps:restart ya16sdb

Notes
=====

Editing the vault
-----------------
::

  EDITOR="$EMACSCLIENT -nw" ansible-vault edit --vault-password-file vault_pass.txt deploy/secrets.yml

Don't commit ``vault_pass.txt``, and make sure the file mode is ``0600``!
