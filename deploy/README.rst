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

  export LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 FLASK_DEBUG=1
  flask run

Or using gunicorn::

  gunicorn app:app.server

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

Create the S3 bucket and IAM role, and install and configure the
application (only needs to be done once, or to update dokku
configuration)::

  deploy/deploy-dash.yml -i deploy/hosts --vault-password-file vault_pass.txt -e TARGET=dokku-stack-prod

After the IAM role is created, create credentials and add them to
``secrets.yml``::

  aws iam create-access-key --user-name Ya16sdbUser > Ya16sdbUser.json

Add a git remote (once for each newly-cloned repo)::

  git remote add dokku-stack-prod dokku@dokku-stack-prod:ya16sdb

Push the application to the Dokku instance::

  git subtree push --prefix dash dokku-stack-prod master

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
