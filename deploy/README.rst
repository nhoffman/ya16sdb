===============================
 Application deployment to AWS
===============================

.. contents::

Requirements
============

* python 3.6+
* Credentials for an AWS IAM user with necessary privileges
* An existing dokku (or Heroku, maybe) instance
* An SSH key (created in the AWS console and saved in ``secrets.yml``)

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

Deployment
==========

Copy the feather file to the S3 bucket::

  aws s3 cp filter_details.feather.gz s3://ya16sdb-data/filter_details.feather.gz

Install and configure the application::

  deploy/deploy-dash.yml -i deploy/hosts --vault-password-file vault_pass.txt -e TARGET=dokku-stack-dev

Push the application to the Dokku instance::

  git subtree push --prefix dash dokku-stack-dev master

Notes
=====

Editing the vault
-----------------
::

  EDITOR="$EMACSCLIENT -nw" ansible-vault edit --vault-password-file vault_pass.txt deploy/secrets.yml

Don't commit ``vault_pass.txt``!
