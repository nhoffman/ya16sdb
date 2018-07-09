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

AWS credentials
---------------

Notes
=====

Editing the vault (don't commit ``vault_pass.txt``)::

  EDITOR="$EMACSCLIENT -nw" ansible-vault edit --vault-password-file vault_pass.txt deploy/secrets.yml


