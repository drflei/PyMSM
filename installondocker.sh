#!/bin/bash -f 
cd /opt
git clone -b ts89c-ext --single-branch https://github.com/drflei/IRBEM.git
cd IRBEM/
make OS=linux64 ENV=gfortran64 all
make install
#
apt update
apt -y install python3-virtualenv
cd /opt
virtualenv -p python3 msmtest
source msmtest/bin/activate
cd PyMSM
pip install -r requirements.txt -e ./
cd ../IRBEM/python/
pip install -r requirements.txt
cd /opt
