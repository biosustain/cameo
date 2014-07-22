# -*- mode: ruby -*-
# vi: set ft=ruby :

$script = <<SCRIPT

sudo apt-get update
sudo apt-get install -y make git libglpk-dev swig glpk-utils swig python-dev python-numpy python-scipy
sudo apt-get install -y python-pip
sudo pip install ply
cd /home/vagrant
mkdir tmp
cd tmp
wget http://www.dcc.fc.up.pt/~jpp/code/python-glpk/python-glpk_0.4.45.orig.tar.gz
tar xzf python-glpk_0.4.45.orig.tar.gz
cd python-glpk-0.4.45/src/
sudo make install
cd /home/vagrant
sudo pip install python-libsbml-experimental cython numpy scipy -f https://dl.dropboxusercontent.com/u/22461024/pypi.drosophi.la/index.html --no-index
sudo pip install nose
cd /vagrant
sudo pip install -r requirements.txt
# sudo python setup.py develop

SCRIPT

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

  config.vm.box = "ubuntu/trusty64"

  # Disable automatic box update checking. If you disable this, then
  # boxes will only be checked for updates when the user runs
  # `vagrant box outdated`. This is not recommended.
  # config.vm.box_check_update = false

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine. In the example below,
  # accessing "localhost:8080" will access port 80 on the guest machine.
  config.vm.network "forwarded_port", guest: 8888, host: 8888

  config.vm.provision :shell, :inline => $script

end
