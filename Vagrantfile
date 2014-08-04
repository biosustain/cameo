# -*- mode: ruby -*-
# vi: set ft=ruby :

$script = <<SCRIPT

sudo apt-get update
sudo apt-get install -y make git libglpk-dev swig glpk-utils swig python-dev python-numpy python-scipy redis-server
sudo apt-get install -y python-pip
sudo pip install python-libsbml-experimental cython numpy scipy -f https://dl.dropboxusercontent.com/u/22461024/pypi.drosophi.la/index.html --no-index
sudo apt-get install curl

mkdir tmp
cd tmp
curl -O ftp://ftp.gnu.org/gnu/glpk/glpk-4.45.tar.gz
tar xzfv glpk-4.45.tar.gz
cd glpk-4.45/
./configure
make
sudo make install

cd /vagrant
sudo pip install -r requirements.txt
# get rid of the stupid openpyxl error message
pip uninstall openpyxl
pip install openpyxl==1.8.6
pip install nose rednose coverage
sudo pip install redis
sudo python setup.py develop

SCRIPT

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

  config.vm.box = "ubuntu/trusty64"

  config.vm.provider :virtualbox do |vb|
    vb.customize ["modifyvm", :id, "--memory", "1024"]
  end

  # Disable automatic box update checking. If you disable this, then
  # boxes will only be checked for updates when the user runs
  # `vagrant box outdated`. This is not recommended.
  # config.vm.box_check_update = false

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine. In the example below,
  # accessing "localhost:8080" will access port 80 on the guest machine.
  config.vm.network "forwarded_port", guest: 6666, host: 6666

  config.vm.provision :shell, :inline => $script


end
