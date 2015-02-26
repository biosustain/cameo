# -*- mode: ruby -*-
# vi: set ft=ruby :

$script = <<SCRIPT

sudo apt-get update
sudo apt-get install -y make git libglpk-dev swig glpk-utils swig python-dev python-numpy python-scipy python-matplotlib
sudo apt-get install -y python-pip
sudo pip install python-libsbml-experimental cython -f https://dl.dropboxusercontent.com/u/22461024/pypi.drosophi.la/index.html --no-index
sudo pip install ipython[notebook]
sudo pip install nose
cd /vagrant
sudo pip install -r requirements.txt
sudo pip install escher --pre
sudo python setup.py develop
# get rid of the stupid openpyxl error message
pip uninstall openpyxl
pip install openpyxl==1.8.6
pip install nose rednose coverage
sudo pip install redis
sudo python setup.py develop
ipython notebook --ip=0.0.0.0 --port=8888

SCRIPT

# Vagrantfile API/syntax version. Don't touch unless you know what you're doing!
VAGRANTFILE_API_VERSION = "2"

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

  config.vm.box = "ubuntu/trusty64"

  config.vm.provider "virtualbox" do |v|
  host = RbConfig::CONFIG['host_os']

  # Give VM 1/4 system memory & access to all cpu cores on the host
  if host =~ /darwin/
    cpus = `sysctl -n hw.ncpu`.to_i
    # sysctl returns Bytes and we need to convert to MB
    mem = `sysctl -n hw.memsize`.to_i / 1024 / 1024 / 4
  elsif host =~ /linux/
    cpus = `nproc`.to_i
    # meminfo shows KB and we need to convert to MB
    mem = `grep 'MemTotal' /proc/meminfo | sed -e 's/MemTotal://' -e 's/ kB//'`.to_i / 1024 / 4
  else # sorry Windows folks, I can't help you
    cpus = 2
    mem = 1024
  end

  v.customize ["modifyvm", :id, "--memory", mem]
  v.customize ["modifyvm", :id, "--cpus", cpus]
end

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
