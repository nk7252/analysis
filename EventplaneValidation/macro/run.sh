#!/usr/bin/bash
export USER="$(id -u -n)"
export LOGNAOME=${USER}
export HOME=/sphenix/u/${USER}

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/eumaka/sPHENIX/install

echo "running root.exe -q -b Fun4All_EP1.C"
echo "producing recentering histograms"
root.exe -q -b 'Fun4All_EP1.C(1,"dst_calo.list")';

echo "running root.exe -q -b Fun4All_EP2.C"
echo "producing shifting histograms"
root.exe -q -b 'Fun4All_EP2.C(1,"dst_calo.list")';

echo "running root.exe -q -b Fun4All_EP3.C"
echo "reconstruction pass"
root.exe -q -b 'Fun4All_EP3.C(1,"dst_calo.list")';

