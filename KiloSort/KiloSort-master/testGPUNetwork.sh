
me="$(whoami)"
horus=$(ssh $me@horus  "bash /home/buzadmin/scripts/testGPUsforMATLAB.sh")
isi=$(ssh buzadmin@hyperion  "bash testGPUsforMATLAB.sh")
buzmonster=$(bash /home/buzadmin/Scripts/testGPUsforMATLAB.sh)

id -u -n
