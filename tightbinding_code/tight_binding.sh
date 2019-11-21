echo "############################################################################ TIGHT_BINDING_CALCULATIONS_FOR ############################################################################"
echo
now_start=$(date)
echo "_____________Program for 4 parametrizations has started at : $now_start"
echo

python3 "./main.py" parametrization9
#sleep 60
echo
python3 "./main.py" parametrization10
#sleep 60
echo
python3 "./main.py" parametrization11
#sleep 60
echo
python3 "./main.py" parametrization12
#sleep 60
#echo
#python3 "./main.py" parametrization5
#sleep 60
#echo
#python3 "./main.py" parametrization6
#sleep 60
#echo
#python3 "./main.py" parametrization7
#sleep 60
#echo
#python3 "./main.py" parametrization8
#echo
now_end=$(date)
echo "_____________Program for 8 parametrizations has ended at : $now_end"
