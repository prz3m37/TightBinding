echo "############################################################################ TIGHT_BINDING_CALCULATIONS  ############################################################################"
echo
now_start=$(date)
echo "_____________Program for 5 parametrizations has started at : $now_start"
echo
python3 "./main.py" parametrization
echo
now_end=$(date)
echo "_____________Program for 5 parametrizations has ended at : $now_end"
