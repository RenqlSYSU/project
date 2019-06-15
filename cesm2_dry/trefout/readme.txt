this is An example of source code modifications to output the relaxation temperature profile as a history field

Here, the original modules are named filename-ORIG and the modified modules are named filename. Users can compare these two files to identify the modifications that have been made and should then translate these modifications into the appropriate module versions within their own CESM version. 

The result of these modifications are that an additional history field containing the equilibrium temperature profile (TREF) is written to the output when specified as an output field in $CASEDIR/user_nl_cam